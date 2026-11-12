#!/usr/bin/env python3
"""Matched public-R method benchmark; prepare, timings and analyze separately.

ROOT is an isolated artifact directory. LIBRARY is the freshly validated
cleanup installation, with package-source-manifest.json in its parent.
Existing outputs are never overwritten by prepare or timings.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import random
import re
import statistics
import subprocess
import time

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
METHODS = ('legendre', 'meissel', 'lehmer')
SEED = 2026111101


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, data):
    Path(path).write_text(json.dumps(data, indent=2) + '\n')


def invoke(command, log=None, timeout=300, build=False):
    result = subprocess.run(list(map(str, command)), capture_output=True,
                            text=True, timeout=timeout)
    if log is not None:
        log.write_text(result.stdout + result.stderr)
    # R CMD SHLIB reports compiler/SDK identification on stderr even on a
    # clean build. Preserve it while still rejecting compiler diagnostics.
    bad_stderr = result.stderr and not build
    build_diagnostic = build and re.search(r'\b(?:warning|error):', result.stdout + result.stderr, re.I)
    if result.returncode or bad_stderr or build_diagnostic:
        raise RuntimeError(f'{command}\nexit={result.returncode}\n{result.stdout}\n{result.stderr}')
    return result.stdout


def candidate(n, threads):
    if n < 100000:
        return 'lookup/sieve'
    if n < 10000000000:
        return 'legendre'
    if n >= 1000000000000 and threads > 1:
        return 'lehmer'
    return 'meissel'


def prepare(root, library):
    out = root / 'method-benchmark'
    if out.exists():
        raise RuntimeError('refusing existing method-benchmark directory')
    out.mkdir()
    manifest = json.loads((library.parent / 'package-source-manifest.json').read_text())
    source_hashes = {}
    for name, digest in manifest.items():
        if name.startswith(('src/', 'inst/include/', 'R/')):
            if sha(REPO / name) != digest:
                raise RuntimeError('installed benchmark source differs: ' + name)
            source_hashes[name] = digest
    dll = library / 'RcppAlgos/libs/RcppAlgos.so'
    if not dll.is_file():
        raise RuntimeError('expected isolated installed DLL missing')
    timer_source = out / 'timer.cpp'
    timer_source.write_bytes((HERE / 'timer.cpp').read_bytes())
    timer_dll = out / 'timer.so'
    invoke(['R', 'CMD', 'SHLIB', timer_source, '-o', timer_dll], out / 'timer-build.log', build=True)
    cases = []
    for n in (10**8, 10**9, 3*10**9, 10**10, 3*10**10, 10**11,
              3*10**11, 10**12, 3*10**12, 10**13):
        for threads in (1, 2, 4, 8):
            cases.append(dict(id=f'n{n}-t{threads}', n=n, threads=threads))
    for threads in (1, 4, 8):
        cases.append(dict(id=f'n{10**14}-t{threads}', n=10**14, threads=threads))
    suite = dict(cases=cases, methods=list(METHODS), blocks=3, repeats=4, seed=SEED,
                 candidate='lookup/sieve below 1e5; Legendre below 1e10; Meissel otherwise, '
                           'except Lehmer at n>=1e12 with more than one usable thread.',
                 methodology='One fresh R process per method/case/block; warm wrapper only through small auto sieve; '
                             'one first Phi call and three warmed calls. Complete public R call includes tables and teardown. '
                             'Test-only steady_clock timer with per-process origin; no timer-overhead subtraction.',
                 comparison='Use median of warmed calls per process, then matched method comparisons per block. '
                            'Flag persistent candidate regressions vs historical Legendre only when every block is slower '
                            'and paired median increases exceed both 5% and 0.05ms. Also report all data and fastest-method gaps.',
                 library=str(library), dll=str(dll), dll_sha256=sha(dll),
                 timer_dll=str(timer_dll), timer_sha256=sha(timer_dll),
                 source_sha256=source_hashes,
                 runner_sha256={str(HERE / name): sha(HERE / name) for name in ('call.R', 'benchmark.py', 'timer.cpp')})
    save(out / 'suite.json', suite)
    references = []
    for n in sorted({c['n'] for c in cases}):
        values = []
        for algorithm in ('--gourdon', '--deleglise-rivat'):
            command = ['/opt/homebrew/bin/primecount', str(n), '--threads=1', algorithm]
            text = invoke(command, out / f'reference-{n}-{algorithm[2:]}.log')
            values.append(int(text.strip()))
            references.append(dict(n=n, algorithm=algorithm, result=values[-1], command=command))
        if values[0] != values[1]:
            raise RuntimeError('independent reference algorithms disagree')
    save(out / 'references.json', dict(version=invoke(['/opt/homebrew/bin/primecount', '--version']),
         references=references))
    metadata = invoke(['Rscript', '--vanilla', '-e',
        'cat(R.version.string,"\\n");cat(as.character(packageVersion("microbenchmark")),"\\n");'
        'cat(as.character(packageVersion("cpp11")),"\\n");print(sessionInfo())'])
    (out / 'r-environment.txt').write_text(metadata)
    print(f'Prepared {len(cases)} cases; independent counts agree.', flush=True)


def verify(root):
    out = root / 'method-benchmark'
    suite = json.loads((out / 'suite.json').read_text())
    if sha(suite['dll']) != suite['dll_sha256']:
        raise RuntimeError('installed DLL changed')
    if sha(suite['timer_dll']) != suite['timer_sha256']:
        raise RuntimeError('timer DLL changed')
    for path, digest in suite['runner_sha256'].items():
        if sha(path) != digest:
            raise RuntimeError('benchmark runner changed')
    return out, suite


def timings(root):
    out, suite = verify(root)
    destination = out / 'timings.csv.gz'
    if destination.exists():
        raise RuntimeError('refusing existing timings')
    references = json.loads((out / 'references.json').read_text())['references']
    expected = {r['n']: r['result'] for r in references}
    rng = random.Random(SEED)
    runs, records, writer = [], [], None
    started = time.time()
    with gzip.open(destination, 'wt', newline='') as stream:
        for block in range(suite['blocks']):
            cases = list(suite['cases'])
            rng.shuffle(cases)
            for position, case in enumerate(cases):
                methods = list(METHODS)
                rng.shuffle(methods)
                for order, method in enumerate(methods):
                    command = ['Rscript', '--vanilla', HERE / 'call.R', suite['library'],
                               case['n'], case['threads'], method, expected[case['n']], suite['timer_dll']]
                    result = invoke(command)
                    rows = list(csv.DictReader(io.StringIO(result)))
                    if len(rows) != suite['repeats']:
                        raise RuntimeError('unexpected result count: ' + result)
                    process = len(runs)
                    for repeat, row in enumerate(rows):
                        if (row['method'] != method or int(float(row['n'])) != case['n'] or
                            int(row['threads']) != case['threads'] or int(row['repeat_index']) != repeat or
                            int(float(row['result'])) != expected[case['n']] or float(row['total_ms']) <= 0):
                            raise RuntimeError('incorrect benchmark result: ' + str(row))
                        row.update(case_id=case['id'], block=block, order=order, process=process)
                        if writer is None:
                            writer = csv.DictWriter(stream, fieldnames=list(row))
                            writer.writeheader()
                        writer.writerow(row)
                        records.append(row)
                    runs.append(dict(case_id=case['id'], method=method, block=block,
                                     order=order, process=process, command=list(map(str, command))))
                stream.flush()
                if (position + 1) % 10 == 0 or position + 1 == len(cases):
                    print(f'block {block+1}/{suite["blocks"]} cases {position+1}/{len(cases)} '
                          f'elapsed={time.time()-started:.1f}s', flush=True)
    verify(root)
    save(out / 'timings-metadata.json', dict(rows=len(records), processes=len(runs), runs=runs,
         sha256=sha(destination), suite_sha256=sha(out / 'suite.json'), elapsed_seconds=time.time()-started))
    print(f'PASS {len(records)} matching counts from {len(runs)} fresh processes', flush=True)


def analyze(root):
    out, suite = verify(root)
    metadata = json.loads((out / 'timings-metadata.json').read_text())
    if sha(out / 'timings.csv.gz') != metadata['sha256']:
        raise RuntimeError('timing data changed')
    with gzip.open(out / 'timings.csv.gz', 'rt') as stream:
        rows = list(csv.DictReader(stream))
    cases = {c['id']: c for c in suite['cases']}
    groups, seen = {}, set()
    for row in rows:
        key = row['case_id'], row['method'], int(row['block']), int(row['repeat_index'])
        if key in seen:
            raise RuntimeError('duplicate timing')
        seen.add(key)
        state = 'first' if key[-1] == 0 else 'warm'
        groups.setdefault(key[:3] + (state,), []).append(float(row['total_ms']))
    wanted = {(cid, m, b, r) for cid in cases for m in METHODS for b in range(3) for r in range(4)}
    if seen != wanted:
        raise RuntimeError('incomplete matrix')
    medians = {k: statistics.median(v) for k,v in groups.items()}
    summaries, comparisons, flags = [], [], []
    for cid, case in cases.items():
        selected = candidate(case['n'], case['threads'])
        for state in ('first', 'warm'):
            values = {m: [medians[cid,m,b,state] for b in range(3)] for m in METHODS}
            best = min(METHODS, key=lambda m: statistics.median(values[m]))
            summaries.append(dict(case_id=cid,n=case['n'],threads=case['threads'],state=state,
                candidate=selected, fastest_median_method=best,
                **{m+'_ms':statistics.median(values[m]) for m in METHODS}))
            for method, control in ((selected,'legendre'),(selected,best),('lehmer','meissel'),('meissel','legendre'),('lehmer','legendre')):
                deltas = [a-b for a,b in zip(values[method],values[control])]
                ratios = [100*(a/b-1) for a,b in zip(values[method],values[control])]
                item = dict(case_id=cid,n=case['n'],threads=case['threads'],state=state,
                    method=method,control=control,paired_delta_ms=statistics.median(deltas),
                    paired_percent=statistics.median(ratios),all_slower=all(d>0 for d in deltas),
                    all_faster=all(d<0 for d in deltas),delta_min_ms=min(deltas),delta_max_ms=max(deltas))
                comparisons.append(item)
                if (method==selected and control=='legendre' and item['all_slower'] and
                    item['paired_delta_ms']>.05 and item['paired_percent']>5):
                    if item not in flags: flags.append(item)
    for name, data in (('summary.csv',summaries),('comparisons.csv',comparisons)):
        with (out/name).open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)
    save(out/'analysis.json',dict(rows=len(rows),processes=metadata['processes'],
        candidate_vs_legendre_flags=flags,summary=summaries,comparisons=comparisons,
        note='Three matched blocks; medians of fresh-process first calls or per-process medians of three warmed calls. '
             'First and warm are separate objectives. No sample removal. This is one host, not universal crossover proof.'))
    print('PASS analysis; candidate-vs-Legendre persistent regression flags:',len(flags),flush=True)
    for flag in flags: print(json.dumps(flag),flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=('prepare','timings','analyze'))
    parser.add_argument('root',type=Path)
    parser.add_argument('--library',type=Path)
    args=parser.parse_args()
    root=args.root.resolve()
    if args.mode=='prepare':
        if args.library is None: parser.error('prepare requires --library')
        prepare(root,args.library.resolve())
    elif args.mode=='timings': timings(root)
    else: analyze(root)
