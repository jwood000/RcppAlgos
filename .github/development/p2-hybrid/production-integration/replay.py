#!/usr/bin/env python3
"""Small, frozen parity replay: integrated P2 versus the approved candidate.

Run prepare while builds/tests are allowed, then timings only after all other
CPU-intensive work finishes. ROOT is the production-integration artifact root.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import random
import statistics
import subprocess
import time

HERE = Path(__file__).resolve().parent
EXPERIMENT = HERE.parent
REPO = HERE.parents[3]
CASE_IDS = (
    'p2-zero-80000000000-t8-L6over5',
    'p2-meissel-10000000000-t4-L2over1',
    'p2-lehmer-10000000000-t8-L1over1',
    'p2-meissel-60000000000-t3-L1over1',
    'p2-tail2048-3000000000-t1-L1over1',
    'p2-meissel-1000000000000-t8-L1over1',
    'phases-meissel-1000000000000-t1-L1over1',
    'phases-lehmer-1000000000000-t1-L1over1',
    'phases-meissel-1000000000000-t8-L1over1',
    'phases-lehmer-1000000000000-t8-L1over1',
    'phases-meissel-100000000000000-t8-L1over1',
    'phases-lehmer-100000000000000-t8-L1over1',
)
FIELDS = ('setup_ms', 'phi_ms', 'p2_ms', 'p3_ms', 'teardown_ms', 'total_ms')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n')


def require(value, message):
    if not value:
        raise RuntimeError(message)


def invoke(command, log=None):
    result = subprocess.run(list(map(str, command)), capture_output=True, text=True, timeout=300)
    if log is not None:
        log.write_text(result.stdout + result.stderr)
    if result.returncode or result.stderr:
        raise RuntimeError(f'{command}\nexit={result.returncode}\n{result.stdout}\n{result.stderr}')
    return result.stdout


def prepare(root):
    out = root / 'parity'
    require(not out.exists(), 'refusing existing parity directory')
    out.mkdir()
    source = root / 'native-validation/source'
    old = EXPERIMENT / 'results-consolidated'
    generated = old / 'generated'
    old_suite = json.loads((old / 'suite.json').read_text())
    old_expected = {r['case_id']: r for r in json.loads((old / 'expectations.json').read_text())}
    available = {c['case_id']: c for c in old_suite['cases']}
    cases = []
    for cid in CASE_IDS:
        case = dict(available[cid])
        case.update(blocks=3, repeats=4, variants=['approved', 'integrated'],
                    expected=old_expected[cid]['expected'], provenance=old_expected[cid]['provenance'])
        cases.append(case)
    suite = dict(cases=cases, seed=2026110721, variants=['approved', 'integrated'],
                 design='Review-approved twelve-case integration parity check; no policy tuning.',
                 diagnostic='Both builds use the frozen candidate to describe decisions after timing; '
                            'production chooser equivalence is checked separately against that candidate.',
                 screening='Report primary first/warm >5% AND >0.05ms paired median increases with every block slower; '
                           'also retain all samples and smaller changes. No samples discarded.',
                 runner_sha256=sha(Path(__file__)), source_sha256=sha(source / 'src/P2.cpp'))
    save(out / 'suite.json', suite)
    flags = ['clang++', '-std=c++11', '-O3', '-g', '-fno-omit-frame-pointer', '-pthread',
             '-Wall', '-Wextra', '-Wpedantic', '-I' + str(source / 'inst/include'), '-I' + str(generated)]
    common, commands = [], []
    for name in ('Phi', 'PhiParallel', 'P3', 'Eratosthenes'):
        target = out / (name + '.o')
        command = flags + ['-c', source / ('src/' + name + '.cpp'), '-o', target]
        invoke(command, out / (name + '-build.log'))
        commands.append(list(map(str, command)))
        common.append(target)
    paths = {'approved': generated / 'P2HybridInitial.cpp', 'integrated': source / 'src/P2.cpp'}
    for variant, path in paths.items():
        command = flags + ['-DP2_HYBRID_SELECTOR_CANDIDATE', generated / 'selector_profile.cpp', path]
        command += common + ['-o', out / ('profile-' + variant)]
        invoke(command, out / (variant + '-build.log'))
        commands.append(list(map(str, command)))
    inputs = [p for p in source.rglob('*') if p.is_file()]
    inputs += [p for p in generated.iterdir() if p.suffix in ('.h', '.cpp')]
    inputs += [Path(__file__), out / 'suite.json']
    save(out / 'build.json', dict(commands=commands, inputs={str(p): sha(p) for p in inputs},
         binaries={v: sha(out / ('profile-' + v)) for v in suite['variants']},
         compiler=invoke(['clang++', '--version'])))
    print('PASS parity prepare: 12 cases, 288 planned timed results', flush=True)


def verify(root):
    out = root / 'parity'
    build = json.loads((out / 'build.json').read_text())
    for path, digest in build['inputs'].items():
        require(sha(path) == digest, 'changed replay input: ' + path)
    for variant, digest in build['binaries'].items():
        require(sha(out / ('profile-' + variant)) == digest, 'changed replay binary: ' + variant)
    suite = json.loads((out / 'suite.json').read_text())
    require(sha(REPO / 'src/P2.cpp') == suite['source_sha256'], 'live P2 changed after snapshot')
    return out, suite


def timings(root):
    out, suite = verify(root)
    require((root / 'native-validation/validation.json').is_file(), 'native validation must complete before timings')
    output = out / 'timings.csv.gz'
    require(not output.exists(), 'refusing existing timings')
    rng = random.Random(suite['seed'])
    runs, records, writer = [], [], None
    started = time.time()
    with gzip.open(output, 'wt', newline='') as stream:
        for block in range(3):
            cases = list(suite['cases'])
            rng.shuffle(cases)
            for position, case in enumerate(cases):
                variants = list(suite['variants'])
                rng.shuffle(variants)
                for order, variant in enumerate(variants):
                    command = [out / ('profile-' + variant), case['mode'], case['cutoff'], case['x'],
                               case['threads'], case['repeats'], case['expected'], case['numerator'], case['denominator']]
                    rows = list(csv.DictReader(io.StringIO(invoke(command))))
                    require(len(rows) == case['repeats'], 'unexpected call count')
                    pid = len(runs)
                    for repeat, row in enumerate(rows):
                        require(int(row['result']) == case['expected'] and row['verification'] == 'expected-match', 'wrong count')
                        require(row['build_kind'] == 'candidate' and row['measurement'] == 'timing', 'unexpected harness mode')
                        require((row['mode'], row['method'], int(row['x']), int(row['threads']), int(row['repeat']),
                                 int(row['table_numerator']), int(row['table_denominator'])) ==
                                (case['mode'], case['cutoff'], case['x'], case['threads'], repeat, case['numerator'], case['denominator']),
                                'configuration mismatch')
                        row.update(case_id=case['case_id'], variant=variant, block=block, order=order,
                                   process=pid, reference_provenance=case['provenance'])
                        if writer is None:
                            writer = csv.DictWriter(stream, fieldnames=list(row))
                            writer.writeheader()
                        writer.writerow(row)
                        records.append(row)
                    runs.append(dict(case_id=case['case_id'], variant=variant, block=block,
                                     order=order, process=pid, command=list(map(str, command))))
                stream.flush()
                print(f'block {block + 1} case {position + 1}/12 {case["case_id"]} elapsed={time.time()-started:.1f}s', flush=True)
    verify(root)
    save(out / 'timings-metadata.json', dict(rows=len(records), processes=len(runs), runs=runs,
         sha256=sha(output), suite_sha256=sha(out / 'suite.json'), elapsed_seconds=time.time() - started))
    print('PASS parity timings:', len(records), 'results', flush=True)


def analyze(root):
    out, suite = verify(root)
    metadata = json.loads((out / 'timings-metadata.json').read_text())
    require(sha(out / 'timings.csv.gz') == metadata['sha256'], 'timings digest')
    require(sha(out / 'suite.json') == metadata['suite_sha256'], 'suite digest')
    cases = {c['case_id']: c for c in suite['cases']}
    with gzip.open(out / 'timings.csv.gz', 'rt', newline='') as stream:
        rows = list(csv.DictReader(stream))
    expected = {(cid, v, b, r) for cid, c in cases.items() for v in suite['variants'] for b in range(3) for r in range(4)}
    seen, groups, geometry, processes = set(), {}, {}, {}
    for row in rows:
        cid, variant, block, repeat = row['case_id'], row['variant'], int(row['block']), int(row['repeat'])
        key = cid, variant, block, repeat
        require(key in expected and key not in seen, 'unexpected/duplicate timed row')
        seen.add(key)
        require(int(row['result']) == cases[cid]['expected'], 'incorrect timed result')
        identity = cid, variant, block
        require(processes.setdefault(int(row['process']), identity) == identity, 'process reused')
        layout = tuple(row[k] for k in ('a', 'selected_q', 'selected_k', 'selection_reason', 'phi_queries',
                                       'stream_queries', 'phi_worker_budget', 'stream_worker_budget'))
        require(geometry.setdefault(cid, layout) == layout, 'partition differs across variants or calls')
        require(abs(sum(float(row[f]) for f in FIELDS[:-1]) - float(row['total_ms'])) < 1e-7, 'phase sum')
        state = 'first' if repeat == 0 else 'warm'
        groups.setdefault(identity + (state,), []).append(row)
    require(seen == expected and len(rows) == metadata['rows'] == 288, 'incomplete timing matrix')
    require(len(processes) == metadata['processes'] == 72, 'process count')
    summaries, flags = [], []
    for cid, case in cases.items():
        for state in ('first', 'warm'):
            for field in FIELDS:
                values = {}
                for variant in suite['variants']:
                    values[variant] = [statistics.median(float(row[field]) for row in groups[cid, variant, block, state]) for block in range(3)]
                deltas = [p - c for p, c in zip(values['integrated'], values['approved'])]
                ratios = [100 * (p / c - 1) for p, c in zip(values['integrated'], values['approved']) if c > 0]
                item = dict(case_id=cid, state=state, field=field,
                            approved_median_ms=statistics.median(values['approved']),
                            integrated_median_ms=statistics.median(values['integrated']),
                            paired_delta_ms=statistics.median(deltas), delta_min_ms=min(deltas), delta_max_ms=max(deltas),
                            paired_percent=statistics.median(ratios) if ratios else None,
                            percent_min=min(ratios) if ratios else None, percent_max=max(ratios) if ratios else None,
                            all_slower=all(d > 0 for d in deltas), all_faster=all(d < 0 for d in deltas))
                summaries.append(item)
                primary = 'p2_ms' if case['mode'] == 'p2' else 'total_ms'
                if field == primary and item['all_slower'] and item['paired_delta_ms'] > 0.05 and item['paired_percent'] > 5:
                    flags.append(item)
    with (out / 'summary.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)
    save(out / 'analysis.json', dict(rows=len(rows), processes=len(processes), source_sha256=suite['source_sha256'],
         primary_flags=flags, comparisons=summaries,
         note='Paired medians of fresh-process first calls and median of three warmed calls; three randomized blocks. '
              'Integration parity evidence, not a new selector fit or a statistical equivalence proof.'))
    print('PASS parity analysis;', len(flags), 'primary regression flags', flush=True)
    for item in summaries:
        primary = 'p2_ms' if cases[item['case_id']]['mode'] == 'p2' else 'total_ms'
        if item['field'] == primary:
            print(item['case_id'], item['state'], f'{item["paired_percent"]:+.2f}%', f'{item["paired_delta_ms"]:+.6f} ms',
                  'all slower' if item['all_slower'] else 'all faster' if item['all_faster'] else 'mixed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'timings', 'analyze'))
    parser.add_argument('root', type=Path)
    args = parser.parse_args()
    globals()[args.action](args.root.resolve())


if __name__ == '__main__':
    main()
