#!/usr/bin/env python3
"""Frozen, temporary common-K1 selector comparison; no production edits.

The plan, candidate, native inputs, runner and executables are hashed before
measurement. Focus and held-out cases are frozen together. Each measurement
process supplies one first call followed by warmed calls; no diagnostic build
or reference calculation runs concurrently with the timing loops.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
import platform
import random
import shutil
import subprocess
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
VARIANTS = ['baseline', 'phi1', 'v1', 'v2', 'v3',
            'fixed8', 'fixed16', 'fixed32', 'fixed64', 'candidate']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n')


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def invoke(command, timeout=900):
    process = subprocess.run(list(map(str, command)), capture_output=True,
                             text=True, timeout=timeout)
    if process.returncode or process.stderr:
        raise RuntimeError(f'{command}\nexit={process.returncode}\n'
                           f'{process.stdout}\n{process.stderr}')
    return process.stdout


def kind(variant):
    if variant == 'baseline':
        return 'baseline'
    if variant == 'phi1':
        return 'phi-initial-1'
    if variant.startswith('fixed'):
        return 'fixed'
    return 'candidate' if variant == 'candidate' else 'selector-' + variant


def definitions(variant):
    if variant == 'baseline':
        return ['-DP2_SELECTOR_BASELINE']
    if variant == 'phi1':
        return ['-DP2_SELECTOR_BASELINE', '-DP2_INITIAL_CHUNK=1']
    if variant.startswith('fixed'):
        return ['-DP2_HYBRID_MULTIPLIER=' + variant[5:]]
    suffix = {'v1': '', 'v2': '_V2', 'v3': '_V3',
              'candidate': '_CANDIDATE'}[variant]
    return ['-DP2_HYBRID_SELECTOR' + suffix]


def freeze(root, plan):
    require(plan is not None, 'freeze requires --plan')
    require(not (root / 'suite.json').exists(), 'suite already frozen')
    request = json.loads(plan.read_text())
    cases = []
    for original in request['cases']:
        case = dict(original)
        for key, default in [('numerator', 1), ('denominator', 1),
                             ('blocks', 3 if case['x'] >= 10**13 else 5),
                             ('repeats', 4), ('variants', VARIANTS)]:
            case.setdefault(key, default)
        case.setdefault('case_id', '{mode}-{cutoff}-{x}-t{threads}-L{numerator}over{denominator}'.format(**case))
        require(case['cohort'] in ('focus', 'holdout'), 'invalid cohort')
        require(case['mode'] in ('p2', 'phases'), 'invalid mode')
        require(case['blocks'] >= 3 and case['repeats'] >= 2, 'insufficient replication')
        require(0 < case['denominator'] <= case['numerator'], 'invalid table ratio')
        require(len(case['variants']) == len(set(case['variants'])), 'duplicate variants')
        require(set(case['variants']) <= set(VARIANTS), 'unknown variants')
        require({'baseline', 'phi1', 'candidate'} <= set(case['variants']), 'missing controls')
        cases.append(case)
    require(len(cases) == len({c['case_id'] for c in cases}), 'duplicate case id')
    require({c['cohort'] for c in cases} == {'focus', 'holdout'}, 'both cohorts required')
    variants = [v for v in VARIANTS if any(v in c['variants'] for c in cases)]
    seeds = request.get('seeds', dict(focus=2026110711, holdout=2026110712))
    require(set(seeds) == {'focus', 'holdout'}, 'one seed per cohort required')
    suite = dict(variants=variants, cases=cases, seeds=seeds,
                 design=request.get('design', 'One common K1 queue for Phi-only and hybrid prefixes; frozen selector controls.'),
                 evaluation=request.get('evaluation', {}),
                 runner_sha256=sha(Path(__file__)), plan_sha256=sha(plan),
                 candidate_sha256=sha(HERE / 'P2SelectorCandidate.h'))
    save(root / 'suite.json', suite)
    for cohort in ('focus', 'holdout'):
        group = [c for c in cases if c['cohort'] == cohort]
        print(f'froze {cohort}: {len(group)} cases, '
              f'{sum(c["blocks"] * c["repeats"] * len(c["variants"]) for c in group)} timed rows', flush=True)


def build(root):
    directory, source = root / 'build', root / 'source'
    suite = json.loads((root / 'suite.json').read_text())
    require(suite['runner_sha256'] == sha(Path(__file__)), 'runner changed after freeze')
    require(suite['candidate_sha256'] == sha(HERE / 'P2SelectorCandidate.h'), 'candidate changed after freeze')
    require(not (root / 'build-metadata.json').exists(), 'build already frozen')
    if not (directory / 'generation.json').exists():
        invoke(['python3', HERE / 'generate.py', source, directory])
    candidate = directory / 'P2SelectorCandidate.h'
    if candidate.exists():
        require(sha(candidate) == suite['candidate_sha256'], 'generated candidate mismatch')
    else:
        shutil.copy2(HERE / 'P2SelectorCandidate.h', candidate)
    flags = ['clang++', '-std=c++11', '-O3', '-g', '-fno-omit-frame-pointer',
             '-pthread', '-Wall', '-Wextra', '-Wpedantic',
             '-I' + str(source / 'inst/include'), '-I' + str(directory)]
    commands, common = [], []
    for name in ('Phi', 'PhiParallel', 'P3', 'Eratosthenes'):
        target = directory / (name + '.o')
        command = flags + ['-c', source / ('src/' + name + '.cpp'), '-o', target]
        invoke(command)
        commands.append(list(map(str, command)))
        common.append(target)
    for variant in suite['variants']:
        native = source / 'src/P2.cpp' if variant == 'baseline' else directory / (
            'P2Initial.cpp' if variant == 'phi1' else 'P2HybridInitial.cpp')
        command = flags + definitions(variant) + [directory / 'selector_profile.cpp', native] + common + ['-o', root / ('profile-' + variant)]
        invoke(command)
        commands.append(list(map(str, command)))
        print('built ' + variant, flush=True)
    inputs = [p for p in source.rglob('*') if p.is_file()]
    inputs += [p for p in directory.iterdir() if p.is_file() and p.suffix in ('.cpp', '.h', '.json')]
    inputs += [root / 'suite.json', root / 'commit.txt']
    save(root / 'build-commands.json', commands)
    save(root / 'build-metadata.json', dict(
        compiler=invoke(['clang++', '--version']), platform=platform.platform(),
        commit=(root / 'commit.txt').read_text().strip(),
        inputs={str(p.relative_to(root)): sha(p) for p in inputs},
        binaries={v: sha(root / ('profile-' + v)) for v in suite['variants']},
        runner_sha256=sha(Path(__file__)), candidate_sha256=sha(candidate)))


def verify(root):
    metadata = json.loads((root / 'build-metadata.json').read_text())
    require(sha(Path(__file__)) == metadata['runner_sha256'], 'runner changed after build')
    require(sha(HERE / 'P2SelectorCandidate.h') == metadata['candidate_sha256'], 'candidate source changed after build')
    for name, digest in metadata['inputs'].items():
        require(sha(root / name) == digest, 'input changed: ' + name)
    for variant, digest in metadata['binaries'].items():
        require(sha(root / ('profile-' + variant)) == digest, 'binary changed: ' + variant)


def command(root, variant, case, repeats, expected, mode=None):
    return [root / ('profile-' + variant), mode or case['mode'], case['cutoff'],
            case['x'], case['threads'], repeats, expected,
            case['numerator'], case['denominator']]


def references(root):
    verify(root)
    require(not (root / 'expectations.json').exists(), 'references already frozen')
    suite = json.loads((root / 'suite.json').read_text())
    expectations, known, p2, references = [], {}, {}, []
    terms_path = HERE.parent / 'results-matched-profile/reference-terms.csv.gz'
    with gzip.open(terms_path, 'rt') as stream:
        terms = list(csv.DictReader(stream))
    cli = shutil.which('primecount')
    require(cli is not None, 'primecount executable missing')
    for x in sorted({c['x'] for c in suite['cases'] if c['mode'] == 'phases'}):
        counts = []
        for method in ('--gourdon', '--deleglise-rivat'):
            cmd = [cli, x, method, '--threads=1']
            value = int(invoke(cmd).strip())
            counts.append(value)
            references.append(dict(command=list(map(str, cmd)), result=value))
        require(counts[0] == counts[1], 'independent complete counts disagree')
        known[x] = counts[0]
    for case in suite['cases']:
        if case['mode'] == 'phases':
            value, provenance = known[case['x']], 'primecount-gourdon-and-deleglise-rivat'
        else:
            key = (case['x'], case['cutoff'], case['numerator'], case['denominator'])
            if key not in p2:
                serial = dict(case, threads=1)
                row = list(csv.DictReader(io.StringIO(invoke(command(root, 'baseline', serial, 1, -1, 'reference')))))[0]
                value, provenance = int(row['result']), 'production-P2-differential'
                if case['x'] == 8 * 10**10:
                    independent = sum(int(t['term']) for t in terms if int(t['i']) > int(row['a']))
                    require(value == independent, 'matched independent terms disagree')
                    provenance = 'independent-saved-sieve-and-primecount-terms'
                p2[key] = value, provenance
            value, provenance = p2[key]
        expectations.append(dict(case_id=case['case_id'], expected=value, provenance=provenance))
    save(root / 'expectations.json', expectations)
    save(root / 'references-metadata.json', dict(
        independent_counts=references, unique_p2=len(p2),
        matched_term_sha256=sha(terms_path), cli=cli,
        cli_sha256=sha(Path(cli).resolve()), expectations_sha256=sha(root / 'expectations.json')))
    print(f'prepared {len(expectations)} references', flush=True)


def check_row(row, case, expected, repeat, variant):
    require(int(row['result']) == expected and row['verification'] == 'expected-match', 'result mismatch')
    require((row['mode'], row['method'], int(row['x']), int(row['threads']),
             int(row['repeat']), int(row['table_numerator']), int(row['table_denominator'])) == (
                 case['mode'], case['cutoff'], case['x'], case['threads'], repeat,
                 case['numerator'], case['denominator']), 'row configuration mismatch')
    require(row['build_kind'] == kind(variant), 'build kind mismatch')


def preflight(root):
    verify(root)
    suite = json.loads((root / 'suite.json').read_text())
    records = []
    for cutoff in ('zero', 'meissel', 'lehmer'):
        case = dict(mode='p2', cutoff=cutoff, x=10**8, threads=1,
                    numerator=1, denominator=1)
        oracle = list(csv.DictReader(io.StringIO(invoke(command(root, 'baseline', case, 1, -1, 'oracle')))))[0]
        expected = int(oracle['result'])
        for threads in (1, 2, 4, 8):
            parallel = dict(case, threads=threads)
            for variant in suite['variants']:
                rows = list(csv.DictReader(io.StringIO(invoke(command(root, variant, parallel, 2, expected)))))
                require(len(rows) == 2, 'preflight call count mismatch')
                for index, row in enumerate(rows):
                    check_row(row, parallel, expected, index, variant)
                records.append(dict(cutoff=cutoff, threads=threads, variant=variant,
                                    expected=expected, calls=len(rows)))
    save(root / 'harness-preflight.json', dict(records=records,
                                            calls=sum(r['calls'] for r in records)))
    print(f'harness {sum(r["calls"] for r in records)} independent matches', flush=True)


def timings(root, cohort):
    verify(root)
    suite = json.loads((root / 'suite.json').read_text())
    reference_metadata = json.loads((root / 'references-metadata.json').read_text())
    require(sha(root / 'expectations.json') == reference_metadata['expectations_sha256'], 'references changed')
    references = {r['case_id']: r for r in json.loads((root / 'expectations.json').read_text())}
    selected = [c for c in suite['cases'] if c['cohort'] == cohort]
    output = root / (cohort + '.csv.gz')
    require(not output.exists(), 'refusing to overwrite timings')
    if cohort == 'holdout':
        require((root / 'focus-metadata.json').exists(), 'complete focus before held-out timing')
    rng = random.Random(suite['seeds'][cohort])
    runs, writer, row_count, process = [], None, 0, 0
    started = time.time()
    with gzip.open(output, 'wt', newline='') as stream:
        for block in range(max(c['blocks'] for c in selected)):
            cases = [c for c in selected if block < c['blocks']]
            rng.shuffle(cases)
            for position, case in enumerate(cases):
                variants = list(case['variants'])
                rng.shuffle(variants)
                reference = references[case['case_id']]
                for order, variant in enumerate(variants):
                    cmd = command(root, variant, case, case['repeats'], reference['expected'])
                    rows = list(csv.DictReader(io.StringIO(invoke(cmd))))
                    require(len(rows) == case['repeats'], 'timing call count mismatch')
                    for index, row in enumerate(rows):
                        check_row(row, case, reference['expected'], index, variant)
                        row.update(case_id=case['case_id'], cohort=cohort, variant=variant,
                                   block=block, order=order, process=process,
                                   initial_k=32 if variant == 'baseline' else 1,
                                   reference_provenance=reference['provenance'])
                        if writer is None:
                            writer = csv.DictWriter(stream, fieldnames=list(row))
                            writer.writeheader()
                        writer.writerow(row)
                        row_count += 1
                    runs.append(dict(block=block, case_id=case['case_id'], variant=variant,
                                     process=process, command=list(map(str, cmd))))
                    process += 1
                stream.flush()
                print(f'{cohort} block {block + 1} case {position + 1}/{len(cases)} '
                      f'{case["case_id"]} elapsed={time.time() - started:.1f}s', flush=True)
    verify(root)
    expected_rows = sum(c['blocks'] * c['repeats'] * len(c['variants']) for c in selected)
    require(row_count == expected_rows, 'incomplete timing matrix')
    save(root / (cohort + '-metadata.json'), dict(
        rows=row_count, processes=process, seed=suite['seeds'][cohort], runs=runs,
        sha256=sha(output), elapsed_seconds=time.time() - started,
        suite_sha256=sha(root / 'suite.json'),
        candidate_sha256=sha(HERE / 'P2SelectorCandidate.h')))
    print(f'PASS {cohort}: {row_count} timed results', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('freeze', 'build', 'references', 'preflight',
                                          'timings-focus', 'timings-holdout', 'verify'))
    parser.add_argument('root', type=Path)
    parser.add_argument('--plan', type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    if args.action == 'freeze':
        freeze(root, args.plan)
    elif args.action.startswith('timings-'):
        timings(root, args.action.split('-', 1)[1])
    else:
        globals()[args.action](root)


if __name__ == '__main__':
    main()
