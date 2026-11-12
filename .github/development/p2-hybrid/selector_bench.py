#!/usr/bin/env python3
"""Frozen selector-v1 experiment; all artifacts and executables stay outside the repo.

Build, validate, generate references, and collect separate memory probes before
running either timing matrix. Do not retune the policy using these results.
"""
import argparse
import csv
import hashlib
import io
import json
import platform
import random
import shutil
import statistics
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
FIELDS = ('setup_ms', 'phi_ms', 'p2_ms', 'p3_ms', 'teardown_ms', 'total_ms')
TRACE = ('a', 'table_limit', 'selected_q', 'selected_k', 'selection_reason',
         'largest_quotient', 'non_table_queries', 'phi_queries', 'stream_queries',
         'direct_queries', 'phi_worker_budget', 'phi_workers', 'stream_worker_budget',
         'stream_workers', 'seed_queries', 'phi_seed_queries')


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def save_csv(path, rows):
    if not rows:
        raise RuntimeError('no rows to save')
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def invoke(command, timeout=900):
    result = subprocess.run(command, text=True, capture_output=True, timeout=timeout)
    if result.returncode or result.stderr:
        raise RuntimeError(f'{command}\nexit={result.returncode}\n{result.stdout}\n{result.stderr}')
    return result.stdout


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_suite(root):
    suite = json.loads((root / 'suite.json').read_text())
    if digest(HERE / 'P2Selector.h') != suite['policy_sha256']:
        raise RuntimeError('frozen policy changed; create a separate version/experiment')
    return suite


def definitions(variant):
    if variant == 'baseline':
        return ['-DP2_SELECTOR_BASELINE']
    if variant == 'selector':
        return ['-DP2_HYBRID_SELECTOR']
    return ['-DP2_HYBRID_MULTIPLIER=' + variant[len('hybrid'):]]


def build(root, suite):
    source = root / 'source'
    directory = root / 'timing-build'
    directory.mkdir(exist_ok=True)
    local = ('selector_profile.cpp', 'selector_cases.h', 'P2Hybrid.cpp',
             'P2Selector.h', 'hybrid_metrics.h')
    for name in local:
        shutil.copy2(HERE / name, directory / name)
    flags = ['clang++', '-std=c++11', '-O3', '-g', '-fno-omit-frame-pointer',
             '-pthread', '-Wall', '-Wextra', '-Wpedantic',
             '-I' + str(source / 'inst/include'), '-I' + str(directory)]
    common, commands = [], []
    for name in ('Phi.cpp', 'PhiParallel.cpp', 'P3.cpp', 'Eratosthenes.cpp'):
        target = directory / (Path(name).stem + '.o')
        command = flags + ['-c', str(source / 'src' / name), '-o', str(target)]
        invoke(command)
        commands.append(command)
        common.append(target)
    variants = suite['variants'] + ['selector-memory']
    for variant in variants:
        memory = variant == 'selector-memory'
        defines = definitions('selector' if memory else variant)
        if memory:
            defines += ['-DP2_HYBRID_MEMORY']
        objects = []
        units = [directory / 'selector_profile.cpp',
                 source / 'src/P2.cpp' if variant == 'baseline' else directory / 'P2Hybrid.cpp']
        for path in units:
            target = directory / (variant + '-' + path.stem + '.o')
            command = flags + defines + ['-c', str(path), '-o', str(target)]
            invoke(command)
            commands.append(command)
            objects.append(target)
        command = ['clang++', '-pthread'] + list(map(str, common + objects)) + [
            '-o', str(root / ('profile-' + variant))]
        invoke(command)
        commands.append(command)
        print('Built ' + variant, flush=True)
    paths = [p for p in source.rglob('*') if p.is_file()]
    paths += [directory / name for name in local] + [root / 'suite.json']
    save_json(root / 'build-commands.json', commands)
    save_json(root / 'build-metadata.json', dict(
        commit=(root / 'commit.txt').read_text().strip(), platform=platform.platform(),
        compiler=invoke(['clang++', '--version']), variants=variants,
        input_sha256={str(p.relative_to(root)): digest(p) for p in paths},
        binary_sha256={v: digest(root / ('profile-' + v)) for v in variants},
        timing='phase-sum-excludes-diagnostics; selection inside P2, trace after P2',
        policy_sha256=suite['policy_sha256']))


def verify_build(root):
    metadata = json.loads((root / 'build-metadata.json').read_text())
    for name, checksum in metadata['input_sha256'].items():
        if digest(root / name) != checksum:
            raise RuntimeError('build input changed: ' + name)
    for variant, checksum in metadata['binary_sha256'].items():
        if digest(root / ('profile-' + variant)) != checksum:
            raise RuntimeError('build executable changed: ' + variant)


def command(root, variant, case, repeats, expected, mode=None):
    return [str(root / ('profile-' + variant)), mode or case['mode'], case['cutoff'],
            str(case['x']), str(case['threads']), str(repeats), str(expected),
            str(case['table_multiplier'])]


def parse(output, case, repeats, mode, verification):
    rows = list(csv.DictReader(io.StringIO(output)))
    if len(rows) != repeats:
        raise RuntimeError('unexpected row count')
    for repeat, row in enumerate(rows):
        if (row['mode'], row['method'], int(row['x']), int(row['threads']),
            int(row['repeat']), int(row['table_multiplier']), row['verification']) != (
                mode, case['cutoff'], case['x'], case['threads'], repeat,
                case['table_multiplier'], verification):
            raise RuntimeError('harness configuration/provenance mismatch')
    return rows


def references(root, suite):
    verify_build(root)
    executable = shutil.which('primecount')
    if executable is None:
        raise RuntimeError('existing primecount CLI required for independent complete references')
    records, known = [], {}
    anchors = {10**12: 37607912018, 10**13: 346065536839, 10**14: 3204941750802}
    for x in sorted({c['x'] for c in suite['cases'] if c['mode'] == 'phases'}):
        values = []
        for algorithm in ('--gourdon', '--deleglise-rivat'):
            cmd = [executable, str(x), algorithm, '--threads=1']
            output = invoke(cmd).strip()
            values.append(int(output))
            records.append(dict(command=cmd, result=int(output)))
        if values[0] != values[1] or (x in anchors and values[0] != anchors[x]):
            raise RuntimeError('independent prime counts disagree')
        known[x] = values[0]
    raw = root / 'reference-raw'
    raw.mkdir(exist_ok=True)
    expectations, partials = [], {}
    for case in suite['cases']:
        if case['mode'] == 'phases':
            expectations.append(dict(id=case['id'], expected=known[case['x']],
                                     provenance='primecount-gourdon-and-deleglise-rivat'))
            continue
        key = (case['x'], case['cutoff'], case['table_multiplier'])
        if key not in partials:
            refcase = dict(case, threads=1)
            output = invoke(command(root, 'baseline', refcase, 1, -1, 'reference'))
            (raw / (case['id'] + '-reference.csv')).write_text(output)
            row = parse(output, refcase, 1, 'reference', 'unchecked-production-reference')[0]
            expected = int(row['result'])
            provenance = 'production-P2-differential'
            if case['x'] <= 10**8:
                oracle = invoke(command(root, 'baseline', refcase, 1, -1, 'oracle'))
                (raw / (case['id'] + '-oracle.csv')).write_text(oracle)
                oracle_row = parse(oracle, refcase, 1, 'oracle', 'independent-dense-pi')[0]
                if int(oracle_row['result']) != expected:
                    raise RuntimeError('independent P2 oracle mismatch')
                provenance = 'independent-dense-pi-and-production-P2'
            partials[key] = (expected, provenance, int(row['a']))
        expected, provenance, a = partials[key]
        expectations.append(dict(id=case['id'], expected=expected, provenance=provenance, a=a))
    save_json(root / 'expectations.json', expectations)
    save_json(root / 'references-metadata.json', dict(
        primecount_path=executable, primecount_sha256=digest(Path(executable).resolve()),
        primecount_version=invoke([executable, '--version']), independent_counts=records,
        unique_p2_references=len(partials), cases=len(expectations)))
    print(f'Prepared {len(expectations)} expectations; {len(partials)} unique P2 references', flush=True)


def memory(root, suite):
    verify_build(root)
    expected = {r['id']: r['expected'] for r in json.loads((root / 'expectations.json').read_text())}
    rows = []
    raw = root / 'memory-raw'
    raw.mkdir(exist_ok=True)
    for i, case in enumerate(suite['cases']):
        output = invoke(command(root, 'selector-memory', case, 2, expected[case['id']]))
        (raw / (case['id'] + '.csv')).write_text(output)
        parsed = parse(output, case, 2, case['mode'], 'expected-match')
        for row in parsed:
            if row['measurement'] != 'memory-probe' or int(row['result']) != expected[case['id']]:
                raise RuntimeError('invalid memory probe result')
            row['case_id'] = case['id']
            rows.append(row)
        for field in TRACE + ('sieve_capacity_bytes', 'next_capacity_bytes',
                              'phi_sum_bytes', 'stream_sum_bytes', 'conservative_scratch_bytes'):
            if parsed[0][field] != parsed[1][field]:
                raise RuntimeError('first/warm deterministic memory/work fields differ')
        if (i + 1) % 20 == 0:
            print(f'Memory cases {i + 1}/{len(suite["cases"])}', flush=True)
    save_csv(root / 'memory.csv', rows)
    save_json(root / 'memory-metadata.json', dict(verified_rows=len(rows),
        repeats=2, instrumented='selector only; separate executable, no timing',
        scratch='conservative capacity payload; excludes tables/cache/stacks/ThreadBatch/allocator/RSS'))
    print(f'All {len(rows)} memory/work results matched', flush=True)


def summarize(rows, variants):
    result = []
    for case_id in sorted({r['case_id'] for r in rows}):
        group = [r for r in rows if r['case_id'] == case_id]
        processes = sorted({r['process'] for r in group})
        for state in ('first', 'warm'):
            medians = {}
            for variant in variants:
                for process in processes:
                    selected = [r for r in group if r['variant'] == variant and r['process'] == process
                                and ((int(r['repeat']) == 0) == (state == 'first'))]
                    medians[(variant, process)] = {
                        f: statistics.median(float(r[f]) for r in selected) for f in FIELDS}
            metric = 'total_ms' if group[0]['mode'] == 'phases' else 'p2_ms'
            fixed = [v for v in variants if v != 'selector']
            best = min(fixed, key=lambda v: statistics.median(medians[(v, p)][metric] for p in processes))
            trace = next(r for r in group if r['variant'] == 'selector')
            chosen = 'hybrid' + trace['selected_k'] if int(trace['selected_k']) else 'baseline'
            for variant in variants:
                original = next(r for r in group if r['variant'] == variant)
                row = {k: original[k] for k in ('case_id', 'mode', 'method', 'x', 'threads', 'table_multiplier') + TRACE}
                row.update(state=state, variant=variant, processes=len(processes),
                           metric=metric, best_fixed=best, selector_chosen=chosen)
                for field in FIELDS:
                    values = [medians[(variant, p)][field] for p in processes]
                    row[field] = statistics.median(values)
                    row[field.replace('_ms', '_min_ms')] = min(values)
                    row[field.replace('_ms', '_max_ms')] = max(values)
                for label, reference in (('baseline', 'baseline'), ('best_fixed', best), ('chosen', chosen)):
                    denominators = [medians[(reference, p)][metric] for p in processes]
                    ratios = [medians[(variant, p)][metric] / medians[(reference, p)][metric]
                              for p in processes] if min(denominators) > 0 else []
                    row[label + '_paired_change_percent'] = 100 * (statistics.median(ratios) - 1) if ratios else ''
                    row[label + '_paired_min_change_percent'] = 100 * (min(ratios) - 1) if ratios else ''
                    row[label + '_paired_max_change_percent'] = 100 * (max(ratios) - 1) if ratios else ''
                    deltas = [medians[(variant, p)][metric] - medians[(reference, p)][metric] for p in processes]
                    row[label + '_paired_delta_ms'] = statistics.median(deltas)
                result.append(row)
    return result


def benchmark(root, suite, mode):
    verify_build(root)
    if (root / (mode + '.csv')).exists():
        raise RuntimeError('timing rows already exist; use a fresh artifact to preserve this run')
    expected = {r['id']: r for r in json.loads((root / 'expectations.json').read_text())}
    cases = [c for c in suite['cases'] if c['mode'] == mode]
    blocks = [(c, p) for c in cases for p in range(suite['processes'])]
    rng = random.Random(suite['seeds'][mode])
    rng.shuffle(blocks)
    raw = root / 'timing-raw'
    raw.mkdir(exist_ok=True)
    rows, traces, partials = [], {}, {}
    for block, (case, process) in enumerate(blocks):
        block_rows = []
        order = list(suite['variants'])
        rng.shuffle(order)
        for position, variant in enumerate(order):
            ref = expected[case['id']]
            output = invoke(command(root, variant, case, suite['repeats'], ref['expected']))
            (raw / f'{case["id"]}-p{process}-{variant}.csv').write_text(output)
            parsed = parse(output, case, suite['repeats'], mode, 'expected-match')
            for row in parsed:
                if row['measurement'] != 'timing' or int(row['result']) != ref['expected']:
                    raise RuntimeError('invalid timed result')
                if any(not row[f] or float(row[f]) < 0 for f in FIELDS):
                    raise RuntimeError('invalid timing field')
                key = (case['id'], variant)
                trace = tuple(row[f] for f in TRACE)
                if traces.setdefault(key, trace) != trace:
                    raise RuntimeError('deterministic decision changed across processes/repeats')
                pkey = (case['x'], case['cutoff'], case['table_multiplier'])
                if partials.setdefault(pkey, int(row['p2'])) != int(row['p2']):
                    raise RuntimeError('P2 values differ across variants/budgets/repeats')
                row.update(case_id=case['id'], variant=variant, process=process,
                           block=block, order=position, reference_provenance=ref['provenance'])
                rows.append(row)
                block_rows.append(row)
        with (root / (mode + '.csv')).open('w' if block == 0 else 'a', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(block_rows[0]))
            if block == 0:
                writer.writeheader()
            writer.writerows(block_rows)
        print(f'{mode}: block {block + 1}/{len(blocks)}; {case["id"]}', flush=True)
    save_csv(root / (mode + '-summary.csv'), summarize(rows, suite['variants']))
    save_json(root / (mode + '-metadata.json'), dict(
        mode=mode, seed=suite['seeds'][mode], configurations=len(cases), blocks=len(blocks),
        processes=suite['processes'], repeats=suite['repeats'], variants=suite['variants'],
        verified_rows=len(rows), policy_sha256=suite['policy_sha256'],
        first=('fresh process; complete Phi phase precedes P2' if mode == 'phases'
               else 'fresh process; isolated P2 has no preceding Phi phase'),
        warm='median of three subsequent calls, then median across processes',
        best_fixed='hindsight minimum of baseline Phi-only and fixed8/16/32/64 medians per configuration/state; optimistic comparison',
        comparisons='paired process ratios/deltas vs baseline, hindsight best fixed, and selector-chosen fixed'))
    print(f'{mode}: all {len(rows)} results matched', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('build', 'references', 'memory', 'phases', 'p2'))
    parser.add_argument('artifact_dir', type=Path)
    args = parser.parse_args()
    root = args.artifact_dir.resolve()
    if root == REPO or REPO in root.parents:
        parser.error('artifacts must be outside the repository')
    suite = load_suite(root)
    if args.action == 'build': build(root, suite)
    elif args.action == 'references': references(root, suite)
    elif args.action == 'memory': memory(root, suite)
    else: benchmark(root, suite, args.action)


if __name__ == '__main__':
    main()
