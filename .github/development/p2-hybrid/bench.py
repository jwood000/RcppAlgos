#!/usr/bin/env python3
"""Build and benchmark isolated current-cache and bounded-streaming P2 copies.

Build, validate, collect diagnostics, then time separately. No production file
is changed. Each timing block runs the selected variants in fresh processes,
in randomized order. The default retains the original 2L/4L/8L experiment.
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
PROFILE = HERE.parent / 'p2-profile'
REPO = HERE.parents[2]
VARIANTS = ('baseline', 'hybrid2', 'hybrid4', 'hybrid8')
AVAILABLE_VARIANTS = VARIANTS + ('hybrid16', 'hybrid32', 'hybrid64', 'hybrid128')
COUNTS = {10**10: 455052511, 10**12: 37607912018,
          10**13: 346065536839, 10**14: 3204941750802}
FIELDS = ('setup_ms', 'phi_ms', 'p2_ms', 'p3_ms', 'teardown_ms', 'total_ms')


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def save_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def execute(command, timeout=600):
    result = subprocess.run(command, text=True, capture_output=True, timeout=timeout)
    if result.returncode:
        raise RuntimeError(f'{command}\n{result.stdout}\n{result.stderr}')
    return result


def build(root, variants=VARIANTS):
    source = root / 'source'
    if not source.is_dir():
        raise RuntimeError('create a committed production snapshot in artifact/source first')
    phi = (source / 'src/Phi.cpp').read_text()
    if 'CACHE_MAX_X = 131071' not in phi:
        raise RuntimeError('this experiment requires the adopted doubled cache')
    directory = root / 'timing-build'
    directory.mkdir(exist_ok=True)
    shutil.copy2(PROFILE / 'profile.cpp', directory / 'profile.cpp')
    shutil.copy2(HERE / 'P2Hybrid.cpp', directory / 'P2Hybrid.cpp')
    flags = ['clang++', '-std=c++11', '-O3', '-g', '-fno-omit-frame-pointer',
             '-pthread', '-Wall', '-Wextra', '-Wpedantic',
             '-I' + str(source / 'inst/include')]
    common = []
    commands = []
    inputs = [directory / 'profile.cpp'] + [source / 'src' / name for name in
              ('Phi.cpp', 'PhiParallel.cpp', 'P3.cpp', 'Eratosthenes.cpp')]
    for path in inputs:
        output = directory / (path.stem + '.o')
        command = flags + ['-c', str(path), '-o', str(output)]
        result = execute(command)
        (directory / (path.stem + '.log')).write_text(result.stdout + result.stderr)
        if result.stderr:
            raise RuntimeError('compiler diagnostics: ' + result.stderr)
        common.append(output)
        commands.append(command)
    for variant in variants:
        path = source / 'src/P2.cpp' if variant == 'baseline' else directory / 'P2Hybrid.cpp'
        defines = [] if variant == 'baseline' else [
            '-DP2_HYBRID_MULTIPLIER=' + variant[len('hybrid'):]]
        output = directory / (variant + '.o')
        command = flags + defines + ['-c', str(path), '-o', str(output)]
        result = execute(command)
        (directory / (variant + '.log')).write_text(result.stdout + result.stderr)
        if result.stderr:
            raise RuntimeError('compiler diagnostics: ' + result.stderr)
        commands.append(command)
        command = ['clang++', '-pthread'] + list(map(str, common + [output])) + [
            '-o', str(root / ('profile-' + variant))]
        result = execute(command)
        if result.stderr:
            raise RuntimeError('linker diagnostics: ' + result.stderr)
        commands.append(command)
    save_json(root / 'timing-build-commands.json', commands)
    save_json(root / 'timing-build-metadata.json', dict(
        commit=(root / 'commit.txt').read_text().strip(), platform=platform.platform(),
        compiler=execute(['clang++', '--version']).stdout,
        input_sha256={str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in inputs + [source / 'src/P2.cpp', directory / 'P2Hybrid.cpp']},
        binary_sha256={v: hashlib.sha256((root / ('profile-' + v)).read_bytes()).hexdigest()
                       for v in variants},
        variants=list(variants), cache_max_x=131071,
        timing='uninstrumented kernels, complete-driver phase order; no R adapter'))
    print(f'{len(variants)} timing binaries built without compiler diagnostics', flush=True)


def summarize(rows, variants=VARIANTS):
    summaries = []
    cases = sorted({(r['mode'], r['method'], int(r['x']), int(r['threads'])) for r in rows})
    for mode, method, x, threads in cases:
        selected = [r for r in rows if (r['mode'], r['method'], int(r['x']), int(r['threads']))
                    == (mode, method, x, threads)]
        processes = sorted({r['process'] for r in selected})
        for state in ('first', 'warm'):
            medians = {}
            for variant in variants:
                for process in processes:
                    group = [r for r in selected if r['variant'] == variant and
                             r['process'] == process and
                             ((int(r['repeat']) == 0) == (state == 'first'))]
                    medians[(variant, process)] = {
                        field: statistics.median(float(r[field]) for r in group) for field in FIELDS}
            for variant in variants:
                row = dict(mode=mode, method=method, x=x, threads=threads,
                           state=state, variant=variant, processes=len(processes))
                for field in FIELDS:
                    values = [medians[(variant, p)][field] for p in processes]
                    row[field] = statistics.median(values)
                    row[field.replace('_ms', '_min_ms')] = min(values)
                    row[field.replace('_ms', '_max_ms')] = max(values)
                for field in ('p2_ms', 'total_ms'):
                    ratios = [medians[(variant, p)][field] / medians[('baseline', p)][field]
                              for p in processes]
                    row[field.replace('_ms', '_paired_change_percent')] = 100 * (statistics.median(ratios) - 1)
                    row[field.replace('_ms', '_paired_min_change_percent')] = 100 * (min(ratios) - 1)
                    row[field.replace('_ms', '_paired_max_change_percent')] = 100 * (max(ratios) - 1)
                    for reference in ('hybrid8', 'hybrid32', 'hybrid64'):
                        if reference not in variants:
                            continue
                        ratios_ref = [medians[(variant, p)][field] / medians[(reference, p)][field]
                                      for p in processes]
                        prefix = '_vs' + reference[len('hybrid'):] + '_paired'
                        row[field.replace('_ms', prefix + '_change_percent')] = 100 * (statistics.median(ratios_ref) - 1)
                        row[field.replace('_ms', prefix + '_min_change_percent')] = 100 * (min(ratios_ref) - 1)
                        row[field.replace('_ms', prefix + '_max_change_percent')] = 100 * (max(ratios_ref) - 1)
                row['p2_share_percent'] = 100 * row['p2_ms'] / row['total_ms']
                summaries.append(row)
    return summaries


def benchmark(root, mode, processes, repeats, variants=VARIANTS, seed=None):
    metadata = json.loads((root / 'timing-build-metadata.json').read_text())
    if list(variants) != metadata['variants']:
        raise RuntimeError('selected variants do not match timing build metadata')
    for name, digest in metadata['input_sha256'].items():
        if hashlib.sha256((root / name).read_bytes()).hexdigest() != digest:
            raise RuntimeError('timing source changed after build: ' + name)
    for variant, digest in metadata['binary_sha256'].items():
        if hashlib.sha256((root / ('profile-' + variant)).read_bytes()).hexdigest() != digest:
            raise RuntimeError('timing binary changed after build: ' + variant)
    if seed is None:
        seed = 202611091 if mode == 'phases' else 202611092
    rng = random.Random(seed)
    cases = [(method, x, threads, process) for method in ('meissel', 'lehmer')
             for x in COUNTS if mode == 'phases' or x >= 10**12
             for threads in (1, 8) for process in range(processes)]
    rng.shuffle(cases)
    partials = {(r['method'], r['x']): r['count'] for r in json.loads(
        (root / 'expected-p2.json').read_text())} if mode == 'p2' else {}
    rows = []
    raw = root / 'timing-raw'
    raw.mkdir(exist_ok=True)
    for block, (method, x, threads, process) in enumerate(cases):
        order_variants = list(variants)
        rng.shuffle(order_variants)
        for order, variant in enumerate(order_variants):
            expected = COUNTS[x] if mode == 'phases' else partials[(method, x)]
            command = [str(root / ('profile-' + variant)), mode, method, str(x),
                       str(threads), str(repeats), str(expected)]
            result = execute(command)
            label = f'{mode}-{method}-{x}-t{threads}-p{process}-{variant}'
            (raw / (label + '.csv')).write_text(result.stdout)
            if result.stderr:
                raise RuntimeError(result.stderr)
            parsed = list(csv.DictReader(io.StringIO(result.stdout)))
            if len(parsed) != repeats:
                raise RuntimeError('unexpected profiler row count')
            for repeat, row in enumerate(parsed):
                if (row['mode'], row['method'], int(row['x']), int(row['threads']), int(row['repeat'])) != (
                        mode, method, x, threads, repeat):
                    raise RuntimeError('unexpected profiler configuration')
                if mode == 'phases' and int(row['prime_count']) != COUNTS[x]:
                    raise RuntimeError('incorrect complete count')
                count = int(row['p2'])
                if partials.setdefault((method, x), count) != count:
                    raise RuntimeError('P2 differs between variants/repetitions')
                row.update(variant=variant, process=process, block=block, order=order)
                rows.append(row)
        save_csv(root / (mode + '.csv'), rows)
        print(f'{mode}: block {block + 1}/{len(cases)}; {method}, {x}, threads={threads}', flush=True)
    save_csv(root / (mode + '-summary.csv'), summarize(rows, variants))
    if mode == 'phases':
        save_json(root / 'expected-p2.json', [dict(method=m, x=x, count=count)
                  for (m, x), count in sorted(partials.items())])
    save_json(root / (mode + '-metadata.json'), dict(
        seed=seed, processes=processes, repeats=repeats, verified_rows=len(rows),
        blocks=len(cases), variants=list(variants), mode=mode,
        first='fresh process per configuration; Phi precedes P2 in phases mode',
        warm='median of subsequent calls within each process; tables rebuilt in phases mode',
        summary='median of process medians; ratios pair each candidate with its block baseline and, when present, 8L/32L/64L'))
    print(f'{mode}: all {len(rows)} counts matched', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('build', 'phases', 'p2'))
    parser.add_argument('artifact_dir', type=Path)
    parser.add_argument('--processes', type=int, default=3)
    parser.add_argument('--repeats', type=int, default=4)
    parser.add_argument('--variants', nargs='+', choices=AVAILABLE_VARIANTS,
                        default=list(VARIANTS))
    parser.add_argument('--seed', type=int, help='randomized block/order seed; defaults preserve original study')
    args = parser.parse_args()
    if args.processes < 1 or args.repeats < 2:
        parser.error('positive process count and at least two repeats required')
    if len(set(args.variants)) != len(args.variants) or 'baseline' not in args.variants:
        parser.error('select unique variants including baseline')
    root = args.artifact_dir.resolve()
    if root == REPO or REPO in root.parents:
        parser.error('artifacts must be outside the repository')
    if args.action == 'build':
        build(root, args.variants)
    else:
        benchmark(root, args.action, args.processes, args.repeats, args.variants, args.seed)


if __name__ == '__main__':
    main()
