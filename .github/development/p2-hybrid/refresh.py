#!/usr/bin/env python3
"""Refresh P2 work counters and diagnostic timing for the 131071 Phi cache.

Build first, then run diagnostics only while other builds/benchmarks are idle.
Production and historical profiling sources are read only. The extra quotient
bands measure Phi time (not total kernel time) and carry clock/probe overhead.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import random
import shutil
import statistics
import subprocess

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / 'p2-profile'
REPO = HERE.parents[2]
COUNTS = {10**12: 37607912018, 10**13: 346065536839}


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def command(args, log, timeout=300):
    result = subprocess.run(args, capture_output=True, text=True, timeout=timeout)
    log.write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f'command failed ({result.returncode}): {args}; see {log}')
    return result


def json_file(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def csv_file(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def build(root):
    source = root / 'source'
    destination = root / 'diagnostics'
    if destination.exists():
        raise RuntimeError('diagnostic directory already exists; use a fresh artifact root')
    if root == REPO or REPO in root.parents:
        raise RuntimeError('artifact root must be outside the repository')
    if 'constexpr std::uint64_t CACHE_MAX_X = 131071;' not in (source / 'src/Phi.cpp').read_text():
        raise RuntimeError('this profile requires the committed 131071 cache')
    instrument = load('p2_instrument', OLD / 'instrument.py')
    replace = instrument.replace_once
    p2 = instrument.instrument_p2((source / 'src/P2.cpp').read_text())
    phi = instrument.instrument_phi((source / 'src/Phi.cpp').read_text())
    # Only the band build starts the per-query clock after the quotient has
    # been formed. Its Phi duration remains measured after cutoff/root lookup,
    # as in the original timing probe. Do not compare its root-lookup clock
    # directly with the original probe, which also includes quotient division.
    p2 = replace(p2,
        '            P2ProfileMetrics::Query profileQuery(i);\n'
        '            const std::int64_t quotient = x / primes[i];\n',
        '#ifndef P2_PROFILE_QUOTIENT_BANDS\n'
        '            P2ProfileMetrics::Query profileQuery(i);\n'
        '#endif\n'
        '            const std::int64_t quotient = x / primes[i];\n'
        '#ifdef P2_PROFILE_QUOTIENT_BANDS\n'
        '            P2ProfileMetrics::Query profileQuery(quotient, tables.limit());\n'
        '#endif\n')
    header = (OLD / 'profile_metrics.h').read_text()
    header = replace(header, '    constexpr std::size_t BINS = 8;\n',
        '#ifdef P2_PROFILE_QUOTIENT_BANDS\n'
        '    constexpr std::size_t BINS = 4;\n'
        '#else\n'
        '    constexpr std::size_t BINS = 8;\n'
        '#endif\n')
    header = replace(header,
        '        explicit Query(std::int64_t index) : active(local().active) {\n'
        '            if (!active) return;\n'
        '            local().bin = static_cast<std::size_t>(\n'
        '                (index - local().first) * BINS / (local().last - local().first));\n',
        '#ifdef P2_PROFILE_QUOTIENT_BANDS\n'
        '        Query(std::int64_t quotient, std::int64_t limit) : active(local().active) {\n'
        '            if (!active) return;\n'
        '            local().bin = quotient <= 2 * limit ? 0 :\n'
        '                (quotient <= 4 * limit ? 1 : (quotient <= 8 * limit ? 2 : 3));\n'
        '#else\n'
        '        explicit Query(std::int64_t index) : active(local().active) {\n'
        '            if (!active) return;\n'
        '            local().bin = static_cast<std::size_t>(\n'
        '                (index - local().first) * BINS / (local().last - local().first));\n'
        '#endif\n')
    destination.mkdir()
    (destination / 'P2.cpp').write_text(p2)
    (destination / 'Phi.cpp').write_text(phi)
    (destination / 'profile_metrics.h').write_text(header)
    shutil.copy2(OLD / 'profile.cpp', destination / 'profile.cpp')
    commands = []
    for kind in ('counts', 'timing', 'bands'):
        args = ['clang++', '-std=c++11', '-O3', '-g', '-fno-omit-frame-pointer',
                '-pthread', '-Wall', '-Wextra', '-Wpedantic', '-DP2_PROFILE_METRICS']
        if kind != 'counts':
            args.append('-DP2_PROFILE_TIMING')
        if kind == 'bands':
            args.append('-DP2_PROFILE_QUOTIENT_BANDS')
        args += ['-I' + str(source / 'inst/include'), '-I' + str(destination),
                 str(destination / 'profile.cpp'), str(destination / 'P2.cpp'),
                 str(destination / 'Phi.cpp')]
        args += [str(source / 'src' / name) for name in ('P3.cpp', 'PhiParallel.cpp', 'Eratosthenes.cpp')]
        args += ['-o', str(destination / ('profile-' + kind))]
        command(args, destination / ('build-' + kind + '.log'))
        commands.append(args)
    files = [destination / name for name in ('profile.cpp', 'P2.cpp', 'Phi.cpp', 'profile_metrics.h')]
    files += [source / 'src' / name for name in ('P2.cpp', 'Phi.cpp', 'P3.cpp', 'PhiParallel.cpp', 'Eratosthenes.cpp')]
    json_file(destination / 'build-metadata.json', dict(
        source_commit=(root / 'commit.txt').read_text().strip(), cache_max_x=131071,
        platform=platform.platform(), compiler=subprocess.check_output(['clang++', '--version'], text=True),
        commands=commands,
        hashes={str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        band_mapping={'0': '(L,2L]', '1': '(2L,4L]', '2': '(4L,8L]', '3': '>8L'},
        band_clock='Phi clock excludes outer quotient formation; per-query clock overhead remains'))
    print(f'Built three diagnostic binaries in {destination}', flush=True)


def collect(root):
    destination = root / 'diagnostics'
    raw = destination / 'raw'
    raw.mkdir(exist_ok=True)
    parser = load('p2_diagnostic_parser', OLD / 'diagnostics.py')
    runs, partials = [], {}
    seed = 202611092
    rng = random.Random(seed)
    for mode in ('phases', 'p2'):
        cases = [(kind, method, x, threads) for kind in ('counts', 'timing')
                 for method in ('meissel', 'lehmer') for x in COUNTS for threads in (1, 8)]
        if mode == 'phases':
            cases += [('bands', method, x, threads)
                      for method in ('meissel', 'lehmer') for x in COUNTS for threads in (1, 8)]
        rng.shuffle(cases)
        for index, (kind, method, x, threads) in enumerate(cases):
            expected = COUNTS[x] if mode == 'phases' else partials[(method, x)]
            args = [str(destination / ('profile-' + kind)), mode, method, str(x),
                    str(threads), '3', str(expected)]
            label = f'{kind}-{mode}-{method}-{x}-t{threads}'
            result = command(args, raw / (label + '.txt'))
            if result.stderr:
                raise RuntimeError('unexpected diagnostic stderr: ' + result.stderr)
            groups = parser.parse(result.stdout, kind == 'counts')
            if len(groups) != 3:
                raise RuntimeError('expected one first and two warm calls')
            for repeat, group in enumerate(groups):
                row = group['row']
                if (row['mode'], row['method'], int(row['x']), int(row['threads']), int(row['repeat'])) != (
                        mode, method, x, threads, repeat):
                    raise RuntimeError('diagnostic configuration mismatch')
                if mode == 'phases' and int(row['prime_count']) != COUNTS[x]:
                    raise RuntimeError('known pi(x) mismatch')
                value = int(row['p2'])
                previous = partials.setdefault((method, x), value)
                if value != previous:
                    raise RuntimeError('P2 mismatch across contexts/instrumentation')
            runs.append(dict(kind=kind, mode=mode, method=method, x=x, threads=threads, groups=groups))
            json_file(destination / 'diagnostics.json', runs)
            print(f'{mode} diagnostic {index + 1}/{len(cases)}: {label}', flush=True)
    # Warm counts must agree across isolated/driver contexts and thread budgets.
    for method in ('meissel', 'lehmer'):
        for x in COUNTS:
            counts = [next(m for m in g['metrics'] if m['kind'] == 'summary')['counters']
                      for run in runs if run['kind'] == 'counts' and run['method'] == method and run['x'] == x
                      for g in run['groups'][1:]]
            if any(c != counts[0] for c in counts):
                raise RuntimeError('warmed work counts differ across contexts/workers')
    json_file(destination / 'run-metadata.json', dict(seed=seed, verified_rows=sum(len(r['groups']) for r in runs),
        repeats=3, inputs=COUNTS, thread_budgets=[1, 8], kinds=['counts', 'timing', 'bands'],
        isolated_p2_reference='derived from known complete identities; small independent tests are separate',
        timings='diagnostic clocks are not uninstrumented benchmark results'))
    summarize(root)


def summarize(root):
    destination = root / 'diagnostics'
    runs = json.loads((destination / 'diagnostics.json').read_text())
    rows, bands = [], []
    for run in runs:
        for group in run['groups']:
            metrics = group['metrics']
            summary = next(m for m in metrics if m['kind'] == 'summary')
            workers = [m for m in metrics if m['kind'] == 'worker']
            bins = [m for m in metrics if m['kind'] == 'quantile']
            row = {k: run[k] for k in ('kind', 'mode', 'method', 'x', 'threads')}
            row.update(repeat=int(group['row']['repeat']), p2_ms=float(group['row']['p2_ms']),
                       workers=summary['worker_slots'], direct_terms=summary['direct_terms'],
                       root_lookup_ns=summary['root_lookup_ns'], phi_ns=summary['phi_ns'])
            row.update(summary['counters'])
            row.update(worker_min_ns=min(w['wall_ns'] for w in workers),
                       worker_max_ns=max(w['wall_ns'] for w in workers),
                       max_chunk_ns=max(w['max_chunk_ns'] for w in workers))
            if run['kind'] != 'bands':
                row.update(first_eighth_phi_ns=sum(b['phi_ns'] for b in bins if b['bin'] == 0),
                           first_quarter_phi_ns=sum(b['phi_ns'] for b in bins if b['bin'] < 2))
                rows.append(row)
            else:
                for end, multiple in enumerate((2, 4, 8)):
                    band = dict(row)
                    band.update(multiple=multiple,
                                tail_queries=sum(b['counters']['queries'] for b in bins if b['bin'] <= end),
                                tail_phi_ns=sum(b['phi_ns'] for b in bins if b['bin'] <= end))
                    band['tail_query_percent'] = 100 * band['tail_queries'] / row['queries']
                    band['tail_phi_percent'] = 100 * band['tail_phi_ns'] / row['phi_ns']
                    bands.append(band)
    csv_file(destination / 'diagnostics-summary.csv', rows)
    csv_file(destination / 'quotient-bands.csv', bands)
    print(f'Exported {len(rows)} diagnostic rows and {len(bands)} cumulative quotient bands.', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=('build', 'run', 'summarize'))
    parser.add_argument('artifact_dir', type=Path)
    args = parser.parse_args()
    root = args.artifact_dir.resolve()
    {'build': build, 'run': collect, 'summarize': summarize}[args.phase](root)


if __name__ == '__main__':
    main()
