#!/usr/bin/env python3
"""Separate untimed candidate capacity/selection probes, never timed binaries.

The kernel instrumentation records vector capacities and task ranges. The
existing harness independently describes the selector's choice and checks it
against those actual ranges. This is allocation-capacity accounting, not peak
resident memory, and excludes allocator overhead and the process-wide Phi cache.
"""
import argparse
import csv
import io
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
import bench


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    bench.verify(root)
    suite = json.loads((root / 'suite.json').read_text())
    references = {r['case_id']: r for r in json.loads((root / 'expectations.json').read_text())}
    output = root / 'memory.csv'
    bench.require(not output.exists(), 'refusing memory probe overwrite')
    directory, source = root / 'build', root / 'source'
    flags = ['clang++', '-std=c++11', '-O3', '-g', '-fno-omit-frame-pointer',
             '-pthread', '-Wall', '-Wextra', '-Wpedantic',
             '-DP2_HYBRID_MEMORY', '-DP2_HYBRID_SELECTOR_CANDIDATE',
             '-I' + str(source / 'inst/include'), '-I' + str(directory),
             '-I' + str(bench.HERE.parent)]
    binary = root / 'memory-candidate'
    command = flags + [directory / 'selector_profile.cpp', directory / 'P2HybridInitial.cpp']
    command += [directory / (n + '.o') for n in ('Phi', 'PhiParallel', 'P3', 'Eratosthenes')]
    command += ['-o', binary]
    bench.invoke(command)
    records, runs = [], []
    for case in suite['cases']:
        expected = references[case['case_id']]['expected']
        command_run = [binary, case['mode'], case['cutoff'], case['x'], case['threads'],
                       1, expected, case['numerator'], case['denominator']]
        rows = list(csv.DictReader(io.StringIO(bench.invoke(command_run))))
        bench.require(len(rows) == 1, 'memory probe row count')
        row = rows[0]
        bench.check_row(row, case, expected, 0, 'candidate')
        bench.require(row['measurement'] == 'memory-probe', 'probe misclassified')
        bench.require(row['p2_ms'] == row['total_ms'] == '', 'probe contains timing')
        for field in ('sieve_capacity_bytes', 'next_capacity_bytes', 'phi_sum_bytes',
                      'stream_sum_bytes', 'conservative_scratch_bytes'):
            bench.require(row[field] != '' and int(row[field]) >= 0, 'missing capacity: ' + field)
        row.update(case_id=case['case_id'], cohort=case['cohort'], variant='candidate',
                   reference_provenance=references[case['case_id']]['provenance'])
        records.append(row)
        runs.append(dict(case_id=case['case_id'], command=list(map(str, command_run))))
    with output.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    bench.verify(root)
    metadata = dict(rows=len(records), instrumented_build=list(map(str, command)),
        binary_sha256=bench.sha(binary), raw_sha256=bench.sha(output), runs=runs,
        instrumentation_sha256=bench.sha(bench.HERE.parent / 'hybrid_metrics.h'),
        runner_sha256=bench.sha(Path(__file__)),
        scope='Vector capacities and actual selected ranges, checked against diagnostic selector description. '
              'Not RSS; excludes allocator bookkeeping, stack/thread storage and shared Phi cache. '
              'One first process per frozen case; no timings reported.')
    bench.save(root / 'memory-metadata.json', metadata)
    print(f'PASS {len(records)} untimed memory/selection probes', flush=True)


if __name__ == '__main__':
    main()
