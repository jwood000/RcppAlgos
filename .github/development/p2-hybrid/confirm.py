#!/usr/bin/env python3
"""Separate confirmation of the isolated Lehmer 1e12 serial timing outlier.

Retains the primary matrix unchanged. Run only after all other CPU work ends.
"""
import argparse
import csv
import hashlib
import io
import json
import random
import statistics
from pathlib import Path
from bench import execute, save_csv, save_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('artifact_dir', type=Path)
    args = parser.parse_args()
    root = args.artifact_dir.resolve()
    metadata = json.loads((root / 'timing-build-metadata.json').read_text())
    for variant in ('baseline', 'hybrid8'):
        actual = hashlib.sha256((root / ('profile-' + variant)).read_bytes()).hexdigest()
        if actual != metadata['binary_sha256'][variant]:
            raise RuntimeError('timing executable changed')
    expected = next(r['count'] for r in json.loads((root / 'expected-p2.json').read_text())
                    if r['method'] == 'lehmer' and r['x'] == 10**12)
    seed = 202611093
    rng = random.Random(seed)
    rows = []
    for process in range(5):
        variants = ['baseline', 'hybrid8']
        rng.shuffle(variants)
        for order, variant in enumerate(variants):
            result = execute([str(root / ('profile-' + variant)), 'p2', 'lehmer',
                              str(10**12), '1', '4', str(expected)])
            if result.stderr:
                raise RuntimeError(result.stderr)
            parsed = list(csv.DictReader(io.StringIO(result.stdout)))
            if len(parsed) != 4:
                raise RuntimeError('unexpected row count')
            for repeat, row in enumerate(parsed):
                if (row['mode'], row['method'], int(row['x']), int(row['threads']),
                    int(row['repeat']), int(row['p2'])) != ('p2', 'lehmer', 10**12, 1, repeat, expected):
                    raise RuntimeError('incorrect configuration or count')
                row.update(variant=variant, process=process, order=order)
                rows.append(row)
    save_csv(root / 'confirmation.csv', rows)
    summary = []
    for state in ('first', 'warm'):
        medians = {(v, p): statistics.median(float(r['p2_ms']) for r in rows
                   if r['variant'] == v and r['process'] == p and
                   ((int(r['repeat']) == 0) == (state == 'first')))
                   for v in ('baseline', 'hybrid8') for p in range(5)}
        changes = [100 * (medians[('hybrid8', p)] / medians[('baseline', p)] - 1)
                   for p in range(5)]
        summary.append(dict(state=state,
            baseline_ms=statistics.median(medians[('baseline', p)] for p in range(5)),
            hybrid8_ms=statistics.median(medians[('hybrid8', p)] for p in range(5)),
            paired_change_percent=statistics.median(changes),
            paired_min_change_percent=min(changes), paired_max_change_percent=max(changes)))
    save_csv(root / 'confirmation-summary.csv', summary)
    save_json(root / 'confirmation-metadata.json', dict(seed=seed, processes=5, repeats=4,
        verified_rows=len(rows), mode='p2', method='lehmer', x=10**12, threads=1,
        reason='one slow hybrid8 process in the primary isolated matrix; original data retained'))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
