#!/usr/bin/env python3
"""Check retained final P2 data and print first/warmed candidate comparisons.

This does not rerun benchmarks or regenerate the old exploratory report trees.
Each warmed process contributes its median; changes pair variants by block.
"""
import collections
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys

HERE = Path(__file__).resolve().parent
DATA = HERE / 'results-consolidated'


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    provenance = json.loads((DATA / 'provenance.json').read_text())
    for name, expected in provenance['raw_files'].items():
        require(hashlib.sha256((DATA / name).read_bytes()).hexdigest() == expected,
                'changed retained data: ' + name)
    suite = json.loads((DATA / 'suite.json').read_text())
    expectations = {row['case_id']: row for row in json.loads((DATA / 'expectations.json').read_text())}
    output = csv.writer(sys.stdout)
    output.writerow(('cohort', 'case_id', 'state', 'field', 'baseline_ms', 'candidate_ms', 'paired_change_percent'))
    for cohort in ('focus', 'holdout'):
        cases = {row['case_id']: row for row in suite['cases'] if row['cohort'] == cohort}
        expected_keys = {(cid, variant, block, repeat)
                         for cid, case in cases.items() for variant in case['variants']
                         for block in range(case['blocks']) for repeat in range(case['repeats'])}
        seen, groups, processes = set(), collections.defaultdict(list), set()
        with gzip.open(DATA / (cohort + '.csv.gz'), 'rt', newline='') as stream:
            for row in csv.DictReader(stream):
                cid, variant = row['case_id'], row['variant']
                block, repeat = int(row['block']), int(row['repeat'])
                key = cid, variant, block, repeat
                require(key in expected_keys and key not in seen, 'unexpected or repeated observation')
                seen.add(key)
                case = cases[cid]
                require(row['cohort'] == cohort and row['mode'] == case['mode'], 'case mismatch')
                require(int(row['result']) == int(row['expected']) == expectations[cid]['expected'], 'count mismatch')
                require(row['reference_provenance'] == expectations[cid]['provenance'], 'reference mismatch')
                require(row['measurement'] == 'timing' and row['verification'] == 'expected-match', 'non-timing row')
                state = 'first' if repeat == 0 else 'warm'
                require(row['cache_state'] == state, 'cache-state mismatch')
                field = 'p2_ms' if case['mode'] == 'p2' else 'total_ms'
                elapsed = float(row[field])
                require(math.isfinite(elapsed) and elapsed > 0, 'invalid elapsed time')
                groups[cid, variant, block, state].append(elapsed)
                processes.add((cid, variant, block, row['process']))
        require(seen == expected_keys, 'incomplete timing matrix')
        metadata = provenance['cohorts'][cohort]
        require(len(seen) == metadata['rows'] and len(processes) == metadata['processes'], 'cohort size mismatch')
        for cid, case in cases.items():
            field = 'p2_ms' if case['mode'] == 'p2' else 'total_ms'
            for state in ('first', 'warm'):
                values = {}
                for variant in ('baseline', 'candidate'):
                    values[variant] = []
                    for block in range(case['blocks']):
                        samples = groups[cid, variant, block, state]
                        require(len(samples) == (1 if state == 'first' else case['repeats'] - 1), 'repeat count mismatch')
                        values[variant].append(statistics.median(samples))
                changes = [100 * (new / old - 1) for old, new in zip(values['baseline'], values['candidate'])]
                output.writerow((cohort, cid, state, field, statistics.median(values['baseline']),
                                 statistics.median(values['candidate']), statistics.median(changes)))
        print('PASS ' + cohort + ': ' + str(len(seen)) + ' checked timed results', file=sys.stderr)


if __name__ == '__main__':
    main()
