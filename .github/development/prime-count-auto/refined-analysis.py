#!/usr/bin/env python3
"""Audit both frozen timing matrices and evaluate the final six-thread policy.

The initial candidate's suite/analysis is never rewritten. This independently
reconstructs process medians and paired differences from raw R timing rows.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import statistics

METHODS = ('legendre', 'meissel', 'lehmer')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def selection(n, threads):
    if n < 10000000000:
        return 'legendre'
    return 'lehmer' if n >= 1000000000000 and threads >= 6 else 'meissel'


def audit(path, label):
    suite = json.loads((path / 'suite.json').read_text())
    metadata = json.loads((path / 'timings-metadata.json').read_text())
    references = json.loads((path / 'references.json').read_text())['references']
    assert sha(path / 'timings.csv.gz') == metadata['sha256']
    assert sha(path / 'suite.json') == metadata['suite_sha256']
    expected = {}
    for reference in references:
        values = expected.setdefault(reference['n'], set())
        values.add(reference['result'])
    assert all(len(values) == 1 for values in expected.values())
    expected = {n: values.pop() for n, values in expected.items()}
    cases = {case['id']: case for case in suite['cases']}
    rows = list(csv.DictReader(gzip.open(path / 'timings.csv.gz', 'rt')))
    observed, processes, samples = set(), {}, {}
    for row in rows:
        case = cases[row['case_id']]
        method, block, repeat = row['method'], int(row['block']), int(row['repeat_index'])
        key = row['case_id'], method, block, repeat
        assert key not in observed
        observed.add(key)
        assert int(float(row['n'])) == case['n']
        assert int(row['threads']) == case['threads'] <= int(row['maxThreads'])
        assert int(float(row['result'])) == expected[case['n']]
        assert math.isfinite(float(row['total_ms'])) and float(row['total_ms']) > 0
        identity = row['case_id'], method, block
        assert processes.setdefault(int(row['process']), identity) == identity
        state = 'first' if repeat == 0 else 'warm'
        samples.setdefault(identity + (state,), []).append(float(row['total_ms']))
    wanted = {(cid, m, b, r) for cid in cases for m in METHODS
              for b in range(suite['blocks']) for r in range(suite['repeats'])}
    assert observed == wanted
    assert len(rows) == metadata['rows'] and len(processes) == metadata['processes']
    medians = {key: statistics.median(values) for key, values in samples.items()}
    summaries, flags = [], []
    for cid, case in cases.items():
        selected = selection(case['n'], case['threads'])
        for state in ('first', 'warm'):
            values = {m: [medians[cid, m, b, state] for b in range(suite['blocks'])] for m in METHODS}
            fastest = min(METHODS, key=lambda m: statistics.median(values[m]))
            for control in dict.fromkeys(('legendre', fastest)):
                deltas = [chosen - baseline for chosen, baseline in zip(values[selected], values[control])]
                ratios = [100 * (chosen / baseline - 1) for chosen, baseline in zip(values[selected], values[control])]
                record = dict(study=label, case_id=cid, n=case['n'], threads=case['threads'],
                    state=state, selected=selected, control=control, fastest_median_method=fastest,
                    selected_median_ms=statistics.median(values[selected]),
                    control_median_ms=statistics.median(values[control]),
                    paired_delta_ms=statistics.median(deltas), paired_percent=statistics.median(ratios),
                    all_slower=all(d > 0 for d in deltas), all_faster=all(d < 0 for d in deltas))
                summaries.append(record)
                if control == 'legendre' and record['all_slower'] and record['paired_delta_ms'] > .05 and record['paired_percent'] > 5:
                    flags.append(record)
    return dict(study=label, rows=len(rows), processes=len(processes), cases=len(cases),
                raw_sha256=sha(path / 'timings.csv.gz'), suite_sha256=sha(path / 'suite.json'),
                comparisons=summaries, candidate_vs_legendre_flags=flags)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    output = root / 'refined-analysis.json'
    if output.exists():
        raise RuntimeError('refusing existing refined analysis')
    records = [audit(root / 'method-benchmark', 'initial'),
               audit(root / 'confirmation/method-benchmark', 'confirmation')]
    result = dict(policy='Historical lookup/sieve below1e5; Legendre below1e10; Meissel thereafter, '
                  'except Lehmer from1e12 with at least6usable threads.',
                  total_rows=sum(r['rows'] for r in records), total_processes=sum(r['processes'] for r in records),
                  studies=records, script_sha256=sha(Path(__file__)),
                  note='Post-study final-policy evaluation; preserves initial candidate analysis. '
                       'The repeated3e10/four-thread case remains in both studies. '
                       'No claim that all input sizes/thread budgets or cache histories have been measured.')
    output.write_text(json.dumps(result, indent=2) + '\n')
    print('PASS independent raw-data audit; counts:', result['total_rows'],
          'processes:', result['total_processes'])
    for record in records:
        print(record['study'], 'persistent first/warm flags:', len(record['candidate_vs_legendre_flags']))
