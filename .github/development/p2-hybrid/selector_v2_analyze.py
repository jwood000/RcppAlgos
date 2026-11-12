#!/usr/bin/env python3
"""Summarize the frozen v2 experiment without discarding any raw samples.

Run against an artifact or retained results directory. Percent summaries use
the predeclared 0.01 ms fixed-control floor for isolated P2. All values come
from the runner's block-paired process-median comparisons.
"""
import argparse
import csv
import json
import statistics
from pathlib import Path


def distribution(rows, label):
    values = [float(r[label + '_paired_change_percent']) for r in rows]
    return dict(
        median_percent=statistics.median(values),
        min_percent=min(values), max_percent=max(values),
        improved=sum(v < 0 for v in values),
        slower_by_more_than_5_percent=sum(v > 5 for v in values),
        within_percent={str(b): sum(v <= b for v in values) for b in (5, 10, 20)})


def analyze(path):
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    lookup = {(r['case_id'], r['state'], r['variant']): r for r in rows}
    if len(lookup) != len(rows):
        raise RuntimeError('duplicate summary row')
    output = {}
    for state in ('first', 'warm'):
        all_v2 = [r for r in rows if r['variant'] == 'selector-v2' and r['state'] == state]

        def same_partition(row):
            other = lookup[row['case_id'], state, 'selector-v1']
            return (row['phi_queries'], row['stream_queries']) == (
                other['phi_queries'], other['stream_queries'])

        relative = [r for r in all_v2 if r['mode'] == 'phases' or
                    float(lookup[r['case_id'], state, r['best_fixed']]['p2_ms']) >= 0.01]
        tiny = [r for r in all_v2 if r not in relative]
        groups = {
            'all': relative,
            'regression': [r for r in relative if r['cohort'] == 'regression'],
            'new-evaluation': [r for r in relative if r['cohort'] == 'new-evaluation'],
            'same-partition': [r for r in relative if same_partition(r)],
            'changed-partition': [r for r in relative if not same_partition(r)],
        }
        record = dict(configurations=len(all_v2), relative_cases=len(relative),
                      tiny_cases=len(tiny), groups={})
        for name, group in groups.items():
            if not group:
                continue
            entry = dict(cases=len(group), selectors={})
            for variant in ('selector-v1', 'selector-v2'):
                selected = [lookup[r['case_id'], state, variant] for r in group]
                entry['selectors'][variant] = {
                    label: distribution(selected, label)
                    for label in ('baseline', 'best_fixed', 'chosen', 'v1')}
            record['groups'][name] = entry
        if tiny:
            values = [float(r['v1_paired_delta_ms']) for r in tiny]
            record['tiny_v2_vs_v1'] = dict(
                min_delta_ms=min(values), max_delta_ms=max(values),
                median_delta_ms=statistics.median(values),
                improved=sum(v < 0 for v in values), regressed=sum(v > 0 for v in values))
        fields = ('case_id', 'cohort', 'selection_reason', 'selected_k', 'best_fixed',
                  'v1_paired_change_percent', 'v1_paired_delta_ms',
                  'best_fixed_paired_change_percent', 'best_fixed_paired_delta_ms')
        for label in ('v1', 'best_fixed'):
            worst = sorted(relative, key=lambda r: float(r[label + '_paired_change_percent']),
                           reverse=True)[:10]
            record['worst_vs_' + label] = [{k: r[k] for k in fields} for r in worst]
        output[state] = record
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results_dir', type=Path)
    parser.add_argument('--mode', choices=('p2', 'phases', 'all'), default='all')
    args = parser.parse_args()
    modes = ('p2', 'phases') if args.mode == 'all' else (args.mode,)
    output = {'definitions': {
        'source': 'runner summary CSVs; each ratio pairs process medians within blocks',
        'case_aggregation': 'unweighted distribution, not an overall user-workload speedup',
        'relative_p2_floor_ms': 0.01,
        'negative_change': 'less elapsed time',
        'same_partition': 'equal Phi/stream query counts for the same configuration',
        'best_fixed': 'hindsight best of baseline and fixed8/16/32/64; excludes both selectors'}}
    for mode in modes:
        output[mode] = analyze(args.results_dir / (mode + '-summary.csv'))
    destination = args.results_dir / ('evaluation-summary.json' if args.mode == 'all'
                                    else args.mode + '-evaluation.json')
    destination.write_text(json.dumps(output, indent=2) + '\n')
    print(destination)


if __name__ == '__main__':
    main()
