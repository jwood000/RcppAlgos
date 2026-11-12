#!/usr/bin/env python3
"""Audit and summarize the frozen common-K1 P2 comparison, without running it.

First calls remain separate from warmed calls. Warm repetitions are reduced
within their process, then variants are paired by randomized block. No samples
are discarded, and the diagnostic flags do not tune the candidate.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import math
import statistics
from pathlib import Path

FIELDS = ('setup_ms', 'phi_ms', 'p2_ms', 'p3_ms', 'teardown_ms', 'total_ms')
GEOMETRY = ('a', 'table_limit', 'root', 'b', 'non_table_queries', 'phi_queries',
            'stream_queries', 'direct_queries', 'phi_worker_budget', 'phi_workers',
            'stream_worker_budget', 'stream_workers', 'selected_k', 'selected_q',
            'candidate_q', 'candidate_stream_first', 'candidate_stream_queries',
            'candidate_stream_span', 'candidate_stream_worker_budget')


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n')


def write_csv(path, records):
    require(bool(records), 'no records for ' + str(path))
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def kind(variant):
    if variant == 'baseline':
        return 'baseline'
    if variant == 'phi1':
        return 'phi-initial-1'
    if variant.startswith('fixed'):
        return 'fixed'
    return 'candidate' if variant == 'candidate' else 'selector-' + variant


def dist(values):
    require(bool(values), 'empty distribution')
    return dict(n=len(values), median=statistics.median(values),
                minimum=min(values), maximum=max(values))


def audit(root, cohort, suite, expectations):
    raw = root / (cohort + '.csv.gz')
    metadata = json.loads((root / (cohort + '-metadata.json')).read_text())
    require(sha(raw) == metadata['sha256'], 'raw digest mismatch')
    require(sha(root / 'suite.json') == metadata['suite_sha256'], 'suite digest mismatch')
    require(suite['candidate_sha256'] == metadata['candidate_sha256'], 'candidate not frozen')
    cases = {c['case_id']: c for c in suite['cases'] if c['cohort'] == cohort}
    expected_keys = {(cid, variant, block, repeat)
                     for cid, case in cases.items() for variant in case['variants']
                     for block in range(case['blocks']) for repeat in range(case['repeats'])}
    with gzip.open(raw, 'rt', newline='') as stream:
        rows = list(csv.DictReader(stream))
    seen, process_ids, process_orders, geometries = set(), {}, {}, {}
    groups = collections.defaultdict(list)
    numeric_counts = {}
    for row in rows:
        cid, variant = row['case_id'], row['variant']
        block, repeat, pid = (int(row[k]) for k in ('block', 'repeat', 'process'))
        key = cid, variant, block, repeat
        require(key in expected_keys, 'unexpected row ' + str(key))
        require(key not in seen, 'duplicate row ' + str(key))
        seen.add(key)
        case = cases[cid]
        require(row['cohort'] == cohort, 'cohort mismatch')
        require((row['mode'], row['method'], int(row['x']), int(row['threads']),
                 int(row['table_numerator']), int(row['table_denominator'])) == (
                    case['mode'], case['cutoff'], case['x'], case['threads'],
                    case['numerator'], case['denominator']), 'case mismatch ' + cid)
        require(row['build_kind'] == kind(variant), 'build kind mismatch')
        require(int(row['initial_k']) == (32 if variant == 'baseline' else 1), 'queue label mismatch')
        require(row['measurement'] == 'timing' and row['verification'] == 'expected-match', 'unchecked or diagnostic row')
        require(row['cache_state'] == ('first' if repeat == 0 else 'warm'), 'cache state mismatch')
        require(row['timing_definition'] == 'phase-sum-excludes-diagnostics', 'timing definition mismatch')
        expected = expectations[cid]
        require(int(row['result']) == int(row['expected']) == expected['expected'], 'count mismatch')
        require(row['reference_provenance'] == expected['provenance'], 'reference mismatch')
        identity = cid, variant, block
        require(process_ids.setdefault(pid, identity) == identity, 'process identifier reused')
        require(process_orders.setdefault(identity, (pid, int(row['order']))) == (pid, int(row['order'])), 'process/order changed across repeats')
        require(0 <= int(row['order']) < len(case['variants']), 'invalid variant order')
        geometry = {name: row[name] for name in GEOMETRY}
        geometry['selection_reason'] = row['selection_reason']
        require(geometries.setdefault((cid, variant), geometry) == geometry, 'selection changed between calls')
        # A cutoff has the same count for every table size and thread budget.
        numeric_key = row['mode'], row['x'], row['a']
        require(numeric_counts.setdefault(numeric_key, row['result']) == row['result'], 'same numeric count differs')
        values = {field: float(row[field]) for field in FIELDS}
        require(all(math.isfinite(v) and v >= 0 for v in values.values()), 'invalid duration')
        phase_sum = sum(values[f] for f in FIELDS if f != 'total_ms')
        require(abs(phase_sum - values['total_ms']) <= max(1e-8, values['total_ms'] * 1e-9), 'phase sum mismatch')
        groups[identity + (row['cache_state'],)].append(values)
    require(seen == expected_keys, 'incomplete row matrix')
    require(len(rows) == metadata['rows'], 'metadata row count mismatch')
    require(len(process_ids) == metadata['processes'], 'metadata process count mismatch')
    require(len(metadata['runs']) == len(process_ids), 'run manifest size mismatch')
    for run in metadata['runs']:
        identity = run['case_id'], run['variant'], run['block']
        require(process_ids.get(run['process']) == identity, 'run manifest identity mismatch')
    for cid, case in cases.items():
        for block in range(case['blocks']):
            orders = [process_orders[cid, v, block][1] for v in case['variants']]
            require(sorted(orders) == list(range(len(case['variants']))), 'variant order not a permutation')
    medians = {}
    process_records = []
    for (cid, variant, block, state), values in sorted(groups.items()):
        expected_repeats = 1 if state == 'first' else cases[cid]['repeats'] - 1
        require(len(values) == expected_repeats, 'wrong first/warm count per process')
        record = dict(cohort=cohort, case_id=cid, variant=variant, block=block,
                      cache_state=state, samples=len(values))
        for field in FIELDS:
            median = statistics.median(v[field] for v in values)
            medians[cid, variant, block, state, field] = median
            record[field] = median
        process_records.append(record)
    return cases, medians, geometries, process_records, dict(
        cohort=cohort, rows=len(rows), processes=len(process_ids), cases=len(cases),
        first_warm_groups=len(groups), raw_sha256=sha(raw),
        checked='Complete matrix, identities, hashes, references, exact results, stable selection, phase sums, process pairing.')


def summarize(root, cohort, suite, expectations):
    cases, medians, geometries, process_records, audit_result = audit(root, cohort, suite, expectations)
    summaries, comparisons, flags, best_records = [], [], [], []
    for cid, case in cases.items():
        for variant in case['variants']:
            geometry = geometries[cid, variant]
            for state in ('first', 'warm'):
                for field in FIELDS:
                    values = [medians[cid, variant, block, state, field] for block in range(case['blocks'])]
                    summaries.append(dict(cohort=cohort, case_id=cid, mode=case['mode'],
                        cutoff=case['cutoff'], x=case['x'], threads=case['threads'],
                        numerator=case['numerator'], denominator=case['denominator'],
                        variant=variant, cache_state=state, field=field,
                        median_ms=statistics.median(values), minimum_ms=min(values), maximum_ms=max(values),
                        blocks=len(values), selection_reason=geometry['selection_reason'],
                        selected_k=geometry['selected_k'], phi_budget=geometry['phi_worker_budget'],
                        phi_queries=geometry['phi_queries'], stream_queries=geometry['stream_queries']))
        pairs = [('candidate', v) for v in case['variants'] if v != 'candidate']
        pairs += [(v, 'phi1') for v in ('v1', 'v2', 'v3') if v in case['variants']]
        pairs += [('phi1', 'baseline')]
        for left, right in pairs:
            for state in ('first', 'warm'):
                for field in FIELDS:
                    lvals = [medians[cid, left, b, state, field] for b in range(case['blocks'])]
                    rvals = [medians[cid, right, b, state, field] for b in range(case['blocks'])]
                    deltas = [l - r for l, r in zip(lvals, rvals)]
                    percentages = [100 * (l / r - 1) for l, r in zip(lvals, rvals) if r > 0]
                    item = dict(cohort=cohort, case_id=cid, mode=case['mode'], cutoff=case['cutoff'],
                        x=case['x'], threads=case['threads'], numerator=case['numerator'],
                        denominator=case['denominator'], left=left, right=right,
                        cache_state=state, field=field, blocks=len(deltas),
                        left_median_ms=statistics.median(lvals), right_median_ms=statistics.median(rvals),
                        delta_median_ms=statistics.median(deltas), delta_min_ms=min(deltas), delta_max_ms=max(deltas),
                        percent_median=statistics.median(percentages) if percentages else None,
                        percent_min=min(percentages) if percentages else None,
                        percent_max=max(percentages) if percentages else None,
                        percent_pairs=len(percentages), improved_blocks=sum(d < 0 for d in deltas),
                        regressed_blocks=sum(d > 0 for d in deltas),
                        all_improve=all(d < 0 for d in deltas), all_regress=all(d > 0 for d in deltas))
                    comparisons.append(item)
                    primary = 'p2_ms' if case['mode'] == 'p2' else 'total_ms'
                    if left == 'candidate' and right in ('phi1', 'baseline') and field == primary:
                        # Prospectively agreed screening threshold; it identifies
                        # a repeatable material tradeoff, not statistical proof.
                        short = state == 'warm' and geometries[cid, 'candidate']['selection_reason'] == 'short-full-range'
                        allowance = 0.02 if short else 0.05
                        relative = item['percent_median'] is not None and item['percent_median'] > 5
                        if item['all_regress'] and item['delta_median_ms'] > allowance and (short or relative):
                            flags.append(dict(item, absolute_allowance_ms=allowance,
                                              classification='short-full-absolute-regression' if short else 'material-regression'))
        controls = [v for v in case['variants'] if v in ('baseline', 'phi1') or v.startswith('fixed')]
        primary = 'p2_ms' if case['mode'] == 'p2' else 'total_ms'
        for state in ('first', 'warm'):
            control_medians = {v: statistics.median(medians[cid, v, b, state, primary]
                                                  for b in range(case['blocks'])) for v in controls}
            best = min(controls, key=control_medians.get)
            candidate = statistics.median(medians[cid, 'candidate', b, state, primary]
                                         for b in range(case['blocks']))
            best_records.append(dict(cohort=cohort, case_id=cid, mode=case['mode'], cache_state=state,
                field=primary, candidate_ms=candidate, best_control=best, best_control_ms=control_medians[best],
                candidate_delta_ms=candidate - control_medians[best],
                candidate_percent=100 * (candidate / control_medians[best] - 1) if control_medians[best] > 0 else None,
                interpretation='Optimistic hindsight control; not a prospective dispatch rule.'))
    write_csv(root / (cohort + '-process-medians.csv'), process_records)
    write_csv(root / (cohort + '-summary.csv'), summaries)
    write_csv(root / (cohort + '-comparisons.csv'), comparisons)
    write_csv(root / (cohort + '-best-controls.csv'), best_records)
    write_json(root / (cohort + '-flags.json'), flags)
    audit_result.update(summary_distributions=len(summaries), paired_distributions=len(comparisons),
                        flagged_comparisons=len(flags), analyzer_sha256=sha(Path(__file__)))
    write_json(root / (cohort + '-analysis-audit.json'), audit_result)
    print(json.dumps(audit_result, indent=2))
    if flags:
        print('Material regression flags:')
        for flag in flags:
            print(f'  {flag["case_id"]} {flag["cache_state"]} vs {flag["right"]}: '
                  f'{flag["percent_median"]:+.2f}% / {flag["delta_median_ms"]:+.6f} ms')
    return audit_result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--cohort', choices=('focus', 'holdout', 'both'), default='both')
    args = parser.parse_args()
    root = args.root.resolve()
    suite = json.loads((root / 'suite.json').read_text())
    expectations = {r['case_id']: r for r in json.loads((root / 'expectations.json').read_text())}
    require(set(expectations) == {c['case_id'] for c in suite['cases']}, 'reference coverage mismatch')
    cohorts = ('focus', 'holdout') if args.cohort == 'both' else (args.cohort,)
    audits = [summarize(root, c, suite, expectations) for c in cohorts]
    if args.cohort == 'both':
        write_json(root / 'analysis-audit.json', dict(cohorts=audits,
                   rows=sum(a['rows'] for a in audits), processes=sum(a['processes'] for a in audits),
                   policy_retuned=False, analyzer_sha256=sha(Path(__file__))))


if __name__ == '__main__':
    main()
