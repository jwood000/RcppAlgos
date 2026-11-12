#!/usr/bin/env python3
"""Analyze the frozen initial-chunk experiment; never run counting kernels.

One fresh process supplies each variant/case/block. Repeat zero is first;
later calls are reduced to a median within that process. Ratios pair process
medians in the same block. Across-case summaries are unweighted descriptions
of this targeted matrix, not predictions for a user's workload.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import statistics

FIELDS = ('setup_ms', 'phi_ms', 'p2_ms', 'p3_ms', 'teardown_ms', 'total_ms')
HEADS = ('head1', 'head2', 'head4', 'head8')
RELATIVE_FLOOR_MS = 0.01


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def case_id(case):
    value = case.get('case_id', case.get('id'))
    require(value is not None, 'suite case lacks id')
    return value


def distribution(values):
    require(bool(values), 'empty distribution')
    return dict(n=len(values), median=statistics.median(values),
                min=min(values), max=max(values), values=values)


def proportional_changes(values):
    result = distribution(values)
    result.update(improved=sum(v < 0 for v in values), regressed=sum(v > 0 for v in values),
                  all_improve=all(v < 0 for v in values), all_regress=all(v > 0 for v in values),
                  slower_by_more_than_5_percent=sum(v > 5 for v in values))
    return result


def state(row):
    return 'first' if int(row['repeat']) == 0 else 'warm'


def geometry(row):
    return {name: int(row[name]) for name in (
        'a', 'table_limit', 'root', 'b', 'non_table_queries', 'phi_queries',
        'stream_queries', 'direct_queries', 'phi_worker_budget', 'phi_workers',
        'stream_worker_budget', 'stream_workers', 'selected_k', 'selected_q')}


def save_csv(path, records):
    require(bool(records), 'empty CSV')
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def load(root):
    suite_path = root / 'suite.json'
    expected_path = root / 'expectations.json'
    raw_path = root / 'timings.csv.gz'
    suite = json.loads(suite_path.read_text())
    cases = {case_id(c): c for c in suite['cases']}
    require(len(cases) == len(suite['cases']), 'duplicate suite case')
    expectations = {case_id(c): c for c in json.loads(expected_path.read_text())}
    require(set(expectations) == set(cases), 'expectation coverage')
    with gzip.open(raw_path, 'rt', newline='') as stream:
        rows = list(csv.DictReader(stream))
    variants = suite['variants']
    seen, processes = set(), collections.defaultdict(list)
    process_identity, process_order = {}, {}
    counts, geometries = {}, {}
    numeric_p2 = {}
    for row in rows:
        identity = (row['case_id'], row['variant'], int(row['block']), int(row['repeat']))
        require(identity not in seen, 'duplicate timed row: ' + str(identity))
        seen.add(identity)
        cid, variant, block, repeat = identity
        require(cid in cases and variant in variants, 'unknown case or variant')
        case = cases[cid]
        process_id = int(row['process'])
        require(process_identity.setdefault(process_id, (cid, variant, block)) == (cid, variant, block),
                'process identifier reused for different configurations')
        require(process_order.setdefault((cid, variant, block), (process_id, int(row['order']))) ==
                (process_id, int(row['order'])), 'process/order changed within repeats')
        numerator = case.get('numerator', case.get('table_numerator', case.get('table_multiplier', 1)))
        denominator = case.get('denominator', case.get('table_denominator', 1))
        cutoff = case.get('cutoff', case.get('method'))
        require((row['mode'], row['method'], int(row['x']), int(row['threads']),
                 int(row['table_numerator']), int(row['table_denominator'])) ==
                (case['mode'], cutoff, case['x'], case['threads'], numerator, denominator),
                'configuration differs from frozen suite: ' + cid)
        require(0 <= block < case['blocks'] and 0 <= repeat < case['repeats'], 'sample bounds')
        require(row['measurement'] == 'timing' and row['verification'] == 'expected-match',
                'diagnostic or unchecked row in timing data')
        kind = 'initial-chunks' if variant.startswith('head') else 'fixed' if variant == 'hybrid16' else 'baseline'
        require(row['build_kind'] == kind, 'variant/build identity mismatch')
        if 'initial_k' in row:
            require(int(row['initial_k']) == (int(variant[4:]) if variant.startswith('head') else 0),
                    'initial chunk label mismatch')
        require(row['reference_provenance'] == expectations[cid]['provenance'], 'reference provenance mismatch')
        require(row['cache_state'] == state(row) and row['timing_definition'] == 'phase-sum-excludes-diagnostics',
                'timing state/definition changed')
        require(int(row['result']) == int(row['expected']) == expectations[cid]['expected'],
                'expected count mismatch: ' + cid)
        values = {f: float(row[f]) for f in FIELDS}
        require(all(math.isfinite(v) and v >= 0 for v in values.values()), 'invalid timing')
        require(abs(sum(values[f] for f in FIELDS if f != 'total_ms') - values['total_ms']) <= 0.00001,
                'phase sum mismatch')
        p2 = int(row['p2'])
        p2_key = (int(row['x']), int(row['a']))
        require(numeric_p2.setdefault(p2_key, p2) == p2,
                'P2 depends on variant, threads, or table coverage')
        if row['mode'] == 'phases':
            require(int(row['phi']) + int(row['a']) - 1 - p2 - int(row['p3']) == int(row['result']),
                    'complete formula result mismatch')
        g = geometry(row)
        require(g['table_limit'] == math.isqrt(int(row['x'])) * numerator // denominator,
                'incorrect rational table limit')
        require(g['phi_queries'] + g['stream_queries'] == g['non_table_queries'], 'query partition')
        require(g['phi_queries'] + g['stream_queries'] + g['direct_queries'] == g['b'] - g['a'],
                'full P2 range coverage')
        require(geometries.setdefault((cid, variant), g) == g, 'changing call geometry')
        require(counts.setdefault(cid, int(row['result'])) == int(row['result']), 'case result drift')
        processes[(cid, variant, block, state(row))].append(values)
    required = {(cid, variant, block, repeat)
                for cid, case in cases.items() for variant in variants
                for block in range(case['blocks']) for repeat in range(case['repeats'])}
    require(seen == required, 'missing or extra timing samples')
    for cid, case in cases.items():
        for block in range(case['blocks']):
            require({process_order[cid, v, block][1] for v in variants} == set(range(len(variants))),
                    'invalid randomized variant order')
    for cid in cases:
        baseline = geometries[(cid, 'baseline')]
        for variant in variants:
            if variant.startswith('head'):
                require(geometries[(cid, variant)] == baseline,
                        'initial-chunk variant changed partition/budgets')
    process_medians = {identity: {f: statistics.median(v[f] for v in values) for f in FIELDS}
                       for identity, values in processes.items()}
    audit = dict(raw_rows=len(rows), distinct_numeric_p2=len(numeric_p2),
                 cases=len(cases), variants=len(variants),
                 process_count=sum(c['blocks'] for c in cases.values()) * len(variants),
                 all_results_expected=True, all_numeric_p2_agree=True,
                 head_geometry_matches_baseline=True, all_rational_limits_exact=True,
                 all_phase_sums_match=True, samples_complete=True,
                 input_sha256={p.name: digest(p) for p in (suite_path, expected_path, raw_path)})
    metadata_path = root / 'timings-metadata.json'
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text())
        require(metadata['rows'] == len(rows) and metadata['processes'] == audit['process_count'],
                'runner metadata dimensions mismatch')
        require(metadata['sha256'] == digest(raw_path), 'runner raw timing hash mismatch')
        audit['input_sha256'][metadata_path.name] = digest(metadata_path)
        audit['runner_metadata_matches'] = True
    return suite, cases, expectations, geometries, process_medians, audit


def analyze(root):
    suite, cases, expectations, geometries, process, audit = load(root)
    variants = suite['variants']
    aggregates, comparisons, flat = [], [], []
    lookup = {}
    for cid, case in cases.items():
        primary = 'p2_ms' if case['mode'] == 'p2' else 'total_ms'
        for variant in variants:
            for s in ('first', 'warm'):
                medians = [process[cid, variant, block, s] for block in range(case['blocks'])]
                values = {f: distribution([m[f] for m in medians]) for f in FIELDS}
                record = dict(case_id=cid, variant=variant, mode=case['mode'],
                    cutoff=case.get('cutoff', case.get('method')), x=case['x'],
                    threads=case['threads'], state=s, primary_field=primary,
                    geometry=geometries[cid, variant], metrics=values)
                aggregates.append(record)
                lookup[cid, variant, s] = record
                flat.append(dict(case_id=cid, variant=variant, mode=record['mode'], cutoff=record['cutoff'],
                    x=record['x'], threads=record['threads'], actual_phi_workers=record['geometry']['phi_workers'],
                    table_limit=record['geometry']['table_limit'], state=s, blocks=case['blocks'],
                    primary_field=primary, primary_ms=values[primary]['median'],
                    primary_min_ms=values[primary]['min'], primary_max_ms=values[primary]['max'],
                    **{f: values[f]['median'] for f in FIELDS}))
        for s in ('first', 'warm'):
            for candidate in variants:
                references = ['baseline']
                if candidate in HEADS:
                    references += ['head32', 'hybrid16']
                for reference in references:
                    if candidate == reference:
                        continue
                    for metric in (('total_ms', 'p2_ms') if case['mode'] == 'phases' else ('p2_ms',)):
                        pairs = [(process[cid, candidate, block, s][metric],
                                  process[cid, reference, block, s][metric])
                                 for block in range(case['blocks'])]
                        require(all(base > 0 for _, base in pairs), 'zero comparison denominator')
                        changes = [100 * (value / base - 1) for value, base in pairs]
                        deltas = [value - base for value, base in pairs]
                        comparisons.append(dict(case_id=cid, mode=case['mode'],
                            cutoff=case.get('cutoff', case.get('method')), x=case['x'], threads=case['threads'],
                            actual_phi_workers=geometries[cid, 'baseline']['phi_workers'],
                            state=s, candidate=candidate, reference=reference, metric=metric,
                            primary_metric=metric == primary,
                            baseline_ms=lookup[cid, 'baseline', s]['metrics'][metric]['median'],
                            paired_percent=proportional_changes(changes), paired_delta_ms=distribution(deltas)))
    groups = []
    predicates = {
        'all-parallel': lambda r: r['actual_phi_workers'] > 1,
        'serial-controls': lambda r: r['actual_phi_workers'] == 1,
        'isolated-zero-parallel': lambda r: r['mode'] == 'p2' and r['cutoff'] == 'zero' and r['actual_phi_workers'] > 1,
        'isolated-meissel-parallel': lambda r: r['mode'] == 'p2' and r['cutoff'] == 'meissel' and r['actual_phi_workers'] > 1,
        'isolated-lehmer-parallel': lambda r: r['mode'] == 'p2' and r['cutoff'] == 'lehmer' and r['actual_phi_workers'] > 1,
        'complete-meissel-parallel': lambda r: r['mode'] == 'phases' and r['cutoff'] == 'meissel' and r['actual_phi_workers'] > 1,
        'complete-lehmer-parallel': lambda r: r['mode'] == 'phases' and r['cutoff'] == 'lehmer' and r['actual_phi_workers'] > 1,
    }
    for name, predicate in predicates.items():
        for s in ('first', 'warm'):
            for candidate in variants:
                if candidate == 'baseline':
                    continue
                records = [r for r in comparisons if r['reference'] == 'baseline' and r['candidate'] == candidate
                           and r['state'] == s and r['primary_metric'] and predicate(r)
                           and r['baseline_ms'] >= RELATIVE_FLOOR_MS]
                if not records:
                    continue
                groups.append(dict(group=name, state=s, variant=candidate, cases=len(records),
                    case_ids=[r['case_id'] for r in records],
                    case_median_changes=proportional_changes([r['paired_percent']['median'] for r in records]),
                    cases_all_blocks_improve=sum(r['paired_percent']['all_improve'] for r in records),
                    cases_all_blocks_regress=sum(r['paired_percent']['all_regress'] for r in records)))
    winners = []
    for cid, case in cases.items():
        primary = 'p2_ms' if case['mode'] == 'p2' else 'total_ms'
        for s in ('first', 'warm'):
            ranked = sorted(HEADS, key=lambda v: lookup[cid, v, s]['metrics'][primary]['median'])
            winners.append(dict(case_id=cid, state=s, metric=primary, best_initial=ranked[0],
                                ranked=ranked, note='Hindsight comparison, not a selected production policy.'))
    result = dict(definitions=dict(
        first='repeat 0 in a fresh process; no preceding Phi in isolated mode, main Phi precedes P2 in phases mode',
        warm='median of repeat>0 within each process; aggregate median/min/max of those process medians',
        primary_metric='isolated P2: p2_ms; complete native formula: total_ms; complete P2 phase also retained separately',
        percentage='100*(candidate/reference-1), paired by case/block/state within this cohort; negative means less time',
        variability='min/max are collected block ranges, not confidence intervals',
        geometry='requested threads and actual P2 Phi workers are separate; head variants retain baseline partition/budget',
        serial='worker==1 bypasses the queue; serial variant differences cannot establish a granularity benefit',
        relative_aggregate_floor_ms=RELATIVE_FLOOR_MS,
        aggregation='unweighted case summaries; not an overall workload speedup',
        hindsight_best='best median among head1/head2/head4/head8, exploratory only'),
        analyzer_sha256=digest(Path(__file__)), audit=audit, aggregates=aggregates,
        comparisons=comparisons, group_summaries=groups, hindsight_initial_winners=winners)
    save_csv(root / 'timings-summary.csv', flat)
    pair_rows = []
    for r in comparisons:
        pair_rows.append({k: v for k, v in r.items() if k not in ('paired_percent', 'paired_delta_ms')} |
                         {'paired_percent_median': r['paired_percent']['median'],
                          'paired_percent_min': r['paired_percent']['min'],
                          'paired_percent_max': r['paired_percent']['max'],
                          'paired_delta_ms_median': r['paired_delta_ms']['median'],
                          'all_blocks_improve': r['paired_percent']['all_improve'],
                          'all_blocks_regress': r['paired_percent']['all_regress']})
    save_csv(root / 'comparisons.csv', pair_rows)
    (root / 'analysis.json').write_text(json.dumps(result, indent=2) + '\n')
    lines = [f"Validated {audit['raw_rows']} timed results across {audit['cases']} configurations and {audit['variants']} variants.",
             'Percentages below are medians of matched process ratios; negative means less time.']
    for mode in ('p2', 'phases'):
        lines.append('\n' + mode + ' parallel summary versus original baseline:')
        for cutoff in ('zero', 'meissel', 'lehmer'):
            group = ('isolated-' if mode == 'p2' else 'complete-') + cutoff + '-parallel'
            for s in ('first', 'warm'):
                for variant in variants:
                    match = [r for r in groups if r['group'] == group and r['state'] == s and r['variant'] == variant]
                    if match:
                        r = match[0]; p = r['case_median_changes']
                        lines.append(f"{cutoff:7} {s:5} {variant:8}: median case {p['median']:+.2f}% "
                                     f"[{p['min']:+.2f},{p['max']:+.2f}], improved {p['improved']}/{p['n']}; "
                                     f"all-block improvement {r['cases_all_blocks_improve']}/{p['n']}")
    lines.append('\nMatched 8e10 cases, two requested threads, primary milliseconds:')
    for r in aggregates:
        if r['x'] == 80000000000 and r['threads'] == 2 and r['mode'] == 'p2':
            lines.append(f"{r['cutoff']:7} {r['state']:5} {r['variant']:8} {r['metrics']['p2_ms']['median']:.6f}")
    text = '\n'.join(lines) + '\n'
    (root / 'analysis.txt').write_text(text)
    print(text, end='')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    args = parser.parse_args()
    analyze(args.root)


if __name__ == '__main__':
    main()
