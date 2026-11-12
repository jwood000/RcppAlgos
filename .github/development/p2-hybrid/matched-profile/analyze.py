#!/usr/bin/env python3
"""Analyze retained matched-cutoff experiments; never run counting kernels.

First is repeat zero after any explicit history. Warmed means the median of
later repeats in the same process. Aggregate medians and paired ratios use
these process-level values, never pooled call samples from different processes.
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

KEYS = ('mode', 'cutoff', 'threads', 'history')
SPLIT, LAST = 1725, 20886
STAGES = ('prefix', 'tail', 'seed')
EPS = 0.00001


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def state(row):
    return 'first' if int(row['repeat']) == 0 else 'warm'


def key(row):
    return tuple(row[k] for k in KEYS)


def expected(row, ref):
    if row['mode'] == 'prefix':
        return ref[row['cutoff'] + '_prefix']
    if row['mode'].startswith('tail-'):
        return ref['sums']['common_tail']
    return ref[row['cutoff'] + '_p2']


def summary(values):
    require(bool(values), 'empty summary')
    return dict(n=len(values), median=statistics.median(values),
                min=min(values), max=max(values), values=values)


def normalized(row, controls):
    if not controls:
        return row
    return dict(mode=row['variant'], cutoff=row['method'],
                threads=int(row['threads']), history='none',
                block=int(row['block']), repeat=int(row['repeat']),
                elapsed_us=float(row['p2_ms']) * 1000,
                result=int(row['result']))


def load(root, name):
    path = root / (name + '.jsonl.gz')
    metadata_path = root / (name + '-metadata.json')
    metadata = json.loads(metadata_path.read_text())
    require(sha(path) == metadata['sha256'], 'raw hash mismatch: ' + name)
    with gzip.open(path, 'rt') as stream:
        rows = [normalized(json.loads(line), name == 'controls') for line in stream]
    require(len(rows) == metadata['rows'], 'raw row count mismatch: ' + name)
    seen = set()
    for row in rows:
        identity = key(row) + (row['block'], row['repeat'])
        require(identity not in seen, 'duplicate call: ' + str(identity))
        seen.add(identity)
        require(math.isfinite(row['elapsed_us']) and row['elapsed_us'] >= 0,
                'invalid elapsed time')
    return rows, dict(rows=len(rows), sha256=sha(path),
                      metadata_sha256=sha(metadata_path),
                      processes=len({key(r) + (r['block'],) for r in rows}))


def check_stats(c):
    require(all(isinstance(v, int) and v >= 0 for v in c.values()),
            'non-integer/negative counter')
    require(c['cache_gets'] == c['cache_eligible'] + c['cache_ineligible'],
            'cache eligibility accounting')
    require(c['cache_eligible'] == c['cache_hits'] + c['cache_misses'],
            'cache hit/miss accounting')
    require(c['cache_misses'] == c['cache_no_page'] + c['cache_empty_cell'],
            'cache miss accounting')
    require(c['expansions'] == c['cache_misses'] + c['cache_ineligible'],
            'expansion accounting')
    require(c['kernel_entries'] == c['kernel_tiny'] + c['cache_gets'],
            'kernel accounting')
    require(c['put_calls'] == c['expansions'], 'put accounting')
    require(c['tiny_calls'] == c['expansions'] + c['kernel_tiny'] + c['public_tiny'],
            'tiny accounting')
    exits = sum(c['public_' + v] for v in ('zero', 'a_zero', 'one', 'tiny', 'pi'))
    require(c['kernel_entries'] == c['recursive_entries'] + c['public_phi'] - exits,
            'recursive/public accounting')
    require(c['top_cache_eligible'] == c['top_cache_hits'] + c['top_cache_misses'],
            'top cache accounting')
    require(c['page_alloc_attempts'] == c['page_alloc_success'] + c['page_alloc_failed'],
            'page allocation accounting')
    require(c['cache_stores'] <= c['put_eligible'] <= c['put_calls'],
            'cache store accounting')


def expected_ranges(row):
    first = 1 if row['cutoff'] == 'zero' else 591
    mode = row['mode']
    if mode == 'baseline':
        return {0: (first, LAST)}
    if mode in ('hybrid', 'prefix'):
        return {0: (first, SPLIT)}
    if mode == 'split-phi':
        return {0: (first, SPLIT), 1: (SPLIT, LAST)}
    if mode == 'tail-phi':
        return {0: (SPLIT, LAST)}
    if mode == 'tail-stream':
        return {}
    raise RuntimeError('unknown probe mode: ' + mode)


def check_trace(row, ref, term_prefix):
    """Check temporal bounds and exact query coverage; return counter totals."""
    ranges = expected_ranges(row)
    groups = collections.defaultdict(list)
    counters = collections.defaultdict(collections.Counter)
    worker_ids = set()
    for worker in row['workers']:
        group, slot = worker['group'], worker['slot']
        require((group, slot) not in worker_ids, 'duplicate worker slot')
        worker_ids.add((group, slot))
        require(0 <= slot < row['threads'] and group in (0, 1), 'worker index')
        pb, pe = row['phases'][group]
        wb, we = worker['begin_us'], worker['end_us']
        require(-EPS <= pb <= wb + EPS <= we + 2 * EPS and
                we <= pe + EPS and pe <= row['elapsed_us'] + EPS, 'phase/worker times')
        previous_end = wb
        for first, last, begin, end in worker['chunks']:
            require(first < last and wb - EPS <= begin <= end + EPS and
                    end <= we + EPS and begin >= previous_end - EPS, 'chunk bounds')
            previous_end = end
            groups[group].append((first, last, begin, end, slot))
        require(all(worker[v] >= 0 for v in ('seed_us', 'setup_us', 'fill_us', 'fills')),
                'negative stream metric')
        require(worker['seed_us'] + worker['setup_us'] + worker['fill_us'] <= we - wb + EPS,
                'stream subtimers exceed worker duration')
        if 'stats' in worker:
            require(len(worker['stats']) == 3, 'stats stage count')
            for stage, stats in zip(STAGES, worker['stats']):
                check_stats(stats)
                counters[(group, stage)].update(stats)
    require(set(groups) == set(ranges), 'unexpected/missing chunk groups')
    for group, (first, last) in ranges.items():
        chunks = sorted(groups[group])
        cursor = first
        for begin, end, _, _, _ in chunks:
            require(begin == cursor, 'query overlap or gap')
            require(end <= last, 'query outside assigned range')
            if row['threads'] == 2:
                require(end == min(begin + 32, last), 'nonstandard dynamic chunk')
            cursor = end
        require(cursor == last, 'incomplete query range')
    # Independent reference terms cover exactly the mathematical pieces used.
    covered = sum(term_prefix[last - 1] - term_prefix[first - 1]
                  for first, last in ranges.values())
    if row['mode'] in ('hybrid', 'tail-stream'):
        covered += ref['sums']['common_tail']
    if row['mode'] in ('baseline', 'hybrid', 'split-phi'):
        covered += ref['sums']['direct']
    require(covered == row['result'], 'independent covered-term sum')
    if counters:
        expected_counts = collections.Counter()
        for group, (first, last) in ranges.items():
            expected_counts[(group, 'prefix')] = max(0, min(last, SPLIT) - first)
            expected_counts[(group, 'tail')] = max(0, last - max(first, SPLIT))
        if row['mode'] in ('hybrid', 'tail-stream'):
            expected_counts[(1, 'seed')] = int(row['threads'] == 2)
        for identity in set(counters) | set(expected_counts):
            require(counters[identity]['public_phi'] == expected_counts[identity],
                    'per-group public Phi coverage: ' + str(identity))
    return counters


def timeline_metrics(row):
    metrics = dict(elapsed_us=row['elapsed_us'])
    chunks = [c for w in row['workers'] for c in w['chunks']]
    prefix = [c for c in chunks if c[0] < SPLIT]
    tail = [c for c in chunks if c[0] >= SPLIT]
    boundary = [c for c in chunks if c[0] < SPLIT < c[1]]
    first = 1 if row['cutoff'] == 'zero' else 591
    initial = [c for c in chunks if c[0] == first]
    metrics['boundary_chunks'] = len(boundary)
    metrics['boundary_queries'] = sum(c[1] - c[0] for c in boundary)
    metrics['boundary_tail_queries'] = sum(c[1] - SPLIT for c in boundary)
    metrics['boundary_work_us'] = sum(c[3] - c[2] for c in boundary)
    metrics['tail_only_chunks'] = len(tail)
    metrics['tail_only_queries'] = sum(c[1] - c[0] for c in tail)
    metrics['tail_only_worker_us'] = sum(c[3] - c[2] for c in tail)
    if initial:
        require(len(initial) == 1, 'multiple initial chunks')
        metrics['initial_chunk_begin_us'] = initial[0][2]
        metrics['initial_chunk_end_us'] = initial[0][3]
        metrics['initial_chunk_work_us'] = initial[0][3] - initial[0][2]
    if prefix:
        finish = max(c[3] for c in prefix)
        metrics['prefix_associated_last_end_us'] = finish
        metrics['tail_only_chunks_completed_before_prefix_finish'] = sum(c[3] <= finish for c in tail)
        metrics['tail_only_queries_completed_before_prefix_finish'] = sum(c[1] - c[0] for c in tail if c[3] <= finish)
        metrics['tail_only_worker_us_before_prefix_finish'] = sum(max(0, min(c[3], finish) - c[2]) for c in tail)
    if tail:
        metrics['tail_only_first_begin_us'] = min(c[2] for c in tail)
        metrics['tail_only_last_end_us'] = max(c[3] for c in tail)
        if prefix:
            metrics['all_tail_only_completed_before_prefix_finish'] = int(max(c[3] for c in tail) <= finish)
            metrics['tail_completion_minus_prefix_finish_us'] = max(c[3] for c in tail) - finish
        if initial:
            metrics['all_tail_only_completed_before_initial_chunk'] = int(max(c[3] for c in tail) <= initial[0][3])
            metrics['tail_completion_minus_initial_end_us'] = max(c[3] for c in tail) - initial[0][3]
    for group, (begin, end) in enumerate(row['phases']):
        workers = [w for w in row['workers'] if w['group'] == group]
        if not workers:
            continue
        prefix = 'phase' + str(group) + '.'
        metrics[prefix + 'begin_us'] = begin
        metrics[prefix + 'end_us'] = end
        metrics[prefix + 'duration_us'] = end - begin
        metrics[prefix + 'first_worker_start_us'] = min(w['begin_us'] for w in workers)
        metrics[prefix + 'last_worker_end_us'] = max(w['end_us'] for w in workers)
        metrics[prefix + 'launch_gap_us'] = min(w['begin_us'] for w in workers) - begin
        metrics[prefix + 'join_gap_us'] = end - max(w['end_us'] for w in workers)
        metrics[prefix + 'finish_gap_us'] = max(w['end_us'] for w in workers) - min(w['end_us'] for w in workers)
        for slot, worker in enumerate(sorted(workers, key=lambda w: w['slot'])):
            wp = prefix + 'slot' + str(slot) + '.'
            for name in ('begin_us', 'end_us', 'seed_us', 'setup_us', 'fill_us', 'fills'):
                metrics[wp + name] = worker[name]
            metrics[wp + 'duration_us'] = worker['end_us'] - worker['begin_us']
            metrics[wp + 'other_worker_us'] = (worker['end_us'] - worker['begin_us'] -
                worker['seed_us'] - worker['setup_us'] - worker['fill_us'])
    if all(end > begin for begin, end in row['phases']):
        metrics['phase_handoff_gap_us'] = row['phases'][1][0] - row['phases'][0][1]
    return metrics


def counter_metrics(counters):
    metrics = {}
    totals = collections.defaultdict(collections.Counter)
    for (group, stage), counts in counters.items():
        totals[stage].update(counts)
        for name, value in counts.items():
            metrics['group' + str(group) + '.' + stage + '.' + name] = value
    for stage, counts in totals.items():
        for name, value in counts.items():
            metrics[stage + '.' + name] = value
    return metrics


def aggregate(rows, observations):
    per_process = collections.defaultdict(lambda: collections.defaultdict(list))
    for row, metrics in zip(rows, observations):
        identity = key(row) + (state(row), row['block'])
        for name, value in metrics.items():
            require(math.isfinite(value), 'nonfinite aggregate metric')
            per_process[identity][name].append(value)
    medians = {k: {name: statistics.median(values) for name, values in v.items()}
               for k, v in per_process.items()}
    groups = collections.defaultdict(lambda: collections.defaultdict(list))
    for identity, values in sorted(medians.items()):
        for name, value in values.items():
            groups[identity[:-1]][name].append(value)
    records = []
    for identity, metrics in sorted(groups.items()):
        record = dict(zip(KEYS + ('state',), identity))
        record['metrics'] = {name: summary(values) for name, values in metrics.items()}
        records.append(record)
    return records, medians


def comparisons(medians, controls=False):
    result = []
    if controls:
        pairs = [('hybrid16', 'baseline'), ('selector-v2', 'baseline'),
                 ('selector-v3', 'baseline'), ('selector-v3', 'selector-v2')]
    else:
        pairs = [('hybrid', 'baseline'), ('split-phi', 'baseline'),
                 ('hybrid', 'split-phi'), ('tail-stream', 'tail-phi')]
    for candidate, reference in pairs:
        ratios, differences = collections.defaultdict(list), collections.defaultdict(list)
        for identity, values in sorted(medians.items()):
            if identity[0] != candidate:
                continue
            other = (reference,) + identity[1:]
            if other not in medians:
                continue
            base = medians[other]['elapsed_us']
            value = values['elapsed_us']
            require(base > 0, 'zero paired reference time')
            config = identity[1:-1]
            ratios[config].append(100 * (value / base - 1))
            differences[config].append(value - base)
        for config in sorted(ratios):
            result.append(dict(candidate=candidate, reference=reference,
                **dict(zip(KEYS[1:] + ('state',), config)),
                paired_percent=summary(ratios[config]),
                paired_difference_us=summary(differences[config])))
    return result


def concise(result):
    lines = ['Matched-cutoff analysis (microseconds; process medians).',
             'Negative paired percentages mean less time for the candidate.']
    for name in ('controls', 'plain'):
        if name not in result['datasets']:
            continue
        lines.append('\n' + name + ' matched comparisons:')
        for r in result['datasets'][name]['comparisons']:
            if r['threads'] != 2 or (r['history'] != 'none' and r['state'] == 'warm'):
                continue
            p = r['paired_percent']
            lines.append(f"{r['cutoff']:7} t{r['threads']} {r['history']:14} {r['state']:5} "
                         f"{r['candidate']}/{r['reference']}: {p['median']:+.2f}% "
                         f"[{p['min']:+.2f}, {p['max']:+.2f}] n={p['n']}")
    if 'timeline' in result['datasets']:
        lines.append('\nTimeline baseline tail/prefix overlap (instrumented):')
        for r in result['datasets']['timeline']['aggregates']:
            if r['mode'] != 'baseline':
                continue
            m = r['metrics']
            lines.append(f"{r['cutoff']} {r['state']}: initial chunk "
                         f"{m['initial_chunk_work_us']['median']:.1f}us; pure tail ends "
                         f"{m['tail_completion_minus_initial_end_us']['median']:+.1f}us "
                         'relative to initial chunk end; fraction per-process median ' +
                         f"all tail completed first={m['all_tail_only_completed_before_initial_chunk']['median']:.2f}")
    if 'counts' in result['datasets']:
        lines.append('\nCounts only (not usable as timing measurements):')
        for r in result['datasets']['counts']['aggregates']:
            if r['threads'] != 2 or r['mode'] not in ('baseline', 'hybrid', 'tail-phi'):
                continue
            m = r['metrics']
            pieces = []
            for stage in STAGES:
                if stage + '.public_phi' not in m or m[stage + '.public_phi']['median'] == 0:
                    continue
                pieces.append(stage + ': ' + ', '.join(
                    n + '=' + str(m[stage + '.' + n]['median'])
                    for n in ('public_phi', 'cache_hits', 'cache_misses', 'expansions')))
            lines.append(f"{r['mode']} {r['cutoff']} {r['history']} {r['state']}: " + '; '.join(pieces))
    return '\n'.join(lines) + '\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--datasets', nargs='+', default=['plain', 'controls', 'timeline', 'counts'])
    parser.add_argument('--output', default='analysis')
    args = parser.parse_args()
    ref_path = args.root / 'references.json'
    ref = json.loads(ref_path.read_text())
    terms_path = args.root / 'reference-terms.csv'
    require(sha(terms_path) == ref['reference_terms_sha256'], 'reference terms hash')
    term_prefix = [0]
    with terms_path.open() as stream:
        for i, term in enumerate(csv.DictReader(stream), 1):
            require(int(term['i']) == i, 'reference term order')
            term_prefix.append(term_prefix[-1] + int(term['term']))
    require(len(term_prefix) == ref['b'] + 1, 'reference term count')
    result = dict(definitions=dict(
        time_unit='microseconds',
        first='repeat zero after the explicit history, if any',
        warm='median of repeat>0 within each process',
        aggregate='median/min/max of process-level medians; values in block order',
        comparison='100*(candidate/reference-1), matched by block within the same dataset',
        uncertainty='min/max span the collected blocks; not confidence intervals',
        timeline='Instrumented separate run, not added to uninstrumented timing results.',
        boundary='Chunks crossing i=1725 are prefix-associated. Pure-tail metrics exclude their tail queries; no boundary chunk time is split.',
        overlap='Worker time sums overlapping chunks; not elapsed time or claimed speedup.',
        counts='Counter-build times are excluded from performance analysis.'),
        references_sha256=sha(ref_path), analyzer_sha256=sha(Path(__file__)), datasets={})
    for name in args.datasets:
        rows, audit = load(args.root, name)
        observations = []
        traced, counted = 0, 0
        overlap_events = collections.defaultdict(collections.Counter)
        for row in rows:
            require(row['result'] == expected(row, ref), 'result differs from independent reference')
            if row.get('workers'):
                counters = check_trace(row, ref, term_prefix)
                traced += 1
                if counters:
                    observations.append(counter_metrics(counters))
                    counted += 1
                else:
                    metrics = timeline_metrics(row)
                    observations.append(metrics)
                    if (row['mode'] == 'baseline' and
                            'all_tail_only_completed_before_initial_chunk' in metrics):
                        event = metrics['all_tail_only_completed_before_initial_chunk']
                        overlap_events[(row['cutoff'], state(row))][bool(event)] += 1
            else:
                observations.append(dict(elapsed_us=row['elapsed_us']))
        audit.update(independent_result_matches=len(rows), trace_coverage_checks=traced,
                     counter_coverage_checks=counted)
        if overlap_events:
            audit['baseline_pure_tail_completed_before_initial_chunk'] = [
                dict(cutoff=cutoff, state=s, yes=events[True], no=events[False])
                for (cutoff, s), events in sorted(overlap_events.items())]
        records, medians = aggregate(rows, observations)
        timed = counted == 0 and name in ('plain', 'controls')
        result['datasets'][name] = dict(audit=audit, aggregates=records,
            comparisons=comparisons(medians, name == 'controls') if timed else [])
    (args.root / (args.output + '.json')).write_text(json.dumps(result, indent=2) + '\n')
    report = concise(result)
    (args.root / (args.output + '.txt')).write_text(report)
    print(report, end='')


if __name__ == '__main__':
    main()
