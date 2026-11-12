#!/usr/bin/env python3
"""Small synthetic checks for data reduction and frozen regression thresholds."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import sys
import tempfile

sys.dont_write_bytecode = True
import analyze


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path, help='artifact root for the check record')
    args = parser.parse_args()
    fixture = Path(tempfile.mkdtemp(prefix='p2-analysis-fixture-'))
    cases = [dict(case_id=label, cohort='focus', mode='p2', cutoff='zero',
                  x=10000 + index, threads=1, numerator=1, denominator=1,
                  blocks=3, repeats=4, variants=['baseline', 'phi1', 'candidate'])
             for index, label in enumerate(('normal', 'short'))]
    suite = dict(cases=cases, candidate_sha256='synthetic-not-a-native-build')
    (fixture / 'suite.json').write_text(json.dumps(suite))
    expected = {c['case_id']: dict(case_id=c['case_id'], expected=100,
                                 provenance='synthetic-audit-fixture') for c in cases}
    rows, runs, process = [], [], 0
    for case in cases:
        for block in range(case['blocks']):
            for order, variant in enumerate(case['variants']):
                runs.append(dict(case_id=case['case_id'], variant=variant,
                                 block=block, process=process))
                values = {'baseline': [10, 2, 4, 6], 'phi1': [5, 1, 2, 3],
                          'candidate': [6, 1.5, 2.5, 3.5]}[variant]
                if case['case_id'] == 'short':
                    values = [10.03, 1.03, 1.03, 1.03] if variant == 'candidate' else [10, 1, 1, 1]
                for repeat, duration in enumerate(values):
                    row = {g: '0' for g in analyze.GEOMETRY}
                    row.update(a='0', table_limit='100', root='100', b='25',
                        selection_reason='short-full-range' if case['case_id'] == 'short' and variant == 'candidate' else 'fixture',
                        case_id=case['case_id'], variant=variant, block=block,
                        repeat=repeat, process=process, order=order, cohort='focus',
                        mode='p2', method='zero', x=case['x'], threads=1,
                        table_numerator=1, table_denominator=1,
                        build_kind=analyze.kind(variant), initial_k=32 if variant == 'baseline' else 1,
                        measurement='timing', verification='expected-match',
                        cache_state='first' if repeat == 0 else 'warm',
                        timing_definition='phase-sum-excludes-diagnostics',
                        result=100, expected=100, reference_provenance='synthetic-audit-fixture')
                    row.update({f: 0 for f in analyze.FIELDS})
                    row.update(p2_ms=duration, total_ms=duration)
                    rows.append(row)
                process += 1

    def write(data):
        raw = fixture / 'focus.csv.gz'
        with gzip.open(raw, 'wt', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
        (fixture / 'focus-metadata.json').write_text(json.dumps(dict(
            rows=len(data), processes=process, runs=runs, sha256=analyze.sha(raw),
            suite_sha256=analyze.sha(fixture / 'suite.json'),
            candidate_sha256=suite['candidate_sha256'])))

    write(rows)
    analyze.summarize(fixture, 'focus', suite, expected)
    with (fixture / 'focus-comparisons.csv').open() as stream:
        comparisons = list(csv.DictReader(stream))
    for state, percentage, delta in [('first', 20, 1), ('warm', 25, 0.5)]:
        result = next(r for r in comparisons if r['case_id'] == 'normal' and
                      r['left'] == 'candidate' and r['right'] == 'phi1' and
                      r['cache_state'] == state and r['field'] == 'p2_ms')
        analyze.require(abs(float(result['percent_median']) - percentage) < 1e-10, 'incorrect paired percentage')
        analyze.require(abs(float(result['delta_median_ms']) - delta) < 1e-10, 'incorrect paired absolute delta')
    flags = json.loads((fixture / 'focus-flags.json').read_text())
    short_flags = [r for r in flags if r['case_id'] == 'short']
    analyze.require(len(short_flags) == 2 and all(r['cache_state'] == 'warm' for r in short_flags),
                    'short-range absolute allowance must apply only to warm calls')
    rejections = []
    for name in ('corrupt count', 'duplicate row'):
        broken = [dict(r) for r in rows]
        if name == 'corrupt count':
            broken[0]['result'] = 101
            message = 'count mismatch'
        else:
            broken[-1] = dict(broken[0])
            message = 'duplicate row'
        write(broken)
        try:
            analyze.audit(fixture, 'focus', suite, expected)
        except RuntimeError as error:
            analyze.require(message in str(error), 'unexpected corruption rejection')
            rejections.append(name)
        else:
            raise RuntimeError('accepted ' + name)
    metadata = dict(scope='Synthetic analysis fixture only; not a native counting test',
        rows=len(rows), processes=process, cases=len(cases),
        first_paired_percent=20, warm_paired_percent=25,
        short_full_absolute_allowance_warm_only=True, rejected_corruptions=rejections,
        analyzer_sha256=analyze.sha(Path(analyze.__file__)), runner_sha256=analyze.sha(Path(__file__)))
    (args.root / 'analysis-self-check.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print('PASS analysis self-check, including warm-only short-range allowance')


if __name__ == '__main__':
    main()
