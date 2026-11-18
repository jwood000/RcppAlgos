#!/usr/bin/env python3
"""Separate follow-up to the frozen initial matrix; no original output changes.

Bracket the usable-thread boundary and independently confirm the first-call
regression around 3e10/four threads. Uses the same frozen R wrapper/timer.
"""
import argparse
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('auto_benchmark', HERE / 'benchmark.py')
bench = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bench)
SEED = 2026111102


def prepare(root):
    original_out, original = bench.verify(root)
    followup = root / 'confirmation'
    out = followup / 'method-benchmark'
    if out.exists():
        raise RuntimeError('refusing existing confirmation matrix')
    out.mkdir(parents=True)
    cases = []
    for n in (2000000000000, 8000000000000, 80000000000000):
        for threads in (4, 6, 8):
            cases.append(dict(id=f'n{n}-t{threads}', n=n, threads=threads))
    for n in (20000000000, 30000000000, 40000000000):
        cases.append(dict(id=f'n{n}-t4', n=n, threads=4))
    suite = dict(original)
    suite.update(cases=cases, seed=SEED,
        purpose='Focused follow-up after reviewing the complete original43-case matrix: '
                'bracket thread admission at4/6/8 using new inputs; independently repeat3e10/four-thread first-call '
                'regression with adjacent2e10/4e10. Keep original observations and candidate analysis unchanged.',
        original_suite_sha256=bench.sha(original_out / 'suite.json'))
    suite['runner_sha256'] = dict(suite['runner_sha256'])
    suite['runner_sha256'][str(HERE / 'confirmation.py')] = bench.sha(HERE / 'confirmation.py')
    bench.save(out / 'suite.json', suite)
    references = []
    for n in sorted({c['n'] for c in cases}):
        values = []
        for algorithm in ('--gourdon', '--deleglise-rivat'):
            command = ['/opt/homebrew/bin/primecount', str(n), '--threads=1', algorithm]
            result = int(bench.invoke(command, out / f'reference-{n}-{algorithm[2:]}.log').strip())
            values.append(result)
            references.append(dict(n=n, algorithm=algorithm, result=result, command=command))
        if values[0] != values[1]:
            raise RuntimeError('independent reference algorithms disagree')
    bench.save(out / 'references.json', dict(references=references,
        version=bench.invoke(['/opt/homebrew/bin/primecount', '--version'])))
    print(f'Prepared{len(cases)} confirmation cases; independent counts agree.', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('prepare', 'timings', 'analyze'))
    parser.add_argument('root', type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    bench.SEED = SEED
    if args.mode == 'prepare':
        prepare(root)
    elif args.mode == 'timings':
        bench.timings(root / 'confirmation')
    else:
        bench.analyze(root / 'confirmation')
