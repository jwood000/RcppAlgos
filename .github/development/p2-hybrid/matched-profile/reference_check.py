#!/usr/bin/env python3
"""Independent references for the fixed 8e10, 6/5-table matched P2 study.

Uses a Python Eratosthenes sieve for every quotient after the Meissel cutoff,
and the already installed primecount C API for the 590 larger quotients.
No RcppAlgos kernels or prepared tables are used to calculate these answers.
This is a reference generator, never a timed benchmark dependency.
"""
import argparse
import csv
import ctypes
import hashlib
import json
import math
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sieve(limit):
    flags = bytearray(b'\x01') * (limit + 1)
    flags[0:2] = b'\x00\x00'
    for p in range(2, math.isqrt(limit) + 1):
        if flags[p]:
            start = p * p
            flags[start::p] = b'\x00' * ((limit - start) // p + 1)
    return flags


def cube_root(x):
    lo, hi = 0, x + 1
    while lo + 1 < hi:
        middle = (lo + hi) // 2
        if middle ** 3 <= x:
            lo = middle
        else:
            hi = middle
    return lo


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('artifact', type=Path)
    parser.add_argument('--library', type=Path,
                        default=Path('/opt/homebrew/lib/libprimecount.dylib'))
    args = parser.parse_args()
    args.artifact.mkdir(parents=True, exist_ok=True)
    output = args.artifact / 'references.json'
    terms_path = args.artifact / 'reference-terms.csv'
    if output.exists() or terms_path.exists():
        raise RuntimeError('refusing to overwrite existing reference data')

    x, threads, numerator, denominator = 80_000_000_000, 2, 6, 5
    root = math.isqrt(x)
    limit = root * numerator // denominator
    flags = sieve(root)
    primes = [0] + [p for p in range(2, root + 1) if flags[p]]
    b = len(primes) - 1
    cbrt = cube_root(x)
    a = sum(p <= cbrt for p in primes[1:])
    # Enumerate the mathematical membership predicates, independently of
    # the production binary-search implementation of either boundary.
    quotients = [0] + [x // p for p in primes[1:]]
    stream_first = next(i for i in range(1, b + 1)
                        if quotients[i] <= 16 * limit)
    last = next(i for i in range(1, b + 1) if quotients[i] <= limit)
    span = quotients[stream_first] - limit
    assert (root, limit, b, a, stream_first, last, span) == (
        282842, 339410, 24668, 590, 1725, 20886, 5089103)

    intervals = []
    lower = limit
    step, extra = divmod(span, threads)
    for slot in range(threads):
        upper = lower + step + (slot < extra)
        indices = [i for i in range(stream_first, last)
                   if lower < quotients[i] <= upper]
        assert indices and indices == list(range(indices[0], indices[-1] + 1))
        intervals.append(dict(slot=slot, lower=lower, upper=upper,
                              first=indices[0], last=indices[-1] + 1,
                              finish=quotients[indices[0]], queries=len(indices)))
        lower = upper
    assert sum(v['queries'] for v in intervals) == last - stream_first

    # All terms from the Meissel cutoff onward fit in this independent sieve.
    sieve_limit = quotients[a + 1]
    flags = sieve(sieve_limit)
    wanted = set(quotients[a + 1:])
    wanted.update(v['lower'] for v in intervals)
    wanted.update((limit, root, cbrt))
    pi = {}
    count = 0
    for n, is_prime in enumerate(flags):
        count += is_prime
        if n in wanted:
            pi[n] = count
    assert len(pi) == len(wanted)

    library = args.library.resolve(strict=True)
    api = ctypes.CDLL(str(library))
    api.primecount_pi.argtypes = [ctypes.c_int64]
    api.primecount_pi.restype = ctypes.c_int64
    api.primecount_set_num_threads.argtypes = [ctypes.c_int]
    api.primecount_set_num_threads.restype = None
    api.primecount_version.argtypes = []
    api.primecount_version.restype = ctypes.c_char_p
    api.primecount_set_num_threads(1)
    for q in quotients[1:a + 1]:
        value = api.primecount_pi(q)
        if value < 0:
            raise RuntimeError('primecount C API failed at ' + str(q))
        pi[q] = value

    cross_targets = sorted(set(
        [limit, root, cbrt, quotients[a + 1], quotients[stream_first],
         quotients[last - 1], quotients[last], quotients[b]] +
        [v['lower'] for v in intervals]))
    for q in cross_targets:
        if api.primecount_pi(q) != pi[q]:
            raise RuntimeError('independent sieve/API disagreement at ' + str(q))

    records, sums = [], dict(extra_zero_prefix=0, common_prefix=0,
                            common_tail=0, direct=0)
    for i in range(1, b + 1):
        part = ('extra_zero_prefix' if i <= a else
                'common_prefix' if i < stream_first else
                'common_tail' if i < last else 'direct')
        q = quotients[i]
        term = pi[q] - (i - 1)
        assert term > 0
        sums[part] += term
        records.append(dict(i=i, prime=primes[i], quotient=q, pi=pi[q],
                            term=term, part=part,
                            source='primecount-C-API' if i <= a else 'Python-sieve'))
    for interval in intervals:
        interval['pi_lower'] = pi[interval['lower']]
        interval['sum'] = sum(records[i - 1]['term']
                              for i in range(interval['first'], interval['last']))
    assert sum(v['sum'] for v in intervals) == sums['common_tail']

    zero = sum(sums.values())
    meissel = zero - sums['extra_zero_prefix']
    old_path = Path(__file__).resolve().parents[1] / 'results-selector-v3/expectations.json'
    old_expected = {v['id']: v['expected'] for v in json.loads(old_path.read_text())}
    assert zero == old_expected['p2-zero-80000000000-t2-L6over5']
    assert meissel == old_expected['p2-meissel-80000000000-t2-L6over5']

    with terms_path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    result = dict(
        x=x, threads=threads, table_numerator=numerator,
        table_denominator=denominator, root=root, table_limit=limit,
        cube_root=cbrt, meissel_a=a, b=b, stream_first=stream_first,
        non_table_last=last, stream_span=span,
        zero_prefix_queries=stream_first - 1,
        meissel_prefix_queries=stream_first - a - 1,
        common_tail_queries=last - stream_first,
        direct_queries=b + 1 - last,
        sums=sums, zero_prefix=sums['extra_zero_prefix'] + sums['common_prefix'],
        meissel_prefix=sums['common_prefix'], zero_p2=zero, meissel_p2=meissel,
        stream_intervals=intervals,
        reference_provenance=dict(
            description='Independent Python Eratosthenes sieve for i>=591; '
                        'installed primecount Gourdon C API for i=1..590. '
                        'No RcppAlgos kernel participates in calculation.',
            sieve_limit=sieve_limit, sieve_terms=b - a, primecount_terms=a,
            primecount_library=str(library), primecount_sha256=digest(library),
            primecount_version=api.primecount_version().decode(),
            primecount_threads=1, sieve_api_cross_targets=cross_targets,
            prior_v3_expectations_match=True,
            prior_v3_expectations_sha256=digest(old_path)),
        reference_terms_sha256=digest(terms_path),
        script_sha256=digest(Path(__file__)),
    )
    output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
