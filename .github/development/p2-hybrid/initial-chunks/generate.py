#!/usr/bin/env python3
"""Generate a temporary P2 copy with only the first 32 queries rechunked."""
import argparse
import hashlib
import json
from pathlib import Path

EXPECTED_P2_SHA256 = '9234a232b99861dc4d447d62ee3c0482aa43174f083dd0f47cd69cf4cf711af7'


def once(source, old, new):
    if source.count(old) != 1:
        raise RuntimeError('expected one source anchor: ' + old)
    return source.replace(old, new)


def generate(source):
    if hashlib.sha256(source.encode()).hexdigest() != EXPECTED_P2_SHA256:
        raise RuntimeError('frozen production P2 source changed')
    result = '''// Generated temporary experiment; never compiled by the R package.
#include "InitialChunks.h"
#ifndef P2_INITIAL_CHUNK
#error "Define P2_INITIAL_CHUNK as 1, 2, 4, 8, or 32"
#endif
#ifndef INITIAL_CHUNK_WORKER
#define INITIAL_CHUNK_WORKER(slot)
#endif
#ifndef INITIAL_CHUNK_RANGE
#define INITIAL_CHUNK_RANGE(begin,end)
#endif
#ifndef INITIAL_CHUNK_QUERY
#define INITIAL_CHUNK_QUERY(i)
#endif
''' + source
    result = once(result,
        '            const std::int64_t quotient = x / primes[i];',
        '            INITIAL_CHUNK_QUERY(i);\n'
        '            const std::int64_t quotient = x / primes[i];')
    result = once(result,
        '''        // Early quotients can be much more expensive than later ones. Workers
        // take small contiguous chunks until the sum is exhausted. The atomic
        // counter assigns unique indices only; it publishes no other data, so
        // relaxed ordering suffices. Launch/join synchronize tables/results.
        constexpr std::int64_t CHUNK_SIZE = 32;
        std::atomic<std::int64_t> next(first);''',
        '''        // One dynamic ticket queue: smaller chunks for the first 32
        // queries, then the original 32-query chunks. No intermediate join.
        InitialChunks::Queue<P2_INITIAL_CHUNK> next(first, last);''')
    result = once(result, '            batch.launch([&, slot] {',
        '            batch.launch([&, slot] {\n                INITIAL_CHUNK_WORKER(slot);')
    result = once(result,
        '''                    const std::int64_t begin = next.fetch_add(CHUNK_SIZE,
                        std::memory_order_relaxed);
                    if (begin >= last) break;
                    subtotal += countP2Range(x, begin,
                        std::min(begin + CHUNK_SIZE, last), tables);''',
        '''                    std::int64_t begin, end;
                    if (!next.claim(begin, end)) break;
                    INITIAL_CHUNK_RANGE(begin, end);
                    subtotal += countP2Range(x, begin, end, tables);''')
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source_root', type=Path)
    p.add_argument('output_directory', type=Path)
    a = p.parse_args()
    original = a.source_root / 'src/P2.cpp'
    header = Path(__file__).with_name('InitialChunks.h')
    destination = a.output_directory
    files = [destination / n for n in ('P2Initial.cpp', 'InitialChunks.h', 'generation.json')]
    if any(p.exists() or p.is_symlink() for p in files):
        raise RuntimeError('refusing to overwrite generated outputs')
    source = original.read_text()
    generated = generate(source)
    destination.mkdir(parents=True, exist_ok=True)
    files[0].write_text(generated)
    files[1].write_bytes(header.read_bytes())
    files[2].write_text(json.dumps(dict(
        production_sha256=EXPECTED_P2_SHA256,
        generated_sha256=hashlib.sha256(generated.encode()).hexdigest(),
        header_sha256=hashlib.sha256(header.read_bytes()).hexdigest(),
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        changes='Only parallel task assignment; first32 use K, remaining use32; no barrier; serial path unchanged.',
        hooks='INITIAL_CHUNK_WORKER(slot), INITIAL_CHUNK_RANGE(begin,end), INITIAL_CHUNK_QUERY(i); defaults are empty.'
    ), indent=2) + '\n')
    print(files[0])


if __name__ == '__main__': main()
