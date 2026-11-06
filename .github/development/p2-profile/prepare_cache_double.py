#!/usr/bin/env python3
"""Prepare temporary Phi cache variants without editing production sources.

Ordinary copies are for timing. Optional memory copies expose a page scan that
must be called only after all evaluations have finished; they are not timed.
"""

import argparse
import hashlib
import itertools
import json
from pathlib import Path


HEADER = """#ifndef RCPPALGOS_CACHE_MEMORY_PROBE_H
#define RCPPALGOS_CACHE_MEMORY_PROBE_H
#include <cstddef>

namespace PrimeCounting {
    struct PhiCacheMemoryStats {
        std::size_t allocatedPages;
        std::size_t pageBytes;
        std::size_t allocatedPageBytes;
        std::size_t cacheObjectBytes;
        std::size_t possiblePages;
        std::size_t possiblePageBytes;
    };

    // Diagnostic copies only. Call after every worker/evaluation has finished.
    // Bytes use sizeof; allocator bookkeeping and tiny Phi tables are excluded.
    PhiCacheMemoryStats PhiCacheStats();
}
#endif
"""


def replace_exact(source, old, new, expected=1):
    count = source.count(old)
    if count != expected:
        raise RuntimeError(f"expected {expected} source matches, found {count}: {old[:100]!r}")
    return source.replace(old, new)


def six_prime_bound(x):
    """Independent finite inclusion-exclusion, not the production Phi kernel."""
    result = 0
    for size in range(7):
        for subset in itertools.combinations((2, 3, 5, 7, 11, 13), size):
            product = 1
            for prime in subset:
                product *= prime
            result += (-1 if size % 2 else 1) * (x // product)
    return result


def doubled(source):
    for literal in (
        '    constexpr int TINY_A = 6;',
        '    constexpr std::size_t CACHE_ROWS = 100;',
        '    constexpr unsigned PAGE_BITS = 12;',
        '    constexpr std::size_t PAGE_SIZE = 1u << PAGE_BITS;',
        '            if (a == TINY_A) return phiTiny(x, TINY_A) * SIGN;\n\n'
        '            const std::int64_t cached = cache.get(x, a);',
        '        if (a <= TINY_A) return phiTiny(x, a);',
    ):
        if source.count(literal) != 1:
            raise RuntimeError('cache proof/source assumptions changed: ' + literal)
    if six_prime_bound(131071) != 25142:
        raise RuntimeError('independent six-prime bound mismatch')
    source = replace_exact(source,
        '    //   PhiCacheRow -> 16 atomic page pointers, indexed by x / 4096',
        '    //   PhiCacheRow -> 32 atomic page pointers, indexed by x / 4096')
    source = replace_exact(source,
        '    // For positive x, 1 <= Phi(x, a) <= x, so cached answers fit in uint16_t.\n',
        '    // Only a >= 7 reaches this cache: a <= 6 uses the tiny tables.\n'
        '    // For positive x <= 131071 and a >= 7, monotonicity gives\n'
        '    // 1 <= Phi(x, a) <= Phi(131071, 6) = 25142 < 65536.\n'
        '    // The bound follows independently by six-prime inclusion-exclusion,\n'
        '    // so doubling the argument range preserves uint16_t answer cells.\n')
    source = replace_exact(source,
        '    constexpr std::size_t PAGE_COUNT = 65536u / PAGE_SIZE;',
        '    constexpr std::uint64_t CACHE_MAX_X = 131071;\n'
        '    constexpr std::size_t PAGE_COUNT = (CACHE_MAX_X + 1) / PAGE_SIZE;\n'
        '    static_assert((CACHE_MAX_X + 1) % PAGE_SIZE == 0,\n'
        '                  "cache pages must cover the complete argument range");')
    return replace_exact(source,
        'if (a >= CACHE_ROWS || x > std::numeric_limits<std::uint16_t>::max())',
        'if (a >= CACHE_ROWS || x > CACHE_MAX_X)', 2)


def memory_probe(source):
    source = '#include "cache_memory_probe.h"\n' + source
    source = replace_exact(source,
        '    public:\n        std::int64_t get(std::uint64_t x, std::uint64_t a) const {',
        '    public:\n'
        '        PhiCacheMemoryStats memoryStats() const {\n'
        '            // The caller has joined all workers before this diagnostic scan.\n'
        '            // No cache entries/pages can change during the snapshot.\n'
        '            std::size_t allocated = 0;\n'
        '            for (const auto& row : rows) {\n'
        '                for (const auto& published : row.pages) {\n'
        '                    allocated += published.value.load(std::memory_order_relaxed)\n'
        '                        != nullptr;\n'
        '                }\n'
        '            }\n'
        '            const std::size_t possible = CACHE_ROWS * PAGE_COUNT;\n'
        '            return {allocated, sizeof(PhiPage), allocated * sizeof(PhiPage),\n'
        '                    sizeof(PhiCache), possible, possible * sizeof(PhiPage)};\n'
        '        }\n\n'
        '        std::int64_t get(std::uint64_t x, std::uint64_t a) const {')
    return replace_exact(source,
        '    PhiTables::PhiTables(std::int64_t limit) {',
        '    PhiCacheMemoryStats PhiCacheStats() {\n'
        '        return sharedCache().memoryStats();\n'
        '    }\n\n'
        '    PhiTables::PhiTables(std::int64_t limit) {')


def validate_output(output, source_root, repository, names, protected):
    for root in (source_root, repository):
        if output == root or root in output.parents:
            raise RuntimeError('use an output directory outside the repository/source snapshot')
    if output.exists() and any(path.is_symlink() for path in output.rglob('*')):
        raise RuntimeError('output directory must not contain symlinks')
    for name in names:
        destination = output / name
        if output not in destination.resolve().parents:
            raise RuntimeError('output destination escapes artifact directory: ' + str(destination))
        if destination.exists() and destination.stat().st_nlink > 1:
            raise RuntimeError('output destination must not be a hard link: ' + str(destination))
        for original in protected:
            if (destination.resolve() == original.resolve() or
                    (destination.exists() and original.exists() and destination.samefile(original))):
                raise RuntimeError('output would overwrite a source file: ' + str(original))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--source-root', type=Path,
                        help='repository or snapshot containing src/Phi.cpp')
    parser.add_argument('--memory-probes', action='store_true')
    args = parser.parse_args()
    repository = Path(__file__).resolve().parents[3]
    source_root = (args.source_root or repository).resolve()
    source_path = source_root / 'src/Phi.cpp'
    original = source_path.read_text()
    candidate = doubled(original)
    copies = {'baseline/Phi.cpp': original, 'double/Phi.cpp': candidate}
    if args.memory_probes:
        copies.update({'baseline-memory/Phi.cpp': memory_probe(original),
                       'double-memory/Phi.cpp': memory_probe(candidate),
                       'cache_memory_probe.h': HEADER})
    output = args.output.resolve()
    manifest_name = 'cache-double-manifest.json'
    validate_output(output, source_root, repository, list(copies) + [manifest_name],
                    [source_path, repository / 'src/Phi.cpp', Path(__file__)])
    manifest = {
        'source_root': str(source_root),
        'source_sha256': hashlib.sha256(original.encode()).hexdigest(),
        'baseline_max_x': 65535, 'candidate_max_x': 131071,
        'cache_rows': 100, 'page_size': 4096,
        'baseline_page_count': 16, 'candidate_page_count': 32,
        'answer_bits': 16, 'six_prime_upper_bound': six_prime_bound(131071),
        'memory_probes': args.memory_probes,
        'generated_sha256': {name: hashlib.sha256(source.encode()).hexdigest()
                             for name, source in copies.items()},
    }
    # Finish all source/guard checks before creating or changing any file.
    for name, source in copies.items():
        path = output / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(source)
    (output / manifest_name).write_text(json.dumps(manifest, indent=2) + '\n')
    print(output)


if __name__ == '__main__':
    main()
