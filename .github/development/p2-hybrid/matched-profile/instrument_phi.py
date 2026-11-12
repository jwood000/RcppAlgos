#!/usr/bin/env python3
"""Generate a counts-only Phi copy; refuse drift, source aliases, and overwrites."""

import argparse
import hashlib
import json
from pathlib import Path


def replace_once(source, old, new):
    found = source.count(old)
    if found != 1:
        raise RuntimeError(f"expected exactly one source match, found {found}: {old[:100]!r}")
    return source.replace(old, new)


def instrument_phi(source):
    for literal, expected in [
        ('    constexpr int TINY_A = 6;', 1),
        ('    constexpr std::size_t CACHE_ROWS = 100;', 1),
        ('    constexpr unsigned PAGE_BITS = 12;', 1),
        ('    constexpr std::uint64_t CACHE_MAX_X = 131071;', 1),
        ('if (a >= CACHE_ROWS || x > CACHE_MAX_X)', 2),
    ]:
        if source.count(literal) != expected:
            raise RuntimeError('Phi cache assumptions differ: ' + literal)
    source = ('#include "phi_metrics.h"\n'
              'namespace MatchedPhi {\n'
              '    thread_local Stats* active = nullptr;\n'
              '    thread_local unsigned depth = 0;\n'
              '}\n' + source)
    source = replace_once(source,
        '    std::int64_t phiTiny(std::int64_t x, std::int64_t a) {\n',
        '    std::int64_t phiTiny(std::int64_t x, std::int64_t a) {\n'
        '        MatchedPhi::add(MatchedPhi::TinyCalls);\n')
    source = replace_once(source,
        '                return 0;\n            }\n\n            // The shift and mask',
        '                MatchedPhi::cacheGet(false, false, false);\n'
        '                return 0;\n            }\n\n            // The shift and mask')
    source = replace_once(source,
        '            return page ? page->cells[x & (PAGE_SIZE - 1)].value.load(\n'
        '                std::memory_order_relaxed\n            ) : 0;\n',
        '            const std::int64_t result = page ? page->cells[x & (PAGE_SIZE - 1)].value.load(\n'
        '                std::memory_order_relaxed\n            ) : 0;\n'
        '            MatchedPhi::cacheGet(true, result != 0, page != nullptr);\n'
        '            return result;\n')
    source = replace_once(source,
        '        void put(std::uint64_t x, std::uint64_t a, std::int64_t value) {\n',
        '        void put(std::uint64_t x, std::uint64_t a, std::int64_t value) {\n'
        '            MatchedPhi::add(MatchedPhi::PutCalls);\n')
    source = replace_once(source,
        '            auto& row = rows[a];\n',
        '            MatchedPhi::add(MatchedPhi::PutEligible);\n'
        '            auto& row = rows[a];\n')
    source = replace_once(source,
        '                std::lock_guard<std::mutex> lock(row.allocationMutex);\n',
        '                MatchedPhi::add(MatchedPhi::LockAttempts);\n'
        '                std::lock_guard<std::mutex> lock(row.allocationMutex);\n')
    source = replace_once(source,
        '                    page = new (std::nothrow) PhiPage;\n'
        '                    if (!page) return; // Caching is optional.\n',
        '                    MatchedPhi::add(MatchedPhi::PageAllocAttempts);\n'
        '                    page = new (std::nothrow) PhiPage;\n'
        '                    if (!page) {\n'
        '                        MatchedPhi::add(MatchedPhi::PageAllocFailed);\n'
        '                        return; // Caching is optional.\n'
        '                    }\n'
        '                    MatchedPhi::add(MatchedPhi::PageAllocSuccess);\n')
    source = replace_once(source,
        '            page->cells[x & (PAGE_SIZE - 1)].value.store(\n',
        '            MatchedPhi::add(MatchedPhi::CacheStores);\n'
        '            page->cells[x & (PAGE_SIZE - 1)].value.store(\n')
    source = replace_once(source,
        '        std::int64_t evaluate(std::int64_t x, std::int64_t a) {\n'
        '            if (a == TINY_A) return phiTiny(x, TINY_A) * SIGN;\n',
        '        std::int64_t evaluate(std::int64_t x, std::int64_t a) {\n'
        '            MatchedPhi::Evaluation profileEvaluation;\n'
        '            if (a == TINY_A) {\n'
        '                MatchedPhi::add(MatchedPhi::KernelTiny);\n'
        '                return phiTiny(x, TINY_A) * SIGN;\n'
        '            }\n')
    source = replace_once(source,
        '            const std::int64_t root = integerSqrt(x);\n',
        '            MatchedPhi::add(MatchedPhi::Expansions);\n'
        '            const std::int64_t root = integerSqrt(x);\n')
    source = replace_once(source,
        '            // Every i < b has a prepared p[i + 1]. Quotients decrease and\n',
        '            MatchedPhi::add(MatchedPhi::PiSuffixTerms, b - i);\n'
        '            // Every i < b has a prepared p[i + 1]. Quotients decrease and\n')
    source = replace_once(source,
        '    std::int64_t Phi(std::int64_t x, std::int64_t a, const PhiTables& tables) {\n',
        '    std::int64_t Phi(std::int64_t x, std::int64_t a, const PhiTables& tables) {\n'
        '        MatchedPhi::add(MatchedPhi::PublicPhi);\n')
    public_start = source.index('    std::int64_t Phi(std::int64_t x, std::int64_t a,')
    public = source[public_start:]
    for condition, expression, counter in [
        ('x == 0', '0', 'PublicZero'), ('a == 0', 'x', 'PublicAZero'),
        ('x <= tables.primes_[a]', '1', 'PublicOne'),
        ('a <= TINY_A', 'phiTiny(x, a)', 'PublicTiny')
    ]:
        public = replace_once(public,
            f'        if ({condition}) return {expression};\n',
            f'        if ({condition}) {{\n'
            f'            MatchedPhi::add(MatchedPhi::{counter});\n'
            f'            return {expression};\n        }}\n')
    source = source[:public_start] + public
    source = replace_once(source,
        '            return static_cast<std::int64_t>(tables.pi_[x]) - a + 1;\n',
        '            MatchedPhi::add(MatchedPhi::PublicPi);\n'
        '            return static_cast<std::int64_t>(tables.pi_[x]) - a + 1;\n')
    return source


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_directory', type=Path)
    parser.add_argument('--source-root', type=Path, required=True,
                        help='frozen source snapshot containing src/Phi.cpp')
    args = parser.parse_args()
    repository = Path(__file__).resolve().parents[4]
    original = args.source_root / 'src' / 'Phi.cpp'
    header = Path(__file__).with_name('phi_metrics.h')
    protected = [original, repository / 'src' / 'Phi.cpp', header]
    output = args.output_directory.resolve()
    destinations = [output / name for name in ('Phi.cpp', 'phi_metrics.h', 'phi-instrumentation.json')]
    for destination in destinations:
        if destination.exists() or destination.is_symlink():
            raise RuntimeError('refusing to overwrite existing output: ' + str(destination))
        if any(destination.resolve() == path.resolve() for path in protected):
            raise RuntimeError('diagnostic output aliases a source: ' + str(destination))
    source = original.read_text()
    generated = instrument_phi(source)  # Complete match checks before any writes.
    header_bytes = header.read_bytes()
    manifest = {
        'source_root': str(args.source_root.resolve()),
        'original_phi_sha256': hashlib.sha256(source.encode()).hexdigest(),
        'instrumented_phi_sha256': hashlib.sha256(generated.encode()).hexdigest(),
        'header_sha256': hashlib.sha256(header_bytes).hexdigest(),
        'generator_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'diagnostics': 'per-thread counts only; no clocks, global counters or diagnostic locks',
        'warning': 'instrumentation changes execution costs and schedules; do not time this build',
        'cache': {'argument_limit': 131071, 'rows': 100, 'page_bits': 12},
    }
    output.mkdir(parents=True, exist_ok=True)
    destinations[0].write_text(generated)
    destinations[1].write_bytes(header_bytes)
    destinations[2].write_text(json.dumps(manifest, indent=2) + '\n')
    print(output)


if __name__ == '__main__':
    main()
