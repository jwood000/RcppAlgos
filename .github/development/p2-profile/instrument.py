#!/usr/bin/env python3
"""Write diagnostic P2/Phi copies with checked, narrowly scoped substitutions."""

import argparse
import hashlib
import json
from pathlib import Path


def replace_once(source, old, new):
    count = source.count(old)
    if count != 1:
        raise RuntimeError(f"expected one source match, found {count}: {old[:100]!r}")
    return source.replace(old, new)


def instrument_p2(source):
    source = '#include "profile_metrics.h"\n' + source
    source = replace_once(source,
        '        const auto& primes = tables.primes();\n'
        '        const auto& pi = tables.pi();\n'
        '        std::int64_t sum = 0;\n',
        '        P2ProfileMetrics::Chunk profileChunk(first, last);\n'
        '        const auto& primes = tables.primes();\n'
        '        const auto& pi = tables.pi();\n'
        '        std::int64_t sum = 0;\n')
    source = replace_once(source,
        '        for (std::int64_t i = first; i < last; ++i) {\n',
        '        for (std::int64_t i = first; i < last; ++i) {\n'
        '            P2ProfileMetrics::Query profileQuery(i);\n')
    source = replace_once(source,
        '            const std::int64_t c = pi[PartialSieveDetail::IntegerSqrt(quotient)];\n',
        '            const std::int64_t c = pi[PartialSieveDetail::IntegerSqrt(quotient)];\n'
        '            profileQuery.root_done();\n')
    source = replace_once(source,
        '            const std::int64_t piQuotient = Phi(quotient, c, tables) + c - 1;\n',
        '            const std::int64_t piQuotient = Phi(quotient, c, tables) + c - 1;\n'
        '            profileQuery.phi_done();\n')
    source = replace_once(source,
        '        PartialSieveDetail::Validate(x, a, tables);\n',
        '        PartialSieveDetail::Validate(x, a, tables);\n'
        '        P2ProfileMetrics::Call profileCall(x, a, nThreads);\n')
    source = replace_once(source,
        '        for (std::int64_t i = last; i <= b; ++i) {\n',
        '        profileCall.range(first, last, b);\n'
        '        profileCall.start_direct();\n'
        '        for (std::int64_t i = last; i <= b; ++i) {\n')
    source = replace_once(source,
        '        if (last == first) return sum;\n',
        '        profileCall.end_direct();\n'
        '        if (last == first) return sum;\n')
    source = replace_once(source,
        '        if (workers == 1) return sum + countP2Range(x, first, last, tables);\n',
        '        profileCall.prepare(workers);\n'
        '        if (workers == 1) {\n'
        '            P2ProfileMetrics::Worker profileWorker(profileCall, 0);\n'
        '            return sum + countP2Range(x, first, last, tables);\n'
        '        }\n')
    source = replace_once(source,
        '            batch.launch([&, slot] {\n',
        '            batch.launch([&, slot] {\n'
        '                P2ProfileMetrics::Worker profileWorker(profileCall, slot);\n')
    return source


def instrument_phi(source):
    source = '#include "profile_metrics.h"\n' + source
    source = replace_once(source,
        '    std::int64_t phiTiny(std::int64_t x, std::int64_t a) {\n',
        '    std::int64_t phiTiny(std::int64_t x, std::int64_t a) {\n'
        '        P2ProfileMetrics::add(P2ProfileMetrics::TinyCalls);\n')
    source = replace_once(source,
        '                return 0;\n            }\n\n            // The shift and mask',
        '                P2ProfileMetrics::cache_get(false, false, false);\n'
        '                return 0;\n            }\n\n            // The shift and mask')
    source = replace_once(source,
        '            return page ? page->cells[x & (PAGE_SIZE - 1)].value.load(\n'
        '                std::memory_order_relaxed\n            ) : 0;\n',
        '            const std::int64_t result = page ? page->cells[x & (PAGE_SIZE - 1)].value.load(\n'
        '                std::memory_order_relaxed\n            ) : 0;\n'
        '            P2ProfileMetrics::cache_get(true, result != 0, page != nullptr);\n'
        '            return result;\n')
    source = replace_once(source,
        '        void put(std::uint64_t x, std::uint64_t a, std::int64_t value) {\n',
        '        void put(std::uint64_t x, std::uint64_t a, std::int64_t value) {\n'
        '            P2ProfileMetrics::add(P2ProfileMetrics::PutCalls);\n')
    source = replace_once(source,
        '            auto& row = rows[a];\n',
        '            P2ProfileMetrics::add(P2ProfileMetrics::PutEligible);\n'
        '            auto& row = rows[a];\n')
    source = replace_once(source,
        '                std::lock_guard<std::mutex> lock(row.allocationMutex);\n',
        '                P2ProfileMetrics::add(P2ProfileMetrics::LockAttempts);\n'
        '                std::lock_guard<std::mutex> lock(row.allocationMutex);\n')
    source = replace_once(source,
        '                    page = new (std::nothrow) PhiPage;\n'
        '                    if (!page) return; // Caching is optional.\n',
        '                    P2ProfileMetrics::add(P2ProfileMetrics::PageAllocAttempts);\n'
        '                    page = new (std::nothrow) PhiPage;\n'
        '                    if (!page) {\n'
        '                        P2ProfileMetrics::add(P2ProfileMetrics::PageAllocFailed);\n'
        '                        return; // Caching is optional.\n'
        '                    }\n'
        '                    P2ProfileMetrics::add(P2ProfileMetrics::PageAllocSuccess);\n')
    source = replace_once(source,
        '            page->cells[x & (PAGE_SIZE - 1)].value.store(\n',
        '            P2ProfileMetrics::add(P2ProfileMetrics::CacheStores);\n'
        '            page->cells[x & (PAGE_SIZE - 1)].value.store(\n')
    source = replace_once(source,
        '        std::int64_t evaluate(std::int64_t x, std::int64_t a) {\n'
        '            if (a == TINY_A) return phiTiny(x, TINY_A) * SIGN;\n',
        '        std::int64_t evaluate(std::int64_t x, std::int64_t a) {\n'
        '            P2ProfileMetrics::Evaluation profileEvaluation;\n'
        '            if (a == TINY_A) {\n'
        '                P2ProfileMetrics::add(P2ProfileMetrics::ExitTiny);\n'
        '                return phiTiny(x, TINY_A) * SIGN;\n'
        '            }\n')
    source = replace_once(source,
        '            const std::int64_t root = integerSqrt(x);\n',
        '            P2ProfileMetrics::add(P2ProfileMetrics::Expansions);\n'
        '            const std::int64_t root = integerSqrt(x);\n')
    source = replace_once(source,
        '            // Every i < b has a prepared p[i + 1]. Quotients decrease and\n',
        '            P2ProfileMetrics::add(P2ProfileMetrics::PiSuffixTerms, b - i);\n'
        '            // Every i < b has a prepared p[i + 1]. Quotients decrease and\n')
    source = replace_once(source,
        '    std::int64_t Phi(std::int64_t x, std::int64_t a, const PhiTables& tables) {\n',
        '    std::int64_t Phi(std::int64_t x, std::int64_t a, const PhiTables& tables) {\n'
        '        P2ProfileMetrics::add(P2ProfileMetrics::PublicPhi);\n')
    public_start = source.index('    std::int64_t Phi(std::int64_t x, std::int64_t a,')
    public = source[public_start:]
    for condition, expression, counter in [
        ('x == 0', '0', 'ExitZero'), ('a == 0', 'x', 'ExitAZero'),
        ('x <= tables.primes_[a]', '1', 'ExitOne'),
        ('a <= TINY_A', 'phiTiny(x, a)', 'ExitTiny')
    ]:
        public = replace_once(public,
            f'        if ({condition}) return {expression};\n',
            f'        if ({condition}) {{\n'
            f'            P2ProfileMetrics::add(P2ProfileMetrics::{counter});\n'
            f'            return {expression};\n        }}\n')
    source = source[:public_start] + public
    source = replace_once(source,
        '            return static_cast<std::int64_t>(tables.pi_[x]) - a + 1;\n',
        '            P2ProfileMetrics::add(P2ProfileMetrics::ExitPi);\n'
        '            return static_cast<std::int64_t>(tables.pi_[x]) - a + 1;\n')
    return source


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_directory', type=Path)
    parser.add_argument('--source-root', type=Path,
                        help='repository or artifact snapshot containing src/')
    args = parser.parse_args()
    repository = Path(__file__).resolve().parents[3]
    source_root = args.source_root or repository
    output = args.output_directory.resolve()
    if output in {(repository / 'src').resolve(), (source_root / 'src').resolve()}:
        parser.error('diagnostic output must not overwrite production or snapshot sources')
    originals = {name: (source_root / 'src' / name).read_text()
                 for name in ('P2.cpp', 'Phi.cpp')}
    copies = {'P2.cpp': instrument_p2(originals['P2.cpp']),
              'Phi.cpp': instrument_phi(originals['Phi.cpp'])}
    # Complete all match checks before writing any diagnostic source.
    args.output_directory.mkdir(parents=True, exist_ok=True)
    for name, source in copies.items():
        (args.output_directory / name).write_text(source)
    header = Path(__file__).with_name('profile_metrics.h').read_bytes()
    (args.output_directory / 'profile_metrics.h').write_bytes(header)
    manifest = {'source_root': str(source_root.resolve()),
                'sources': {name: hashlib.sha256(text.encode()).hexdigest()
                            for name, text in originals.items()},
                'header_sha256': hashlib.sha256(header).hexdigest(),
                'counter_mode': '-DP2_PROFILE_METRICS',
                'timing_mode': '-DP2_PROFILE_METRICS -DP2_PROFILE_TIMING'}
    (args.output_directory / 'instrumentation.json').write_text(
        json.dumps(manifest, indent=2) + '\n')
    print(args.output_directory)


if __name__ == '__main__':
    main()
