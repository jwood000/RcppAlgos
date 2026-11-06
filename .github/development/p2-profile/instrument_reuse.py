#!/usr/bin/env python3
"""Prepare serial, cache-ineligible Phi reuse traces in diagnostic source copies."""

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
    source = '#include "reuse_metrics.h"\n' + source
    source = replace_once(source,
        '        const auto& primes = tables.primes();\n'
        '        const auto& pi = tables.pi();\n'
        '        std::int64_t sum = 0;\n',
        '        P2ReuseMetrics::range(first, last);\n'
        '        const auto& primes = tables.primes();\n'
        '        const auto& pi = tables.pi();\n'
        '        std::int64_t sum = 0;\n')
    source = replace_once(source,
        '        for (std::int64_t i = first; i < last; ++i) {\n',
        '        for (std::int64_t i = first; i < last; ++i) {\n'
        '            P2ReuseMetrics::query(i);\n')
    source = replace_once(source,
        '        PartialSieveDetail::Validate(x, a, tables);\n',
        '        PartialSieveDetail::Validate(x, a, tables);\n'
        '        P2ReuseMetrics::require_serial(nThreads);\n')
    return source


def instrument_phi(source):
    for literal, expected in [
        ('    constexpr int TINY_A = 6;', 1),
        ('    constexpr std::size_t CACHE_ROWS = 100;', 1),
        ('if (a >= CACHE_ROWS || x > std::numeric_limits<std::uint16_t>::max())', 2),
        ('            if (a == TINY_A) return phiTiny(x, TINY_A) * SIGN;', 1)
    ]:
        if source.count(literal) != expected:
            raise RuntimeError('Phi cache/entry assumptions no longer match: ' + literal)
    source = '#include "reuse_metrics.h"\n' + source
    source = replace_once(source,
        '        std::int64_t evaluate(std::int64_t x, std::int64_t a) {\n',
        '        std::int64_t evaluate(std::int64_t x, std::int64_t a) {\n'
        '            P2ReuseMetrics::Evaluation reuseEvaluation(x, a);\n')
    source = replace_once(source,
        '            // Every i < b has a prepared p[i + 1]. Quotients decrease and\n',
        '            P2ReuseMetrics::suffix_terms(b - i);\n'
        '            // Every i < b has a prepared p[i + 1]. Quotients decrease and\n')
    return source


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_directory', type=Path)
    parser.add_argument('--source-root', type=Path,
                        help='repository or snapshot containing src/')
    args = parser.parse_args()
    repository = Path(__file__).resolve().parents[3]
    source_root = args.source_root or repository
    output = args.output_directory.resolve()
    paths = {name: source_root / 'src' / name for name in ('P2.cpp', 'Phi.cpp')}
    header_path = Path(__file__).with_name('reuse_metrics.h')
    protected = list(paths.values()) + [repository / 'src' / name for name in paths]
    protected.append(header_path)
    for name in ('P2.cpp', 'Phi.cpp', 'reuse_metrics.h', 'reuse-instrumentation.json'):
        destination = output / name
        if destination.is_symlink():
            raise RuntimeError('diagnostic output must not follow a file symlink: ' +
                               str(destination))
        for original in protected:
            if (destination.resolve() == original.resolve() or
                    (destination.exists() and original.exists() and
                     destination.samefile(original))):
                raise RuntimeError('diagnostic output would overwrite a source file: ' +
                                   str(original))
    originals = {name: path.read_text() for name, path in paths.items()}
    copies = {'P2.cpp': instrument_p2(originals['P2.cpp']),
              'Phi.cpp': instrument_phi(originals['Phi.cpp'])}
    # Finish source-match checks before writing any file.
    output.mkdir(parents=True, exist_ok=True)
    for name, source in copies.items():
        (output / name).write_text(source)
    header = header_path.read_bytes()
    (output / 'reuse_metrics.h').write_bytes(header)
    manifest = {'source_root': str(source_root.resolve()),
                'sources': {name: hashlib.sha256(text.encode()).hexdigest()
                            for name, text in originals.items()},
                'header_sha256': hashlib.sha256(header).hexdigest(),
                'record_bytes': 32, 'serial_only': True,
                'record_predicate': 'a > 6 && (a >= 100 || x > 65535)'}
    (output / 'reuse-instrumentation.json').write_text(
        json.dumps(manifest, indent=2) + '\n')
    print(output)


if __name__ == '__main__':
    main()
