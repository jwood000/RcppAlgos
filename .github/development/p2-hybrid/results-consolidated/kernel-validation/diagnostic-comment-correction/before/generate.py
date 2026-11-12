#!/usr/bin/env python3
"""Generate common-K=1 temporary P2 kernels and an untimed-choice harness.

The production source and frozen experimental implementations are inputs only.
No generated output belongs in src/ or in the package's normal build.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXPERIMENT = HERE.parent
FROZEN = {
    'P2Hybrid.cpp': '8be6319d634edce0032d1afd41496c2c05c3fb2f0dc2ecc86ea97c47495fce82',
    'selector_v3_profile.cpp': '03d7cc03376facb83924562878a21155b2bb0b35564a107dfd071584acdc8585',
    'P2Selector.h': '2841a2e6e303179b386962bd18281bbbc0f78bed628dcbb7ad03784f400426f6',
    'P2SelectorV2.h': 'e953443f1aa2c9943365bc70cd92dc53996a8c0ab51ec3e675f55c380404cbd8',
    'P2SelectorV3.h': 'b8bf308dc79ac508871887fd9aad223a38b610129da8f8d2f9f09ddd7b233931',
    'selector_v3_cases.h': '3a5f3a9ae4c1adb0b85b53bf9657b119bbd428b7076f0e200e6ab502b12610b5',
    'initial-chunks/generate.py': '964d1eee960c91b703b72208e96f7adc35c963d1edaf059f7f606432583dcb07',
    'initial-chunks/InitialChunks.h': '6a68bef35210a7b2aad798923f03a7d7637c502d43c9f22bf3d480c72154788d',
}


def sha(content):
    if isinstance(content, str):
        content = content.encode()
    return hashlib.sha256(content).hexdigest()


def once(source, old, new):
    if source.count(old) != 1:
        raise RuntimeError('expected one source anchor: ' + old)
    return source.replace(old, new)


def frozen(name, source):
    if sha(source) != FROZEN[name]:
        raise RuntimeError('frozen experiment source changed: ' + name)


def generate_hybrid(source):
    frozen('P2Hybrid.cpp', source)
    result = '// Generated common-K=1 experiment; production remains unchanged.\n' + source
    result = once(result, '#include <atomic>', '#include "InitialChunks.h"')
    result = once(result,
        '// once per call. The four selection flags are mutually exclusive.',
        '// once per call. P2_HYBRID_SELECTOR_CANDIDATE selects the consolidated\n'
        '// policy; these five selection flags are mutually exclusive.')
    result = once(result,
        '     defined(P2_HYBRID_SELECTOR_V3) + \\\n',
        '     defined(P2_HYBRID_SELECTOR_V3) + defined(P2_HYBRID_SELECTOR_CANDIDATE) + \\\n')
    result = once(result,
        '#error "Select only one of P2_HYBRID_SELECTOR, P2_HYBRID_SELECTOR_V2, P2_HYBRID_SELECTOR_V3, or P2_HYBRID_MULTIPLIER"',
        '#error "Select only one P2 hybrid policy"')
    result = once(result,
        '#elif defined(P2_HYBRID_SELECTOR_V3)\n#include "P2SelectorV3.h"',
        '#elif defined(P2_HYBRID_SELECTOR_CANDIDATE)\n#include "P2SelectorCandidate.h"\n'
        '#elif defined(P2_HYBRID_SELECTOR_V3)\n#include "P2SelectorV3.h"')
    result = once(result,
        '#error "Define P2_HYBRID_SELECTOR, P2_HYBRID_SELECTOR_V2, P2_HYBRID_SELECTOR_V3, or P2_HYBRID_MULTIPLIER as 2, 4, 8, 16, 32, 64, or 128"',
        '#error "Define a P2 hybrid selector or fixed multiplier"')
    result = once(result,
        '''        // Retain production's dynamic 32-query chunks. The budget derives
        // from the full original non-table range, so removing cheaper tail
        // queries does not by itself reduce workers on the costly prefix.
        constexpr std::int64_t CHUNK_SIZE = 32;
        std::atomic<std::int64_t> next(first);''',
        '''        // Use the same queue as the Phi-only candidate: individual jobs
        // for the first 32 queries, then 32-query jobs, with no barrier.
        // The budget still derives from the full original non-table range.
        InitialChunks::Queue<1> next(first, last);''')
    result = once(result,
        '''                    const std::int64_t begin = next.fetch_add(CHUNK_SIZE,
                        std::memory_order_relaxed);
                    if (begin >= last) break;
                    subtotal += countP2Range(x, begin,
                        std::min(begin + CHUNK_SIZE, last), tables);''',
        '''                    std::int64_t begin, end;
                    if (!next.claim(begin, end)) break;
                    subtotal += countP2Range(x, begin, end, tables);''')
    positive = '''#if defined(P2_HYBRID_SELECTOR) || defined(P2_HYBRID_SELECTOR_V2) || \\
    defined(P2_HYBRID_SELECTOR_V3)'''
    if result.count(positive) != 2:
        raise RuntimeError('expected both hybrid selector guards')
    result = result.replace(positive, positive + ' || defined(P2_HYBRID_SELECTOR_CANDIDATE)')
    result = once(result,
        '#ifdef P2_HYBRID_SELECTOR_V3\n        const P2SelectorV3::Decision',
        '#ifdef P2_HYBRID_SELECTOR_CANDIDATE\n'
        '        const P2SelectorCandidate::Decision choice = P2SelectorCandidate::choose(\n'
        '            x, tables.limit(), first, last, phiBudget, nThreads, primes);\n'
        '#elif defined(P2_HYBRID_SELECTOR_V3)\n        const P2SelectorV3::Decision')
    result = once(result,
        '''#if !defined(P2_HYBRID_SELECTOR) && !defined(P2_HYBRID_SELECTOR_V2) && \\
    !defined(P2_HYBRID_SELECTOR_V3)''',
        '''#if !defined(P2_HYBRID_SELECTOR) && !defined(P2_HYBRID_SELECTOR_V2) && \\
    !defined(P2_HYBRID_SELECTOR_V3) && !defined(P2_HYBRID_SELECTOR_CANDIDATE)''')
    return result


def generate_profile(source):
    frozen('selector_v3_profile.cpp', source)
    result = '// Generated consolidated comparison; decisions remain outside timers.\n' + source
    result = once(result,
        '#ifdef P2_HYBRID_SELECTOR_V3\n#include "P2SelectorV3.h"\n#endif',
        '#ifdef P2_HYBRID_SELECTOR_V3\n#include "P2SelectorV3.h"\n#endif\n'
        '#ifdef P2_HYBRID_SELECTOR_CANDIDATE\n#include "P2SelectorCandidate.h"\n#endif')
    result = once(result,
        ' + defined(P2_HYBRID_MULTIPLIER)) != 1',
        ' + defined(P2_HYBRID_SELECTOR_CANDIDATE) + defined(P2_HYBRID_MULTIPLIER)) != 1')
    result = once(result,
        '#error "Select exactly one baseline, fixed multiplier, v1, v2, or v3 build"',
        '#error "Select exactly one baseline or hybrid policy build"')
    result = once(result,
        '#if defined(P2_SELECTOR_BASELINE)\n        return "baseline";',
        '#if defined(P2_SELECTOR_BASELINE) && defined(P2_INITIAL_CHUNK)\n'
        '        return "phi-initial-1";\n'
        '#elif defined(P2_SELECTOR_BASELINE)\n        return "baseline";\n'
        '#elif defined(P2_HYBRID_SELECTOR_CANDIDATE)\n        return "candidate";')
    start = '#ifdef P2_HYBRID_SELECTOR_V3\n        const auto choice = P2SelectorV3::choose('
    end = '#elif defined(P2_HYBRID_SELECTOR_V2)'
    position = result.index(start)
    finish = result.index(end, position)
    branch = result[position:finish]
    candidate = branch.replace('P2_HYBRID_SELECTOR_V3', 'P2_HYBRID_SELECTOR_CANDIDATE').replace(
        'P2SelectorV3::', 'P2SelectorCandidate::')
    result = once(result, branch, candidate + branch.replace('#ifdef ', '#elif defined(', 1).replace(
        'P2_HYBRID_SELECTOR_V3\n', 'P2_HYBRID_SELECTOR_V3)\n', 1))
    return result


def generate(source_root, output):
    inputs = {name: (EXPERIMENT / name).read_bytes() for name in FROZEN}
    for name, content in inputs.items():
        frozen(name, content)
    original = source_root / 'src/P2.cpp'
    spec = importlib.util.spec_from_file_location('initial_generate', EXPERIMENT / 'initial-chunks/generate.py')
    initial = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(initial)
    outputs = {
        'P2Initial.cpp': initial.generate(original.read_text()).encode(),
        'P2HybridInitial.cpp': generate_hybrid(inputs['P2Hybrid.cpp'].decode()).encode(),
        'selector_profile.cpp': generate_profile(inputs['selector_v3_profile.cpp'].decode()).encode(),
        'InitialChunks.h': inputs['initial-chunks/InitialChunks.h'],
    }
    for name in ('P2Selector.h', 'P2SelectorV2.h', 'P2SelectorV3.h', 'selector_v3_cases.h'):
        outputs[name] = inputs[name]
    if any((output / name).exists() or (output / name).is_symlink()
           for name in list(outputs) + ['generation.json']):
        raise RuntimeError('refusing to overwrite generated outputs')
    metadata = {
        'production_sha256': sha(original.read_bytes()),
        'frozen_experiment_sha256': FROZEN,
        'generator_sha256': sha(Path(__file__).read_bytes()),
        'output_sha256': {name: sha(content) for name, content in outputs.items()},
        'changes': 'Phi-only and hybrid prefix share K=1 first32/32-query remainder mapping. '
                   'No new barrier, worker-budget change, serial-path change or counting change. '
                   'Hybrid Phi/stream phase barrier remains unchanged. '
                   'Candidate selection and diagnostic branches added; frozen policies unmodified.',
        'candidate': 'P2SelectorCandidate.h supplied separately through include path; not generated.',
    }
    output.mkdir(parents=True, exist_ok=True)
    for name, content in outputs.items():
        (output / name).write_bytes(content)
    (output / 'generation.json').write_text(json.dumps(metadata, indent=2) + '\n')
    return metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source_root', type=Path)
    parser.add_argument('output_directory', type=Path)
    args = parser.parse_args()
    generate(args.source_root.resolve(), args.output_directory.resolve())
    print(args.output_directory / 'generation.json')


if __name__ == '__main__':
    main()
