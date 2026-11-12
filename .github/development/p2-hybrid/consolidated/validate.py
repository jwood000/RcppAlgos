#!/usr/bin/env python3
"""Sanitizer/correctness validation for common-K=1 experimental P2 kernels.

Run before timing; deliberately does not collect performance observations.
The common check remains the original independent/differential experiment.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import subprocess

HERE = Path(__file__).resolve().parent
EXPERIMENT = HERE.parent
SANITIZERS = {'asan-ubsan': 'address,undefined', 'tsan': 'thread'}
VARIANTS = {
    'baseline': ('production', []),
    'phi1': ('initial', ['-DP2_INITIAL_CHUNK=1']),
    'v1': ('hybrid', ['-DP2_HYBRID_SELECTOR']),
    'v2': ('hybrid', ['-DP2_HYBRID_SELECTOR_V2']),
    'v3': ('hybrid', ['-DP2_HYBRID_SELECTOR_V3']),
    'fixed8': ('hybrid', ['-DP2_HYBRID_MULTIPLIER=8']),
    'fixed16': ('hybrid', ['-DP2_HYBRID_MULTIPLIER=16']),
    'fixed32': ('hybrid', ['-DP2_HYBRID_MULTIPLIER=32']),
    'fixed64': ('hybrid', ['-DP2_HYBRID_MULTIPLIER=64']),
    'candidate': ('hybrid', ['-DP2_HYBRID_SELECTOR_CANDIDATE']),
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def invoke(command, log):
    result = subprocess.run(list(map(str, command)), capture_output=True, text=True, timeout=300)
    log.write_text(result.stdout + result.stderr)
    if result.returncode or result.stderr:
        raise RuntimeError(f'failed {command}\n{result.stdout}\n{result.stderr}')
    return result.stdout


def generator_checks(generator, source, generated):
    checks = 0
    for name, method in [('P2Hybrid.cpp', generator.generate_hybrid),
                         ('selector_v3_profile.cpp', generator.generate_profile)]:
        original = (EXPERIMENT / name).read_text()
        for mutated in (original + '\n', original.replace('std::int64_t', 'long long', 1)):
            try:
                method(mutated)
            except RuntimeError:
                checks += 1
            else:
                raise RuntimeError('generator accepted source drift: ' + name)
    before = {p.name: sha(p) for p in generated.iterdir() if p.is_file()}
    try:
        generator.generate(source, generated)
    except RuntimeError as error:
        if 'refusing to overwrite' not in str(error):
            raise
        checks += 1
    else:
        raise RuntimeError('generator accepted existing outputs')
    if before != {p.name: sha(p) for p in generated.iterdir() if p.is_file()}:
        raise RuntimeError('rejected generation changed outputs')
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    source = root / 'source'
    out = root / 'kernel-validation'
    if out.exists():
        raise RuntimeError('refusing existing validation directory')
    candidate = HERE / 'P2SelectorCandidate.h'
    if not candidate.is_file():
        raise RuntimeError('candidate header must be frozen before validation')
    spec = importlib.util.spec_from_file_location('consolidated_generate', HERE / 'generate.py')
    generator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(generator)
    generated = out / 'generated'
    generator.generate(source, generated)
    (generated / candidate.name).write_bytes(candidate.read_bytes())
    rejection_checks = generator_checks(generator, source, generated)
    units = ('Phi.cpp', 'PhiParallel.cpp', 'P3.cpp', 'PiMeissel.cpp', 'PiLehmer.cpp', 'Eratosthenes.cpp')
    commands, checks = [], []
    for name, sanitizer in SANITIZERS.items():
        directory = out / name
        directory.mkdir()
        flags = ['clang++', '-std=c++11', '-O1', '-g', '-pthread', '-Wall', '-Wextra', '-Wpedantic',
                 '-fno-omit-frame-pointer', '-fno-sanitize-recover=all', '-fsanitize=' + sanitizer,
                 '-I' + str(source / 'inst/include'), '-I' + str(generated), '-I' + str(EXPERIMENT)]
        common = []
        for path, defs in ([(source / 'src' / n, []) for n in units] +
                           [(source / 'src/P2.cpp', ['-DP2=P2Reference']), (EXPERIMENT / 'check.cpp', [])]):
            target = directory / ('P2Reference.o' if defs else path.stem + '.o')
            command = flags + defs + ['-c', path, '-o', target]
            invoke(command, directory / (target.name + '.log'))
            commands.append(list(map(str, command)))
            common.append(target)
        sources = {'production': source / 'src/P2.cpp', 'initial': generated / 'P2Initial.cpp',
                   'hybrid': generated / 'P2HybridInitial.cpp'}
        for variant, (kind, defs) in VARIANTS.items():
            target = directory / ('check-' + variant)
            command = flags + defs + [sources[kind]] + common + ['-o', target]
            invoke(command, directory / ('build-' + variant + '.log'))
            commands.append(list(map(str, command)))
            result = invoke([target], directory / ('run-' + variant + '.log'))
            commands.append([str(target)])
            match = re.fullmatch(r'small_checks=(\d+)\nPASS checks=(\d+)\n', result)
            if not match:
                raise RuntimeError('unexpected common checker output: ' + result)
            checks.append(dict(sanitizer=name, variant=variant, checks=int(match[2]),
                               small_checks=int(match[1])))
            print(name, variant, result.strip(), flush=True)
        target = directory / 'check-queue'
        command = flags + [EXPERIMENT / 'initial-chunks/queue-check.cpp', '-o', target]
        invoke(command, directory / 'build-queue.log')
        commands.append(list(map(str, command)))
        result = invoke([target], directory / 'run-queue.log')
        commands.append([str(target)])
        match = re.fullmatch(r'PASS queue_checks=(\d+)\n', result)
        if not match:
            raise RuntimeError('unexpected queue checker output: ' + result)
        checks.append(dict(sanitizer=name, variant='queue', checks=int(match[1])))
        print(name, result.strip(), flush=True)
        target = directory / 'check-selector'
        command = flags + [HERE / 'selector-check.cpp', '-o', target]
        invoke(command, directory / 'build-selector.log')
        commands.append(list(map(str, command)))
        result = invoke([target], directory / 'run-selector.log')
        commands.append([str(target)])
        match = re.fullmatch(r'PASS selector_checks=(\d+)\n', result)
        if not match:
            raise RuntimeError('unexpected selector checker output: ' + result)
        checks.append(dict(sanitizer=name, variant='selector', checks=int(match[1])))
        print(name, result.strip(), flush=True)
    paths = [HERE / n for n in ('generate.py', 'validate.py', 'selector-check.cpp', 'P2SelectorCandidate.h')]
    paths += [EXPERIMENT / 'check.cpp', EXPERIMENT / 'initial-chunks/queue-check.cpp']
    paths += [source / 'src' / n for n in units + ('P2.cpp',)]
    paths += [p for p in generated.iterdir() if p.is_file()]
    summary = dict(checks=checks, generator_rejections=rejection_checks, commands=commands,
                   input_sha256={str(path): sha(path) for path in paths},
                   scope='Unchanged independent factor oracle, quotient/segment boundaries, '
                         'large oracle stream calls, production differential checks, concurrent '
                         'P2 calls, input rejection, and complete Meissel/Lehmer counts. '
                         'Pure queue tests include overflow, exactly-once coverage and no barrier. '
                         'No benchmark observations collected.')
    (out / 'validation.json').write_text(json.dumps(summary, indent=2) + '\n')
    print('saved', out / 'validation.json', flush=True)


if __name__ == '__main__':
    main()
