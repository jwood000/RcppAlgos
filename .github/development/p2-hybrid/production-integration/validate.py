#!/usr/bin/env python3
"""Run the retained P2 checks against current production, in an isolated directory.

No historical checkout, generated variant, or package-source mutation is needed.
The small frozen fixtures supplement independent oracles; see fixtures/README.md.
"""
import argparse
import hashlib
import json
import os
import re
from pathlib import Path
import shutil
import subprocess

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
REPO = HERE.parents[3]
UNITS = ('Phi.cpp', 'PhiParallel.cpp', 'P3.cpp', 'PiMeissel.cpp',
         'PiLehmer.cpp', 'Eratosthenes.cpp')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path, help='Artifact root; native-validation/ must not exist')
    args = parser.parse_args()
    out = args.output.resolve() / 'native-validation'
    if out.exists():
        raise RuntimeError('refusing existing native-validation directory')
    source = out / 'source'
    (source / 'src').mkdir(parents=True)
    shutil.copytree(REPO / 'inst/include', source / 'inst/include')
    for name in UNITS + ('P2.cpp',):
        shutil.copy2(REPO / 'src' / name, source / 'src' / name)
    shutil.copy2(source / 'src/P2.cpp', out / 'ProductionP2.cpp')
    shutil.copytree(STUDY / 'fixtures', out / 'fixtures')
    # The frozen P2 reference used the former header path/namespace. Adapt
    # its temporary build environment without rewriting the historical source.
    (source / 'inst/include/NumbersUtils/ThreadBatch.h').write_text(
        '#pragma once\n#include "../ThreadBatch.h"\n'
        'namespace PrimeCounting { using RcppAlgos::ThreadBatch; }\n')
    shutil.copy2(HERE / 'helpers.cpp', out / 'helpers.cpp')
    shutil.copy2(STUDY / 'check.cpp', out / 'check.cpp')
    shutil.copy2(Path(__file__), out / 'validate.py')
    checker = (out / 'check.cpp').read_text()
    if checker.count('#include "P2Hybrid.cpp"') != 1:
        raise RuntimeError('private stream-check include changed')
    (out / 'stream-check.cpp').write_text(checker.replace(
        '#include "P2Hybrid.cpp"', '#include "ProductionP2.cpp"'))
    inputs = {str(p.relative_to(out)): digest(p) for p in out.rglob('*') if p.is_file()}
    commands, results = [], {}
    environment = os.environ.copy()
    environment.update(ASAN_OPTIONS='halt_on_error=1', UBSAN_OPTIONS='halt_on_error=1',
                       TSAN_OPTIONS='halt_on_error=1')

    def run(command, log):
        command = list(map(str, command))
        commands.append(command)
        result = subprocess.run(command, capture_output=True, text=True,
                                timeout=300, env=environment)
        log.write_text(result.stdout + result.stderr)
        if result.returncode or result.stderr:
            raise RuntimeError('failed command; see ' + str(log))
        return result.stdout.strip()

    for label, sanitizer in (('asan-ubsan', 'address,undefined'), ('tsan', 'thread')):
        directory = out / label
        directory.mkdir()
        flags = ['clang++', '-std=c++11', '-O1', '-g', '-pthread', '-Wall', '-Wextra',
                 '-Wpedantic', '-fno-omit-frame-pointer', '-fno-sanitize-recover=all',
                 '-fsanitize=' + sanitizer, '-I' + str(source / 'inst/include'),
                 '-I' + str(out), '-I' + str(out / 'fixtures')]
        objects = []
        for name in UNITS:
            target = directory / (Path(name).stem + '.o')
            run(flags + ['-c', source / 'src' / name, '-o', target],
                directory / (name + '-build.log'))
            objects.append(target)
        reference = directory / 'P2Reference.o'
        run(flags + ['-DP2=P2Reference', '-c', out / 'fixtures/P2Reference.cpp', '-o', reference],
            directory / 'reference-build.log')
        suites = {
            'common': [out / 'check.cpp', source / 'src/P2.cpp', reference],
            'helpers': [out / 'helpers.cpp'],
            'stream': ['-DP2_HYBRID_STREAM_CHECK', out / 'stream-check.cpp']
        }
        patterns = {
            'common': r'small_checks=(?P<small_checks>\d+)\nPASS checks=(?P<checks>\d+)',
            'helpers': r'PASS queue_checks=(?P<queue_checks>\d+) selector_checks=(?P<selector_checks>\d+) worker_checks=(?P<worker_checks>\d+)',
            'stream': r'PASS stream_checks=(?P<stream_checks>\d+)'
        }
        results[label] = {}
        for name, sources in suites.items():
            executable = directory / ('check-' + name)
            run(flags + sources + objects + ['-o', executable], directory / (name + '-build.log'))
            result = run([executable], directory / (name + '-run.log'))
            match = re.fullmatch(patterns[name], result)
            if not match:
                raise RuntimeError('unexpected checker output: ' + name)
            results[label][name] = {key: int(value) for key, value in match.groupdict().items()}
            print(label, name, result, flush=True)
    for name, expected in inputs.items():
        if digest(out / name) != expected:
            raise RuntimeError('test input changed: ' + name)
    (out / 'validation.json').write_text(json.dumps(dict(
        results=results, commands=commands, inputs=inputs, compiler_warnings=0,
        scope='Current P2 independent and differential counts; signed queue coverage and concurrency; '
              'frozen selector boundaries; actual worker launches; independent private stream windows.'),
        indent=2) + '\n')
    print('PASS current-production P2 validation', flush=True)


if __name__ == '__main__':
    main()
