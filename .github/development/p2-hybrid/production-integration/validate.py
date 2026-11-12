#!/usr/bin/env python3
"""Validate the integrated production P2 against its reviewed temporary source.

ROOT is an isolated artifact directory with baseline/ extracted from commit
8b43ed790435b22eb9ec02627379082fe7deb103. No package source is modified here.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess

HERE = Path(__file__).resolve().parent
EXPERIMENT = HERE.parent
REPO = HERE.parents[3]
UNITS = ('Phi.cpp', 'PhiParallel.cpp', 'P3.cpp', 'PiMeissel.cpp',
         'PiLehmer.cpp', 'Eratosthenes.cpp')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n')


def invoke(command, log, timeout=300):
    result = subprocess.run(list(map(str, command)), capture_output=True, text=True, timeout=timeout)
    log.write_text(result.stdout + result.stderr)
    if result.returncode or result.stderr:
        raise RuntimeError(f'{command}\nexit={result.returncode}\n{result.stdout}\n{result.stderr}')
    return result.stdout


def body(source, name):
    pattern = (r'class\s+' + name + r'\s*\{') if name == 'PrimeCountStream' else r'\b' + name + r'\s*\('
    match = re.search(pattern, source)
    if not match:
        raise RuntimeError('missing body: ' + name)
    start = match.start()
    opening = source.index('{', match.end() - 1)
    depth = 0
    for index in range(opening, len(source)):
        if source[index] == '{':
            depth += 1
        elif source[index] == '}':
            depth -= 1
            if depth == 0:
                return source[start:index + 1]
    raise RuntimeError('unclosed body: ' + name)


def normalized(source):
    source = source.replace('InitialChunks::Queue<1>', 'PhiWorkQueue')
    source = re.sub(r'constexpr\s+std::int64_t\s+MIN_WORK_PER_THREAD\s*=\s*1048576\s*;', '', source)
    source = source.replace('MIN_STREAM_WORK_PER_THREAD', '1048576').replace('MIN_WORK_PER_THREAD', '1048576')
    source = re.sub(r'\s+', '', source)
    # The maintainer added braces around this unchanged single return.
    return source.replace('if(x<=tables.limit()){returntables.pi()[x];}',
                          'if(x<=tables.limit())returntables.pi()[x];')


def audit_source(snapshot, out):
    generated = EXPERIMENT / 'results-consolidated/generated'
    flags = ['clang++', '-std=c++11', '-E', '-P', '-I' + str(snapshot / 'inst/include'), '-I' + str(generated)]
    production = invoke(flags + [snapshot / 'src/P2.cpp'], out / 'production-preprocessed.log')
    candidate = invoke(flags + ['-DP2_HYBRID_SELECTOR_CANDIDATE', generated / 'P2HybridInitial.cpp'],
                       out / 'candidate-preprocessed.log')
    if len(re.findall(r'constexpr\s+std::int64_t\s+MIN_STREAM_WORK_PER_THREAD\s*=\s*1048576\s*;', production)) != 1:
        raise RuntimeError('shared stream worker grain differs from the reviewed literal')
    compared = {}
    for name in ('countP2Range', 'piAt', 'countStreamRange'):
        left, right = normalized(body(production, name)), normalized(body(candidate, name))
        if left != right:
            (out / (name + '-production.txt')).write_text(left)
            (out / (name + '-candidate.txt')).write_text(right)
            raise RuntimeError('source audit differs: ' + name)
        compared[name] = hashlib.sha256(left.encode()).hexdigest()
    # These bodies intentionally differ after the signed-index cleanup and
    # job-budget fix. Record both bodies instead of claiming token identity.
    # Their behavior is checked by the independent/differential common suite,
    # private stream-window oracle, and launch-observing helper suite below.
    revised = {}
    for name, reason in {
        'countPhiPrefix': 'Cap actual launches by fixed-queue job count, tested by observing real ThreadBatch launches.',
        'PrimeCountStream': 'Signed offset/length and renamed crossing vector; checked against independent stream-window prime oracle.',
        'countStreamingTail': 'Shared stream worker budget and bounded int conversion; common differential/parallel coverage.'
    }.items():
        left, right = normalized(body(production, name)), normalized(body(candidate, name))
        (out / (name + '-production.txt')).write_text(left)
        (out / (name + '-candidate.txt')).write_text(right)
        revised[name] = dict(production_sha256=hashlib.sha256(left.encode()).hexdigest(),
                             candidate_sha256=hashlib.sha256(right.encode()).hexdigest(), reason=reason)
    for name in ('chooseStreamFirst', 'countStreamingTail'):
        if len(re.findall(r'\bstreamWorkerBudget\s*\(', body(production, name))) != 1:
            raise RuntimeError(name + ' does not use the shared worker-budget helper exactly once')
    result = dict(compared_bodies=compared, intentionally_revised_bodies=revised,
                  production_sha256=sha(snapshot / 'src/P2.cpp'),
                  candidate_sha256=sha(generated / 'P2HybridInitial.cpp'),
                  shared_worker_grain=1048576,
                  normalization='Preprocess instrumentation away; discard whitespace; rename Queue<1> to PhiWorkQueue; '
                                'replace identical local/shared 1048576 constants with their literal value; '
                                'normalize only the added braces around the piAt direct-table return.',
                  selector='Compared at runtime against the frozen candidate, including cap/admission boundaries.',
                  queue='Signed tickets compared against frozen K=1 mapping, negative/extreme tickets, independent coverage.',
                  worker_budget='Shared helper used once by chooser and launcher; runtime grain/clamp/extreme checks; '
                                'actual Phi launch counts checked against independently generated jobs.')
    save(out / 'source-audit.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    out = root / 'native-validation'
    if out.exists():
        raise RuntimeError('refusing existing native-validation directory')
    snapshot = out / 'source'
    (snapshot / 'src').mkdir(parents=True)
    shutil.copytree(REPO / 'inst/include', snapshot / 'inst/include')
    for name in UNITS + ('P2.cpp',):
        shutil.copy2(REPO / 'src' / name, snapshot / 'src' / name)
    shutil.copy2(root / 'baseline/src/P2.cpp', out / 'BaselineP2.cpp')
    shutil.copy2(snapshot / 'src/P2.cpp', out / 'ProductionP2.cpp')
    # Reuse the retained private stream-window oracle against this production
    # source, without changing the historical checker or adding package hooks.
    stream_check = out / 'stream-check.cpp'
    stream_source = (EXPERIMENT / 'check.cpp').read_text()
    if stream_source.count('#include "P2Hybrid.cpp"') != 1:
        raise RuntimeError('private stream-check include changed')
    stream_check.write_text(stream_source.replace('#include "P2Hybrid.cpp"', '#include "ProductionP2.cpp"'))
    source_audit = audit_source(snapshot, out)
    input_paths = [p for p in snapshot.rglob('*') if p.is_file()]
    input_paths += [out / 'BaselineP2.cpp', out / 'ProductionP2.cpp', stream_check,
                    HERE / 'helpers.cpp', HERE / 'validate.py',
                    EXPERIMENT / 'check.cpp', EXPERIMENT / 'consolidated/P2SelectorCandidate.h',
                    EXPERIMENT / 'initial-chunks/InitialChunks.h']
    before = {str(p): sha(p) for p in input_paths}
    records, commands = [], []
    for label, sanitizer in (('asan-ubsan', 'address,undefined'), ('tsan', 'thread')):
        directory = out / label
        directory.mkdir()
        flags = ['clang++', '-std=c++11', '-O1', '-g', '-pthread', '-Wall', '-Wextra', '-Wpedantic',
                 '-fno-omit-frame-pointer', '-fno-sanitize-recover=all', '-fsanitize=' + sanitizer,
                 '-I' + str(snapshot / 'inst/include'), '-I' + str(out),
                 '-I' + str(EXPERIMENT / 'consolidated'), '-I' + str(EXPERIMENT / 'initial-chunks')]
        common = []
        for name in UNITS:
            target = directory / (Path(name).stem + '.o')
            command = flags + ['-c', snapshot / 'src' / name, '-o', target]
            invoke(command, directory / (name + '-build.log'))
            commands.append(list(map(str, command)))
            common.append(target)
        reference = directory / 'P2Reference.o'
        command = flags + ['-DP2=P2Reference', '-c', out / 'BaselineP2.cpp', '-o', reference]
        invoke(command, directory / 'reference-build.log')
        commands.append(list(map(str, command)))
        executable = directory / 'check-production'
        command = flags + [EXPERIMENT / 'check.cpp', snapshot / 'src/P2.cpp', reference] + common + ['-o', executable]
        invoke(command, directory / 'common-build.log')
        commands.append(list(map(str, command)))
        result = invoke([executable], directory / 'common-run.log')
        commands.append([str(executable)])
        match = re.fullmatch(r'small_checks=(\d+)\nPASS checks=(\d+)\n', result)
        if not match:
            raise RuntimeError('unexpected common checker output: ' + result)
        record = dict(sanitizer=label, small_checks=int(match[1]), common_checks=int(match[2]))
        print(label, result.strip(), flush=True)
        executable = directory / 'check-helpers'
        command = flags + [HERE / 'helpers.cpp'] + common + ['-o', executable]
        invoke(command, directory / 'helpers-build.log')
        commands.append(list(map(str, command)))
        result = invoke([executable], directory / 'helpers-run.log')
        commands.append([str(executable)])
        match = re.fullmatch(r'PASS queue_checks=(\d+) selector_checks=(\d+) worker_checks=(\d+)\n', result)
        if not match:
            raise RuntimeError('unexpected helper output: ' + result)
        record.update(queue_checks=int(match[1]), selector_checks=int(match[2]), worker_checks=int(match[3]))
        print(label, result.strip(), flush=True)
        executable = directory / 'check-stream'
        command = flags + ['-DP2_HYBRID_STREAM_CHECK', stream_check] + common + ['-o', executable]
        invoke(command, directory / 'stream-build.log')
        commands.append(list(map(str, command)))
        result = invoke([executable], directory / 'stream-run.log')
        commands.append([str(executable)])
        match = re.fullmatch(r'PASS stream_checks=(\d+)\n', result)
        if not match:
            raise RuntimeError('unexpected stream checker output: ' + result)
        record.update(stream_checks=int(match[1]))
        records.append(record)
        print(label, result.strip(), flush=True)
    for path, digest in before.items():
        if sha(path) != digest:
            raise RuntimeError('validation input changed: ' + path)
    if sha(REPO / 'src/P2.cpp') != source_audit['production_sha256']:
        raise RuntimeError('live production P2 changed during validation')
    save(out / 'validation.json', dict(checks=records, commands=commands, input_sha256=before,
         source_audit=source_audit, warnings=0,
         scope='One unchanged common independent/differential checker per sanitizer, plus production helper equivalence '
               'and signed queue concurrency/overflow/no-barrier checks, actual Phi launch counts, shared stream budgets, '
               'and independent private stream-window checks. Three unchanged helper bodies retain token equivalence; '
               'three intentionally revised bodies are recorded and behaviorally tested. '
               'No benchmark timings collected.'))
    print('PASS production native validation', flush=True)


if __name__ == '__main__':
    main()
