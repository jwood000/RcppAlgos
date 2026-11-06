#!/usr/bin/env python3
"""Build diagnostic copies and measure exact serial, warm P2 state reuse.

Run this separately from performance benchmarks. The trace/sort/replay costs
are diagnostic overhead, not an implementation of a faster Phi cache.
"""

import argparse
import hashlib
import json
import platform
import shutil
import subprocess
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
CASES = [
    ("meissel", 10**12, 23729370364, 77269, 8022045, 879801, 27515549),
    ("lehmer", 10**12, 39098548800, 78330, 10790592, 1186596, 33647030),
    ("meissel", 10**13, 220518863542, 225230, 68318691, 7501776, 187553851),
    ("lehmer", 10**13, 361639241781, 227372, 93769798, 10689403, 226651637),
]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def checked(command, output, timeout=180):
    result = subprocess.run(command, text=True, capture_output=True, timeout=timeout)
    output.write_text(result.stdout)
    output.with_suffix(output.suffix + ".stderr").write_text(result.stderr)
    if result.returncode:
        raise RuntimeError(f"command failed ({result.returncode}): {command}\n{result.stderr}")
    return result.stdout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--skip-build", action="store_true")
    args = parser.parse_args()
    output = args.output.resolve()
    if output == REPO or REPO in output.parents:
        raise RuntimeError("use a temporary artifact directory outside the repository")
    output.mkdir(parents=True, exist_ok=True)
    # Repeated runs reuse these artifacts. Reject pre-existing aliases before
    # a compiler, copy, or log write can follow one outside the artifact tree.
    if any(path.is_symlink() for path in output.rglob("*")):
        raise RuntimeError("artifact directory must not contain symlinks")
    source = output / "source"
    if not source.exists():
        for directory in ("src", "inst/include"):
            shutil.copytree(REPO / directory, source / directory,
                            ignore=shutil.ignore_patterns("*.o", "*.so", "*.dylib"))
    instrumented = output / "instrumented"
    executable = output / "reuse"
    if not args.skip_build:
        for name in ("run_reuse.py", "instrument_reuse.py"):
            shutil.copy2(HERE / name, output / name)
        checked(["python3", str(HERE / "instrument_reuse.py"), str(instrumented),
                 "--source-root", str(source)], output / "instrument.log")
        shutil.copy2(HERE / "reuse.cpp", output / "reuse.cpp")
        command = ["clang++", "-std=c++11", "-O3", "-g", "-pthread",
                   "-Wall", "-Wextra", "-Wpedantic", "-I" + str(source / "inst/include"),
                   "-I" + str(instrumented), str(output / "reuse.cpp"),
                   str(instrumented / "P2.cpp"), str(instrumented / "Phi.cpp"),
                   str(source / "src/Eratosthenes.cpp"), "-o", str(executable)]
        (output / "build-command.json").write_text(json.dumps(command, indent=2) + "\n")
        checked(command, output / "build.log")
    checked([str(executable), "--self-test"], output / "self-test.json")
    results = []
    for method, x, expected, queries, kernels, records, suffix in CASES:
        command = [str(executable), method, str(x), str(expected), str(records)]
        print(f"Tracing and replaying {method} x={x}", flush=True)
        result = json.loads(checked(command, output / f"{method}-{x}.json"))
        for key, value in (("p2", expected), ("queries", queries),
                           ("total_kernels", kernels), ("records", records),
                           ("total_suffix", suffix)):
            if result[key] != value:
                raise RuntimeError(f"{method} {x}: {key}={result[key]}, expected {value}")
        results.append(result)
    # Fresh-process repetition checks deterministic trace-derived statistics.
    method, x, expected, _, _, records, _ = CASES[-1]
    repeated = json.loads(checked([str(executable), method, str(x), str(expected), str(records)],
                                 output / "repeat-lehmer-10000000000000.json"))
    if repeated != results[-1]:
        raise RuntimeError("largest-case repeat produced different statistics")
    (output / "reuse-results.json").write_text(json.dumps(results, indent=2) + "\n")
    metadata = {
        "platform": platform.platform(),
        "compiler": subprocess.check_output(["clang++", "--version"], text=True),
        "methodology": "Serial P2, one untraced warmup, one traced call per fresh process; "
                       "largest case repeated in a second fresh process. Counts only; no timing claims.",
        "counter_reference": "Initial P2 profiling study: warm isolated serial count probe.",
        "large_p2_reference": "Uninstrumented complete counting identities in initial study.",
        "files": {str(p.relative_to(output)): digest(p) for p in
                  [output / "reuse.cpp", instrumented / "reuse_metrics.h",
                   instrumented / "P2.cpp", instrumented / "Phi.cpp",
                   output / "run_reuse.py", output / "instrument_reuse.py"]},
        "snapshot_hashes": {str(p.relative_to(source)): digest(p) for p in source.rglob("*")
                            if p.is_file()},
        "deterministic_repeat": True,
    }
    (output / "reuse-metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"All four cases and the repeated largest case passed: {output}", flush=True)


if __name__ == "__main__":
    main()
