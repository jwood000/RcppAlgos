#!/usr/bin/env python3
"""Confirm selected sieve sizes without per-kernel instrumentation.

Run after run.R finishes, from any directory:
  python3 confirm-native.py RESULTS_DIRECTORY

The only sieve-source substitution disables instrumented workers while retaining
the segment override, policy, and coarse Main timers. The temporary build is kept
for inspection; production sources are never edited.
"""

import argparse
import csv
import hashlib
import json
import math
import pathlib
import platform
import shlex
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=pathlib.Path)
    parser.add_argument("--max-threads", type=int, default=10)
    parser.add_argument("--max-cores", type=int, default=10)
    parser.add_argument("--repeats", type=int, default=15)
    parser.add_argument("--seed", type=int, default=20261102)
    parser.add_argument("--cxx", default="clang++")
    args = parser.parse_args()
    if min(args.max_threads, args.max_cores, args.repeats) < 1:
        parser.error("thread/core counts and repeats must be positive")
    if not 0 <= args.seed <= 2**32 - 1:
        parser.error("seed must fit an unsigned 32-bit integer")

    here = pathlib.Path(__file__).resolve().parent
    root = here.parents[2]
    results = args.results.resolve()
    output = results / "native-confirmation.csv"
    metadata_path = results / "native-confirmation-metadata.json"
    if output.exists() or metadata_path.exists():
        parser.error("native confirmation output already exists; use a fresh results directory")
    selected_path = results / "selected-configurations.txt"
    runs_path = results / "runs.csv"
    selected = selected_path.read_text().splitlines()
    if len(selected) != 3 or selected[0] != "default" or len(set(selected)) != 3:
        parser.error("expected default and two distinct finalists in selected-configurations.txt")
    if any(not name or "," in name or "\n" in name for name in selected):
        parser.error("configuration names must be nonempty and contain no CSV separators")

    budgets = {name: set() for name in selected[1:]}
    with runs_path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            if row["config"] in budgets:
                value = float(row["budget_bytes"])
                if not math.isfinite(value) or not value.is_integer() or value < 3754:
                    parser.error("finalist byte budgets must be finite whole numbers >= 3754")
                budgets[row["config"]].add(int(value))
    if any(len(values) != 1 for values in budgets.values()):
        parser.error("each finalist must have one unambiguous byte budget in runs.csv")
    configurations = [("default", 0)] + [
        (name, next(iter(budgets[name]))) for name in selected[1:]
    ]
    if (args.max_cores * 6 * 1024**2) ** 2 >= 4_000_000_001_234_567:
        parser.error("the fixed serial workload must naturally use Big under --max-cores")

    original = root / "src/Eratosthenes.cpp"
    source = original.read_text()
    old = "const auto worker = inspector ? PrimeWorker<T, true> : PrimeWorker<T, false>;"
    new = "const auto worker = PrimeWorker<T, false>;"
    if source.count(old) != 1:
        parser.error("worker dispatch changed: expected exactly one documented substitution")
    build = pathlib.Path(tempfile.mkdtemp(prefix="rcppalgos-native-confirmation-"))
    clone = build / "sieve-clone.cpp"
    clone.write_text(source.replace(old, new))
    harness = build / "confirm-native.cpp"
    shutil.copyfile(here / harness.name, harness)
    binary = build / "confirm-native"
    command = [args.cxx, "-std=c++11", "-O3", "-pthread",
               "-I" + str(root / "inst/include"), str(harness),
               str(root / "src/CacheInspector.cpp"), "-o", str(binary)]
    run_command = [str(binary), str(output), str(args.repeats), str(args.seed),
                   str(args.max_threads), str(args.max_cores)]
    for name, budget in configurations:
        run_command.extend([name, str(budget)])

    source_paths = [original, root / "src/CacheInspector.cpp", harness,
                    pathlib.Path(__file__).resolve(), clone,
                    selected_path, runs_path]
    source_paths.extend(sorted((root / "inst/include/NumbersUtils").glob("*.h")))
    metadata = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "host": platform.platform(), "temporary_build_retained": str(build),
        "source_substitution": {"before": old, "after": new},
        "configurations": configurations, "requested_threads": [1, 4],
        "max_threads": args.max_threads, "max_cores": args.max_cores,
        "repeats": args.repeats, "seed": args.seed,
        "compiler": subprocess.check_output([args.cxx, "--version"], text=True),
        "compile_command": shlex.join(command),
        "run_command": shlex.join(run_command),
        "hash_algorithm": "sha256",
        "source_hashes": {str(path): sha256(path) for path in source_paths},
        "timing": "PrimeSieveMain only; excludes output concatenation and comparison, R, "
                  "and report serialization. Kernel instrumentation is disabled; coarse "
                  "Main clocks and inspector setup remain enabled for every configuration.",
        "verification": "Every warmup and timed output is compared exactly, outside the "
                        "timer, with the current null-inspector default-policy output.",
        "ordering": "Each configuration runs once per round, shuffled with std::mt19937. "
                    "The CSV records actual order; shuffle order may vary across C++ libraries.",
        "status": "compiling"
    }
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print("Retained temporary build:", build, flush=True)
    print("Compile:", shlex.join(command), flush=True)
    try:
        subprocess.run(command, check=True)
        metadata["binary_sha256"] = sha256(binary)
        metadata["status"] = "running"
        metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
        print("Run:", shlex.join(run_command), flush=True)
        subprocess.run(run_command, check=True)
    except (OSError, subprocess.CalledProcessError) as error:
        metadata["status"] = "failed"
        metadata["error"] = str(error)
        raise
    else:
        metadata["status"] = "completed"
        metadata["output_sha256"] = sha256(output)
    finally:
        metadata["finished_utc"] = datetime.now(timezone.utc).isoformat()
        metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
