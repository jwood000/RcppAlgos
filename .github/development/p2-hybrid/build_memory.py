#!/usr/bin/env python3
"""Build/run hybrid P2 memory probes separately from timed work."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import shutil
import subprocess

HERE = Path(__file__).resolve().parent
MULTIPLIERS = (2, 4, 8)
SUPPORTED_MULTIPLIERS = (2, 4, 8, 16, 32, 64, 128)
COUNTS = {10**10: 455052511, 10**12: 37607912018,
          10**13: 346065536839, 10**14: 3204941750802}


def invoke(command):
    result = subprocess.run(command, text=True, capture_output=True, timeout=300)
    if result.returncode:
        raise RuntimeError(f"failed command: {command}\n{result.stdout}\n{result.stderr}")
    return result


def build(root, multipliers):
    source = root / "source"
    if not (source / "src/Phi.cpp").is_file():
        raise RuntimeError("artifact must contain the immutable production snapshot")
    directory = root / "memory-source"
    directory.mkdir(exist_ok=True)
    for name in ("memory.cpp", "P2Hybrid.cpp", "hybrid_metrics.h"):
        shutil.copy2(HERE / name, directory / name)
    commands = []
    for multiplier in multipliers:
        executable = root / f"memory-{multiplier}"
        command = ["clang++", "-std=c++11", "-O3", "-pthread", "-Wall", "-Wextra",
                   "-Wpedantic", "-DP2_HYBRID_MEMORY",
                   f"-DP2_HYBRID_MULTIPLIER={multiplier}",
                   "-I" + str(source / "inst/include"), "-I" + str(directory),
                   str(directory / "memory.cpp"), str(directory / "P2Hybrid.cpp")]
        command += [str(source / "src" / name) for name in
                    ("P3.cpp", "PhiParallel.cpp", "Phi.cpp", "Eratosthenes.cpp")]
        command += ["-o", str(executable)]
        result = invoke(command)
        (root / f"memory-{multiplier}-build.log").write_text(result.stdout + result.stderr)
        commands.append(command)
    manifest = {
        "multipliers": list(multipliers),
        "commands": commands,
        "files": {name: hashlib.sha256((directory / name).read_bytes()).hexdigest()
                  for name in ("memory.cpp", "P2Hybrid.cpp", "hybrid_metrics.h")},
        "production_files": {name: hashlib.sha256((source / "src" / name).read_bytes()).hexdigest()
                             for name in ("P3.cpp", "PhiParallel.cpp", "Phi.cpp", "Eratosthenes.cpp")},
    }
    (root / "memory-build.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Built {len(multipliers)} separate memory probes", flush=True)


def measure(root, multipliers):
    rows = []
    for multiplier in multipliers:
        for method in ("meissel", "lehmer"):
            for x, expected in COUNTS.items():
                for threads in (1, 8):
                    result = invoke([str(root / f"memory-{multiplier}"), method, str(x),
                                     str(threads), "2", str(expected)])
                    if result.stderr:
                        raise RuntimeError(result.stderr)
                    parsed = list(csv.DictReader(io.StringIO(result.stdout)))
                    if len(parsed) != 2:
                        raise RuntimeError("unexpected probe row count")
                    for repeat, row in enumerate(parsed):
                        if (int(row["multiplier"]), row["method"], int(row["x"]),
                            int(row["threads"]), int(row["repeat"]), int(row["count"])) != (
                                multiplier, method, x, threads, repeat, expected):
                            raise RuntimeError("probe configuration/count mismatch")
                        rows.append(row)
                    print(f"Memory: multiplier={multiplier} {method} x={x} threads={threads}",
                          flush=True)
    with (root / "hybrid-memory.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (root / "hybrid-memory-metadata.json").write_text(json.dumps({
        "multipliers": list(multipliers),
        "verified_rows": len(rows), "repeats_per_process": 2,
        "timing": "none; separately compiled capacity/work instrumentation",
        "conservative_scratch_bytes":
            "max(phi subtotal capacity, stream subtotal capacity + sum of all active stream "
            "sieve/next capacities); phases do not overlap, worker allocations may not all overlap",
        "exclusions": ["diagnostic snapshot/vector", "ThreadBatch vectors", "thread stacks",
                       "allocator overhead", "Phi cache", "table construction scratch", "RSS"],
        "phase_order": "tables, PhiParallel, P2, optional P3, table destruction",
    }, indent=2) + "\n")
    print(f"All {len(rows)} memory/work rows verified", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("build", "run"))
    parser.add_argument("artifact_dir", type=Path)
    parser.add_argument("--multipliers", type=int, nargs="+", choices=SUPPORTED_MULTIPLIERS,
                        default=MULTIPLIERS,
                        help="variants to build or measure (default: 2 4 8)")
    options = parser.parse_args()
    if len(set(options.multipliers)) != len(options.multipliers):
        parser.error("multipliers must be distinct")
    root = options.artifact_dir.resolve()
    if not root.is_dir():
        parser.error("artifact directory must already exist")
    (build if options.mode == "build" else measure)(root, options.multipliers)


if __name__ == "__main__":
    main()
