#!/usr/bin/env python3
"""Build and compare isolated 65535/131071 Phi-cache variants.

Use phases build, bench, large, memory sequentially. Finish sanitizer runs and
other builds before bench/large. Timing binaries contain no memory probes.
"""
import argparse
import csv
import hashlib
import io
import json
import platform
import random
import shutil
import statistics
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
VARIANTS = ("baseline", "double")
METHODS = ("legendre", "meissel", "lehmer")
COUNTS = {10**6: 78498, 10**8: 5761455, 10**10: 455052511,
          10**12: 37607912018, 10**13: 346065536839,
          10**14: 3204941750802, 10**15: 29844570422669}
FIELDS = ("method", "n", "threads", "repeat", "count", "total_ms")
SOURCES = ("PiLegendre.cpp", "PiMeissel.cpp", "PiLehmer.cpp",
           "PhiParallel.cpp", "P2.cpp", "P3.cpp", "Eratosthenes.cpp")


def run(command, timeout=300):
    result = subprocess.run(command, capture_output=True, text=True, timeout=timeout)
    if result.returncode:
        raise RuntimeError(f"failed: {command}\n{result.stderr}")
    return result


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def save_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def build(root):
    source = root / "source"
    if not source.exists():
        for directory in ("src", "inst/include"):
            shutil.copytree(REPO / directory, source / directory,
                            ignore=shutil.ignore_patterns("*.o", "*.so", "*.dylib"))
    manifest_path = root / "variants/cache-double-manifest.json"
    if not manifest_path.exists():
        result = run(["python3", str(HERE / "prepare_cache_double.py"),
                      str(root / "variants"), "--source-root", str(source), "--memory-probes"])
        (root / "prepare.log").write_text(result.stdout + result.stderr)
    manifest = json.loads(manifest_path.read_text())
    if (not manifest["memory_probes"] or manifest["candidate_max_x"] != 131071 or
            manifest["source_sha256"] != hashlib.sha256((source / "src/Phi.cpp").read_bytes()).hexdigest()):
        raise RuntimeError("prepared variants do not match this experiment")
    for name, digest in manifest["generated_sha256"].items():
        if hashlib.sha256((root / "variants" / name).read_bytes()).hexdigest() != digest:
            raise RuntimeError("prepared variant changed: " + name)
    shutil.copy2(HERE.parent / "prime-count-methods-benchmark.cpp", root / "benchmark.cpp")
    shutil.copy2(HERE / "cache_double_memory.cpp", root / "memory.cpp")
    commands = []
    for memory in (False, True):
        for variant in VARIANTS:
            directory = root / "variants" / (variant + ("-memory" if memory else ""))
            executable = root / (("memory-" if memory else "benchmark-") + variant)
            command = ["clang++", "-std=c++11", "-O3", "-g", "-pthread",
                       "-Wall", "-Wextra", "-Wpedantic", "-I" + str(source / "inst/include"),
                       "-I" + str(root / "variants"),
                       str(root / ("memory.cpp" if memory else "benchmark.cpp"))]
            command += [str(source / "src" / name) for name in SOURCES]
            command += [str(directory / "Phi.cpp"), "-o", str(executable)]
            result = run(command)
            (root / (executable.name + "-build.log")).write_text(result.stdout + result.stderr)
            commands.append(command)
    save_json(root / "build-commands.json", commands)
    files = [root / "benchmark.cpp", root / "memory.cpp"]
    files += list((root / "variants").rglob("*.cpp"))
    files += list((root / "variants").rglob("*.h"))
    save_json(root / "build-metadata.json", {
        "platform": platform.platform(),
        "compiler": run(["clang++", "--version"]).stdout,
        "files": {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        "production_sources": {name: hashlib.sha256((source / "src" / name).read_bytes()).hexdigest()
                               for name in SOURCES},
    })
    print("Four separate timing/memory binaries built", flush=True)


def summarize(rows):
    output = []
    cases = sorted({(r["method"], int(r["n"]), int(r["threads"])) for r in rows})
    for method, n, threads in cases:
        selected = [r for r in rows if (r["method"], int(r["n"]), int(r["threads"])) ==
                    (method, n, threads)]
        for state in ("first", "warm"):
            groups = {}
            for r in selected:
                if (int(r["repeat"]) == 0) == (state == "first"):
                    groups.setdefault((r["variant"], r["process"]), []).append(float(r["total_ms"]))
            medians = {key: statistics.median(value) for key, value in groups.items()}
            processes = sorted({process for _, process in medians})
            ratios = [medians[("double", p)] / medians[("baseline", p)] for p in processes]
            row = dict(method=method, n=n, threads=threads, state=state, processes=len(processes))
            for variant in VARIANTS:
                values = [medians[(variant, p)] for p in processes]
                row[variant + "_ms"] = statistics.median(values)
                row[variant + "_min_ms"] = min(values)
                row[variant + "_max_ms"] = max(values)
            row["paired_change_percent"] = 100 * (statistics.median(ratios) - 1)
            row["paired_min_change_percent"] = 100 * (min(ratios) - 1)
            row["paired_max_change_percent"] = 100 * (max(ratios) - 1)
            output.append(row)
    return output


def benchmark(root, large):
    seed = 202611081 if not large else 202611082
    rng = random.Random(seed)
    if large:
        cases = [(m, 10**14, t, p, 3) for m in METHODS for t in (1, 8) for p in range(3)]
        # A short screen at the user's input, not a three-process confirmation.
        cases += [(m, 10**15, 8, 0, 2) for m in METHODS]
    else:
        cases = [(m, n, t, p, 4) for m in METHODS for n in list(COUNTS)[:5]
                 for t in (1, 8) for p in range(3)]
    rng.shuffle(cases)
    rows = []
    name = "cache-double-large" if large else "cache-double"
    for index, (method, n, threads, process, repeats) in enumerate(cases):
        variants = list(VARIANTS)
        rng.shuffle(variants)
        for order, variant in enumerate(variants):
            command = [str(root / ("benchmark-" + variant)), method, str(n),
                       str(threads), str(repeats), str(COUNTS[n])]
            result = run(command)
            if result.stderr:
                raise RuntimeError(result.stderr)
            parsed = list(csv.DictReader(io.StringIO(result.stdout), fieldnames=FIELDS))
            if len(parsed) != repeats:
                raise RuntimeError("unexpected benchmark row count")
            for repetition, row in enumerate(parsed):
                if (row["method"], int(row["n"]), int(row["threads"]),
                    int(row["repeat"]), int(row["count"])) != (
                        method, n, threads, repetition, COUNTS[n]):
                    raise RuntimeError("benchmark configuration/result mismatch")
                row.update(variant=variant, process=process, pair_order=order)
                rows.append(row)
        save_csv(root / (name + ".csv"), rows)
        print(f"{name}: {index + 1}/{len(cases)} paired cases; {method} {n} threads={threads}", flush=True)
    save_csv(root / (name + "-summary.csv"), summarize(rows))
    save_json(root / (name + "-metadata.json"), {
        "seed": seed, "verified_rows": len(rows), "paired_cases": len(cases),
        "methodology": "Fresh process per variant/configuration/pair, randomized paired order; "
                       "first call separate, warm median per process; complete native drivers "
                       "include table construction/destruction, exclude process/R startup.",
        "summary": "Median of paired ratios; min/max are variability ranges, not confidence intervals.",
        "large_screen": "1e15: one fresh process per variant, first plus one warm call." if large else None,
    })


def memory(root):
    rows = []
    for method in METHODS:
        for n in (10**8, 10**12, 10**13, 10**14):
            for threads in (1, 8):
                for variant in VARIANTS:
                    result = run([str(root / ("memory-" + variant)), method, str(n),
                                  str(threads), "2", str(COUNTS[n])])
                    parsed = [json.loads(line) for line in result.stdout.splitlines()]
                    if len(parsed) != 2:
                        raise RuntimeError("unexpected memory probe output")
                    for row in parsed:
                        row["variant"] = variant
                        rows.append(row)
            print(f"memory: {method} {n}", flush=True)
    save_json(root / "cache-double-memory.json", rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=("build", "bench", "large", "memory"))
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    root = args.output.resolve()
    if root == REPO or REPO in root.parents:
        raise RuntimeError("use an artifact directory outside the repository")
    root.mkdir(parents=True, exist_ok=True)
    if any(p.is_symlink() for p in root.rglob("*")):
        raise RuntimeError("artifact directory must not contain symlinks")
    if args.phase == "build": build(root)
    elif args.phase == "memory": memory(root)
    else: benchmark(root, args.phase == "large")


if __name__ == "__main__":
    main()
