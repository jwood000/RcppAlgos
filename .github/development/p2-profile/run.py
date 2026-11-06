#!/usr/bin/env python3
"""Run the uninstrumented phase/P2 matrix after all builds have finished."""
import argparse
import csv
import io
import json
from pathlib import Path
import random
import statistics
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_dir", type=Path)
    parser.add_argument("--processes", type=int, default=3)
    parser.add_argument("--repeats", type=int, default=4)
    parser.add_argument("--seed", type=int, default=202611071)
    options = parser.parse_args()
    if options.processes < 1 or options.repeats < 2:
        parser.error("positive process count and at least two repeats required")
    root = options.artifact_dir.resolve()
    expected_pi = {10**10: 455052511, 10**12: 37607912018, 10**13: 346065536839}
    methods = ("meissel", "lehmer")
    threads = (1, 4, 8)
    randomizer = random.Random(options.seed)
    rows, expected_p2 = [], {}

    # First validate each complete identity against a known pi(x). Partial
    # counts must then agree across every repetition/configuration and the
    # isolated P2 runs. This is not a new independent oracle for P2 itself.
    for mode in ("phases", "p2"):
        cases = [(x, t, process) for x in expected_pi for t in threads
                 for process in range(options.processes)]
        randomizer.shuffle(cases)
        for index, (x, budget, process) in enumerate(cases):
            order = list(methods)
            randomizer.shuffle(order)
            for method in order:
                expected = expected_pi[x] if mode == "phases" else expected_p2[(method, x)]
                command = [str(root / "profile"), mode, method, str(x), str(budget),
                           str(options.repeats), str(expected)]
                result = subprocess.run(command, text=True, capture_output=True, check=True)
                parsed = list(csv.DictReader(io.StringIO(result.stdout)))
                if len(parsed) != options.repeats:
                    raise RuntimeError("unexpected profiler output")
                for row in parsed:
                    if (row["mode"], row["method"], int(row["x"]), int(row["threads"])) != (
                            mode, method, x, budget):
                        raise RuntimeError("unexpected configuration in profiler output")
                    if mode == "phases" and int(row["prime_count"]) != expected_pi[x]:
                        raise RuntimeError("incorrect prime count")
                    count = int(row["p2"])
                    previous = expected_p2.setdefault((method, x), count)
                    if count != previous:
                        raise RuntimeError("P2 differs between configurations")
                    row["process"] = process
                    rows.append(row)
            if (index + 1) % 3 == 0:
                print(f"{mode}: completed {index + 1}/{len(cases)} paired configurations",
                      flush=True)

    with (root / "baseline.csv").open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    summaries = []
    fields = ("setup_ms", "phi_ms", "p2_ms", "p3_ms", "teardown_ms", "total_ms")
    for mode in ("phases", "p2"):
        for x in expected_pi:
            for budget in threads:
                for method in methods:
                    selected = [r for r in rows if (r["mode"], r["method"], int(r["x"]),
                                int(r["threads"])) == (mode, method, x, budget)]
                    for state in ("first", "warm"):
                        summary = dict(mode=mode, method=method, x=x, threads=budget,
                                       state=state)
                        for field in fields:
                            per_process = []
                            for process in range(options.processes):
                                values = [float(r[field]) for r in selected
                                          if r["process"] == process and
                                          ((int(r["repeat"]) == 0) == (state == "first"))]
                                per_process.append(statistics.median(values))
                            summary[field] = statistics.median(per_process)
                        # Ratio of medians, not a claim about the median of ratios.
                        summary["p2_share_percent"] = (100 * summary["p2_ms"] /
                                                       summary["total_ms"])
                        summaries.append(summary)
    with (root / "baseline-summary.csv").open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)
    (root / "expected-p2.json").write_text(json.dumps([
        dict(method=m, x=x, count=count) for (m, x), count in sorted(expected_p2.items())
    ], indent=2) + "\n")
    (root / "baseline-metadata.json").write_text(json.dumps(dict(
        seed=options.seed, processes=options.processes, repeats=options.repeats,
        verified_rows=len(rows), timing="uninstrumented production kernels",
        compiler_command_file=("build-command.json" if (root / "build-command.json").exists()
                               else None), thread_budgets=list(threads),
        input_counts=expected_pi), indent=2) + "\n")
    print(f"All {len(rows)} measured counts matched; raw and summary CSV written to {root}")


if __name__ == "__main__":
    main()
