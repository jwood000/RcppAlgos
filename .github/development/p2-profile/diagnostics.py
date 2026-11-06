#!/usr/bin/env python3
"""Collect diagnostic counters/timers separately from the baseline timings."""
import argparse
import csv
import json
from pathlib import Path
import random
import subprocess


def parse(output, count_mode):
    lines = output.splitlines()
    header = lines[0]
    groups = []
    for line in lines[1:]:
        if line.startswith("{"):
            groups[-1]["metrics"].append(json.loads(line))
        else:
            groups.append(dict(row=next(csv.DictReader([header, line])), metrics=[]))
    for group in groups:
        summary = next(m for m in group["metrics"] if m["kind"] == "summary")
        counts = summary["counters"]
        if summary["mixed_cases"] or summary["calls"] != 1:
            raise RuntimeError("unexpected diagnostic session scope")
        if counts["queries"] != summary["last"] - summary["first"] or \
                counts["queries"] != counts["completed_queries"]:
            raise RuntimeError("incomplete P2 query instrumentation")
        if count_mode:
            invariants = [
                counts["public_phi"] == counts["queries"],
                counts["kernel_entries"] == counts["queries"] + counts["recursive_entries"],
                counts["cache_gets"] == counts["cache_eligible"] + counts["cache_ineligible"],
                counts["cache_eligible"] == counts["cache_hits"] + counts["cache_misses"],
                counts["cache_misses"] == counts["cache_no_page"] + counts["cache_empty_cell"],
                counts["expansions"] == counts["cache_gets"] - counts["cache_hits"],
                counts["top_cache_gets"] == counts["queries"],
                counts["top_cache_eligible"] == 0,
            ]
            # The kernel-entry and top-cache equalities rely on this matrix's
            # large inputs/cutoffs; they are not general Phi contracts.
            if not all(invariants):
                raise RuntimeError("inconsistent diagnostic counters")
    return groups


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_dir", type=Path)
    options = parser.parse_args()
    root = options.artifact_dir.resolve()
    partials = {(r["method"], r["x"]): r["count"]
                for r in json.loads((root / "expected-p2.json").read_text())}
    primes = {10**12: 37607912018, 10**13: 346065536839}
    cases = [(kind, mode, method, x, threads)
             for kind in ("counts", "timing") for mode in ("phases", "p2")
             for method in ("meissel", "lehmer") for x in primes for threads in (1, 8)]
    random.Random(202611072).shuffle(cases)
    raw = root / "raw-diagnostics"
    raw.mkdir(exist_ok=True)
    collected = []
    for index, (kind, mode, method, x, threads) in enumerate(cases):
        expected = primes[x] if mode == "phases" else partials[(method, x)]
        result = subprocess.run([str(root / ("profile-" + kind)), mode, method,
                                 str(x), str(threads), "3", str(expected)],
                                text=True, capture_output=True, check=True)
        label = f"{kind}-{mode}-{method}-{x}-t{threads}"
        (raw / (label + ".txt")).write_text(result.stdout)
        if result.stderr:
            (raw / (label + ".stderr")).write_text(result.stderr)
        groups = parse(result.stdout, kind == "counts")
        if len(groups) != 3 or any(int(g["row"]["p2"]) != partials[(method, x)] for g in groups):
            raise RuntimeError("diagnostic count mismatch")
        collected.append(dict(kind=kind, mode=mode, method=method, x=x,
                              threads=threads, groups=groups))
        print(f"Diagnostic {index + 1}/{len(cases)}: {label}", flush=True)
    (root / "diagnostics.json").write_text(json.dumps(collected, indent=2) + "\n")
    print(f"All {3 * len(cases)} diagnostic results and counter invariants matched.")


if __name__ == "__main__":
    main()
