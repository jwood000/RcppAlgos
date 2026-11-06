#!/usr/bin/env python3
"""Export compact, reviewable profiling results from the raw artifact directory."""
import argparse
import csv
import json
from pathlib import Path
import re
import shutil


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    options = parser.parse_args()
    root, output = options.artifact_dir, options.output_dir
    output.mkdir(parents=True, exist_ok=True)
    for name in ("baseline.csv", "baseline-summary.csv", "baseline-metadata.json",
                 "expected-p2.json"):
        shutil.copy2(root / name, output / name)
    rows = []
    for run in json.loads((root / "diagnostics.json").read_text()):
        for group in run["groups"]:
            summary = next(m for m in group["metrics"] if m["kind"] == "summary")
            workers = [m for m in group["metrics"] if m["kind"] == "worker"]
            bins = [m for m in group["metrics"] if m["kind"] == "quantile"]
            row = {key: run[key] for key in ("kind", "mode", "method", "x", "threads")}
            row.update(repeat=group["row"]["repeat"], p2_ms=group["row"]["p2_ms"],
                       workers=summary["worker_slots"], direct_terms=summary["direct_terms"],
                       root_lookup_ns=summary["root_lookup_ns"], phi_ns=summary["phi_ns"])
            row.update(summary["counters"])
            row.update(worker_min_ns=min(w["wall_ns"] for w in workers),
                       worker_max_ns=max(w["wall_ns"] for w in workers),
                       max_chunk_ns=max(w["max_chunk_ns"] for w in workers),
                       first_eighth_phi_ns=sum(b["phi_ns"] for b in bins if b["bin"] == 0),
                       first_quarter_phi_ns=sum(b["phi_ns"] for b in bins if b["bin"] < 2))
            rows.append(row)
    with (output / "diagnostics-summary.csv").open("w", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    samples = []
    for file in sorted(root.glob("sample-*.txt")):
        report = file.read_text()
        threads = [int(n) for n in re.findall(r"^    (\d+) Thread_", report, re.M)]
        if not threads:
            continue
        collapsed = report.split("Sort by top of stack, same collapsed (when >= 5):", 1)[1]
        collapsed = collapsed.split("Binary Images:", 1)[0]
        leaves = {}
        for line in collapsed.splitlines():
            match = re.match(r"\s*(.*?)\s{2,}(\d+)\s*$", line)
            if match:
                leaves[match.group(1)] = int(match.group(2))
        samples.append(dict(file=file.name, thread_samples=sum(threads),
                            main_thread_samples=threads[0], leaf_samples=leaves))
    (output / "sampling-summary.json").write_text(json.dumps(samples, indent=2) + "\n")
    print(f"Exported {len(rows)} diagnostic rows and {len(samples)} stack summaries.")


if __name__ == "__main__":
    main()
