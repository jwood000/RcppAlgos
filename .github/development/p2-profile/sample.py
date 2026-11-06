#!/usr/bin/env python3
"""Sample only the temporary P2 profiler process started by this script (macOS)."""
import argparse
import json
from pathlib import Path
import subprocess
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_dir", type=Path)
    parser.add_argument("method", choices=("meissel", "lehmer"))
    parser.add_argument("threads", type=int)
    parser.add_argument("--x", type=int, default=10**13)
    options = parser.parse_args()
    root = options.artifact_dir.resolve()
    expected = next(r["count"] for r in json.loads((root / "expected-p2.json").read_text())
                    if r["method"] == options.method and r["x"] == options.x)
    label = f"sample-{options.method}-{options.x}-t{options.threads}"
    with (root / (label + ".csv")).open("w") as output, \
            (root / (label + ".stderr")).open("w") as errors:
        child = subprocess.Popen([str(root / "profile"), "sample", options.method,
                                  str(options.x), str(options.threads), "8000", str(expected)],
                                 stdout=output, stderr=errors)
        try:
            time.sleep(1)
            if child.poll() is not None:
                raise RuntimeError("profiler exited before sampling")
            sampled = subprocess.run(["sample", str(child.pid), "5", "1", "-file",
                                      str(root / (label + ".txt"))],
                                     text=True, capture_output=True)
            (root / (label + ".sampler-log")).write_text(sampled.stdout + sampled.stderr)
            print(sampled.stdout + sampled.stderr)
            sampled.check_returncode()
            result = child.wait(timeout=30)
            if result:
                raise RuntimeError(f"profile child failed with exit {result}")
        finally:
            if child.poll() is None:
                child.terminate()
                child.wait(timeout=10)
    print(f"Checked sampling run: {root / (label + '.csv')}")


if __name__ == "__main__":
    main()
