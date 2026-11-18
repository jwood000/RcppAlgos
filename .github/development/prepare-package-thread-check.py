#!/usr/bin/env python3
"""Extract current launchers for native control-flow checks, without changing source.

Only mathematical workers and preparation helpers are replaced by the fixture.
The launchers, storage declarations, captures, and real ThreadBatch are retained.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re

FUNCTIONS = (
    ("MotleyPrimes.cpp", "MotleyMain"),
    ("DivNumSieve.cpp", "DivisorMain"),
    ("PollardRhoContainer.cpp", "PollardRhoMain"),
    ("ThreadSafeComb.cpp", "ThreadSafeCombinations"),
    ("ThreadSafePerm.cpp", "ThreadSafePermutations"),
    ("ComboResMain.cpp", "ComboResMain"),
    ("PermuteResMain.cpp", "PermuteResMain"),
    ("ConstraintsSpecial.cpp", "ConstraintsSpecial"),
    ("ThreadSafeParts.cpp", "StandardPartitions"),
    ("ThreadSafeParts.cpp", "GeneralPartitions"),
    ("SampCombPermStd.cpp", "ThreadSafeSample"),
    ("SamplePartitions.cpp", "ThreadSafeSample"),
    ("RankCombPermMain.cpp", "RankCPParallel"),
    ("RankPartitionMain.cpp", "RankPartsParallel"),
    ("GetProduct.cpp", "PureOutputMain"),
    ("GetComboGroups.cpp", "GroupsMain"),
)


def extract(source, name):
    starts = list(re.finditer(r"^(?:template\s*<[^>]+>\s*)?void " + name + r"\s*\(",
                             source, re.MULTILINE))
    if len(starts) != 1:
        raise RuntimeError("expected one definition for " + name)
    start = starts[0].start()
    # Mask strings/comments so their braces cannot change the extracted span.
    masked = re.sub(r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"',
                    lambda match: " " * len(match.group()), source, flags=re.S)
    opening = masked.index("{", starts[0].end())
    depth = 1
    end = opening + 1
    while depth:
        depth += (masked[end] == "{") - (masked[end] == "}")
        end += 1
    return source[start:end]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    args.output.mkdir(parents=True, exist_ok=True)
    batch_path = root / "inst/include/ThreadBatch.h"
    batch = batch_path.read_text()
    storage = "std::vector<std::thread> threads;"
    if batch.count(storage) != 1 or batch.count("work();") != 1:
        raise RuntimeError("ThreadBatch structure changed; review the fixture")
    batch = batch.replace("#pragma once\n", "", 1).replace(
        storage, "std::vector<PackageThreadFaults::Thread> threads;")
    manifest = {str(batch_path.relative_to(root)): hashlib.sha256(batch_path.read_bytes()).hexdigest()}
    pieces = ["// Generated current-source launcher fixture. See manifest.json.\n", batch]
    for index, (filename, name) in enumerate(FUNCTIONS):
        path = root / "src" / filename
        source = path.read_text()
        function = extract(source, name)
        if function.count("RcppAlgos::ThreadBatch batch(nThreads);") != 1:
            raise RuntimeError("expected one migrated batch in " + name)
        if function.count("batch.finish();") != 1 or "std::thread" in function:
            raise RuntimeError("unexpected launch/finish structure in " + name)
        manifest[str(path.relative_to(root))] = hashlib.sha256(path.read_bytes()).hexdigest()
        pieces.append("namespace Launcher%d {\n%s\n}\n" % (index, function))
    (args.output / "package-thread-core.inc").write_text("\n".join(pieces))
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(args.output / "package-thread-core.inc")


if __name__ == "__main__":
    main()
