#!/usr/bin/env python3
"""Prepare a temporary, instrumented parallel-Phi/Legendre copy for tests."""

import argparse
import hashlib
from pathlib import Path
import re


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("--source", type=Path,
                        help="optional before-change PiLegendre.cpp snapshot")
    args = parser.parse_args()
    repository = Path(__file__).resolve().parents[2]
    driver_path = args.source or repository / "src/PiLegendre.cpp"
    source = driver_path.read_text()
    if "std::vector<std::thread> threads;" not in source:
        # Include the actual helper and scheduler in the temporary test unit.
        # The production header must not be included a second time, since its
        # thread type is replaced below only in this copied definition.
        root = driver_path.resolve().parents[1]
        batch = (root / "inst/include/NumbersUtils/ThreadBatch.h").read_text()
        parallel = (root / "src/PhiParallel.cpp").read_text()
        include = '#include "NumbersUtils/ThreadBatch.h"\n'
        if parallel.count(include) != 1:
            raise RuntimeError("expected one ThreadBatch include in PhiParallel.cpp")
        parallel = parallel.replace(include, "")
        source = batch.replace("#pragma once\n", "", 1) + "\n" + parallel + "\n" + source
    source_hash = hashlib.sha256(source.encode()).hexdigest()

    storage = "std::vector<std::thread> threads;"
    if source.count(storage) != 1:
        raise RuntimeError("expected exactly one batch thread vector")
    source = source.replace(storage, "std::vector<ThreadFaults::Thread> threads;")

    foreman = re.compile(
        r"    void phiForeman\(std::int64_t &mySum, std::int64_t lowerBound,\s+"
        r"std::int64_t upperBound, std::int64_t x,\s+"
        r"const PhiTables& tables\) \{\n"
    )
    if len(list(foreman.finditer(source))) != 1:
        raise RuntimeError("phiForeman entry no longer matches the test fixture")
    source = foreman.sub(lambda match: match.group(0) +
                        "        ThreadFaults::BeforeForeman();\n", source)

    construction = re.compile(
        r"^([ \t]*)(?:PhiThreadBatch batch\(nThreads, tables\)|ThreadBatch batch\(nThreads\));$",
        re.MULTILINE
    )
    sites = list(construction.finditer(source))
    if len(sites) != 4:
        raise RuntimeError("expected first, middle, final, and small-input batch sites")
    regions = iter(("First", "Middle", "Final", "Small"))

    def tag_batch(match):
        indent = match.group(1)
        return (indent + "ThreadFaults::BeginBatch(ThreadFaults::Region::" + next(regions) + ");\n"
                + match.group(0))

    source = construction.sub(tag_batch, source)
    args.output_directory.mkdir(parents=True, exist_ok=True)
    output = args.output_directory / "prime-count-thread-core.cpp"
    output.write_text("// Test-only copy; production source SHA-256: " + source_hash + "\n" + source)
    print(output)


if __name__ == "__main__":
    main()
