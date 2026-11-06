#!/usr/bin/env python3
"""Prepare a synchronous, bounds-checked copy of the production scheduler."""

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
    extracted = "std::vector<std::thread> threads;" not in source
    if extracted:
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
    source = source.replace(storage, "std::vector<ScheduleCheck::Thread> threads;")

    foreman = re.compile(
        r"    void phiForeman\(std::int64_t &mySum, std::int64_t lowerBound,\s+"
        r"std::int64_t upperBound, std::int64_t x,\s+"
        r"const PhiTables& tables\) \{\n"
    )
    if len(list(foreman.finditer(source))) != 1:
        raise RuntimeError("phiForeman entry no longer matches the test fixture")
    source = foreman.sub(lambda match: match.group(0) +
                        "        if (ScheduleCheck::BeforeForeman(mySum, lowerBound,\n"
                        "                upperBound, tables.primeCount())) return;\n", source)

    # This location follows the production cap, when present, and precedes
    # every allocation sized by nThreads. It makes an unfixed INT_MAX request
    # fail safely instead of attempting a multi-gigabyte test allocation.
    bound = "a" if extracted else "piSqrtx"
    allocation_guard = "            std::int64_t myRange = (" + bound + " - strt) + 1;"
    if source.count(allocation_guard) != 1:
        raise RuntimeError("parallel schedule entry no longer matches the fixture")
    source = source.replace(allocation_guard,
                            "            ScheduleCheck::BeforeSchedule(nThreads);\n" +
                            allocation_guard)

    construction = re.compile(
        r"^([ \t]*)(?:PhiThreadBatch batch\(nThreads, tables\)|ThreadBatch batch\(nThreads\));$",
        re.MULTILINE
    )
    if len(list(construction.finditer(source))) != 4:
        raise RuntimeError("expected first, middle, final, and small-input batch sites")
    sites = iter((("first", "intermediate"), ("middle", "intermediate"),
                  ("final", "intermediate"), ("small", "mySums")))

    def tag_batch(match):
        region, sums = next(sites)
        return (match.group(1) + 'ScheduleCheck::BeginBatch("' + region + '", ' +
                sums + ");\n" + match.group(0))

    source = construction.sub(tag_batch, source)
    args.output_directory.mkdir(parents=True, exist_ok=True)
    output = args.output_directory / "prime-count-schedule-core.cpp"
    output.write_text("// Test-only copy; production source SHA-256: " + source_hash +
                      "\n#define SCHEDULE_HAS_PHI_PARALLEL " + str(int(extracted)) +
                      "\n" + source)
    print(output)


if __name__ == "__main__":
    main()
