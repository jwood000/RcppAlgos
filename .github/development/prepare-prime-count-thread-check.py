#!/usr/bin/env python3
"""Prepare a temporary, narrowly instrumented Legendre copy for native tests."""

import argparse
import hashlib
from pathlib import Path
import re


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_directory", type=Path)
    args = parser.parse_args()
    repository = Path(__file__).resolve().parents[2]
    source = (repository / "src/PiLegendre.cpp").read_text()
    source_hash = hashlib.sha256(source.encode()).hexdigest()

    storage = "std::vector<std::thread> threads;"
    if source.count(storage) != 1:
        raise RuntimeError("expected exactly one PhiThreadBatch thread vector")
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

    construction = re.compile(r"^([ \t]*)PhiThreadBatch batch\(nThreads, tables\);$", re.MULTILINE)
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
