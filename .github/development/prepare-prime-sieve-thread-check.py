#!/usr/bin/env python3
"""Prepare an instrumented copy of the real sieve launcher and ThreadBatch."""

import argparse
import hashlib
from pathlib import Path
import re


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_directory", type=Path)
    args = parser.parse_args()
    repository = Path(__file__).resolve().parents[2]
    source = (repository / "src/Eratosthenes.cpp").read_text()
    batch = (repository / "inst/include/NumbersUtils/ThreadBatch.h").read_text()
    hashes = (hashlib.sha256(source.encode()).hexdigest(),
              hashlib.sha256(batch.encode()).hexdigest())
    include = '#include "NumbersUtils/ThreadBatch.h"\n'
    if source.count(include) != 1:
        raise RuntimeError("expected one ThreadBatch include in Eratosthenes.cpp")
    source = source.replace(include, "")
    storage = "std::vector<std::thread> threads;"
    if batch.count(storage) != 1:
        raise RuntimeError("expected one ThreadBatch thread vector")
    batch = batch.replace("#pragma once\n", "", 1).replace(
        storage, "std::vector<SieveThreadFaults::Thread> threads;")
    worker = re.compile(r"(    void PrimeWorker\([^{}]+?\) \{\n)", re.DOTALL)
    if len(list(worker.finditer(source))) != 1:
        raise RuntimeError("expected one PrimeWorker definition")
    source = worker.sub(lambda m: m.group(1) +
                        "        SieveThreadFaults::BeforeWorker();\n", source)
    # Keep the altered test class distinct if a bridge links real kernels too.
    batch = re.sub(r"\bThreadBatch\b", "SieveTestThreadBatch", batch)
    source = re.sub(r"\bThreadBatch\b", "SieveTestThreadBatch", source)
    args.output_directory.mkdir(parents=True, exist_ok=True)
    output = args.output_directory / "prime-sieve-thread-core.cpp"
    output.write_text("// Test-only copy. Eratosthenes.cpp SHA256: " + hashes[0] +
                      "\n// ThreadBatch.h SHA256: " + hashes[1] + "\n" + batch + "\n" + source)
    print(output)


if __name__ == "__main__":
    main()
