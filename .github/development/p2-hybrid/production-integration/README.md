# Current P2 validation

From the repository root, run:

```sh
python3 .github/development/p2-hybrid/production-integration/validate.py /tmp/p2-check
```

The output root's `native-validation/` directory must not exist. The runner
copies current C++ sources and headers, its checkers, and three small frozen
[fixtures](../fixtures/README.md) into that directory. It requires Python 3,
Clang, and working Address/UndefinedBehavior and Thread sanitizers. It modifies
no package files and does not need the old experimental checkout.

Each sanitizer runs the common independent/differential P2 checks, private
stream windows, and helper checks for queue coverage, selector boundaries,
useful worker budgets and actual thread launches. Compiler diagnostics,
sanitizer diagnostics, nonzero exits, or unrecognized checker output fail the
run. Logs, commands, input hashes and parsed check counts are saved in the
output directory, not committed to the development tree.

The fixtures supplement independent trial division and sieve references; they
are not all independent mathematical oracles. See the
[retained findings](../CONSOLIDATED.md) for their purpose and the benchmark
history. Historical variant replay scripts remain recoverable from the commit
identified there; this directory's supported runner tests current production.
