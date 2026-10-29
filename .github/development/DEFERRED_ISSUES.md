# Deferred Investigations

## Unsafe floating-point to integer conversions

Status: Needs investigation
Priority: High
Discovered during: PartitionType refactoring

### Observation

SetTolerance(), CheckIsInteger(), and CheckPartition()
may cast out-of-range, NaN, or infinite double values
to std::int64_t, resulting in undefined behavior.

### Proposed investigation

- Reproduce the issue with targeted tests.
- Review all affected conversion paths.
- Evaluate std::isfinite() and std::round() approaches.
- Agree on the fix before implementation.

### Decision

Defer to a dedicated investigation and PR.
