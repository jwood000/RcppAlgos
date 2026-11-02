test_that("primeSieveInspect preserves results and reports each sieve kernel", {
    # Fixed dispatch inputs make these kernel choices independent of the host.
    inspect <- function(lower, upper) {
        .Call(RcppAlgos:::`_RcppAlgos_PrimeSieveInspectCpp`,
              lower, upper, 1L, 1L, 1L, NULL)
    }
    lowerBounds <- c(Small = 1e6, Medium = 1.5e9, Big = 1e14)

    for (kernel in names(lowerBounds)) {
        lower <- unname(lowerBounds[[kernel]])
        values <- seq.int(lower, lower + 1000)
        expected <- values[isPrimeRcpp(values)]
        if (tail(values, 1L) <= .Machine$integer.max) {
            expected <- as.integer(expected)
        }
        result <- inspect(lower, lower + 1000)

        expect_identical(result$primes, expected, info = kernel)
        expect_identical(unique(result$kernels$kernel), kernel)
        expect_equal(result$summary$requested_threads, 1)
        expect_equal(result$summary$effective_threads, 1)
        expect_null(result$summary$requested_segment_bytes)
        expect_equal(result$summary$effective_segment_span, 0)

        expect_s3_class(result$phases, "data.frame")
        expect_named(result$phases, c("worker", "kernel", "phase", "seconds", "calls"))
        expect_true(nrow(result$phases) > 0)
        expect_true(all(is.finite(result$phases$seconds)))
        expect_true(all(result$phases$seconds >= 0))
        expect_true(all(result$phases$calls >= 1))
        expect_true(all(result$phases$worker == 1))

        expect_s3_class(result$kernels, "data.frame")
        expect_named(result$kernels, c(
            "worker", "kernel", "lower", "upper", "segment_span",
            "sieve_bytes", "offset_bytes", "bucket_entries",
            "bucket_capacity_bytes", "buffer_bytes", "prime_count",
            "output_capacity_bytes", "segments"
        ))
        expect_equal(sum(result$kernels$prime_count), length(expected))
        expect_true(all(result$kernels$segments >= 1))
        memoryFields <- c("sieve_bytes", "offset_bytes", "bucket_entries",
                          "bucket_capacity_bytes", "buffer_bytes",
                          "output_capacity_bytes")
        expect_true(all(as.matrix(result$kernels[memoryFields]) >= 0))
        expect_named(result$shared, c("base_prime_bytes", "dense_prime_bytes",
                                      "bucket_prime_bytes", "wheel_bytes"))
        expect_true(all(unlist(result$shared) >= 0))
        elapsed <- unlist(result$summary[c("base_seconds", "policy_seconds",
                                          "sieve_seconds", "result_seconds")])
        expect_true(all(is.finite(elapsed) & elapsed >= 0))
    }

    # A single worker can call more than one kernel at a dispatch boundary.
    cutoff <- (32340 - 100)^2
    values <- seq.int(cutoff - 100, cutoff + 100)
    result <- inspect(cutoff - 100, cutoff + 100)
    expect_identical(result$primes, as.integer(values[isPrimeRcpp(values)]))
    expect_identical(result$kernels$kernel, c("Small", "Medium"))
    expect_equal(sum(result$kernels$prime_count), length(result$primes))
})

test_that("primeSieveInspect keeps public bound handling and host context", {
    inspect <- RcppAlgos:::primeSieveInspect
    for (bounds in list(c(1, 1), c(100, 100), c(101, 101),
                       c(10.2, 10.8), c(1e6 + 101, 1e6))) {
        result <- inspect(bounds[1], bounds[2], nThreads = 1L)
        expect_identical(result$primes,
                         primeSieve(bounds[1], bounds[2], nThreads = 1L))
    }
    result <- inspect(1000)
    expect_identical(result$primes, primeSieve(1000))
    expect_named(result$host, c("system", "R", "physical_cores", "max_threads"))
    expect_identical(result$host$R, R.version.string)
})

test_that("Big segment budgets preserve boundaries and repeated bucket carry", {
    inspect <- function(lower, upper, bytes) {
        .Call(RcppAlgos:::`_RcppAlgos_PrimeSieveInspectCpp`,
              lower, upper, 1L, 1L, 1L, bytes)
    }
    segmentSize <- 30030
    boundary <- floor(1e14 / segmentSize) * segmentSize
    intervals <- list(
        aligned = c(boundary, boundary + 100),
        crossing = c(boundary - 51, boundary + 51),
        several_segments = c(boundary + 17, boundary + 4 * segmentSize + 101)
    )

    for (caseName in names(intervals)) {
        bounds <- intervals[[caseName]]
        values <- seq.int(bounds[1], bounds[2])
        expected <- values[isPrimeRcpp(values)]
        result <- inspect(bounds[1], bounds[2], 3754)
        expect_identical(result$primes, expected, info = caseName)
        expect_equal(result$summary$requested_segment_bytes, 3754)
        expect_equal(result$summary$effective_segment_span, segmentSize)
        expect_true(all(result$kernels$segment_span == segmentSize))
    }

    # A budget need not be an exact wheel payload; round down, never up.
    result <- inspect(boundary - 51, boundary + 51, 8000)
    expect_identical(result$primes,
                     inspect(boundary - 51, boundary + 51, 3754)$primes)
    expect_equal(result$summary$effective_segment_span, 2 * segmentSize)
})

test_that("primeSieveInspect rejects invalid or inapplicable segment budgets", {
    inspect <- function(bytes, lower = 1e14) {
        .Call(RcppAlgos:::`_RcppAlgos_PrimeSieveInspectCpp`,
              lower, lower + 100, 1L, 1L, 1L, bytes)
    }
    for (bytes in list(NA_real_, NaN, Inf, 0, -1, 3753, 3754.5,
                       1e20, c(3754, 8000), numeric(), "8000", TRUE)) {
        expect_error(inspect(bytes), "segmentBytes")
    }
    expect_error(inspect(8000, 1e6), "Big")
    expect_error(inspect(8000, 1.5e9), "Big")

    # Bucket preparation may advance twice after its last admitted strike.
    # Reject this span before allocating buckets or evaluating int offsets.
    expect_error(.Call(RcppAlgos:::`_RcppAlgos_PrimeSieveInspectCpp`,
                       1e15, 1e15 + 2.05e9, 1L, 1L, 1L, 3754),
                 "segmentBytes")
})

test_that("primeSieveInspect records actual parallel work and serial fallback", {
    maxThreads <- get("pkgEnv", asNamespace("RcppAlgos"))$maxThreads
    skip_if(maxThreads < 2L, "At least two supported threads are required")
    inspect <- RcppAlgos:::primeSieveInspect

    lower <- 32340 + 17
    upper <- lower + 17 * 32340 + 101
    result <- inspect(lower, upper, nThreads = 2L)
    expect_identical(result$primes, primeSieve(lower, upper, nThreads = 1L))
    expect_equal(result$summary$requested_threads, 2)
    expect_equal(result$summary$effective_threads, 2)
    expect_equal(sort(unique(result$kernels$worker)), 1:2)
    expect_equal(sum(result$kernels$prime_count), length(result$primes))

    result <- inspect(1e6, 1e6 + 100, nThreads = 2L)
    expect_equal(result$summary$requested_threads, 2)
    expect_equal(result$summary$effective_threads, 1)
    expect_equal(unique(result$kernels$worker), 1)
})

test_that("Big inspection preserves output when requested workers are clamped", {
    maxThreads <- get("pkgEnv", asNamespace("RcppAlgos"))$maxThreads
    skip_if(maxThreads < 3L, "At least three supported threads are required")
    segmentSize <- 30030
    lower <- floor(1e14 / segmentSize) * segmentSize + 17
    upper <- lower + 13 * segmentSize + 101
    values <- seq.int(lower, upper)
    expected <- values[isPrimeRcpp(values)]

    # Ten requested output chunks, but only three actual workers: retain the
    # final worker's output even when unused chunk slots separate it from others.
    result <- .Call(RcppAlgos:::`_RcppAlgos_PrimeSieveInspectCpp`,
                    lower, upper, 10L, 1L, 3L, 3754)
    expect_identical(result$primes, expected)
    expect_equal(result$summary$requested_threads, 10)
    expect_equal(result$summary$effective_threads, 3)
    expect_equal(sort(unique(result$kernels$worker)), 1:3)
    expect_equal(sum(result$kernels$prime_count), length(expected))
    expect_true(all(result$kernels$segment_span == segmentSize))
})
