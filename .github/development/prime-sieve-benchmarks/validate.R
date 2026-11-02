# Source after loading the development RcppAlgos library. All reference work
# happens outside the performance measurements in the benchmark driver.
validate_inspector <- function(inspect, output_dir, big_lower, source_root) {
    testthat::test_file(file.path(source_root, "tests/testthat/testCacheInspector.R"),
                        reporter = "summary", stop_on_failure = TRUE)
    testthat::test_file(file.path(source_root, "tests/testthat/testPrimeSieve.R"),
                        reporter = "summary", stop_on_failure = TRUE)

    # Select a prime square above the supplied naturally-Big lower bound.
    # seq.int can return integers here; square in double precision, where
    # these integer products remain exact within the supported 2^53 domain.
    candidates <- as.double(seq.int(ceiling(sqrt(big_lower + 100)), length.out = 2000L))
    candidates <- candidates[RcppAlgos::isPrimeRcpp(candidates)]
    candidates <- candidates[candidates * candidates - 100 >= big_lower]
    if (!length(candidates)) {
        stop("Could not select a prime square above big_lower")
    }
    square <- candidates[1L] * candidates[1L]
    maximum <- 2^53 - 1
    if (square + 101 > maximum) {
        stop("The selected prime square is too close to the supported maximum")
    }

    budgets <- c(3754, 8000, 64 * 1024, 1024 * 1024, 4 * 1024 * 1024)
    results <- vector("list", length(budgets) * 6L)
    index <- 0L

    for (budget in budgets) {
        span <- floor(8 * budget / 30030) * 30030
        # Move forward so even the lower crossing endpoint remains above
        # big_lower. These budget-specific bounds test geometry, not speed.
        boundary <- (floor(big_lower / span) + 1) * span
        if (boundary - 51 < big_lower) boundary <- boundary + span
        intervals <- list(
            aligned_start = c(boundary, boundary + 100),
            crossing_boundary = c(boundary - 51, boundary + 51),
            square_at_start = c(square, square + 100),
            square_at_end = c(square - 100, square),
            after_square = c(square + 1, square + 101),
            near_supported_maximum = c(maximum - 1002, maximum - 2)
        )

        for (case_name in names(intervals)) {
            bounds <- intervals[[case_name]]
            values <- seq.int(bounds[1L], bounds[2L])
            expected <- values[RcppAlgos::isPrimeRcpp(values)]
            actual <- inspect(bounds[1L], bounds[2L], nThreads = 1L,
                              segmentBytes = budget)
            index <- index + 1L
            results[[index]] <- data.frame(
                case = case_name,
                lower = bounds[1L],
                upper = bounds[2L],
                budget = budget,
                prime_count = length(actual$primes),
                pass = identical(actual$primes, expected),
                stringsAsFactors = FALSE
            )
        }
    }

    results <- do.call(rbind, results)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(results, file.path(output_dir, "correctness.csv"), row.names = FALSE)
    if (!all(results$pass)) {
        stop("Inspector output disagreed with independent primality checks; see correctness.csv")
    }
    results
}
