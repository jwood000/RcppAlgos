test_that("scalar sequence detection preserves finite input semantics", {
    cases <- list(
        list(v = 0, expected = 0L),
        list(v = 3, expected = 1:3),
        list(v = -3, expected = -3:-1),
        list(v = 3.5, expected = 3.5),
        list(v = -3.5, expected = -3.5)
    )

    for (case in cases) {
        expect_identical(comboGeneral(case$v, 1),
                         matrix(case$expected, ncol = 1L))
        expect_identical(comboCount(case$v, 1), length(case$expected))
    }
})

test_that("invalid numeric scalars reach the existing input errors", {
    for (fun in list(comboGeneral, comboCount, comboIter)) {
        for (value in c(NA_real_, NaN)) {
            expect_error(fun(value, 1), "cannot be NA or NaN")
        }
        for (value in c(Inf, -Inf, 2^53, -2^53, 2^63, -2^63)) {
            expect_error(fun(value, 1), "must be less than 2\\^53")
        }
    }
})

test_that("constraints preserve large and infinite vector elements", {
    cases <- list(
        list(v = c(-1e100, 0, 1e100), comparison = "<=",
             expected = c(-1e100, 0)),
        list(v = c(0, Inf), comparison = ">=", expected = Inf),
        list(v = c(-Inf, 0), comparison = "<=", expected = c(-Inf, 0)),
        list(v = c(1, 1e20), comparison = "==", expected = 1),
        # Exactly 2^63 must remain outside the int64 partition path.
        list(v = c(1, 2^63), comparison = "==", expected = 1),
        list(v = c(1, Inf), comparison = "==", expected = 1)
    )

    for (case in cases) {
        for (tolerance in list(NULL, 0)) {
            args <- list(case$v, 1, constraintFun = "sum",
                         comparisonFun = case$comparison,
                         limitConstraints = 1, keepResults = FALSE,
                         tolerance = tolerance)
            expected <- matrix(case$expected, ncol = 1L)
            expect_identical(do.call(comboGeneral, args), expected)
            iter <- do.call(comboIter, args)
            expect_identical(iter@nextNIter(nrow(expected)), expected)
        }
    }
})

test_that("constraint targets preserve the integer output boundary", {
    expected <- rbind(c(1L, 1L, 2L), c(1L, 2L, 3L), c(2L, 2L, 4L))
    for (target in c(2^31 - 1, 2^31)) {
        actual <- comboGeneral(1:2, 2, repetition = TRUE,
                               constraintFun = "sum", comparisonFun = "<=",
                               limitConstraints = target, keepResults = TRUE)
        if (target == 2^31) {
            storage.mode(expected) <- "double"
        }
        expect_identical(actual, expected)
    }
})

test_that("mean targets can scale beyond the signed 64-bit range", {
    # The finite mean limit is accepted even when conversion to a sum limit
    # crosses 2^63. The first candidate is independently known to be all ones.
    for (width in c(1024L, 1025L)) {
        for (tolerance in list(NULL, 0)) {
            actual <- comboGeneral(1:2, width, repetition = TRUE,
                                   constraintFun = "mean", comparisonFun = "<=",
                                   limitConstraints = 2^53 - 1, upper = 1,
                                   keepResults = FALSE, tolerance = tolerance)
            expect_identical(actual, matrix(1, nrow = 1L, ncol = width))
        }
    }
})

test_that("large mapped partitions retain an exact slope", {
    # Both endpoints and the target are exact doubles, but their difference
    # 2^53 + 1 is not. Of the three repeated pairs, only the mixed pair matches.
    v <- c(-(2^53 - 1), 2)
    target <- -(2^53 - 3)
    expected <- matrix(v, nrow = 1L)

    for (tolerance in list(NULL, 0)) {
        expect_identical(
            partitionsGeneral(v, 2, repetition = TRUE, target = target,
                              tolerance = tolerance),
            expected
        )
        expect_identical(
            comboGeneral(v, 2, repetition = TRUE, constraintFun = "sum",
                         comparisonFun = "==", limitConstraints = target,
                         keepResults = FALSE, tolerance = tolerance),
            expected
        )
    }
    expect_identical(partitionsCount(v, 2, repetition = TRUE, target = target), 1L)
    expect_identical(
        partitionsSample(v, 2, repetition = TRUE, target = target, sampleVec = 1),
        expected
    )
    expect_identical(
        partitionsRank(expected, v = v, repetition = TRUE, target = target), 1L
    )

    iter <- partitionsIter(v, 2, repetition = TRUE, target = target)
    expect_identical(iter@nextIter(), v)
    iter@startOver()
    expect_identical(iter@nextNIter(1), expected)
    expect_identical(iter[[1]], v)
})

test_that("mapped searches reject aggregates outside the native range", {
    # Individual values fit int64_t, whereas even the smallest pair sums to
    # 2^63. Reject the search before evaluating those aggregates.
    v <- c(2^62, 2^62 + 1024)
    for (tolerance in list(NULL, 0)) {
        expect_error(
            partitionsGeneral(v, 2, repetition = TRUE, target = 1,
                              tolerance = tolerance),
            "partition mapping search exceeds the supported integer range"
        )
        expect_error(
            comboGeneral(v, 2, repetition = TRUE, constraintFun = "sum",
                         comparisonFun = "==", limitConstraints = 1,
                         keepResults = FALSE, tolerance = tolerance),
            "partition mapping search exceeds the supported integer range"
        )
    }
})

test_that("mapped searches retain the supported aggregate boundary", {
    # Each value is below 2^53. The largest width-1024 aggregate fits int64_t;
    # every width-1025 aggregate exceeds it. Both have no sum equal to one.
    v <- c(2^53 - 2, 2^53 - 1)
    expect_identical(
        partitionsGeneral(v, 1024, repetition = TRUE, target = 1),
        matrix(numeric(), nrow = 0L, ncol = 1024L)
    )
    expect_error(
        partitionsGeneral(v, 1025, repetition = TRUE, target = 1),
        "partition mapping search exceeds the supported integer range"
    )
})

test_that("unspecified widths are checked before searching smaller widths", {
    a <- 2^60
    v <- c(-3 * a, 0, 3 * a, 6 * a)
    # All four elements sum to 6*a. Of the four triples, only the first
    # three sum to zero, but the initial search exceeds the native range.
    expect_error(
        partitionsGeneral(v, target = 0),
        "partition mapping search exceeds the supported integer range"
    )
    expect_error(
        partitionsCount(v, target = 0),
        "partition mapping search exceeds the supported integer range"
    )
})

test_that("mapping range checks handle the minimum signed integer", {
    # abs(INT64_MIN) cannot be represented in int64_t. Reject its source
    # range without first attempting to take that absolute value.
    v <- c(-2^63, -2^63 + 2048)
    expect_error(
        partitionsGeneral(v, 1, target = 0),
        "partition mapping search exceeds the supported integer range"
    )
})

test_that("mapping arithmetic rejects overflowing intermediates", {
    v <- c(-2^61, 0, 2^61)
    # The search fits, but slope * mapped_target is 6 * 2^61. The supported
    # arithmetic range excludes this intermediate even though its quotient fits.
    expect_error(
        partitionsGeneral(v, 3, target = 0),
        "partition mapping arithmetic exceeds the supported integer range"
    )
    expect_error(
        partitionsDesign(v, 3, target = 0),
        "partition mapping arithmetic exceeds the supported integer range"
    )
    # Ranking can preserve descending source order and a negative slope.
    descending <- rev(v)
    expect_error(
        partitionsRank(matrix(descending, nrow = 1L),
                       v = descending, target = 0),
        "partition mapping arithmetic exceeds the supported integer range"
    )
    # The product itself fits here: slope * mapped_target is 2^63 - 2048.
    # Subtracting the negative target would then produce 2^63.
    v <- c(-2^61 - 512, -1024, 2^61 - 1536)
    expect_error(
        partitionsGeneral(v, 2, repetition = TRUE, target = -2048),
        "partition mapping arithmetic exceeds the supported integer range"
    )
})

test_that("large mapped multisets retain their multiplicities and order", {
    a <- 2^60
    v <- c(-a, 0, a, 2 * a)
    expected <- rbind(c(-a, -a, 2 * a), c(-a, 0, a))
    expect_identical(
        partitionsGeneral(v, 3, freqs = c(2, 1, 1, 1), target = 0),
        expected
    )
})

test_that("unrepresentable partition spacing retains generic constraint results", {
    v <- c(-2^62, 2^62)
    # Both values fit int64_t, but their positive spacing of 2^63 does not.
    expect_identical(
        comboGeneral(v, 2, constraintFun = "sum", comparisonFun = "==",
                     limitConstraints = 0, keepResults = FALSE),
        matrix(v, nrow = 1L)
    )
})

test_that("partition mapping reports unsupported search and target ranges", {
    expect_error(
        partitionsGeneral(c(-2^62, 2^61), 3, repetition = TRUE, target = 0),
        "partition mapping search exceeds the supported integer range"
    )
    # The mapped target is 50000 * 50001, exceeding INT_MAX. Preparation
    # reaches the range check without attempting to enumerate these results.
    expect_error(
        partitionsGeneral(seq.int(-50000L, 50000L), 50000L,
                          repetition = TRUE, target = 0),
        "mapped partition target exceeds the supported integer range"
    )
    # Standard weak-composition mapping adds the width to the target. Use a
    # double source so validation reaches the guard without an integer row.
    expect_error(
        compositionsCount(c(0, 1), .Machine$integer.max,
                          repetition = TRUE, weak = TRUE),
        "mapped partition target exceeds the supported integer range"
    )
})
