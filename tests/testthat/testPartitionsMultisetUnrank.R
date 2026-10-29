test_that("A width-one multiset with maxZeros == 0 maps correctly", {
    expect_identical(
        partitionsSample(1:2, 1, freqs = c(2, 1),
                         target = 2, sampleVec = 1),
        matrix(2L, nrow = 1L)
    )
})

test_that("A multiset with maxZeros == 0 retries as distinct", {
    expect_identical(
        partitionsSample(1:3, 2, freqs = c(2, 1, 1),
                         target = 5, sampleVec = 1),
        matrix(c(2L, 3L), nrow = 1L)
    )
})

test_that("Zero-valued parts are accepted with positive frequencies", {
    freqs <- c(3, rep(2:1, 5))
    expected <- partitionsGeneral(0:10, 5, freqs = freqs, target = 20)
    expect_identical(nrow(expected), 65L)
    expect_true(any(expected == 0L))
    expect_identical(
        partitionsSample(0:10, 5, freqs = freqs, target = 20,
                         sampleVec = seq_len(nrow(expected))),
        expected
    )
})

test_that("Multiset sampling rejects leading and consecutive zeros in freqs", {
    for (freqs in list(c(0, 2, 1, 1), c(2, 0, 0, 1))) {
        expect_error(
            partitionsSample(1:4, 2, freqs = freqs,
                             target = 5, sampleVec = 1),
            "Each element in freqs must be a positive number"
        )
    }
})

test_that("Multiset sampling rejects zero width before unranking", {
    expect_error(
        partitionsSample(1:4, 0, freqs = c(1, 1, 2, 1),
                         target = 0, sampleVec = 1),
        "m must be a positive whole number"
    )
})

test_that("Multiset sampling rejects invalid double ranks", {
    freqs <- c(1, 1, 2, 1)
    total <- partitionsCount(1:4, 3, freqs = freqs, target = 8)
    for (index in c(-1, 0)) {
        expect_error(
            partitionsSample(1:4, 3, freqs = freqs,
                             target = 8, sampleVec = index),
            "sampleVec must be a positive number"
        )
    }
    expect_error(
        partitionsSample(1:4, 3, freqs = freqs,
                         target = 8, sampleVec = total + 1),
        "exceeds the maximum number of possible results"
    )
})

test_that("Multiset sampling rejects invalid GMP ranks including total + 1", {
    freqs <- rep(2:4, length.out = 100)
    total <- partitionsCount(1:100, 20, freqs = freqs, target = 500)
    expect_s3_class(total, "bigz")
    for (index in list(gmp::as.bigz(-1), gmp::as.bigz(0))) {
        expect_error(
            partitionsSample(1:100, 20, freqs = freqs,
                             target = 500, sampleVec = index),
            "sampleVec must be a positive"
        )
    }
    expect_error(
        partitionsSample(1:100, 20, freqs = freqs,
                         target = 500, sampleVec = total + 1),
        "exceeds the maximum number of possible results"
    )
})

test_that("Multiset unranking reserves candidates after skipping a block", {
    expected <- rbind(c(1L, 3L, 4L), c(2L, 3L, 3L))

    for (threads in 1:2) {
        expect_identical(
            partitionsSample(
                1:4, 3, freqs = c(1, 1, 2, 1), target = 8,
                sampleVec = 1:2, nThreads = threads
            ),
            expected
        )
    }
})

test_that("Multiset unranking reserves candidates after bucket exhaustion", {
    expect_identical(
        partitionsSample(
            1:4, 4, freqs = c(1, 1, 2, 1), target = 9, sampleVec = 1
        ),
        matrix(c(1L, 2L, 3L, 3L), nrow = 1L)
    )
})

test_that("Exhausting a bucket preserves the value just selected", {
    expect_identical(
        partitionsSample(
            1:3, 2, freqs = c(1, 2, 1), target = 4, sampleVec = 1:2
        ),
        rbind(c(1L, 3L), c(2L, 2L))
    )
})

test_that("GMP multiset unranking agrees with unrestricted multiplicities", {
    v <- 1:100
    width <- 20L
    freqs <- rep(width, length(v))
    total <- partitionsCount(v, width, freqs = freqs, target = 500)
    expect_s3_class(total, "bigz")
    indices <- c(gmp::as.bigz(1), gmp::as.bigz("9007199254740993"),
                 total %/% 2, total - 1, total)
    expected <- partitionsSample(
        v, width, repetition = TRUE, target = 500, sampleVec = indices
    )
    for (threads in c(1L, 2L)) {
        expect_identical(
            partitionsSample(v, width, freqs = freqs, target = 500,
                             sampleVec = indices, nThreads = threads),
            expected
        )
    }
})

test_that("GMP multiset unranking respects binding multiplicities", {
    v <- 1:100
    width <- 20L
    freqs <- rep(2:4, length.out = length(v))
    total <- partitionsCount(v, width, freqs = freqs, target = 500)
    expect_s3_class(total, "bigz")
    indices <- c(gmp::as.bigz(1), gmp::as.bigz("9007199254740993"),
                 total %/% 2, total - 1, total)
    actual <- partitionsSample(
        v, width, freqs = freqs, target = 500, sampleVec = indices
    )
    expect_identical(
        partitionsSample(v, width, freqs = freqs, target = 500,
                         sampleVec = indices, nThreads = 2),
        actual
    )
    expect_true(all(rowSums(actual) == 500))
    expect_true(all(apply(actual, 1, function(row) all(diff(row) >= 0))))
    expect_true(all(apply(actual, 1, function(row) {
        all(tabulate(row, nbins = length(v)) <= freqs)
    })))

    # Recover each exact rank by counting preceding suffixes through the
    # public counting API, independently of the unranking implementation.
    rank_row <- function(row) {
        rank <- gmp::as.bigz(1)
        reps <- freqs
        first <- 1L
        remaining_sum <- 500L
        for (i in seq_len(width - 1L)) {
            if (first < row[i]) {
                for (candidate in seq.int(first, row[i] - 1L)) {
                    if (reps[candidate] == 0L) next
                    available <- reps
                    available[candidate] <- available[candidate] - 1L
                    active <- which(v >= candidate & available > 0L)
                    if (sum(available[active]) < width - i) next
                    count <- partitionsCount(
                        v[active], width - i, freqs = available[active],
                        target = remaining_sum - candidate
                    )
                    rank <- rank + gmp::as.bigz(count)
                }
            }
            reps[row[i]] <- reps[row[i]] - 1L
            remaining_sum <- remaining_sum - row[i]
            first <- row[i]
        }
        rank
    }
    for (i in seq_len(nrow(actual))) {
        expect_identical(rank_row(actual[i, ]), indices[i])
    }

    # The first rows come from sequential enumeration, without a lower bound.
    expect_identical(
        partitionsSample(v, width, freqs = freqs, target = 500,
                         sampleVec = 1:12),
        partitionsGeneral(v, width, freqs = freqs, target = 500, upper = 12)
    )
    expect_identical(
        partitionsSample(v - 1L, width, freqs = freqs, target = 480,
                         sampleVec = indices),
        actual - 1L
    )
    expect_identical(
        partitionsSample(v * 2e9, width, freqs = freqs, target = 1e12,
                         sampleVec = indices),
        actual * 2e9
    )
})


test_that("Multiset bounds index the same partitions through both APIs", {
    freqs <- rep(1:5, 2)
    for (mapping in list(c(1L, 0L), c(1L, -1L), c(2e9, 7))) {
        v <- (1:10) * mapping[1L] + mapping[2L]
        target <- 25L * mapping[1L] + 5L * mapping[2L]
        candidates <- comboGeneral(v, 5, freqs = freqs)
        bench <- candidates[rowSums(candidates) == target, , drop = FALSE]
        total <- nrow(bench)
        expect_identical(total, 74L)
        expect_identical(
            partitionsGeneral(v, 5, freqs = freqs, target = target, upper = 5),
            bench[1:5, , drop = FALSE]
        )
        for (bounds in list(c(1L, 5L), c(2L, 6L),
                            c(total - 1L, total), c(total, total))) {
            lower <- bounds[1L]
            upper <- bounds[2L]
            expected <- bench[lower:upper, , drop = FALSE]
            expect_identical(
                partitionsGeneral(v, 5, freqs = freqs, target = target,
                                  lower = lower, upper = upper),
                expected
            )
            expect_identical(
                comboGeneral(v, 5, freqs = freqs, constraintFun = "sum",
                             comparisonFun = "==", limitConstraints = target,
                             lower = lower, upper = upper, keepResults = FALSE),
                expected
            )
            expect_identical(
                comboGeneral(v, 5, freqs = freqs, constraintFun = "sum",
                             comparisonFun = "==", limitConstraints = target,
                             lower = lower, upper = upper, keepResults = TRUE),
                unname(cbind(expected, target))
            )
        }
        expect_identical(
            partitionsGeneral(v, 5, freqs = freqs, target = target,
                              lower = total - 1L),
            bench[(total - 1L):total, , drop = FALSE]
        )
        expect_error(
            partitionsGeneral(v, 5, freqs = freqs, target = target,
                              lower = 500),
            "bounds cannot exceed"
        )
        expect_error(
            comboGeneral(v, 5, freqs = freqs, constraintFun = "sum",
                         comparisonFun = "==", limitConstraints = target,
                         lower = 500),
            "bounds cannot exceed"
        )
    }
    tab <- table(rep(1:10, freqs))
    expected <- partitionsSample(10, 5, freqs = freqs, target = 25,
                                 sampleVec = 2:6)
    expect_identical(
        partitionsGeneral(tab, 5, target = 25, lower = 2, upper = 6), expected
    )
    expect_identical(
        comboGeneral(tab, 5, constraintFun = "sum", comparisonFun = "==",
                     limitConstraints = 25, lower = 2, upper = 6,
                     keepResults = FALSE),
        expected
    )
})

test_that("Large multiset bounds preserve double and GMP partition ranks", {
    for (case in list(list(cap = 50L, width = 12L, target = 240L),
                      list(cap = 100L, width = 20L, target = 500L))) {
        v <- seq_len(case$cap)
        freqs <- rep(2:4, length.out = case$cap)
        total <- partitionsCount(v, case$width, freqs = freqs,
                                 target = case$target)
        expect_identical(inherits(total, "bigz"), case$cap == 100L)
        expect_gt(total, .Machine$integer.max)
        lower <- total - 2L
        expected <- partitionsSample(v, case$width, freqs = freqs,
                                     target = case$target,
                                     sampleVec = lower + 0:2)
        for (bound in list(lower, as.character(lower))) {
            expect_identical(
                partitionsGeneral(v, case$width, freqs = freqs,
                                  target = case$target, lower = bound),
                expected
            )
            expect_identical(
                comboGeneral(v, case$width, freqs = freqs,
                             constraintFun = "sum", comparisonFun = "==",
                             limitConstraints = case$target, lower = bound,
                             upper = total, keepResults = FALSE),
                expected
            )
        }
    }
})
