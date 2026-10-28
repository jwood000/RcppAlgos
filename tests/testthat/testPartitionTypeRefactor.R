# Regression tests for PartitionType discovery, mapped-zero semantics,
# capability guards, and iterator reconstruction.
#
# Focus of this PR:
#   * Positive-domain multiset compositions with a repeatable minimum can
#     use the specialized CompDistinct*MZWeak generation algorithms without
#     requiring weak = TRUE.
#   * The mapped minimum is a genuine part, not non-weak zero padding.
#   * General multiset compositions remain countable even when their
#     generation algorithms are unsupported.
#   * Iterator reconstruction preserves the mapped-zero budget and complement.
#
# Run with the package loaded using your usual testthat workflow.

# ---- Helpers ---------------------------------------------------------------

sort_rows <- function(x) {
    if (nrow(x) < 2L) return(x)
    x[do.call(order, as.data.frame(x)), , drop = FALSE]
}

check_against_permutations <- function(v, m, fr, tar, expected_type = NULL,
                                       expected_n = NULL, weak = FALSE) {
    actual <- compositionsGeneral(v, m, freqs = fr, target = tar, weak = weak)
    oracle <- permuteGeneral(
        v, m, freqs = fr,
        constraintFun = "sum", comparisonFun = "==",
        limitConstraints = tar
    )
    design <- compositionsDesign(v, m, freqs = fr, target = tar, weak = weak)

    expect_identical(actual, sort_rows(oracle))
    expect_identical(nrow(unique(actual)), nrow(actual))
    expect_true(all(rowSums(actual) == tar))

    for (i in seq_along(v)) {
        expect_true(all(rowSums(actual == v[i]) <= fr[i]))
    }

    expect_equal(as.numeric(design$num_partitions), nrow(actual))
    expect_false(design$count_unknown)
    expect_true(design$is_partition)

    if (!is.null(expected_type)) {
        expect_identical(design$partition_type, expected_type)
    }
    if (!is.null(expected_n)) {
        expect_identical(nrow(actual), as.integer(expected_n))
    }

    invisible(actual)
}

# ---- 1. The original positive-domain regression ----------------------------

test_that("positive-domain repeated minima use MZWeak without weak=TRUE", {
    v <- 2:5
    fr <- c(2L, 1L, 1L, 1L)

    x <- check_against_permutations(
        v, 3L, fr, 9L,
        expected_type = "CompDistinctMZWeak", expected_n = 9L
    )
    expect_identical(x[1L, ], c(2L, 2L, 5L))
    expect_identical(x[9L, ], c(5L, 2L, 2L))

    d <- compositionsDesign(v, 3L, freqs = fr, target = 9L)
    expect_identical(d$constraint_type, "PartMapping")
    expect_true(d$eqn_check)
})

# ---- 2. Capped version: the prior 6-versus-60 failure -----------------------

test_that("capped mapped minima are permuted as actual parts", {
    v <- 2:5
    fr <- c(3L, 1L, 1L, 1L)

    x <- check_against_permutations(
        v, 5L, fr, 16L,
        expected_type = "CompDistinctCapMZWeak", expected_n = 60L
    )

    # Every solution uses two 2s and one each of 3, 4, 5.
    expect_true(all(rowSums(x == 2L) == 2L))
    expect_identical(x[1L, ], c(2L, 2L, 3L, 4L, 5L))

    d <- compositionsDesign(v, 5L, freqs = fr, target = 16L)
    expect_identical(d$mapped_vector, 0:3)
    expect_identical(d$mapped_target, 6L)
    expect_identical(d$first_index_vector, c(1L, 1L, 2L, 3L, 4L))
})

# ---- 3. Positive domains: weak flag must not change mathematics ------------

test_that("requesting weak does not change positive-domain compositions", {
    for (case in list(
        list(v = 2:5, m = 3L, fr = c(2L, 1L, 1L, 1L), tar = 9L),
        list(v = 2:5, m = 5L, fr = c(3L, 1L, 1L, 1L), tar = 16L)
    )) {
        args <- list(case$v, case$m, freqs = case$fr, target = case$tar)
        normal <- do.call(compositionsGeneral, args)
        weak <- do.call(compositionsGeneral, c(args, list(weak = TRUE)))
        expect_identical(normal, weak)

        d_normal <- do.call(compositionsDesign, args)
        d_weak <- do.call(compositionsDesign, c(args, list(weak = TRUE)))
        expect_identical(d_normal$partition_type, d_weak$partition_type)
        expect_identical(d_normal$num_partitions, d_weak$num_partitions)
    }
})

# ---- 4. Additional small mapped-domain cases -------------------------------

test_that("mapped multiset composition output agrees with permutations", {
    cases <- list(
        list(v = 2:7, m = 3L, fr = c(2L, rep(1L, 5)), tar = 11L),
        list(v = 1:8, m = 4L, fr = c(3L, rep(1L, 7)), tar = 12L),
        list(v = 3:8, m = 4L, fr = c(2L, rep(1L, 5)), tar = 18L)
    )

    for (case in cases) {
        check_against_permutations(case$v, case$m, case$fr, case$tar)
    }
})

# ---- 5. Translation of the original positive domain -----------------------

test_that("translation preserves composition order and cardinality", {
    fr <- c(2L, 1L, 1L, 1L)
    x <- compositionsGeneral(2:5, 3L, freqs = fr, target = 9L)
    y <- compositionsGeneral(7:10, 3L, freqs = fr, target = 24L)

    expect_identical(y, x + 5L)
    expect_identical(
        compositionsDesign(2:5, 3L, freqs = fr, target = 9L)$partition_type,
        compositionsDesign(7:10, 3L, freqs = fr, target = 24L)$partition_type
    )
})

# ---- 6. Iterator traversal: single results and batches ---------------------

test_that("mapped-minimum iterator matches full lexicographic generation", {
    v <- 1:10
    fr <- c(3L, rep(1L, 9))
    x <- compositionsGeneral(v, 4L, freqs = fr, target = 14L)
    expect_identical(nrow(x), 168L)

    a <- compositionsIter(v, 4L, freqs = fr, target = 14L)
    single <- matrix(0L, nrow = nrow(x), ncol = ncol(x))
    for (i in seq_len(nrow(x))) {
        single[i, ] <- a@nextIter()
    }
    expect_identical(single, x)

    a@startOver()
    chunks <- matrix(0L, nrow = nrow(x), ncol = ncol(x))
    for (i in seq_len(3L)) {
        rows <- seq.int((i - 1L) * 56L + 1L, i * 56L)
        chunks[rows, ] <- a@nextNIter(56L)
    }
    expect_identical(chunks, x)
})

# ---- 7. Random access and reconstructing the zero budget -------------------

test_that("random-access continuation can reintroduce the mapped minimum", {
    v <- 1:10
    fr <- c(3L, rep(1L, 9))
    x <- compositionsGeneral(v, 4L, freqs = fr, target = 14L)
    a <- compositionsIter(v, 4L, freqs = fr, target = 14L)

    # Locate a zero-free mapped result whose successor contains the minimum.
    minimum_count <- rowSums(x == min(v))
    transitions <- which(
        head(minimum_count, -1L) == 0L &
        tail(minimum_count, -1L) > 0L
    )
    expect_gt(length(transitions), 0L)

    idx <- transitions[1L]
    expect_identical(as.integer(a[[idx]]), as.integer(x[idx, ]))
    expect_identical(a@nextIter(), x[idx + 1L, ])

    # Rebuild again, and cross that same transition using batch iteration.
    a[[idx - 1L]]
    rows <- seq.int(idx, min(idx + 9L, nrow(x)))
    expect_identical(a@nextNIter(length(rows)), x[rows, , drop = FALSE])
})

# ---- 8. General multiset: count is supported; generation is not ------------

test_that("general multiset composition allows counts but rejects generation", {
    v <- 1:7
    fr <- c(2L, 2L, rep(1L, 5))
    d <- compositionsDesign(v, 4L, freqs = fr, target = 14L)

    expect_identical(d$partition_type, "CompMultiset")
    expect_false(d$count_unknown)
    expect_gt(as.numeric(d$num_partitions), 0)
    expect_equal(
        as.numeric(compositionsCount(v, 4L, freqs = fr, target = 14L)),
        as.numeric(d$num_partitions)
    )

    # Independent enumeration; count() and design() share the DP.
    oracle <- permuteGeneral(
        v, 4L, freqs = fr, constraintFun = "sum",
        comparisonFun = "==", limitConstraints = 14L
    )
    expect_equal(as.numeric(d$num_partitions), nrow(oracle))

    expect_error(
        compositionsGeneral(v, 4L, freqs = fr, target = 14L),
        "permuteGeneral"
    )
    expect_error(
        compositionsIter(v, 4L, freqs = fr, target = 14L),
        "permuteIter"
    )
})

# ---- 9. Literal-zero semantics remain distinct ----------------------------

test_that("actual zero distinguishes weak and zero-padded multisets", {
    v <- 0:6
    fr <- c(4L, 2L, 2L, 1L, 1L, 1L, 1L)
    args <- list(v, 4L, freqs = fr, target = 9L)

    padded <- do.call(compositionsDesign, args)
    weak <- do.call(compositionsDesign, c(args, list(weak = TRUE)))

    expect_identical(padded$partition_type, "CompMultisetZero")
    expect_identical(weak$partition_type, "CompMultisetWeak")
    expect_false(padded$count_unknown)
    expect_false(weak$count_unknown)
    expect_equal(
        as.numeric(do.call(compositionsCount, args)),
        as.numeric(padded$num_partitions)
    )
    expect_equal(
        as.numeric(do.call(compositionsCount, c(args, list(weak = TRUE)))),
        as.numeric(weak$num_partitions)
    )

    # Weak case: literal zero is a genuine part. Fixed-width permutation
    # enumeration is an independent oracle for the DP count.
    weak_oracle <- permuteGeneral(
        v, 4L, freqs = fr, constraintFun = "sum",
        comparisonFun = "==", limitConstraints = 9L
    )
    expect_equal(as.numeric(weak$num_partitions), nrow(weak_oracle))

    # Nonweak case: zero is padding, not a permuted part. Enumerate
    # positive-only compositions at every permitted positive length.
    positive <- v[-1L]
    positive_fr <- fr[-1L]
    oracle_count <- 0L
    for (k in seq_len(4L)) {
        if (4L - k > fr[1L]) next  # Insufficient zero padding.
        if (k * max(positive) < 9L) next  # Impossible target at this length.
        sub <- permuteGeneral(
            positive, k, freqs = positive_fr,
            constraintFun = "sum", comparisonFun = "==",
            limitConstraints = 9L
        )
        oracle_count <- oracle_count + nrow(sub)
    }
    expect_equal(as.numeric(padded$num_partitions), oracle_count)
})

# ---- 10. Partition and ordinary-distinct controls --------------------------

test_that("partition generation remains consistent with partition counts", {
    v <- 1:10
    fr <- c(3L, rep(1L, 9))
    x <- partitionsGeneral(v, 4L, freqs = fr, target = 14L)
    d <- partitionsDesign(v, 4L, freqs = fr, target = 14L)

    expect_equal(as.numeric(d$num_partitions), nrow(x))
    expect_equal(
        as.numeric(partitionsCount(v, 4L, freqs = fr, target = 14L)),
        nrow(x)
    )
})

test_that("ordinary distinct compositions are unchanged by weak on positive v", {
    v <- 1:10
    fr <- rep(1L, length(v))
    args <- list(v, 4L, freqs = fr, target = 14L)

    x <- do.call(compositionsGeneral, args)
    y <- do.call(compositionsGeneral, c(args, list(weak = TRUE)))
    expect_identical(x, y)
    expect_match(do.call(compositionsDesign, args)$partition_type, "^CompDistinct")
    expect_identical(nrow(unique(x)), nrow(x))
})

# ---- 11. Capped iterator: step, batches, and repositioning -------------------

test_that("capped MZWeak iterator supports all traversal modes", {
    v <- 2:5
    fr <- c(3L, 1L, 1L, 1L)
    x <- compositionsGeneral(v, 5L, freqs = fr, target = 16L)
    expect_identical(nrow(x), 60L)

    a <- compositionsIter(v, 5L, freqs = fr, target = 16L)
    single <- matrix(0L, nrow(x), ncol(x))
    for (i in seq_len(nrow(x))) {
        single[i, ] <- a@nextIter()
    }
    expect_identical(single, x)

    a@startOver()
    chunks <- matrix(0L, nrow(x), ncol(x))
    for (rows in list(1:17, 18:40, 41:60)) {
        chunks[rows, ] <- a@nextNIter(length(rows))
    }
    expect_identical(chunks, x)

    expect_identical(as.integer(a[[29L]]), as.integer(x[29L, ]))
    expect_identical(a@nextIter(), x[30L, ])

    a[[17L]]
    expect_identical(a@nextNIter(12L), x[18:29, , drop = FALSE])
})

# ---- 12. Rank / unrank both mapped-MZWeak families --------------------------

test_that("both mapped MZWeak types rank and unrank with weak=FALSE", {
    cases <- list(
        list(v = 2:5, m = 3L, fr = c(2L, 1L, 1L, 1L), tar = 9L,
             ptype = "CompDistinctMZWeak"),
        list(v = 2:5, m = 5L, fr = c(3L, 1L, 1L, 1L), tar = 16L,
             ptype = "CompDistinctCapMZWeak")
    )
    for (case in cases) {
        x <- compositionsGeneral(
            case$v, case$m, freqs = case$fr, target = case$tar,
            weak = FALSE
        )
        d <- compositionsDesign(
            case$v, case$m, freqs = case$fr, target = case$tar,
            weak = FALSE
        )
        expect_identical(d$partition_type, case$ptype)

        ranks <- compositionsRank(
            x, v = case$v, freqs = case$fr, target = case$tar,
            weak = FALSE
        )
        expect_identical(as.integer(ranks), seq_len(nrow(x)))

        indices <- unique(c(1L, 2L, nrow(x) %/% 2L, nrow(x)))
        recovered <- compositionsSample(
            case$v, case$m, freqs = case$fr, target = case$tar,
            weak = FALSE, sampleVec = indices
        )
        expect_identical(recovered, x[indices, , drop = FALSE])
        expect_identical(
            as.integer(compositionsRank(
                recovered, v = case$v, freqs = case$fr,
                target = case$tar, weak = FALSE
            )),
            indices
        )
    }
})

# ---- 13. lower / upper bounded generation ----------------------------------

test_that("bounded mapped compositions equal full lexicographic slices", {
    for (case in list(
        list(v = 2:5, m = 3L, fr = c(2L, 1L, 1L, 1L), tar = 9L),
        list(v = 2:5, m = 5L, fr = c(3L, 1L, 1L, 1L), tar = 16L),
        list(v = 1:10, m = 4L, fr = c(3L, rep(1L, 9)), tar = 14L)
    )) {
        x <- compositionsGeneral(
            case$v, case$m, freqs = case$fr, target = case$tar
        )
        n <- nrow(x)
        lo <- max(2L, n %/% 3L)
        hi <- min(n - 1L, lo + 6L)
        actual <- compositionsGeneral(
            case$v, case$m, freqs = case$fr, target = case$tar,
            lower = lo, upper = hi
        )
        expect_identical(actual, x[lo:hi, , drop = FALSE])
        expect_identical(
            compositionsGeneral(
                case$v, case$m, freqs = case$fr,
                target = case$tar, upper = 1L
            ),
            x[1L, , drop = FALSE]
        )
    }
})

# ---- 14. Parallel generation with a sufficiently large result set -----------

test_that("threaded mapped-composition generation matches serial output", {
    # Much larger than the 9/60/168-row fixtures. Whether threads are
    # actually launched is implementation/platform dependent; equality
    # alone does not prove a particular worker count was used.
    v <- 1:20
    fr <- c(3L, rep(1L, 19))
    m <- 6L
    tar <- 50L

    serial <- compositionsGeneral(
        v, m, freqs = fr, target = tar, nThreads = 1L
    )
    # Independently: choose distinct values from 2:20 with sum 50 - t,
    # then permute them alongside t copies of 1 (t = 0, 1, 2, 3).
    # 371*720 + 362*720 + 139*360 + 14*120 = 579480
    expect_identical(nrow(serial), 579480L)
    parallel <- compositionsGeneral(
        v, m, freqs = fr, target = tar, nThreads = 2L
    )
    expect_identical(parallel, serial)

    lo <- 5001L
    hi <- 450000L
    expect_identical(
        compositionsGeneral(
            v, m, freqs = fr, target = tar,
            lower = lo, upper = hi, nThreads = 2L
        ),
        serial[lo:hi, , drop = FALSE]
    )
})

# ------------ 15. Predicates and Factories have matching Enums ---------------

test_that("Ranking capabilities match dispatch factories", {
    expect_true(testRankCoverage())
})

test_that("Counting capabilities match dispatch factories", {
    expect_true(testCountCoverage())
})

# ---- 16. Unsupported generic composition fallbacks -------------------------

test_that("generic composition fallbacks reject generation and iteration", {
    for (case in list(
        list(v = c(1L, 2L, 4L, 7L), target = 9L),
        list(v = c(0L, 1L, 3L, 7L), target = 4L)
    )) {
        for (width in list(3L, NULL)) {
            for (weak in c(FALSE, TRUE)) {
                args <- c(case, list(m = width, repetition = TRUE, weak = weak))
                design <- do.call(compositionsDesign, args)
                expect_identical(
                    design$partition_type,
                    if (is.null(width)) "NotPartition" else "CoarseGrained"
                )
                for (bounds in list(list(), list(upper = 3L),
                                    list(lower = 1L, upper = 20L),
                                    list(lower = 2L, upper = 20L))) {
                    expect_error(
                        do.call(compositionsGeneral, c(args, bounds)),
                        "Generation is not supported for this composition.",
                        fixed = TRUE
                    )
                }
                expect_error(
                    do.call(compositionsIter, args),
                    "Generation is not supported for this composition.",
                    fixed = TRUE
                )
            }
        }
    }
})

test_that("generic combination and permutation fallbacks retain their bounds", {
    v <- c(1L, 2L, 4L, 7L)
    for (generate in list(comboGeneral, permuteGeneral)) {
        for (width in list(3L, NULL)) {
            args <- list(v = v, m = width, repetition = TRUE)
            candidates <- do.call(generate, args)
            matches <- candidates[rowSums(candidates) == 9L, , drop = FALSE]
            constrained <- c(args, list(constraintFun = "sum",
                                       comparisonFun = "==",
                                       limitConstraints = 9L))
            actual <- do.call(generate, constrained)
            expect_identical(sort_rows(actual), matches)
            expect_identical(
                do.call(generate, c(constrained, list(upper = 3L))),
                actual[seq_len(min(3L, nrow(actual))), , drop = FALSE]
            )
            for (lower in 1:2) {
                window <- candidates[lower:20L, , drop = FALSE]
                expect_identical(
                    do.call(generate, c(constrained,
                                        list(lower = lower, upper = 20L))),
                    window[rowSums(window) == 9L, , drop = FALSE]
                )
            }
        }
    }
})

# ---- 17. Known-empty results retain their bounds and shape -----------------

test_that("NoSolution uses matching-result bounds for partitions", {
    for (case in list(
        # Mapping establishes infeasibility; counting can also discover zero.
        list(v = 1:5, m = 3L, target = 100L),
        list(v = 1:5, m = 3L, target = 5L),
        list(v = 1:50, m = 10L, freqs = rep(2L, 50L), target = 10000L)
    )) {
        expect_identical(do.call(partitionsDesign, case)$partition_type,
                         "NoSolution")
        expect_identical(do.call(compositionsDesign, case)$partition_type,
                         "NoSolution")
        for (generate in list(partitionsGeneral, compositionsGeneral)) {
            expect_identical(dim(do.call(generate, case)), c(0L, case$m))
            for (bounds in list(list(lower = 1L), list(upper = 1L),
                                list(lower = 1L, upper = 3L),
                                list(lower = gmp::as.bigz(1)))) {
                expect_error(do.call(generate, c(case, bounds)),
                             "bounds cannot exceed")
            }
        }
        args <- list(v = case$v, m = case$m, freqs = case$freqs,
                     constraintFun = "sum", comparisonFun = "==",
                     limitConstraints = case$target)
        expect_identical(dim(do.call(comboGeneral, args)), c(0L, case$m))
        for (bounds in list(list(lower = 1L), list(upper = 1L),
                            list(lower = 1L, upper = 3L))) {
            expect_error(do.call(comboGeneral, c(args, bounds)),
                         "bounds cannot exceed")
        }
    }
})

test_that("empty constrained matrices retain the result column and type", {
    for (v in list(1:5, as.double(1:5))) {
        for (target in c(5L, 100L)) {
            for (generate in list(comboGeneral, permuteGeneral)) {
                for (keep in c(FALSE, TRUE)) {
                    actual <- generate(
                        v, 3L, constraintFun = "sum", comparisonFun = "==",
                        limitConstraints = target, keepResults = keep
                    )
                    expect_identical(actual, matrix(v[FALSE], 0L, 3L + keep))
                }
            }
        }
    }
})

test_that("NoSolution preserves permutation candidate bounds", {
    for (case in list(list(), list(repetition = TRUE),
                      list(freqs = c(2L, 1L, 2L, 1L, 1L)))) {
        total <- do.call(permuteCount, c(list(v = 1:5, m = 3L), case))
        args <- c(list(v = 1:5, m = 3L, constraintFun = "sum",
                       comparisonFun = "==", limitConstraints = 100L), case)
        for (bounds in list(list(lower = 1L), list(lower = 1L, upper = total),
                            list(lower = total, upper = total))) {
            expect_identical(
                do.call(permuteGeneral,
                        c(args, bounds, list(keepResults = TRUE))),
                matrix(integer(), 0L, 4L)
            )
        }
        for (bounds in list(list(lower = total + 1),
                            list(lower = 1L, upper = total + 1),
                            list(upper = 1L))) {
            expect_error(do.call(permuteGeneral, c(args, bounds)),
                         "bounds cannot exceed")
        }
        for (bounds in list(list(lower = 2L, upper = 1L),
                            list(lower = 3L, upper = 1L))) {
            expect_error(do.call(permuteGeneral, c(args, bounds)),
                         "lowerBound.*upperBound")
        }
        for (bounds in list(list(lower = 0L), list(lower = -1L),
                            list(lower = gmp::as.bigz(0)),
                            list(lower = 1L, upper = 0L))) {
            expect_error(do.call(permuteGeneral, c(args, bounds)))
        }
    }
})

test_that("known-empty GMP permutation ranges do not size or scan candidates", {
    args <- list(v = 1:50, m = 10L, constraintFun = "sum",
                 comparisonFun = "==", limitConstraints = 10000L)
    total <- gmp::as.bigz(permuteCount(1:50, 10L))
    expect_gt(total, gmp::as.bigz(2)^53)
    for (bounds in list(list(lower = 1L), list(lower = 1L, upper = total),
                        list(lower = total, upper = total))) {
        expect_identical(do.call(permuteGeneral, c(args, bounds)),
                         matrix(integer(), 0L, 10L))
    }
    for (bounds in list(list(lower = total + 1L),
                        list(lower = 1L, upper = total + 1L))) {
        expect_error(do.call(permuteGeneral, c(args, bounds)),
                     "bounds cannot exceed")
    }
    expect_error(
        do.call(permuteGeneral,
                c(args, list(lower = total, upper = total - 1L))),
        "lowerBound.*upperBound"
    )
})

test_that("known-empty iterators report zero without searching", {
    check_empty_iterator <- function(a) {
        expect_equal(a@summary()$totalResults, 0)
        for (get_next in list(a@nextIter, function() a@nextNIter(2L),
                              a@nextRemaining)) {
            a@startOver()
            expect_equal(a@summary()$currentIndex, 0)
            expect_equal(a@summary()$totalRemaining, 0)
            invisible(capture.output(value <- a@currIter()))
            expect_null(value)
            invisible(capture.output(value <- get_next()))
            expect_null(value)
            invisible(capture.output(value <- a@currIter()))
            expect_null(value)
        }
        a@startOver()
        expect_equal(a@summary()$currentIndex, 0)
        if (is(a, "Partitions")) {
            expect_error(a@front(), "No random access available")
            expect_error(a@back(), "No random access available")
        }
    }

    # The last case enters setup with a GMP-sized ordinary permutation count.
    expect_gt(gmp::as.bigz(permuteCount(1:50, 10L)), gmp::as.bigz(2)^53)
    for (case in list(list(v = 1:5, m = 3L, target = 5L),
                      list(v = 1:5, m = 3L, target = 100L),
                      list(v = 1:5, target = 100L),
                      list(v = 1:50, m = 10L, target = 10000L))) {
        for (make_iter in list(partitionsIter, compositionsIter)) {
            check_empty_iterator(do.call(make_iter, case))
        }
        for (make_iter in list(comboIter, permuteIter)) {
            a <- make_iter(case$v, case$m, constraintFun = "sum",
                           comparisonFun = "==", limitConstraints = case$target,
                           keepResults = TRUE)
            expect_match(a@summary()$description,
                         if (identical(make_iter, permuteIter)) {
                             "^Permutations of partitions "
                         } else "^Partitions ")
            check_empty_iterator(a)
        }
    }
})

test_that("known-empty generation preserves inferred and candidate widths", {
    for (generate in list(partitionsGeneral, compositionsGeneral)) {
        expect_identical(generate(1:5, target = 100L),
                         matrix(integer(), 0L, 1L))
    }
    args <- list(v = 1:5, constraintFun = "sum", comparisonFun = "==",
                 limitConstraints = 100L)
    for (keep in c(FALSE, TRUE)) {
        for (generate in list(comboGeneral, permuteGeneral)) {
            expect_identical(
                do.call(generate, c(args, list(keepResults = keep))),
                matrix(integer(), 0L, 1L + keep)
            )
        }
        # With lower, omitted m means the ordinary permutation width of five.
        expect_identical(
            do.call(permuteGeneral, c(args, list(lower = 1L, upper = 3L,
                                                keepResults = keep))),
            matrix(integer(), 0L, 5L + keep)
        )
    }
})

# ---- Follow-up fixture to add when confirmed --------------------------------
# * Find a CURRENT, reproducible PartMapping + NotMapped composition after
#   the updated DiscoverPType classification. Do not reuse the former
#   2:5, m=3, freqs=c(2,1,1,1), target=9 case: it now resolves as
#   CompDistinctMZWeak and correctly supports generation.
# * When identified, assert design() reports NotMapped/count_unknown,
#   while count(), general(), and iter() reject it before dispatch.
# * If a bounded search finds no reachable NotMapped composition, treat
#   that as evidence to audit the classification invariant rather than
#   inventing an artificial regression fixture.

# * Any additional mapped distinct types reached by the retry-to-distinct path.
