zero_partition_cases <- list(
    list(name = "partitions", general = partitionsGeneral,
         count = partitionsCount, sample = partitionsSample,
         rank = partitionsRank, iter = partitionsIter, extra = list()),
    list(name = "compositions", general = compositionsGeneral,
         count = compositionsCount, sample = compositionsSample,
         rank = compositionsRank, iter = compositionsIter,
         extra = list(weak = FALSE)),
    list(name = "weak compositions", general = compositionsGeneral,
         count = compositionsCount, sample = compositionsSample,
         rank = compositionsRank, iter = compositionsIter,
         extra = list(weak = TRUE))
)

for (case in zero_partition_cases) {
    test_that(paste("zero-only", case$name, "retain their requested width"), {
        common <- c(list(v = 0L, repetition = TRUE, target = 0L), case$extra)
        args <- c(common, list(m = 3L))
        expected <- matrix(0L, nrow = 1L, ncol = 3L)

        expect_identical(do.call(case$general, args), expected)
        expect_identical(do.call(case$count, args), 1L)
        expect_identical(
            do.call(case$sample, c(args, list(sampleVec = 1L))), expected
        )
        expect_identical(do.call(case$rank, c(list(expected), common)), 1L)

        iter <- do.call(case$iter, args)
        expect_identical(iter@nextIter(), expected[1L, ])
        invisible(capture.output(exhausted <- iter@nextIter()))
        expect_null(exhausted)

        iter@startOver()
        invisible(capture.output(batch <- iter@nextNIter(2L)))
        expect_identical(batch, expected)
        invisible(capture.output(exhausted <- iter@nextIter()))
        expect_null(exhausted)

        iter@startOver()
        expect_identical(iter[[1L]], expected[1L, ])
        invisible(capture.output(exhausted <- iter@nextIter()))
        expect_null(exhausted)
    })

    test_that(paste("zero-only", case$name, "preserve defaults and output type"), {
        args <- c(list(v = 0L, repetition = TRUE, target = 0L), case$extra)
        expect_identical(do.call(case$general, args), matrix(0L, 1L, 1L))
        expect_identical(do.call(case$count, args), 1L)

        # Duplicate input values collapse to one repeatable zero while the
        # numeric input type and explicit width remain part of the result.
        args$v <- c(0, 0)
        args$m <- 2L
        expected <- matrix(0, nrow = 1L, ncol = 2L)
        expect_identical(do.call(case$general, args), expected)
        expect_identical(
            do.call(case$sample, c(args, list(sampleVec = 1L))), expected
        )
        iter <- do.call(case$iter, args)
        expect_identical(iter@nextIter(), expected[1L, ])

        args$target <- 1L
        expect_identical(
            do.call(case$general, args), matrix(numeric(), 0L, 2L)
        )
        expect_identical(do.call(case$count, args), 0L)
    })
}

test_that("adjacent zero domains retain their existing semantics", {
    expected <- matrix(0L, nrow = 1L, ncol = 3L)
    expect_identical(partitionsGeneral(0L, 3L, freqs = 3L, target = 0L),
                     expected)
    expect_identical(partitionsCount(0L, 3L, freqs = 3L, target = 0L), 1L)
    expect_identical(partitionsGeneral(0:1, 3L, repetition = TRUE, target = 0L),
                     expected)
    expect_identical(
        compositionsGeneral(0:1, 3L, repetition = TRUE,
                            target = 0L, weak = TRUE),
        expected
    )

    # The non-weak width is capped at the positive target. Weak compositions
    # retain the requested width, with zero allowed in either position.
    expect_identical(
        compositionsGeneral(0:1, 2L, repetition = TRUE,
                            target = 1L, weak = FALSE),
        matrix(1L, nrow = 1L, ncol = 1L)
    )
    expect_identical(
        compositionsGeneral(0:1, 2L, repetition = TRUE,
                            target = 1L, weak = TRUE),
        rbind(c(0L, 1L), c(1L, 0L))
    )
})

test_that("constrained combinations and permutations retain repeated zeros", {
    expected <- matrix(0L, nrow = 1L, ncol = 3L)
    for (fun in list(comboGeneral, permuteGeneral)) {
        expect_identical(
            fun(0L, 3L, repetition = TRUE, constraintFun = "sum",
                comparisonFun = "==", limitConstraints = 0L,
                keepResults = FALSE),
            expected
        )
    }
})
