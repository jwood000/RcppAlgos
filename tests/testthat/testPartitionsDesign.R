context("testing partitionsDesign")

test_that("irregular inputs retain general constraint classification", {
    v <- c(1L, 2L, 4L, 7L)
    freqs <- c(2L, 2L, 2L, 2L)

    # Omitting m leaves NotPartition and skips SetPartitionDesign. Both native
    # count entry points must initialize the constraint type before this path.
    for (width in list(NULL, 3L)) {
        args <- list(v = v, m = width, freqs = freqs, target = 9L)
        design <- do.call(partitionsDesign, args)
        expect_identical(
            design$partition_type,
            if (is.null(width)) "NotPartition" else "CoarseGrained"
        )
        expect_identical(
            design$constraint_type,
            if (is.null(width)) "General" else "PartitionEsque"
        )
        expect_false(design$is_partition)
        expect_true(design$count_unknown)

        for (count in list(partitionsCount, compositionsCount,
                           partitionsMultisetCount, compositionsMultisetCount)) {
            expect_error(
                do.call(count, args),
                "The count is unknown for this case.", fixed = TRUE
            )
        }
    }
})

test_that("partitionsDesign produces correct results", {
    expect_equal(partitionsDesign(0:10, repetition = TRUE)$partition_type,
                 "RepStdAll")
    expect_equal(compositionsDesign(0:10, repetition = TRUE)$partition_type,
                 "CompRepZero")
    expect_equal(partitionsDesign(
            0:10, freqs = c(10, rep(1, 10))
        )$partition_type, "DistinctStdAll")
    expect_equal(partitionsDesign(0:10, 1)$partition_type,
                 "LengthOne")
    expect_equal(compositionsDesign(0:10, 1)$partition_type,
                 "LengthOne")
    expect_equal(partitionsDesign(0:15, 4, target = 20)$partition_type,
                 "DistinctCapped")
    expect_equal(partitionsDesign(0:15, 4, freqs = c(2, rep(1, 15)),
                                  target = 20)$partition_type,
                 "DistinctCappedMZ")
    expect_equal(partitionsDesign(15, 4, target = 20)$partition_type,
                 "DistinctNoZero")
    expect_equal(partitionsDesign(0:15, 4)$partition_type,
                 "DistinctOneZero")
    expect_equal(partitionsDesign(0:15, 4,
                                  freqs = c(3, rep(1, 15)))$partition_type,
                 "DistinctMZ")
    expect_equal(partitionsDesign(40, 8, freqs = rep(1:5, 8))$partition_type,
                 "Multiset")
    expect_equal(partitionsDesign(40, 8, TRUE, target = 70)$partition_type,
                 "RepCapped")
    expect_equal(partitionsDesign(0:40, 8, repetition = TRUE)$partition_type,
                 "RepShort")
    expect_equal(partitionsDesign(40, 8, repetition = TRUE)$partition_type,
                 "RepNoZero")
    expect_equal(compositionsDesign(0:40, 8, repetition = TRUE)$partition_type,
                 "CompRepZero")
    expect_equal(compositionsDesign(40, 8, repetition = TRUE)$partition_type,
                 "CompRepNoZero")

    verbose <- capture.output(partitionsDesign(15, 4, target = 20,
                                               showDesign = TRUE))
    expect_true(any(grepl("Partition Design Overview", verbose)))

    # https://www.wolframalpha.com/input?i=number+of+distinct+partitions+of+1000
    expect_equal(partitionsDesign(
            0:1000, freqs = c(1000, rep(1, 1000))
        )$num_partitions, gmp::as.bigz("8635565795744155161506"))

    ## N.B. In the partition case, since all values will be less than 2^31 - 1,
    ## the integerness will be preserved. Compare this to the obtaining the sum
    ## of every combination (i.e. the values will exceed 2^31 - 1).
    expect_equal(class(partitionsGeneral(as.integer((1:100) * 2e7),
                                         10, upper = 10)[1, ]), "integer")
    expect_equal(class(comboGeneral(as.integer((1:100) * 2e7), 10, upper = 10,
                                    constraintFun = "sum")[1,]), "numeric")
})
