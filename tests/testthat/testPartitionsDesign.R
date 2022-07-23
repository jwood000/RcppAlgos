context("testing partitionsDesign")

test_that("partitionsDesign produces correct results", {
    expect_equal(partitionsDesign(0:10, repetition = TRUE)$partition_type,
                 "RepStdAll")
    expect_equal(compositionsDesign(0:10, repetition = TRUE)$partition_type,
                 "RepStdAll")
    expect_equal(partitionsDesign(
            0:10, freqs = c(10, rep(1, 10))
        )$partition_type, "DstctStdAll")
    expect_equal(partitionsDesign(0:10, 1)$partition_type,
                 "LengthOne")
    expect_equal(compositionsDesign(0:10, 1)$partition_type,
                 "LengthOne")
    expect_equal(partitionsDesign(0:15, 4, target = 20)$partition_type,
                 "DistCapped")
    expect_equal(partitionsDesign(0:15, 4, freqs = c(2, rep(1, 15)),
                                  target = 20)$partition_type,
                 "DstctCappedMZ")
    expect_equal(partitionsDesign(15, 4, target = 20)$partition_type,
                 "DstctNoZero")
    expect_equal(partitionsDesign(0:15, 4)$partition_type,
                 "DstctOneZero")
    expect_equal(partitionsDesign(0:15, 4,
                                  freqs = c(3, rep(1, 15)))$partition_type,
                 "DstctMultiZero")
    expect_equal(partitionsDesign(40, 8, freqs = rep(1:5, 8))$partition_type,
                 "Multiset")
    expect_equal(partitionsDesign(40, 8, TRUE, target = 70)$partition_type,
                 "RepCapped")
    expect_equal(partitionsDesign(0:40, 8, repetition = TRUE)$partition_type,
                 "RepShort")
    expect_equal(partitionsDesign(40, 8, repetition = TRUE)$partition_type,
                 "RepNoZero")
    expect_equal(compositionsDesign(0:40, 8, repetition = TRUE)$partition_type,
                 "RepShort")
    expect_equal(compositionsDesign(40, 8, repetition = TRUE)$partition_type,
                 "RepNoZero")

    verbose <- capture.output(partitionsDesign(15, 4, target = 20,
                                               showDesign = TRUE))
    expect_true(any(grepl("Partition Design Overview", verbose)))

    # https://www.wolframalpha.com/input?i=number+of+distinct+partitions+of+1000
    expect_equal(partitionsDesign(
            0:1000, freqs = c(1000, rep(1, 1000))
        )$num_partitions, gmp::as.bigz("8635565795744155161506"))

    ## N.B. In the partition case, since all values will be less than 2^31 - 1,
    ## the integerness will be preserverd. Compare this to the obtaining the sum
    ## of every combination (i.e. the values will exceed 2^31 - 1).
    expect_equal(class(partitionsGeneral(as.integer((1:100) * 2e7),
                                         10, upper = 10)[1, ]), "integer")
    expect_equal(class(comboGeneral(as.integer((1:100) * 2e7), 10, upper = 10,
                                    constraintFun = "sum")[1,]), "numeric")
})
