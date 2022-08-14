context("testing partitionsSample and compositionsSample")

test_that("parttionsSample and compositionsSample produces correct results", {

    expect_identical(partitionsRank(v = 100), integer(0))
    expect_identical(compositionsRank(v = 100), integer(0))

    mySamp = partitionsSample((1:10) * 1e13, 1, n = 1, namedSample = TRUE)
    expect_equal(unname(mySamp), matrix(1e14))
    expect_equal(as.character(partitionsRank(mySamp, v = (1:10) * 1e13)),
                 rownames(mySamp))
    mySamp = compositionsSample((1:10) * 1e13, 1, n = 1, namedSample = TRUE)
    expect_equal(unname(mySamp), matrix(1e14))
    expect_equal(as.character(compositionsRank(mySamp, v = (1:10) * 1e13)),
                 rownames(mySamp))

    mySamp = partitionsSample(100, 1, n = 1, namedSample = TRUE)
    expect_equal(unname(mySamp), matrix(100L))
    expect_equal(as.character(partitionsRank(mySamp, v = 100)),
                 rownames(mySamp))

    mySamp = partitionsSample(100, 1, TRUE, n = 1, namedSample = TRUE)
    expect_equal(unname(mySamp), matrix(100L))
    expect_equal(as.character(partitionsRank(mySamp, v = 100, repetition = TRUE)),
                 rownames(mySamp))
    mySamp = compositionsSample(100, 1, TRUE, n = 1, namedSample = TRUE)
    expect_equal(unname(mySamp), matrix(100L))
    expect_equal(as.character(compositionsRank(mySamp, v = 100, repetition = TRUE)),
                 rownames(mySamp))

    mySamp = partitionsSample(100, 1, n = 1, target = 40L, namedSample = TRUE)
    expect_equal(unname(mySamp), matrix(40L))
    expect_equal(as.character(partitionsRank(mySamp, v = 100, target = 40L)),
                 rownames(mySamp))

    mySamp = partitionsSample(100, 1, TRUE, n = 1, target = 40L, namedSample = TRUE)
    expect_equal(unname(mySamp), matrix(40L))
    expect_equal(as.character(partitionsRank(mySamp, v = 100, repetition = TRUE, target = 40L)),
                 rownames(mySamp))

    mySamp = partitionsSample(15, 3, sampleVec = 1:12, namedSample = TRUE)
    expect_equal(unname(mySamp), partitionsGeneral(15, 3))
    expect_equal(as.character(partitionsRank(mySamp, v = 15)),
                 rownames(mySamp))

    mySamp = partitionsSample(15, 3, n = 5, seed = 2, namedSample = TRUE)
    expect_equal(nrow(mySamp), 5)
    expect_equal(as.character(partitionsRank(mySamp, v = 15)),
                 rownames(mySamp))

    mySamp = partitionsSample(1, 1, n = 1, namedSample = TRUE)
    expect_equal(as.vector(mySamp), 1)
    expect_equal(as.character(partitionsRank(mySamp, v = 1)),
                 rownames(mySamp))

    mySamp = partitionsSample(1, 1, n = 1, namedSample = TRUE)
    expect_equal(as.vector(mySamp), 1)
    expect_equal(as.character(partitionsRank(mySamp, v = 1)),
                 rownames(mySamp))

    mySamp = partitionsSample(c(1:15, NA), 3, n = 3, namedSample = TRUE)
    expect_true(is.integer(mySamp[1, ]))
    expect_equal(as.character(partitionsRank(mySamp, v = 15)),
                 rownames(mySamp))

    totalParts = partitionsCount(10, 15, TRUE, target = 27)
    mySamp = partitionsSample(10, 15, TRUE, target = 27,
                              sampleVec = 1:totalParts, namedSample = TRUE)
    expect_equal(unname(mySamp), partitionsGeneral(10, 15, TRUE, target = 27))
    expect_equal(as.character(partitionsRank(mySamp, v = 10,
                                             repetition = TRUE,
                                             target = 27)), rownames(mySamp))

    expect_equal(partitionsCount(500, 25), partitionsCount(500, 25, freqs = rep(1, 500)))

    mySamp = partitionsSample(0:10, repetition = TRUE,
                              sampleVec = 1:42, namedSample = TRUE)
    expect_equal(unname(mySamp), partitionsGeneral(0:10, repetition = TRUE))
    expect_equal(as.character(partitionsRank(mySamp, v = 0:10,
                                             repetition = TRUE)), rownames(mySamp))

    ## Compositions
    mySamp = compositionsSample(0:10, repetition = TRUE,
                                sampleVec = 1:512, namedSample = TRUE)
    expect_equal(unname(mySamp), compositionsGeneral(0:10, repetition = TRUE))
    expect_equal(as.character(compositionsRank(mySamp, v = 0:10,
                                               repetition = TRUE)), rownames(mySamp))

    ## N.B. Weak = TRUE
    mySamp = compositionsSample(0:7, repetition = TRUE, weak = TRUE,
                                sampleVec = 1:1716, namedSample = TRUE)
    expect_equal(unname(mySamp), compositionsGeneral(0:7, repetition = TRUE,
                                                     weak = TRUE))
    expect_equal(compositionsRank(mySamp, v = 0:7, repetition = TRUE,
                                  weak = TRUE),
                 as.integer(rownames(mySamp)))

    mySamp = partitionsSample(500, 20, n = 10, seed = 42, namedSample = TRUE)
    expect_equal(unname(mySamp), partitionsSample(500, 20, n = 10, seed = 42, nThreads = 2))
    expect_equal(as.character(partitionsRank(mySamp, v = 500)), rownames(mySamp))

    mySamp = partitionsSample(500, 20, TRUE, n = 10, seed = 97, namedSample = TRUE)
    expect_equal(unname(mySamp), partitionsSample(500, 20, TRUE, n = 10, seed = 97, nThreads = 2))
    expect_equal(as.character(partitionsRank(mySamp, v = 500, repetition = TRUE)), rownames(mySamp))

    ## Compositions
    mySamp = compositionsSample(500, 20, TRUE, n = 10, seed = 97, namedSample = TRUE)
    expect_equal(unname(mySamp), compositionsSample(500, 20, TRUE, n = 10, seed = 97, nThreads = 2))
    expect_equal(as.character(compositionsRank(mySamp, v = 500, repetition = TRUE)), rownames(mySamp))

    ## N.B. Weak = TRUE
    mySamp = compositionsSample(0:500, 20, TRUE, n = 10, weak = TRUE,
                                seed = 97, namedSample = TRUE)
    expect_equal(unname(mySamp), compositionsSample(0:500, 20, TRUE, n = 10, weak = TRUE,
                                                    seed = 97, nThreads = 2))
    expect_equal(as.character(compositionsRank(mySamp, v = 0:500, weak = TRUE,
                                               repetition = TRUE)), rownames(mySamp))

    mySamp = partitionsSample(1e9 * (1:100), 10, n = 20, seed = 123, namedSample = TRUE)
    expect_equal(unname(mySamp), partitionsSample(1e9 * (1:100), 10, n = 20,
                                                  seed = 123, nThreads = 2))
    expect_equal(as.character(partitionsRank(mySamp, v = 1e9 * (1:100))),
                 rownames(mySamp))

    mySamp = partitionsSample(1e9 * (1:1000), 11, TRUE,
                              n = 20, seed = 8128, namedSample = TRUE)
    expect_equal(unname(mySamp), partitionsSample(1e9 * (1:1000), 11, TRUE,
                                                  n = 20, seed = 8128, nThreads = 2))
    expect_equal(as.character(partitionsRank(mySamp, v = 1e9 * (1:1000), repetition = TRUE)),
                 rownames(mySamp))
    mySamp = compositionsSample(1e9 * (1:1000), 11, TRUE,
                                n = 20, seed = 8128, namedSample = TRUE)
    expect_equal(unname(mySamp), compositionsSample(1e9 * (1:1000), 11, TRUE,
                                                    n = 20, seed = 8128, nThreads = 2))
    expect_equal(as.character(compositionsRank(mySamp, v = 1e9 * (1:1000), repetition = TRUE)),
                 rownames(mySamp))

    mySamp = partitionsSample(100, 10, sampleVec = seq(1, 10000, 1000),
                              namedSample = TRUE)
    expect_equal(as.numeric(rownames(mySamp)), seq(1, 10000, 1000))
    expect_equal(partitionsRank(mySamp, v = 100), seq(1, 10000, 1000))

    mySamp = partitionsSample(40, 3, repetition = TRUE, target = 110,
                              n = 14, namedSample = TRUE)
    expect_equal(sort(as.integer(rownames(mySamp))), 1:14)
    expect_equal(partitionsRank(mySamp, v = 40,
                                repetition = TRUE,
                                target = 110),
                 as.integer(rownames(mySamp)))
})
