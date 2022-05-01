context("testing partitionsSample")

test_that("parttionsSample produces correct results", {

    expect_equal(partitionsSample((1:10) * 1e13, 1, n = 1), matrix(1e14))
    expect_equal(partitionsSample(100, 1, n = 1), matrix(100L))
    expect_equal(partitionsSample(100, 1, TRUE, n = 1), matrix(100L))
    expect_equal(partitionsSample(100, 1, n = 1, target = 40L), matrix(40L))
    expect_equal(partitionsSample(100, 1, TRUE, n = 1, target = 40L), matrix(40L))

    expect_equal(partitionsSample(15, 3, sampleVec = 1:12), partitionsGeneral(15, 3))
    expect_equal(nrow(partitionsSample(15, 3, n = 5, seed = 2)), 5)
    expect_equal(as.vector(partitionsSample(1, 1, n = 1)), 1)
    expect_equal(class(partitionsSample(c(1:15, NA), 3, n = 3)[1, ]), "integer")

    totalParts = partitionsCount(10, 15, TRUE, target = 27)
    expect_equal(partitionsSample(10, 15, TRUE, target = 27, sampleVec = 1:totalParts),
                 partitionsGeneral(10, 15, TRUE, target = 27))

    expect_equal(partitionsCount(500, 25), partitionsCount(500, 25, freqs = rep(1, 500)))

    expect_equal(partitionsSample(0:10, repetition = TRUE,
                                  sampleVec = 1:42),
                 partitionsGeneral(0:10, repetition = TRUE))

    expect_equal(partitionsSample(500, 20, n = 10, seed = 42),
                 partitionsSample(500, 20, n = 10, seed = 42, nThreads = 2))

    expect_equal(partitionsSample(500, 20, TRUE, n = 10, seed = 97),
                 partitionsSample(500, 20, TRUE, n = 10, seed = 97, nThreads = 2))

    expect_equal(partitionsSample(1e9 * (1:100), 10, n = 20, seed = 123),
                 partitionsSample(1e9 * (1:100), 10, n = 20,
                                  seed = 123, nThreads = 2))

    expect_equal(partitionsSample(1e9 * (1:1000), 11, TRUE,
                                  n = 20, seed = 8128),
                 partitionsSample(1e9 * (1:1000), 11, TRUE,
                                  n = 20, seed = 8128, nThreads = 2))

    expect_equal(as.numeric(rownames(
        partitionsSample(100, 10, sampleVec = seq(1, 10000, 1000),
                         namedSample = TRUE)
        )), seq(1, 10000, 1000))
})
