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
    expect_equal(as.character(partitionsRank(mySamp, v = 500, nThreads = 2)), rownames(mySamp))

    mySamp = partitionsSample(500, 20, TRUE, n = 10, seed = 97, namedSample = TRUE)
    expect_equal(unname(mySamp), partitionsSample(500, 20, TRUE, n = 10, seed = 97, nThreads = 2))
    expect_equal(as.character(partitionsRank(mySamp, v = 500, repetition = TRUE, nThreads = 2)),
                 rownames(mySamp))

    ## Compositions
    mySamp = compositionsSample(500, 20, TRUE, n = 10, seed = 97, namedSample = TRUE)
    expect_equal(unname(mySamp), compositionsSample(500, 20, TRUE, n = 10, seed = 97, nThreads = 2))
    expect_equal(as.character(compositionsRank(mySamp, v = 500, repetition = TRUE, nThreads = 2)), rownames(mySamp))

    ## N.B. Weak = TRUE
    mySamp = compositionsSample(0:500, 20, TRUE, n = 10, weak = TRUE,
                                seed = 97, namedSample = TRUE)
    expect_equal(unname(mySamp), compositionsSample(0:500, 20, TRUE, n = 10, weak = TRUE,
                                                    seed = 97, nThreads = 2))
    expect_equal(as.character(compositionsRank(mySamp, v = 0:500, weak = TRUE,
                                               repetition = TRUE, nThreads = 2)), rownames(mySamp))

    mySamp = partitionsSample(1e9 * (1:100), 10, n = 20, seed = 123, namedSample = TRUE)
    expect_equal(unname(mySamp), partitionsSample(1e9 * (1:100), 10, n = 20,
                                                  seed = 123, nThreads = 2))
    expect_equal(as.character(partitionsRank(mySamp, v = 1e9 * (1:100), nThreads = 2)),
                 rownames(mySamp))

    mySamp = partitionsSample(1e9 * (1:1000), 11, TRUE,
                              n = 20, seed = 8128, namedSample = TRUE)
    expect_equal(unname(mySamp), partitionsSample(1e9 * (1:1000), 11, TRUE,
                                                  n = 20, seed = 8128, nThreads = 2))
    expect_equal(as.character(partitionsRank(mySamp, v = 1e9 * (1:1000),
                                             repetition = TRUE, nThreads = 2)),
                 rownames(mySamp))
    mySamp = compositionsSample(1e9 * (1:1000), 11, TRUE,
                                n = 20, seed = 8128, namedSample = TRUE)
    expect_equal(unname(mySamp), compositionsSample(1e9 * (1:1000), 11, TRUE,
                                                    n = 20, seed = 8128, nThreads = 2))
    expect_equal(as.character(compositionsRank(mySamp, v = 1e9 * (1:1000),
                                               repetition = TRUE, nThreads = 2)),
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

    ##### *********** Torture Test < 2^.Machine$double.digits *********** #####
    ## Partitions
    ## "DstctCapped"
    bench = partitionsGeneral(15, 5, target = 40)
    expect_identical(partitionsSample(15, 5, target = 40, nThreads = 2,
                                      sampleVec = seq_len(nrow(bench))), bench)
    expect_identical(partitionsRank(bench, v = 15, target = 40, nThreads = 2),
                     seq_len(nrow(bench)))

    ## "DstctCappedMZ"
    bench = partitionsGeneral(0:15, 5, freqs = c(3, rep(1, 15)), target = 40)
    expect_identical(partitionsSample(0:15, 5, target = 40, nThreads = 2,
                                      freqs = c(3, rep(1, 15)),
                                      sampleVec = seq_len(nrow(bench))), bench)
    expect_identical(partitionsRank(bench, v = 0:15, target = 40, nThreads = 2,
                                    freqs = c(3, rep(1, 15))),
                     seq_len(nrow(bench)))

    ## "DstctOneZero"
    bench = partitionsGeneral(0:70, 8)
    expect_identical(partitionsSample(0:70, 8, nThreads = 2,
                                      sampleVec = seq_len(nrow(bench))), bench)
    expect_identical(partitionsRank(bench, v = 0:70, nThreads = 2),
                     seq_len(nrow(bench)))

    ## "DstctMultiZero"
    bench = partitionsGeneral(0:60, 8, freqs = c(4, rep(1, 60)))
    expect_identical(partitionsSample(0:60, 8, freqs = c(4, rep(1, 60)),
                                      nThreads = 2,
                                      sampleVec = seq_len(nrow(bench))), bench)
    expect_identical(partitionsRank(bench, v = 0:60,
                                    freqs = c(4, rep(1, 60)), nThreads = 2),
                     seq_len(nrow(bench)))

    ## "DstctStdAll"
    bench = partitionsGeneral(0:50, freqs = c(40, rep(1, 50)))
    expect_identical(partitionsSample(0:50, freqs = c(40, rep(1, 50)),
                                      nThreads = 2,
                                      sampleVec = seq_len(nrow(bench))), bench)
    expect_identical(partitionsRank(bench, v = 0:50, freqs = c(40, rep(1, 50)),
                                    nThreads = 2),
                     seq_len(nrow(bench)))

    ## "RepShort"
    bench = partitionsGeneral(0:40, 5, TRUE)
    expect_identical(partitionsSample(0:40, 5, TRUE, nThreads = 2,
                                      sampleVec = seq_len(nrow(bench))), bench)
    expect_identical(partitionsRank(bench, v = 0:40, repetition = TRUE,
                                    nThreads = 2),
                     seq_len(nrow(bench)))

    ## "RepNoZero"
    bench = partitionsGeneral(40, 6, TRUE)
    expect_identical(partitionsSample(40, 6, TRUE, nThreads = 2,
                                      sampleVec = seq_len(nrow(bench))), bench)
    expect_identical(partitionsRank(bench, v = 40, repetition = TRUE,
                                    nThreads = 2),
                     seq_len(nrow(bench)))

    ## Compositions
    ## "CompRepNoZero"
    bench = compositionsGeneral(25, 5, TRUE)
    expect_identical(compositionsSample(25, 5, TRUE, nThreads = 2,
                                        sampleVec = seq_len(nrow(bench))), bench)
    expect_identical(compositionsRank(bench, v = 25, repetition = TRUE,
                                      nThreads = 2),
                     seq_len(nrow(bench)))

    ## "CmpDstctNoZero"
    bench = compositionsGeneral(28, 5)
    expect_identical(compositionsSample(28, 5, nThreads = 2,
                                        sampleVec = seq_len(nrow(bench))), bench)
    expect_identical(compositionsRank(bench, v = 28, nThreads = 2),
                     seq_len(nrow(bench)))

    ## "CmpDstctCapped"
    bench = compositionsGeneral(15, 5, target = 40)
    expect_identical(compositionsSample(15, 5, target = 40, nThreads = 2,
                                      sampleVec = seq_len(nrow(bench))), bench)
    expect_identical(compositionsRank(bench, v = 15, target = 40, nThreads = 2),
                     seq_len(nrow(bench)))

    ## "CmpDstctWeak"
    bench = compositionsGeneral(0:25, 5, weak = TRUE)
    expect_identical(compositionsSample(0:25, 5, weak = TRUE, nThreads = 2,
                                        sampleVec = seq_len(nrow(bench))), bench)
    expect_identical(compositionsRank(bench, v = 0:25,
                                      weak = TRUE, nThreads = 2),
                     seq_len(nrow(bench)))

    ## "CmpDstCapWeak"
    bench = compositionsGeneral(0:13, 5, weak = TRUE, target = 30)
    expect_identical(compositionsSample(0:13, 5, weak = TRUE,
                                        target = 30, nThreads = 2,
                                        sampleVec = seq_len(nrow(bench))), bench)
    expect_identical(compositionsRank(bench, v = 0:13, weak = TRUE,
                                      target = 30, nThreads = 2),
                     seq_len(nrow(bench)))

    ## "CmpDstctZNotWk"
    bench = compositionsGeneral(0:25, 5, freqs = c(3, rep(1, 25)))
    expect_identical(compositionsSample(0:25, 5, freqs = c(3, rep(1, 25)),
                                        nThreads = 2,
                                        sampleVec = seq_len(nrow(bench))), bench)
    expect_identical(compositionsRank(bench, v = 0:25,
                                      freqs = c(3, rep(1, 25)), nThreads = 2),
                     seq_len(nrow(bench)))

    ## "CmpDstCapMZNotWk"
    bench = compositionsGeneral(0:15, 5, freqs = c(3, rep(1, 15)), target = 40)
    expect_identical(compositionsSample(0:15, 5, target = 40, nThreads = 2,
                                        freqs = c(3, rep(1, 15)),
                                        sampleVec = seq_len(nrow(bench))),
                     bench)
    expect_identical(compositionsRank(bench, v = 0:15, target = 40, nThreads = 2,
                                      freqs = c(3, rep(1, 15))),
                     seq_len(nrow(bench)))

    ## "CmpDstctZNotWk"
    bench = compositionsGeneral(0:30, 5)
    expect_identical(compositionsSample(0:30, 5, nThreads = 2,
                                        sampleVec = seq_len(nrow(bench))), bench)
    expect_identical(compositionsRank(bench, v = 0:30,nThreads = 2),
                     seq_len(nrow(bench)))

    ##### *********** Torture Test > 2^.Machine$double.digits *********** #####

    mySamp = partitionsSample(400, 20, target = 1000,
                              n = 4, namedSample = TRUE, nThreads = 2)
    expect_equal(partitionsRank(mySamp, v = 400,
                                target = 1000, nThreads = 2),
                 rownames(mySamp))

    mySamp = partitionsSample(0:200, 15, freqs = c(7, rep(1, 200)),
                              target = 800, n = 4, namedSample = TRUE,
                              nThreads = 2)
    expect_equal(partitionsRank(mySamp, v = 0:200, freqs = c(7, rep(1, 200)),
                                target = 800, nThreads = 2),
                 rownames(mySamp))

    mySamp = partitionsSample(800, 15, n = 4, namedSample = TRUE, nThreads = 2)
    expect_equal(partitionsRank(mySamp, v = 800, nThreads = 2),
                 rownames(mySamp))

    mySamp = partitionsSample(0:700, 15, n = 4, namedSample = TRUE, nThreads = 2)
    expect_equal(partitionsRank(mySamp, v = 0:700, nThreads = 2),
                 rownames(mySamp))

    mySamp = partitionsSample(0:700, 15, freqs = c(10, rep(1, 700)),
                              n = 4, namedSample = TRUE, nThreads = 2)
    expect_equal(partitionsRank(mySamp, v = 0:700, nThreads = 2,
                                freqs = c(10, rep(1, 700))),
                 rownames(mySamp))

    mySamp = partitionsSample(0:600, freqs = c(600, rep(1, 600)),
                              n = 4, namedSample = TRUE, nThreads = 2)
    expect_equal(partitionsRank(mySamp, v = 0:600, nThreads = 2,
                                freqs = c(600, rep(1, 600))),
                 rownames(mySamp))

    mySamp = partitionsSample(150, 15, TRUE, target = 600,
                              n = 4, namedSample = TRUE, nThreads = 2)
    expect_equal(partitionsRank(mySamp, v = 150, repetition = TRUE,
                                target = 600, nThreads = 2),
                 rownames(mySamp))

    mySamp = partitionsSample(0:600, 15, TRUE,
                              n = 4, namedSample = TRUE, nThreads = 2)
    expect_equal(partitionsRank(mySamp, v = 0:600, repetition = TRUE,
                                nThreads = 2),
                 rownames(mySamp))

    mySamp = partitionsSample(0:300, repetition = TRUE,
                              n = 4, namedSample = TRUE, nThreads = 2)
    expect_equal(partitionsRank(mySamp, v = 0:300, repetition = TRUE,
                                nThreads = 2),
                 rownames(mySamp))

    mySamp = compositionsSample(80, 18, repetition = TRUE,
                                n = 4, namedSample = TRUE, nThreads = 2)
    expect_equal(compositionsRank(mySamp, v = 80, repetition = TRUE,
                                  nThreads = 2),
                 rownames(mySamp))

    mySamp = compositionsSample(0:80, 18, repetition = TRUE,
                                n = 4, namedSample = TRUE, nThreads = 2)
    expect_equal(compositionsRank(mySamp, v = 0:80, repetition = TRUE,
                                  nThreads = 2),
                 rownames(mySamp))

    mySamp = compositionsSample(180, 15,
                                n = 4, namedSample = TRUE, nThreads = 2)
    expect_equal(compositionsRank(mySamp, v = 180, nThreads = 2),
                 rownames(mySamp))

    mySamp = compositionsSample(180, 12, target = 300,
                                n = 4, namedSample = TRUE, nThreads = 2)
    expect_equal(compositionsRank(mySamp, v = 180, target = 300, nThreads = 2),
                 rownames(mySamp))

    mySamp = compositionsSample(0:180, 12, weak = TRUE,
                                n = 4, namedSample = TRUE, nThreads = 2)
    expect_equal(compositionsRank(mySamp, v = 0:180, weak = TRUE, nThreads = 2),
                 rownames(mySamp))

    mySamp = compositionsSample(0:100, 12, weak = TRUE, target = 220,
                                n = 4, namedSample = TRUE, nThreads = 2)
    expect_equal(compositionsRank(mySamp, v = 0:100, target = 220,
                                  weak = TRUE, nThreads = 2),
                 rownames(mySamp))

    mySamp = compositionsSample(0:180, 12, freqs = c(8, rep(1, 180)),
                                n = 4, namedSample = TRUE, nThreads = 2)
    expect_equal(compositionsRank(mySamp, v = 0:180, freqs = c(8, rep(1, 180)),
                                  nThreads = 2),
                 rownames(mySamp))

    mySamp = compositionsSample(0:180, 12, freqs = c(8, rep(1, 180)),
                                n = 4, namedSample = TRUE,
                                nThreads = 2, target = 300)
    expect_equal(compositionsRank(mySamp, v = 0:180, freqs = c(8, rep(1, 180)),
                                  nThreads = 2, target = 300),
                 rownames(mySamp))
})
