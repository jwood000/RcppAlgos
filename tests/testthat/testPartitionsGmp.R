context("testing partitionGeneral GMP")

test_that("partitionGeneral Distinct Parallel Lower GMP", {

    #######################################################################
    ## See commentary above
    ########********************** Lower Only *******************##########
    ## partitionsDesign(1000, 15)[c("num_partitions", "partition_type")]
    ## $num_partitions
    ## Big Integer ('bigz') :
    ## [1] 3649675516801903698
    ##
    ## $partition_type
    ## [1] "DstctNoZero"
    bench <- partitionsGeneral(1000, 15, lower = "3649675516801803698")
    expect_identical(partitionsGeneral(1000, 15, lower = "3649675516801803698",
                                       nThreads = 2), bench)
    expect_identical(partitionsRank(bench[1, ], v = 1000),
                     gmp::as.bigz("3649675516801803698"))
    expect_equal(gmp::sub.bigz("3649675516801903698",
                               "3649675516801803698") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 1000))

    ## partitionsDesign((1:1000) * 2e9, 15)[c("num_partitions",
    ##                                        "partition_type",
    ##                                        "mapped_target")]
    ## $num_partitions
    ## Big Integer ('bigz') :
    ## [1] 3649675516801903698
    ##
    ## $partition_type
    ## [1] "DstctNoZero"
    ##
    ## $mapped_target
    ## [1] 1000
    bench <- partitionsGeneral((1:1000) * 2e9, 15,
                               lower = "3649675516801803698")
    expect_identical(partitionsGeneral((1:1000) * 2e9, 15,
                                       lower = "3649675516801803698",
                                       nThreads = 2), bench)
    expect_identical(partitionsRank(bench[1, ], v = (1:1000) * 2e9),
                     gmp::as.bigz("3649675516801803698"))
    expect_equal(gmp::sub.bigz("3649675516801903698",
                               "3649675516801803698") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 1000 * 2e9))

    ## partitionsDesign(0:1000, 15)[c("num_partitions",
    ##                                "partition_type",
    ##                                "mapped_target")]
    ## $num_partitions
    ## Big Integer ('bigz') :
    ## [1] 4556757507869210155
    ##
    ## $partition_type
    ## [1] "DstctOneZero"
    ##
    ## $mapped_target
    ## [1] 1015
    bench <- partitionsGeneral(0:1000, 15,
                               lower = "4556757507869110155")
    expect_identical(partitionsGeneral(0:1000, 15,
                                       lower = "4556757507869110155",
                                       nThreads = 2), bench)
    expect_identical(partitionsRank(bench[1, ], v = 0:1000),
                     gmp::as.bigz("4556757507869110155"))
    expect_equal(gmp::sub.bigz("4556757507869210155",
                               "4556757507869110155") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 1000))

    ## partitionsDesign(0:1000 * 2e9, 15)[c("num_partitions",
    ##                                      "partition_type",
    ##                                      "mapped_target")]
    ## $num_partitions
    ## Big Integer ('bigz') :
    ## [1] 4556757507869210155
    ##
    ## $partition_type
    ## [1] "DstctOneZero"
    ##
    ## $mapped_target
    ## [1] 1015
    bench <- partitionsGeneral(0:1000 * 2e9, 15,
                               lower = "4556757507869110155")
    expect_identical(partitionsGeneral(0:1000 * 2e9, 15,
                                       lower = "4556757507869110155",
                                       nThreads = 2), bench)
    expect_identical(partitionsRank(bench[1, ], v = (0:1000) * 2e9),
                     gmp::as.bigz("4556757507869110155"))
    expect_equal(gmp::sub.bigz("4556757507869210155",
                               "4556757507869110155") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 1000 * 2e9))

    ## partitionsDesign(0:1000, 12,
    ##                  freqs = c(3, rep(1, 1000)))[c("num_partitions",
    ##                                                "mapped_target",
    ##                                                "partition_type")]
    ## $num_partitions
    ## Big Integer ('bigz') :
    ## [1] 39228755152043560
    ##
    ## $mapped_target
    ## [1] 1000
    ##
    ## $partition_type
    ## [1] "DstctMultiZero"
    bench <- partitionsGeneral(0:1000, 12,
                               lower = "39228755151943560",
                               freqs = c(3, rep(1, 1000)))
    expect_identical(partitionsGeneral(0:1000, 12,
                                       lower = "39228755151943560",
                                       freqs = c(3, rep(1, 1000)),
                                       nThreads = 2), bench)
    expect_identical(partitionsRank(bench[1, ], v = 0:1000,
                                    freqs = c(3, rep(1, 1000))),
                     gmp::as.bigz("39228755151943560"))
    expect_equal(gmp::sub.bigz("39228755152043560",
                               "39228755151943560") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 1000))

    #### Mapped version
    ##
    ## 15 * 12 + 1000 * 3 = 3180
    ##
    ## partitionsDesign(15 + 0:1000 * 3, 12,
    ##                  freqs = c(3, rep(1, 1000)),
    ##                  target = 3180)[c("num_partitions",
    ##                                   "mapped_target",
    ##                                   "partition_type")]
    ## $num_partitions
    ## Big Integer ('bigz') :
    ## [1] 39228755152043560
    ##
    ## $mapped_target
    ## [1] 1000
    ##
    ## $partition_type
    ## [1] "DstctMultiZero"
    bench <- partitionsGeneral(15 + 0:1000 * 3, 12,
                               lower = "39228755151943560",
                               freqs = c(3, rep(1, 1000)),
                               target = 3180)
    expect_identical(partitionsGeneral(15 + 0:1000 * 3, 12,
                                       lower = "39228755151943560",
                                       freqs = c(3, rep(1, 1000)),
                                       target = 3180, nThreads = 2), bench)
    expect_identical(partitionsRank(bench[1, ], v = 15 + 0:1000 * 3,
                                    freqs = c(3, rep(1, 1000)),
                                    target = 3180),
                     gmp::as.bigz("39228755151943560"))
    expect_equal(gmp::sub.bigz("39228755152043560",
                               "39228755151943560") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 3180))

    ## partitionsDesign(0:1000, 12,
    ##                  freqs = c(13, rep(1, 1000)))[c("num_partitions",
    ##                                                 "mapped_target",
    ##                                                 "partition_type")]
    ## $num_partitions
    ## Big Integer ('bigz') :
    ## [1] 39233351450439724
    ##
    ## $mapped_target
    ## [1] 1000
    ##
    ## $partition_type
    ## [1] "DstctMultiZero"
    bench <- partitionsGeneral(0:1000, 12,
                               lower = "39233351450339724",
                               freqs = c(13, rep(1, 1000)))
    expect_identical(partitionsGeneral(0:1000, 12,
                                       freqs = c(13, rep(1, 1000)),
                                       lower = "39233351450339724",
                                       nThreads = 2), bench)
    expect_identical(partitionsRank(bench[1, ], v = 0:1000,
                                    freqs = c(13, rep(1, 1000))),
                     gmp::as.bigz("39233351450339724"))
    expect_equal(gmp::sub.bigz("39233351450439724",
                               "39233351450339724") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 1000))

    #### Mapped version
    ##
    ## 19 * 12 + 1000 * 2 = 2228
    ##
    ## partitionsDesign(19 + 0:1000 * 2, 12, target = 2228,
    ##                  freqs = c(13, rep(1, 1000)))[c("num_partitions",
    ##                                                 "mapped_target",
    ##                                                 "partition_type")]
    ## $num_partitions
    ## Big Integer ('bigz') :
    ## [1] 39233351450439724
    ##
    ## $mapped_target
    ## [1] 1000
    ##
    ## $partition_type
    ## [1] "DstctMultiZero"
    bench <- partitionsGeneral(19 + 0:1000 * 2, 12, target = 2228,
                               lower = "39233351450339724",
                               freqs = c(13, rep(1, 1000)))
    expect_identical(partitionsGeneral(19 + 0:1000 * 2, 12, target = 2228,
                                       freqs = c(13, rep(1, 1000)),
                                       lower = "39233351450339724",
                                       nThreads = 2), bench)
    expect_identical(partitionsRank(bench[1, ], v = 19 + 0:1000 * 2,
                                    target = 2228,
                                    freqs = c(13, rep(1, 1000))),
                     gmp::as.bigz("39233351450339724"))
    expect_equal(gmp::sub.bigz("39233351450439724",
                               "39233351450339724") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 2228))

    ########************* These take a few seconds ***************#######
    ## partitionsDesign(500, 10, target = 1380)[c("num_partitions",
    ##                                            "mapped_target",
    ##                                            "partition_type")]
    ## $num_partitions
    ## Big Integer ('bigz') :
    ## [1] 9605186196368891
    ##
    ## $mapped_target
    ## [1] 1380
    ##
    ## $partition_type
    ## [1] "DistCapped"
    bench <- partitionsGeneral(500, 10, lower = "9605186196218891",
                               target = 1380)
    expect_identical(partitionsGeneral(500, 10, lower = "9605186196218891",
                                       target = 1380, nThreads = 2), bench)
    expect_identical(partitionsRank(bench[1, ], v = 500, target = 1380),
                     gmp::as.bigz("9605186196218891"))
    expect_equal(gmp::sub.bigz("9605186196368891",
                               "9605186196218891") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 1380))

    ## partitionsDesign(0:500, 10, target = 1380,
    ##                  freqs = c(3, rep(1, 500)))[c("num_partitions",
    ##                                               "mapped_target",
    ##                                               "partition_type")]
    ## $num_partitions
    ## Big Integer ('bigz') :
    ## [1] 10236925075643716
    ##
    ## $mapped_target
    ## [1] 1380
    ##
    ## $partition_type
    ## [1] "DstctCappedMZ"
    bench <- partitionsGeneral(0:500, 10, freqs = c(3, rep(1, 500)),
                               lower = "10236925075443716", target = 1380)
    expect_identical(partitionsGeneral(0:500, 10, freqs = c(3, rep(1, 500)),
                                       lower = "10236925075443716", target = 1380,
                                       nThreads = 2), bench)
    expect_identical(partitionsRank(bench[1, ], v = 0:500,
                                    freqs = c(3, rep(1, 500)),
                                    target = 1380),
                     gmp::as.bigz("10236925075443716"))
    expect_equal(gmp::sub.bigz("10236925075643716",
                               "10236925075443716") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 1380))
})

test_that("partition/compositionsGeneral and Repetition Parallel Lower GMP", {

    #######################################################################
    ## See commentary above
    ########********************** Lower Only *******************##########
    ## partitionsDesign(0:300, repetition = T)[c("num_partitions",
    ##                                           "mapped_target",
    ##                                           "partition_type")]
    ## $num_partitions
    ## Big Integer ('bigz') :
    ## [1] 9253082936723602
    ##
    ## $mapped_target
    ## [1] 300
    ##
    ## $partition_type
    ## [1] "RepStdAll"
    bench <- partitionsGeneral(0:300, repetition = TRUE,
                               lower = "9253082936523602")
    expect_identical(partitionsGeneral(0:300, repetition = TRUE,
                                       lower = "9253082936523602",
                                       nThreads = 2), bench)
    expect_identical(partitionsRank(bench[1, ], v = 0:300,
                                    repetition = TRUE),
                     gmp::as.bigz("9253082936523602"))
    expect_equal(gmp::sub.bigz("9253082936723602",
                               "9253082936523602") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 300))

    #### Mapped version
    ##
    ## 3 * 300 + 17 * 300 = 6000
    ##
    ## partitionsDesign(3L + (0:300) * 17L, 300, repetition = TRUE,
    ##                  target = 6000L)[c("num_partitions",
    ##                                    "mapped_target",
    ##                                    "partition_type")]
    ## $num_partitions
    ## Big Integer ('bigz') :
    ## [1] 9253082936723602
    ##
    ## $mapped_target
    ## [1] 600
    ##
    ## $partition_type
    ## [1] "RepNoZero"
    bench <- partitionsGeneral(3L + (0:300) * 17L, 300, repetition = TRUE,
                               lower = "9253082936523602", target = 6000L)
    expect_identical(partitionsGeneral(3L + (0:300) * 17L, 300, repetition = TRUE,
                                       lower = "9253082936523602", target = 6000L,
                                       nThreads = 2), bench)
    expect_identical(partitionsRank(bench[1, ], v = 3L + (0:300) * 17L,
                                    target = 6000L, repetition = TRUE),
                     gmp::as.bigz("9253082936523602"))
    expect_equal(gmp::sub.bigz("9253082936723602",
                               "9253082936523602") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 6000L))

    ########********************** Lower Only *******************##########
    ## compositionsDesign(0:60, repetition = T)[c("num_partitions",
    ##                                            "mapped_target",
    ##                                            "partition_type")]
    ## $num_partitions
    ## Big Integer ('bigz') :
    ## [1] 576460752303423488
    ##
    ## $mapped_target
    ## [1] 60
    ##
    ## $partition_type
    ## [1] "RepStdAll"
    bench <- compositionsGeneral(0:60, repetition = TRUE,
                                 lower = "576460752303223488")
    expect_identical(compositionsGeneral(0:60, repetition = TRUE,
                                         lower = "576460752303223488",
                                         nThreads = 2), bench)
    expect_identical(compositionsRank(bench[1:5, ], v = 0:60,
                                      repetition = TRUE),
                     gmp::as.bigz("576460752303223488") + 0:4)
    expect_equal(gmp::sub.bigz("576460752303423488",
                               "576460752303223488") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 60))

    #### Mapped version
    ##
    ## compositionsDesign((0:60) * 19, repetition = T)[c("num_partitions",
    ##                                                   "mapped_target",
    ##                                                   "partition_type")]
    ## $num_partitions
    ## Big Integer ('bigz') :
    ## [1] 576460752303423488
    ##
    ## $mapped_target
    ## [1] 60
    ##
    ## $partition_type
    ## [1] "RepShort"
    bench <- compositionsGeneral((0:60) * 19, repetition = TRUE,
                                 lower = "576460752303223488")
    expect_identical(compositionsGeneral((0:60) * 19, repetition = TRUE,
                                         lower = "576460752303223488",
                                         nThreads = 2), bench)
    expect_identical(compositionsRank(bench[1:5, ], v = (0:60) * 19,
                                      repetition = TRUE),
                     gmp::as.bigz("576460752303223488") + 0:4)
    expect_equal(gmp::sub.bigz("576460752303423488",
                               "576460752303223488") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 60 * 19))

    ## Weak = TRUE
    ##
    ## compositionsDesign(0:60, repetition = T, weak = T)[c("num_partitions",
    ##                                                      "mapped_target",
    ##                                                     "partition_type")]
    ## $num_partitions
    ## Big Integer ('bigz') :
    ## [1] 48307454420181661301946569760686328
    ##
    ## $mapped_target
    ## [1] 120
    ##
    ## $partition_type
    ## [1] "RepNoZero"
    bench <- compositionsGeneral(0:60, repetition = TRUE, weak = TRUE,
                                 lower = "48307454420181661301946569760486328")
    expect_identical(compositionsGeneral(0:60, repetition = TRUE, weak = TRUE,
                                         lower = "48307454420181661301946569760486328",
                                         nThreads = 2), bench)
    expect_identical(compositionsRank(bench[1:5, ], v = 0:60,
                                      repetition = TRUE, weak = TRUE),
                     gmp::as.bigz("48307454420181661301946569760486328") + 0:4)
    expect_equal(gmp::sub.bigz("48307454420181661301946569760686328",
                               "48307454420181661301946569760486328") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 60))

    ## Mapped Version; Weak = TRUE
    ##
    ## compositionsDesign((0:60) * 19, 60, repetition = T,
    ##                    weak = T)[c("num_partitions", "mapped_target",
    ##                                "partition_type")]
    ## $num_partitions
    ## Big Integer ('bigz') :
    ## [1] 48307454420181661301946569760686328
    ##
    ## $mapped_target
    ## [1] 120
    ##
    ## $partition_type
    ## [1] "RepNoZero"
    bench <- compositionsGeneral((0:60) * 19, 60, repetition = TRUE, weak = TRUE,
                                 lower = "48307454420181661301946569760486328")
    expect_identical(compositionsGeneral((0:60) * 19, 60, repetition = TRUE, weak = TRUE,
                                         lower = "48307454420181661301946569760486328",
                                         nThreads = 2), bench)
    expect_identical(compositionsRank(bench[1:5, ], v = (0:60) * 19,
                                      repetition = TRUE, weak = TRUE),
                     gmp::as.bigz("48307454420181661301946569760486328") + 0:4)
    expect_equal(gmp::sub.bigz("48307454420181661301946569760686328",
                               "48307454420181661301946569760486328") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 60 * 19))

    ## partitionsDesign(6000, 10, TRUE)[c("num_partitions",
    ##                                    "mapped_target",
    ##                                    "partition_type")]
    ## $num_partitions
    ## Big Integer ('bigz') :
    ## [1] 7856063261819197461639
    ##
    ## $mapped_target
    ## [1] 6000
    ##
    ## $partition_type
    ## [1] "RepNoZero"
    bench <- partitionsGeneral(6000, 10, TRUE, lower = "7856063261819197261639")
    expect_identical(partitionsGeneral(6000, 10, TRUE,
                                       lower = "7856063261819197261639",
                                       nThreads = 2), bench)
    expect_identical(partitionsRank(bench[1, ], v = 6000,
                                    repetition = TRUE),
                     gmp::as.bigz("7856063261819197261639"))
    expect_equal(gmp::sub.bigz("7856063261819197461639",
                               "7856063261819197261639") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 6000L))

    #### Mapped version
    ##
    ## partitionsDesign(2e9 * 1:6000, 10, TRUE)[c("num_partitions",
    ##                                            "mapped_target",
    ##                                            "partition_type")]
    ## $num_partitions
    ## Big Integer ('bigz') :
    ## [1] 7856063261819197461639
    ##
    ## $mapped_target
    ## [1] 6000
    ##
    ## $partition_type
    ## [1] "RepNoZero"
    bench <- partitionsGeneral(2e9 * 1:6000, 10, TRUE,
                               lower = "7856063261819197261639")
    expect_identical(partitionsGeneral(2e9 * 1:6000, 10, TRUE,
                                       lower = "7856063261819197261639",
                                       nThreads = 2), bench)
    expect_identical(partitionsRank(bench[1, ], v = 2e9 * 1:6000,
                                    repetition = TRUE),
                     gmp::as.bigz("7856063261819197261639"))
    expect_equal(gmp::sub.bigz("7856063261819197461639",
                               "7856063261819197261639") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 1.2e+13))

    ## compositionsDesign(1000, 10, TRUE)[c("num_partitions",
    ##                                      "mapped_target",
    ##                                      "partition_type")]
    ## $num_partitions
    ## Big Integer ('bigz') :
    ## [1] 2634095604619702128324
    ##
    ## $mapped_target
    ## [1] 1000
    ##
    ## $partition_type
    ## [1] "RepNoZero"
    bench <- compositionsGeneral(1000, 10, TRUE, lower = "2634095604619701928324")
    expect_identical(compositionsGeneral(1000, 10, TRUE,
                                       lower = "2634095604619701928324",
                                       nThreads = 2), bench)
    expect_identical(compositionsRank(bench[1:11, ], v = 1000,
                                    repetition = TRUE),
                     gmp::as.bigz("2634095604619701928324") + 0:10)
    expect_equal(gmp::sub.bigz("2634095604619702128324",
                               "2634095604619701928324") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 1000L))

    #### Mapped version
    ##
    ## compositionsDesign(23 + (1:1000) * 123e8, 10, TRUE
    ##                    target = 12300000000230)[c("num_partitions",
    ##                                               "mapped_target",
    ##                                               "partition_type")]
    ## $num_partitions
    ## Big Integer ('bigz') :
    ## [1] 2634095604619702128324
    ##
    ## $mapped_target
    ## [1] 1000
    ##
    ## $partition_type
    ## [1] "RepNoZero"
    bench <- compositionsGeneral(23 + (1:1000) * 123e8, 10, TRUE,
                                 lower = "2634095604619701928324",
                                 target = 12300000000230)
    expect_identical(compositionsGeneral(23 + (1:1000) * 123e8, 10, TRUE,
                                         lower = "2634095604619701928324",
                                         target = 12300000000230,
                                         nThreads = 2), bench)
    expect_identical(compositionsRank(bench[1:11, ], v = 23 + (1:1000) * 123e8,
                                      target = 12300000000230,
                                      repetition = TRUE),
                     gmp::as.bigz("2634095604619701928324") + 0:10)
    expect_equal(gmp::sub.bigz("2634095604619702128324",
                               "2634095604619701928324") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 12300000000230))

    ## partitionsDesign(0:6000, 10, TRUE)[c("num_partitions",
    ##                                      "mapped_target",
    ##                                      "partition_type")]
    ## $num_partitions
    ## Big Integer ('bigz') :
    ## [1] 7974346430545135003177
    ##
    ## $mapped_target
    ## [1] 6010
    ##
    ## $partition_type
    ## [1] "RepShort"
    bench <- partitionsGeneral(0:6000, 10, TRUE, lower = "7974346430545134803177")
    expect_identical(partitionsGeneral(0:6000, 10, TRUE,
                                       lower = "7974346430545134803177",
                                       nThreads = 2), bench)
    expect_identical(partitionsRank(bench[1, ], v = 0:6000,
                                    repetition = TRUE),
                     gmp::as.bigz("7974346430545134803177"))
    expect_equal(gmp::sub.bigz("7974346430545135003177",
                               "7974346430545134803177") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 6000L))

    #### Mapped version
    ##
    ## partitionsDesign(2e9 * 0:6000, 10, TRUE)[c("num_partitions",
    ##                                            "mapped_target",
    ##                                            "partition_type")]
    ## $num_partitions
    ## Big Integer ('bigz') :
    ## [1] 7974346430545135003177
    ##
    ## $mapped_target
    ## [1] 6010
    ##
    ## $partition_type
    ## [1] "RepNoZero"
    bench <- partitionsGeneral(2e9 * 0:6000, 10, TRUE,
                               lower = "7974346430545134803177")
    expect_identical(partitionsGeneral(2e9 * 0:6000, 10, TRUE,
                                       lower = "7974346430545134803177",
                                       nThreads = 2), bench)
    expect_identical(partitionsRank(bench[1, ], v = 2e9 * 0:6000,
                                    repetition = TRUE),
                     gmp::as.bigz("7974346430545134803177"))
    expect_equal(gmp::sub.bigz("7974346430545135003177",
                               "7974346430545134803177") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 1.2e+13))

    ## partitionsDesign(200, 12, TRUE,
    ##                  target = 1000)[c("num_partitions",
    ##                                   "mapped_target",
    ##                                   "partition_type")]
    ## $num_partitions
    ## Big Integer ('bigz') :
    ## [1] 14001484381527012
    ##
    ## $mapped_target
    ## [1] 1000
    ##
    ## $partition_type
    ## [1] "RepCapped"
    bench <- partitionsGeneral(200, 12, TRUE, target = 1000,
                               lower = "14001484381327012")
    expect_identical(partitionsGeneral(200, 12, TRUE, lower = "14001484381327012",
                                       target = 1000, nThreads = 2), bench)
    expect_identical(partitionsRank(bench[1, ], v = 200,
                                    target = 1000,
                                    repetition = TRUE),
                     gmp::as.bigz("14001484381327012"))
    expect_equal(gmp::sub.bigz("14001484381527012",
                               "14001484381327012") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 1000))
})
