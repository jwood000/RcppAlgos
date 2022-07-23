context("testing partitionGeneral Parallel")

test_that("partitionGeneral Distinct Parallel", {

    ## N.B. Parallel has no effect when number of results is less than 40000
    ## partitionsCount(20)
    ## [1] 7
    expect_identical(partitionsGeneral(20, nThreads = 2),
                     partitionsGeneral(20))

    ## For both of the usages below, only 2 threads will be spawned
    ## partitionsCount(0:100, repetition = T)
    ## [1] 190569292
    expect_identical(partitionsGeneral(0:100, repetition = TRUE,
                                       nThreads = 3, upper = 50000),
                     partitionsGeneral(0:100, repetition = TRUE,
                                       nThreads = 8, upper = 50000))

    ######****************** All Results Distinct **************#########
    #### Distinct; Length determined internally; No zero;
    ##
    ## sum(1:41)
    ## [1] 861
    ## sum(1:42)
    ## [1] 903
    ##
    ## partitionsDesign(902)[c("num_partitions", "partition_type")]
    ## $num_partitions
    ## [1] 44583
    ##
    ## $partition_type
    ## [1] "DstctNoZero"
    myParts = partitionsGeneral(902, nThreads = 2)
    expect_identical(myParts, partitionsGeneral(902))
    expect_identical(10000L, partitionsRank(myParts[10000, ], v = 902))

    #### Mapped version
    ##
    ## 902 * 17 + 3 * 41 = 15457
    ##
    ## partitionsDesign(3 + (1:902) * 17, 41,
    ##                  target = 15457)[c("num_partitions",
    ##                                    "mapped_target", "partition_type")]
    ## $num_partitions
    ## [1] 44583
    ##
    ## $mapped_target
    ## [1] 902
    ##
    ## $partition_type
    ## [1] "DstctNoZero"
    myParts = partitionsGeneral(3 + (1:902) * 17, 41,
                               target = 15457, nThreads = 2)
    expect_identical(myParts,
                     partitionsGeneral(3 + (1:902) * 17, 41, target = 15457))
    expect_identical(partitionsRank(myParts[10000, ],
                                    v = 3 + (1:902) * 17,
                                    target = 15457), 10000L)

    #### Distinct; Specific Length; No zero
    ##
    ## partitionsDesign(105, 10)[c("num_partitions", "partition_type")]
    ## $num_partitions
    ## [1] 62740
    ##
    ## $partition_type
    ## [1] "DstctNoZero"
    myParts = partitionsGeneral(105, 10, nThreads = 2)
    expect_identical(myParts, partitionsGeneral(105, 10))
    expect_identical(partitionsRank(myParts[1000, ],
                                    myParts[10000, ], v = 105),
                     list(1000L, 10000L))

    #### Mapped version
    ##
    ## 105 * 3 + 6 * 10 = 375
    ##
    ## partitionsDesign(6 + (1:105) * 3, 10,
    ##                  target = 375)[c("num_partitions",
    ##                                  "mapped_target", "partition_type")]
    ## $num_partitions
    ## [1] 62740
    ##
    ## $mapped_target
    ## [1] 105
    ##
    ## $partition_type
    ## [1] "DstctNoZero"
    myParts = partitionsGeneral(6 + (1:105) * 3, 10, target = 375)
    expect_identical(partitionsGeneral(6 + (1:105) * 3, 10,
                                       target = 375, nThreads = 2),
                     myParts)
    expect_identical(partitionsRank(myParts[seq(1, 60000, 10000), ],
                                    v = 6 + (1:105) * 3, target = 375),
                     seq(1L, 60000L, 10000L))

    #### Distinct; Specific Length; One zero
    ##
    ## partitionsDesign(0:95, 10)[c("num_partitions", "partition_type")]
    ## $num_partitions
    ## [1] 62740
    ##
    ## $partition_type
    ## [1] "DstctOneZero"
    myParts = partitionsGeneral(0:95, 10, nThreads = 2)
    expect_identical(myParts, partitionsGeneral(0:95, 10))
    expect_identical(partitionsRank(myParts[seq(1, 60000, 10000), ],
                                    v = 0:95), seq(1L, 60000L, 10000L))

    #### Mapped version. N.B. partition_type is different since we are in
    #### a mapped case.
    ##
    ## 95 * 7 = 665
    ##
    ## partitionsDesign((0:95) * 7, 10,
    ##                  target = 665)[c("num_partitions",
    ##                                  "mapped_target", "partition_type")]
    ## $num_partitions
    ## [1] 62740
    ##
    ## $mapped_target
    ## [1] 105
    ##
    ## $partition_type
    ## [1] "DstctNoZero"
    myParts = partitionsGeneral((0:95) * 7, 10, target = 665)
    expect_identical(partitionsGeneral((0:95) * 7, 10,
                                       target = 665, nThreads = 2), myParts)
    expect_identical(partitionsRank(myParts[1000, ], v = (0:95) * 7,
                                    target = 665), 1000L)

    #### Distinct; Specific Length; Multiple Zeros; Not enough to maximize
    ##
    ## partitionsDesign(0:77, 8,
    ##                  freqs = c(3, rep(1, 77)))[c("num_partitions",
    ##                                              "partition_type")]
    ## $num_partitions
    ## [1] 50349
    ##
    ## $partition_type
    ## [1] "DstctMultiZero"
    myParts = partitionsGeneral(0:77, 8, freqs = c(3, rep(1, 77)))
    expect_identical(partitionsGeneral(0:77, 8, freqs = c(3, rep(1, 77)),
                                       nThreads = 2), myParts)
    expect_identical(partitionsRank(myParts[50349, ], v = 0:77,
                                    freqs = c(3, rep(1, 77))), 50349L)

    #### Mapped version
    ##
    ## 77 * 3 + 15 * 8 = 351
    ##
    ## partitionsDesign(15 + 0:70 * 3, 8, target = 351,
    ##                  freqs = c(3, rep(1, 70)))[c("num_partitions",
    ##                                              "mapped_target",
    ##                                              "partition_type")]
    ## $num_partitions
    ## [1] 50349
    ##
    ## $mapped_target
    ## [1] 77
    ##
    ## $partition_type
    ## [1] "DstctMultiZero"
    myParts = partitionsGeneral(15 + 0:77 * 3, 8, target = 351,
                               freqs = c(3, rep(1, 77)), nThreads = 2)
    expect_identical(partitionsGeneral(15 + 0:77 * 3, 8, target = 351,
                                       freqs = c(3, rep(1, 77))), myParts)
    expect_identical(partitionsRank(myParts[c(1, 50349), ], v = 15 + 0:77 * 3,
                                    target = 351, freqs = c(3, rep(1, 77))),
                     c(1L, 50349L))

    #### Distinct; Specific Length; Multiple Zeros; Enough to maximize;
    #### Length is restrictive
    ##
    ## partitionsDesign(0:110, 5,
    ##                  freqs = c(8, rep(1, 110)))[c("num_partitions",
    ##                                               "partition_type")]
    ## $num_partitions
    ## [1] 47271
    ##
    ## $partition_type
    ## [1] "DstctMultiZero"
    myParts = partitionsGeneral(0:110, 5, freqs = c(8, rep(1, 110)))
    expect_identical(partitionsGeneral(0:110, 5, freqs = c(8, rep(1, 110)),
                                       nThreads = 2), myParts)
    expect_identical(partitionsRank(myParts[c(1, 20000, 47271), ],
                                    v = 0:110, freqs = c(8, rep(1, 110))),
                     c(1L, 20000L, 47271L))

    #### Mapped version
    ##
    ## 110 * 2 + 19 * 5 = 315
    ##
    ## partitionsDesign(19 + (0:110) * 2, 5, target = 315,
    ##                  freqs = c(8, rep(1, 110)))[c("num_partitions",
    ##                                               "mapped_target",
    ##                                               "partition_type")]
    ## $num_partitions
    ## [1] 47271
    ##
    ## $mapped_target
    ## [1] 110
    ##
    ## $partition_type
    ## [1] "DstctMultiZero"
    myParts = partitionsGeneral(19 + (0:110) * 2, 5, target = 315,
                                freqs = c(8, rep(1, 110)), nThreads = 2)
    expect_identical(partitionsGeneral(19 + (0:110) * 2, 5, target = 315,
                                       freqs = c(8, rep(1, 110))), myParts)
    expect_identical(partitionsRank(myParts[seq(15757, 47271, 15757), ],
                                    v = 19 + (0:110) * 2, target = 315,
                                    freqs = c(8, rep(1, 110))),
                     seq(15757L, 47271L, 15757L))

    #### Distinct; Length determined internally; Multiple Zeros;
    #### Enough to maximize; N.B. There is no mapped version of this case
    ##
    ## partitionsDesign(0:80,
    ##                  freqs = c(80, rep(1, 80)))[c("num_partitions",
    ##                                                 "partition_type")]
    ## $num_partitions
    ## [1] 77312
    ##
    ## $partition_type
    ## [1] "DstctStdAll"
    myParts = partitionsGeneral(0:80, freqs = c(80, rep(1, 80)))
    expect_identical(partitionsGeneral(0:80, freqs = c(80, rep(1, 80)),
                                       nThreads = 2), myParts)
    expect_identical(partitionsRank(myParts[c(1, 12345, 77312), ],
                                    v = 0:80, freqs = c(80, rep(1, 80))),
                     c(1L, 12345L, 77312L))

    #### Distinct; Specific Length; No Zeros; Specific Target; N.B.
    #### Technically this is already a mapped case, however we will
    #### still provide an additional example
    ##
    ## partitionsDesign(40, 10, target = 115)[c("num_partitions",
    ##                                          "partition_type")]
    ## $num_partitions
    ## [1] 180436
    ##
    ## $partition_type
    ## [1] "DistCapped"
    myParts = partitionsGeneral(40, 10, target = 115, nThreads = 2)
    expect_identical(partitionsGeneral(40, 10, target = 115), myParts)
    expect_identical(partitionsRank(myParts[seq(45109, 180436, 45109), ],
                                    v = 40, target = 115),
                     seq(45109L, 180436L, 45109L))

    #### Mapped version
    ##
    ## 1001 * 10 + 115 * 107 = 22315
    ##
    #3 partitionsDesign(1001 + (1:40) * 107, 10,
    ##                  target = 22315)[c("num_partitions",
    ##                                    "mapped_target",
    ##                                    "partition_type")]
    ## $num_partitions
    ## [1] 180436
    ##
    ## $mapped_target
    ## [1] 115
    ##
    ## $partition_type
    ## [1] "DistCapped"
    myParts = partitionsGeneral(1001 + (1:40) * 107, 10,
                                target = 22315, nThreads = 2)
    expect_identical(partitionsGeneral(1001 + (1:40) * 107, 10,
                                       target = 22315), myParts)
    expect_identical(partitionsRank(a = myParts[123456, ], b = myParts[54321, ],
                                    v = 1001 + (1:40) * 107, target = 22315),
                     list(a = 123456L, b = 54321L))

    #### Distinct; Specific Length; Multi Zeros; Specific Target; N.B.
    #### Technically this is already a mapped case, however we will
    #### still provide an additional example
    ##
    ## partitionsDesign(0:30, 9, target = 115,
    ##                  freqs = c(3, rep(1, 30)))[c("num_partitions",
    ##                                              "partition_type")]
    ## $num_partitions
    ## [1] 284705
    ##
    ## $partition_type
    ## [1] "DstctCappedMZ"
    myParts = partitionsGeneral(0:30, 9, freqs = c(3, rep(1, 30)),
                                target = 115, nThreads = 2)
    expect_identical(partitionsGeneral(0:30, 9, freqs = c(3, rep(1, 30)),
                                       target = 115), myParts)
    expect_identical(partitionsRank(myParts[seq(56941, 284705, 56941), ],
                                    v = 0:30, freqs = c(3, rep(1, 30)),
                                    target = 115), seq(56941L, 284705L, 56941L))

    #### Mapped version
    ##
    ## 1001 * 9 + 115 * 107 = 22315
    ##
    ## partitionsDesign(1001 + (0:30) * 107, 9, target = 21314,
    ##                  freqs = c(3, rep(1, 30)))[c("num_partitions",
    ##                                              "partition_type")]
    ## $num_partitions
    ## [1] 284705
    ##
    ## $partition_type
    ## [1] "DstctCappedMZ"
    myParts = partitionsGeneral(1001 + (0:30) * 107, 9,
                                freqs = c(3, rep(1, 30)),
                                target = 21314, nThreads = 2)
    expect_identical(partitionsGeneral(1001 + (0:30) * 107, 9,
                                       freqs = c(3, rep(1, 30)),
                                       target = 21314), myParts)
    expect_identical(partitionsRank(myParts[c(1, 123456, 284705), ],
                                    v = 1001 + (0:30) * 107,
                                    freqs = c(3, rep(1, 30)),
                                    target = 21314),
                     c(1L, 123456L, 284705L))
})

test_that("partitionGeneral Distinct Parallel Lower", {

    #######################################################################
    ## See commentary above
    ########********************** Lower Only *******************##########
    expect_identical(partitionsGeneral(902, lower = 17, nThreads = 2),
                     partitionsGeneral(902, lower = 17))

    expect_identical(partitionsGeneral(3 + (1:902) * 17, 41, lower = 100,
                                       target = 15457, nThreads = 2),
                     partitionsGeneral(3 + (1:902) * 17, 41, lower = 100,
                                       target = 15457))

    expect_identical(partitionsGeneral(105, 10, lower = 501, nThreads = 2),
                     partitionsGeneral(105, 10, lower = 501))

    expect_identical(partitionsGeneral(6 + (1:105) * 3, 10, lower = 13,
                                       target = 375, nThreads = 2),
                     partitionsGeneral(6 + (1:105) * 3, 10, lower = 13,
                                       target = 375))

    expect_identical(partitionsGeneral(0:95, 10, lower = 313, nThreads = 2),
                     partitionsGeneral(0:95, 10, lower = 313))

    expect_identical(partitionsGeneral((0:95) * 7, 10, lower = 123,
                                       target = 665, nThreads = 2),
                     partitionsGeneral((0:95) * 7, 10, lower = 123,
                                       target = 665))

    expect_identical(partitionsGeneral(0:77, 8, lower = 4123,
                                       freqs = c(3, rep(1, 77)),
                                       nThreads = 2),
                     partitionsGeneral(0:77, 8, lower = 4123,
                                       freqs = c(3, rep(1, 77))))

    expect_identical(partitionsGeneral(15 + 0:77 * 3, 8,
                                       lower = 321, target = 351,
                                       freqs = c(3, rep(1, 77)), nThreads = 2),
                     partitionsGeneral(15 + 0:77 * 3, 8,
                                       lower = 321, target = 351,
                                       freqs = c(3, rep(1, 77))))

    expect_identical(partitionsGeneral(0:110, 5, freqs = c(8, rep(1, 110)),
                                       lower = 5000, nThreads = 2),
                     partitionsGeneral(0:110, 5, lower = 5000,
                                       freqs = c(8, rep(1, 110))))

    expect_identical(partitionsGeneral(19 + (0:110) * 2, 5,
                                       lower = 4321, target = 315,
                                       freqs = c(8, rep(1, 110)),
                                       nThreads = 2),
                     partitionsGeneral(19 + (0:110) * 2, 5,
                                       lower = 4321, target = 315,
                                       freqs = c(8, rep(1, 110))))

    expect_identical(partitionsGeneral(0:80, freqs = c(80, rep(1, 80)),
                                       lower = 11111, nThreads = 2),
                     partitionsGeneral(0:80, lower = 11111,
                                       freqs = c(80, rep(1, 80))))

    expect_identical(partitionsGeneral(40, 10, lower = 12345,
                                       target = 115, nThreads = 2),
                     partitionsGeneral(40, 10, lower = 12345, target = 115))

    expect_identical(partitionsGeneral(1001 + (1:40) * 107, 10, lower = 12321,
                                       target = 22315, nThreads = 2),
                     partitionsGeneral(1001 + (1:40) * 107, 10, lower = 12321,
                                       target = 22315))

    expect_identical(partitionsGeneral(0:30, 9, freqs = c(3, rep(1, 30)),
                                       target = 115, lower = 2222,
                                       nThreads = 2),
                     partitionsGeneral(0:30, 9, freqs = c(3, rep(1, 30)),
                                       lower = 2222, target = 115))

    expect_identical(partitionsGeneral(1001 + (0:30) * 107, 9,
                                       freqs = c(3, rep(1, 30)),
                                       lower = 33333,
                                       target = 21314, nThreads = 2),
                     partitionsGeneral(1001 + (0:30) * 107, 9,
                                       lower = 33333,
                                       freqs = c(3, rep(1, 30)),
                                       target = 21314))
})

test_that("partitionGeneral Repetition Parallel", {

    ## N.B. Parallel has no effect when number of results is less than 40000
    ## partitionsCount(0:20, repetition = TRUE)
    ## [1] 627
    myParts = partitionsGeneral(0:20, repetition = TRUE, nThreads = 2)
    expect_identical(partitionsGeneral(0:20, repetition = TRUE), myParts)
    expect_identical(partitionsRank(myParts[c(1L, 321L, 627L), ],
                                    v = 0:20, repetition = TRUE),
                     c(1L, 321L, 627L))
    ## compositionsCount(0:10, repetition = TRUE)
    ## [1] 512
    myComps = compositionsGeneral(0:10, repetition = TRUE, nThreads = 2)
    expect_identical(compositionsGeneral(0:10, repetition = TRUE), myComps)
    expect_identical(compositionsRank(myComps[c(1L, 321L, 512L), ],
                                      v = 0:10, repetition = TRUE),
                     c(1L, 321L, 512L))

    ## For both of the usages below, only 2 threads will be spawned
    ## partitionsCount(0:100, repetition = T)
    ## [1] 190569292
    expect_identical(partitionsGeneral(0:100, repetition = TRUE,
                                       nThreads = 3, upper = 50000),
                     partitionsGeneral(0:100, repetition = TRUE,
                                       nThreads = 8, upper = 50000))
    expect_identical(compositionsGeneral(0:30, repetition = TRUE,
                                         nThreads = 3, upper = 50000),
                     compositionsGeneral(0:30, repetition = TRUE,
                                         nThreads = 8, upper = 50000))

    ######****************** All Results Repetition **************#########
    #### Repetition; Length determined internally; Multiple Zero;
    ##
    ## partitionsDesign(0:45, repetition = TRUE)[c("num_partitions",
    ##                                             "partition_type")]
    ## $num_partitions
    ## [1] 89134
    ##
    ## $partition_type
    ## [1] "RepStdAll"
    expect_identical(partitionsGeneral(0:45, repetition = TRUE,
                                       nThreads = 2),
                     partitionsGeneral(0:45, repetition = TRUE))

    #### Mapped version
    ##
    ## 45 * 3 + 45 * 17 = 900
    ##
    ## partitionsDesign(3L + (0:45) * 17L, 46, repetition = TRUE,
    ##                  target = 900)[c("num_partitions",
    ##                                    "mapped_target", "partition_type")]
    ## $num_partitions
    ## [1] 89134
    ##
    ## $mapped_target
    ## [1] 90
    ##
    ## $partition_type
    ## [1] "RepNoZero"
    myParts = partitionsGeneral(3L + (0:45) * 17L, 45,
                                repetition = TRUE,
                                target = 900L, nThreads = 2)
    expect_identical(myParts,
                     partitionsGeneral(3L + (0:45) * 17L, 45,
                                       repetition = TRUE, target = 900L))
    expect_identical(partitionsRank(myParts[c(1L, 12345L, 89134L), ],
                                    v = 3L + (0:45) * 17L, target = 900L,
                                    repetition = TRUE),
                     c(1L, 12345L, 89134L))

    #### Repetition; Length determined internally; Multiple Zero; Composition
    ##
    ## compositionsDesign(0:17, repetition = TRUE)[c("num_partitions",
    ##                                               "partition_type")]
    ## $num_partitions
    ## [1] 65536
    ##
    ## $partition_type
    ## [1] "RepStdAll"
    expect_identical(compositionsGeneral(0:17, repetition = TRUE,
                                       nThreads = 2),
                     compositionsGeneral(0:17, repetition = TRUE))

    #### Mapped version
    ##
    ## 17 * 17 = 289  // No need to supply below. See GetTarget.R
    ##
    ## compositionsDesign((0:17) * 17, repetition = TRUE)[c(
    ##       "num_partitions", "mapped_target", "partition_type"
    ## )]
    ## $num_partitions
    ## [1] 65536
    ##
    ## $mapped_target
    ## [1] 17
    ##
    ## $partition_type
    ## [1] "RepShort"
    expect_identical(compositionsGeneral((0:17) * 17L, repetition = TRUE,
                                         nThreads = 2),
                     compositionsGeneral((0:17) * 17L, repetition = TRUE))

    #### Repetition; Specific Length; No zero
    ##
    ## partitionsDesign(60, 10, TRUE)[c("num_partitions", "partition_type")]
    ## $num_partitions
    ## [1] 62740
    ##
    ## $partition_type
    ## [1] "RepNoZero"
    myParts = partitionsGeneral(60, 10, TRUE, nThreads = 2)
    expect_identical(partitionsGeneral(60, 10, TRUE), myParts)
    expect_identical(partitionsRank(myParts[c(1L, 12345L, 62740L), ],
                                            v = 60, repetition = TRUE),
                     c(1L, 12345L, 62740L))

    #### Mapped version
    ##
    ## 60 * 3 + 6 * 10 = 240
    ##
    ## partitionsDesign(6 + (1:60) * 3, 10, repetition = TRUE,
    ##                  target = 240)[c("num_partitions",
    ##                                  "mapped_target", "partition_type")]
    ## $num_partitions
    ## [1] 62740
    ##
    ## $mapped_target
    ## [1] 60
    ##
    ## $partition_type
    ## [1] "RepNoZero"
    myParts = partitionsGeneral(6 + (1:60) * 3, 10, TRUE,
                                target = 240, nThreads = 2)
    expect_identical(partitionsGeneral(6 + (1:60) * 3, 10,
                                       TRUE, target = 240), myParts)
    expect_identical(partitionsRank(myParts[c(1L, 12345L, 62740L), ],
                                    v = 6 + (1:60) * 3,
                                    target = 240, repetition = TRUE),
                     c(1L, 12345L, 62740L))

    #### Repetition; Specific Length; No zero; Composition
    ##
    ## compositionsDesign(28, 6, TRUE)[c("num_partitions", "partition_type")]
    ## $num_partitions
    ## [1] 80730
    ##
    ## $partition_type
    ## [1] "RepNoZero"
    myComps = compositionsGeneral(28, 6, TRUE, nThreads = 2)
    expect_identical(compositionsGeneral(28, 6, TRUE), myComps)
    expect_identical(compositionsRank(myComps[c(1L, 12345L, 80730L), ],
                                    v = 28, repetition = TRUE),
                     c(1L, 12345L, 80730L))

    #### Mapped version
    ##
    ## 28 * 3 + 6 * 6 = 120
    ##
    ## compositionsDesign(6 + (1:28) * 3, 6, repetition = TRUE,
    ##                    target = 120)[c("num_partitions",
    ##                                    "mapped_target", "partition_type")]
    ## $num_partitions
    ## [1] 80730
    ##
    ## $mapped_target
    ## [1] 28
    ##
    ## $partition_type
    ## [1] "RepNoZero"
    myComps = compositionsGeneral(6 + (1:28) * 3, 10, TRUE,
                                  target = 120, nThreads = 2)
    expect_identical(compositionsGeneral(6 + (1:28) * 3, 10,
                                         TRUE, target = 120), myComps)
    expect_identical(compositionsRank(myComps[c(1L, 12345L, 80730L), ],
                                      v = 6 + (1:28) * 3,
                                      target = 120, repetition = TRUE),
                     c(1L, 12345L, 80730L))

    #### Repetition; Specific Length; Zero included
    ##
    ## partitionsDesign(0:60, 10, TRUE)[c("num_partitions", "partition_type")]
    ## $num_partitions
    ## [1] 195491
    ##
    ## $partition_type
    ## [1] "RepShort"
    myParts = partitionsGeneral(0:60, 10, TRUE, nThreads = 2)
    expect_identical(myParts, partitionsGeneral(0:60, 10, TRUE))
    expect_identical(partitionsRank(myParts[c(1L, 123456L, 195491L), ],
                                    v = 0:60, repetition = TRUE),
                     c(1L, 123456L, 195491L))

    #### Mapped version
    ##
    ## 60 * 3 + 6 * 10 = 240
    ##
    ## partitionsDesign(6 + (0:60) * 3, 10, repetition = TRUE,
    ##                  target = 240)[c("num_partitions",
    ##                                  "mapped_target", "partition_type")]
    ## $num_partitions
    ## [1] 195491
    ##
    ## $mapped_target
    ## [1] 70
    ##
    ## $partition_type
    ## [1] "RepNoZero"
    myParts = partitionsGeneral(6 + (0:60) * 3, 10, TRUE,
                                target = 240, nThreads = 2)
    expect_identical(partitionsGeneral(6 + (0:60) * 3, 10,
                                       TRUE, target = 240), myParts)
    expect_identical(partitionsRank(myParts[c(1L, 123456L, 195491L), ],
                                    v = 6 + (0:60) * 3,
                                    target = 240, repetition = TRUE),
                     c(1L, 123456L, 195491L))

    #### Repetition; Specific Length; Zero included; Composition
    ##
    ## compositionsDesign(0:28, 6, TRUE)[c("num_partitions", "partition_type")]
    ## $num_partitions
    ## [1] 101584
    ##
    ## $partition_type
    ## [1] "RepShort"
    myComps = compositionsGeneral(0:28, 6, TRUE, nThreads = 2)
    expect_identical(myComps, compositionsGeneral(0:28, 6, TRUE))
    expect_identical(compositionsRank(myComps[c(1L, 12345L, 101584L), ],
                                      v = 0:28, repetition = TRUE),
                     c(1L, 12345L, 101584L))

    #### Mapped version
    ##
    ## 28 * 3 = 84  // No need to supply below. See GetTarget.R
    ##
    ## compositionsDesign((0:28) * 3, 6, repetition = TRUE,
    ##                    target = 84)[c("num_partitions",
    ##                                   "mapped_target", "partition_type")]
    ## $num_partitions
    ## [1] 101584
    ##
    ## $mapped_target
    ## [1] 28
    ##
    ## $partition_type
    ## [1] "RepShort"
    myComps = compositionsGeneral((0:28) * 3, 6, TRUE, nThreads = 2)
    expect_identical(myComps, compositionsGeneral((0:28) * 3, 6, TRUE))
    expect_identical(compositionsRank(myComps[c(1L, 12345L, 101584L), ],
                                      v = (0:28) * 3, repetition = TRUE),
                     c(1L, 12345L, 101584L))

    #### Repetition; Specific Length; No Zeros; Specific Target; N.B.
    #### Technically this is already a mapped case, however we will
    #### still provide an additional example
    ##
    ## partitionsDesign(18, 10, TRUE, target = 100)[c("num_partitions",
    ##                                                "partition_type")]
    ## $num_partitions
    ## [1] 161073
    ##
    ## $partition_type
    ## [1] "RepCapped"
    myParts = partitionsGeneral(18, 10, TRUE,
                                target = 100, nThreads = 2)
    expect_identical(myParts, partitionsGeneral(18, 10, TRUE, target = 100))
    expect_identical(partitionsRank(myParts[c(1, 10000, 161073), ],
                                    v = 18, target = 100,
                                    repetition = TRUE),
                     c(1L, 10000L, 161073L))

    #### Mapped version
    ##
    ## 1001 * 10 + 100 * 107 = 20710
    ##
    ## partitionsDesign(1001 + (1:18) * 107, 10, TRUE,
    ##                  target = 20710)[c("num_partitions",
    ##                                    "mapped_target",
    ##                                    "partition_type")]
    ## $num_partitions
    ## [1] 161073
    ##
    ## $mapped_target
    ## [1] 100
    ##
    ## $partition_type
    ## [1] "RepCapped"
    myParts = partitionsGeneral(1001 + (1:18) * 107, 10, TRUE,
                                target = 20710, nThreads = 2)
    expect_identical(partitionsGeneral(1001 + (1:18) * 107, 10, TRUE,
                                       target = 20710), myParts)
    expect_identical(partitionsRank(myParts[c(1, 10000, 161073), ],
                                    v = 1001 + (1:18) * 107,
                                    target = 20710, repetition = TRUE),
                     c(1L, 10000L, 161073L))
})

test_that("partitionGeneral Repetition Parallel Lower", {

    #######################################################################
    ## See commentary above
    ########********************** Lower Only *******************##########
    expect_identical(partitionsGeneral(0:45, repetition = TRUE,
                                       nThreads = 2, lower = 101),
                     partitionsGeneral(0:45, repetition = TRUE,
                                       lower = 101))

    expect_identical(partitionsGeneral(3L + (0:45) * 17L, 46,
                                       repetition = TRUE, lower = 3333,
                                       target = 903L, nThreads = 2),
                     partitionsGeneral(3L + (0:45) * 17L, 46, lower = 3333,
                                       repetition = TRUE, target = 903L))

    expect_identical(partitionsGeneral(60, 10, TRUE,
                                       nThreads = 2, lower = 22741),
                     partitionsGeneral(60, 10, TRUE, lower = 22741))

    expect_identical(partitionsGeneral(6 + (1:60) * 3, 10, TRUE,
                                       target = 240, lower = 12321,
                                       nThreads = 2),
                     partitionsGeneral(6 + (1:60) * 3, 10,
                                       TRUE, target = 240, lower = 12321))

    expect_identical(partitionsGeneral(0:60, 10, TRUE,
                                       lower = 43212, nThreads = 2),
                     partitionsGeneral(0:60, 10, TRUE, lower = 43212))

    expect_identical(partitionsGeneral(6 + (0:60) * 3, 10, TRUE,
                                       target = 240, nThreads = 2,
                                       lower = 54321),
                     partitionsGeneral(6 + (0:60) * 3, 10,
                                       TRUE, target = 240,
                                       lower = 54321))

    expect_identical(partitionsGeneral(18, 10, TRUE, lower = 2,
                                       target = 100, nThreads = 2),
                     partitionsGeneral(18, 10, TRUE, lower = 2,
                                       target = 100))

    expect_identical(partitionsGeneral(1001 + (1:18) * 107, 10, TRUE,
                                       target = 20710, lower = 2121,
                                       nThreads = 2),
                     partitionsGeneral(1001 + (1:18) * 107, 10, TRUE,
                                       target = 20710, lower = 2121))
})
