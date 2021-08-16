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
    expect_identical(partitionsGeneral(902, nThreads = 2),
                     partitionsGeneral(902))
    
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
    expect_identical(partitionsGeneral(3 + (1:902) * 17, 41,
                                       target = 15457, nThreads = 2),
                     partitionsGeneral(3 + (1:902) * 17, 41, target = 15457))
    
    #### Distinct; Specific Length; No zero
    ##
    ## partitionsDesign(105, 10)[c("num_partitions", "partition_type")]
    ## $num_partitions
    ## [1] 62740
    ## 
    ## $partition_type
    ## [1] "DstctNoZero"
    expect_identical(partitionsGeneral(105, 10, nThreads = 2),
                     partitionsGeneral(105, 10))
    
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
    expect_identical(partitionsGeneral(6 + (1:105) * 3, 10,
                                       target = 375, nThreads = 2),
                     partitionsGeneral(6 + (1:105) * 3, 10, target = 375))
    
    #### Distinct; Specific Length; One zero
    ##
    ## partitionsDesign(0:95, 10)[c("num_partitions", "partition_type")]
    ## $num_partitions
    ## [1] 62740
    ## 
    ## $partition_type
    ## [1] "DstctOneZero"
    expect_identical(partitionsGeneral(0:95, 10, nThreads = 2),
                     partitionsGeneral(0:95, 10))
    
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
    expect_identical(partitionsGeneral((0:95) * 7, 10,
                                       target = 665, nThreads = 2),
                     partitionsGeneral((0:95) * 7, 10, target = 665))
    
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
    expect_identical(partitionsGeneral(0:77, 8, freqs = c(3, rep(1, 77)),
                                       nThreads = 2), 
                     partitionsGeneral(0:77, 8, freqs = c(3, rep(1, 77))))
    
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
    expect_identical(partitionsGeneral(15 + 0:77 * 3, 8, target = 351,
                                       freqs = c(3, rep(1, 77)), nThreads = 2), 
                     partitionsGeneral(15 + 0:77 * 3, 8, target = 351,
                                       freqs = c(3, rep(1, 77))))
    
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
    expect_identical(partitionsGeneral(0:110, 5, freqs = c(8, rep(1, 110)),
                                       nThreads = 2), 
                     partitionsGeneral(0:110, 5, freqs = c(8, rep(1, 110))))
    
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
    expect_identical(partitionsGeneral(19 + (0:110) * 2, 5, target = 315,
                                       freqs = c(8, rep(1, 110)), nThreads = 2), 
                     partitionsGeneral(19 + (0:110) * 2, 5, target = 315,
                                       freqs = c(8, rep(1, 110))))
    
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
    expect_identical(partitionsGeneral(0:80, freqs = c(80, rep(1, 80)),
                                       nThreads = 2),
                     partitionsGeneral(0:80, freqs = c(80, rep(1, 80))))
    
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
    expect_identical(partitionsGeneral(40, 10, target = 115, nThreads = 2),
                     partitionsGeneral(40, 10, target = 115))

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
    expect_identical(partitionsGeneral(1001 + (1:40) * 107, 10,
                                       target = 22315, nThreads = 2),
                     partitionsGeneral(1001 + (1:40) * 107, 10,
                                       target = 22315))
    
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
    expect_identical(partitionsGeneral(0:30, 9, freqs = c(3, rep(1, 30)),
                                       target = 115, nThreads = 2),
                     partitionsGeneral(0:30, 9, freqs = c(3, rep(1, 30)),
                                       target = 115))
    
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
    expect_identical(partitionsGeneral(1001 + (0:30) * 107, 9,
                                       freqs = c(3, rep(1, 30)),
                                       target = 21314, nThreads = 2),
                     partitionsGeneral(1001 + (0:30) * 107, 9,
                                       freqs = c(3, rep(1, 30)),
                                       target = 21314))
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
    expect_equal(gmp::sub.bigz("39228755152043560",
                               "39228755151943560") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 1000))
    
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
    expect_equal(gmp::sub.bigz("39233351450439724",
                               "39233351450339724") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 1000))
    
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
    expect_equal(gmp::sub.bigz("10236925075643716",
                               "10236925075443716") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 1380))
})

test_that("partitionGeneral Repetition Parallel", {
    
    ## N.B. Parallel has no effect when number of results is less than 40000
    ## partitionsCount(20)
    ## [1] 7
    expect_identical(partitionsGeneral(0:20, repetition = TRUE, nThreads = 2), 
                     partitionsGeneral(0:20, repetition = TRUE))
    
    ## For both of the usages below, only 2 threads will be spawned
    ## partitionsCount(0:100, repetition = T)
    ## [1] 190569292
    expect_identical(partitionsGeneral(0:100, repetition = TRUE,
                                       nThreads = 3, upper = 50000), 
                     partitionsGeneral(0:100, repetition = TRUE,
                                       nThreads = 8, upper = 50000))
    
    ######****************** All Results Repetition **************#########
    #### Repetition; Length determined internally; No zero;
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
    ## 46 * 3 + 45 * 17 = 903
    ##
    ## partitionsDesign(3L + (0:45) * 17L, 46, repetition = TRUE,
    ##                  target = 903)[c("num_partitions",
    ##                                    "mapped_target", "partition_type")]
    ## $num_partitions
    ## [1] 89134
    ## 
    ## $mapped_target
    ## [1] 91
    ## 
    ## $partition_type
    ## [1] "RepNoZero"
    expect_identical(partitionsGeneral(3L + (0:45) * 17L, 46,
                                       repetition = TRUE,
                                       target = 903L, nThreads = 2),
                     partitionsGeneral(3L + (0:45) * 17L, 46,
                                       repetition = TRUE, target = 903L))
    
    #### Repetition; Specific Length; No zero
    ##
    ## partitionsDesign(60, 10, TRUE)[c("num_partitions", "partition_type")]
    ## $num_partitions
    ## [1] 62740
    ## 
    ## $partition_type
    ## [1] "RepNoZero"
    expect_identical(partitionsGeneral(60, 10, TRUE, nThreads = 2),
                     partitionsGeneral(60, 10, TRUE))
    
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
    expect_identical(partitionsGeneral(6 + (1:60) * 3, 10, TRUE,
                                       target = 240, nThreads = 2),
                     partitionsGeneral(6 + (1:60) * 3, 10,
                                       TRUE, target = 240))

    #### Repetition; Specific Length; Zero included
    ##
    ## partitionsDesign(0:60, 10, TRUE)[c("num_partitions", "partition_type")]
    ## $num_partitions
    ## [1] 195491
    ## 
    ## $partition_type
    ## [1] "RepShort"
    expect_identical(partitionsGeneral(0:60, 10, TRUE, nThreads = 2),
                     partitionsGeneral(0:60, 10, TRUE))
    
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
    expect_identical(partitionsGeneral(6 + (0:60) * 3, 10, TRUE,
                                       target = 240, nThreads = 2),
                     partitionsGeneral(6 + (0:60) * 3, 10,
                                       TRUE, target = 240))
    
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
    expect_identical(partitionsGeneral(18, 10, TRUE,
                                       target = 100, nThreads = 2),
                     partitionsGeneral(18, 10, TRUE, target = 100))
    
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
    expect_identical(partitionsGeneral(1001 + (1:18) * 107, 10, TRUE,
                                       target = 20710, nThreads = 2),
                     partitionsGeneral(1001 + (1:18) * 107, 10, TRUE,
                                       target = 20710))
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


test_that("partitionGeneral Repetition Parallel Lower GMP", {
    
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
    expect_equal(gmp::sub.bigz("9253082936723602",
                               "9253082936523602") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 300))
    
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
    expect_equal(gmp::sub.bigz("9253082936723602",
                               "9253082936523602") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 6000L))
    
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
    expect_equal(gmp::sub.bigz("7856063261819197461639",
                               "7856063261819197261639") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 6000L))
    
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
    expect_equal(gmp::sub.bigz("7856063261819197461639",
                               "7856063261819197261639") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 1.2e+13))
    
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
    expect_equal(gmp::sub.bigz("7974346430545135003177",
                               "7974346430545134803177") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 6000L))
    
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
    expect_equal(gmp::sub.bigz("14001484381527012",
                               "14001484381327012") + 1, nrow(bench))
    expect_true(all(rowSums(bench) == 1000))
})
