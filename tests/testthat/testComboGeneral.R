context("testing comboGeneral")

test_that("comboGeneral produces correct results with no constraints", {

    expect_equal(comboGeneral(5, 3), t(combn(5, 3)))
    expect_equal(comboGeneral(factor(1:5, ordered = TRUE), 3),
                 t(combn(factor(1:5, ordered = TRUE), 3)))

    expect_equal(comboGeneral(as.raw(1:5), 3), t(combn(as.raw(1:5), 3)))

    expect_equal(comboGeneral(factor(1:5, ordered = TRUE), 5, freqs = rep(3, 5)),
                 comboSample(factor(1:5, ordered = TRUE), 5, freqs = rep(3, 5),
                             sampleVec = 1:comboCount(5, 5, freqs = rep(3, 5))))

    expect_equal(nrow(comboGeneral(6, 3)), choose(6, 3))
    expect_equal(as.vector(comboGeneral(1,1)), 1)
    expect_equal(as.vector(comboGeneral(1,1,TRUE)), 1)

    expect_equal(comboGeneral(15, 8)[500:600, ], comboGeneral(15, 8,
                                                              lower = 500,
                                                              upper = 600))
    expect_equal(comboGeneral(5, 5),
                 comboGeneral(5, 5, freqs = rep(1, 5)))

    expect_equal(comboGeneral(as.complex(1:5), 3),
                 t(combn(as.complex(1:5), 3)))

    expect_equal(comboGeneral(as.raw(1:5), 3),
                 t(combn(as.raw(1:5), 3)))

    set.seed(103)
    myNums = rnorm(5)
    expect_equal(comboGeneral(myNums, 3), t(combn(myNums, 3)))

    myNums2 = 1:15 / 3
    expect_equal(comboGeneral(myNums2, 8, freqs = rep(2, 15))[150000:157000, ],
                 comboGeneral(myNums2, 8, freqs = rep(2, 15), lower = 150000, upper = 157000))

    expect_equal(comboGeneral(letters[1:15], 8, TRUE)[319000:comboCount(15, 8, TRUE), ],
                 comboGeneral(letters[1:15], 8, TRUE, lower = 319000))

    expect_equal(nrow(comboGeneral(myNums, 3, TRUE)), 35)
    expect_equal(comboGeneral(myNums, 3, freqs = rep(3, 5)),
                 comboGeneral(myNums, 3, TRUE))

    expect_equal(comboGeneral(LETTERS[1:5], 3, freqs = rep(3, 5)),
                 comboGeneral(LETTERS[1:5], 3, TRUE))

    expect_equal(comboGeneral(LETTERS[1:5], 3), t(combn(LETTERS[1:5], 3)))

    expect_equal(sum(comboGeneral(3, 3, freqs = c(1, 1, 1))),
                 sum(comboGeneral(3, 3)))

    expect_equal(as.vector(comboGeneral(1, 2, freqs = 2)), c(1, 1))
    expect_equal(ncol(comboGeneral(5, 3)), 3)
    expect_equal(ncol(comboGeneral(5, 3, TRUE)), 3)
    expect_equal(ncol(comboGeneral(5, 3, FALSE, constraintFun = "prod")), 4)
    expect_equal(ncol(comboGeneral(5, 3, TRUE, constraintFun = "prod", keepResults = TRUE)), 4)

    expect_equal(ncol(comboGeneral(5, 3, FALSE,
                                    constraintFun = "prod", freqs = c(1,2,1,2,4),
                                     keepResults = TRUE)), 4)

    expect_equal(nrow(comboGeneral(10, 3, TRUE, upper = 20)), 20)
    expect_equal(nrow(comboGeneral(10, 3, upper = 10)), 10)

    expect_equal(nrow(comboGeneral(1:10 + .01, 3, FALSE, constraintFun = "prod",
                                     keepResults = TRUE, upper = 10)), 10)

    expect_equal(nrow(comboGeneral(5, 5, freqs = 1:5, upper = 10)), 10)

    expect_equal(nrow(comboGeneral(10, 3, TRUE, constraintFun = "prod",
                                     keepResults = TRUE, upper = 10)), 10)

    ##******** BIG TESTS *********##

    ## NO REPETITION
    numR = comboCount(1000, 10)
    n1 = gmp::sub.bigz(numR, 99)

    ## accepts raw values
    expect_equal(nrow(comboGeneral(1000, 10, lower = n1)), 100)
    ## accepts characters
    expect_equal(nrow(comboGeneral(1000, 10, lower = as.character(n1))), 100)
    expect_equal(as.vector(comboGeneral(1000, 10, lower = numR)), 991:1000)

    ## WITH REPETITION
    numR = comboCount(1000, 10, TRUE)
    n1 = gmp::sub.bigz(numR, 99)

    expect_equal(nrow(comboGeneral(1000, 10, TRUE, lower = n1)), 100)
    expect_equal(nrow(comboGeneral(1000, 10, TRUE, lower = as.character(n1))), 100)
    expect_equal(as.vector(comboGeneral(1000, 10, TRUE, lower = numR)), rep(1000, 10))

    ## MULTISETS
    numR = comboCount(1000, 10, freqs = rep(1:4, 250))
    n1 = gmp::sub.bigz(numR, 99)

    expect_equal(nrow(comboGeneral(1000, 10, freqs = rep(1:4, 250), lower = n1)), 100)
    expect_equal(nrow(comboGeneral(1000, 10, freqs = rep(1:4, 250), lower = as.character(n1))), 100)
    expect_equal(as.vector(comboGeneral(1000, 10, freqs = rep(1:4, 250), lower = numR)),
                 rep(997:1000, times = 1:4))
})

test_that("comboGeneral produces correct results with constraints", {
    tinyTol = nrow(comboGeneral(1:5 + 0.00000000001, 3,
                      constraintFun = "mean",
                      comparisonFun = "==",
                      limitConstraints = 3,
                      tolerance = .Machine$double.eps))

    ## The default tolerance is sqrt(.Machine$double.eps)
    defaultTol = nrow(comboGeneral(1:5 + 0.00000000001, 3,
                                     constraintFun = "mean",
                                     comparisonFun = "==", limitConstraints = 3))

    expect_false(tinyTol == defaultTol)

    ## check that classes behave properly N.B. limitContraint > INT_MAX
    expect_equal(class(comboGeneral(10, 5, constraintFun = "prod",
                                    comparisonFun = "<",
                                    limitConstraints = 2^32)[,1]), "numeric")

    ## the greatest product is prod(100:96) > INT_MAX
    expect_equal(class(comboGeneral(100:90, 5, constraintFun = "prod",
                                    comparisonFun = "<",
                                    limitConstraints = 5)[,1]), "numeric")

    expect_equal(class(comboGeneral(5, 5, TRUE, constraintFun = "prod",
                                    comparisonFun = "<",
                                    limitConstraints = 5.5)[,1]), "numeric")

    expect_equal(nrow(comboGeneral(-5:5, 4, FALSE, constraintFun = "sum",
                                   comparisonFun = "==", limitConstraints = 6)),
                 length(which(apply(combn(-5:5, 4), 2, sum) == 6)))

    expect_equal(unique(comboGeneral(5, 5, TRUE,
                                       constraintFun = "sum", comparisonFun = "==",
                                     limitConstraints = 9,
                                       keepResults = TRUE)[,6]), 9)

    expect_true(all(comboGeneral(5, 5, TRUE,
                                 constraintFun = "min", comparisonFun = "<",
                                 limitConstraints = 3,
                                   keepResults = TRUE)[,6] < 3))
    expect_equal(as.vector(comboGeneral(5, 5, freqs = 1:5,
                                        constraintFun = "sum",
                                        comparisonFun = "==",
                                        limitConstraints = 25)), rep(5, 5))

    expect_equal(as.vector(comboGeneral(5, 5, constraintFun = "sum",
                                        comparisonFun = "==",
                                        limitConstraints = 15)), 1:5)

    expect_true(all(comboGeneral(5, 5, TRUE,
                                 constraintFun = "prod", comparisonFun = ">",
                                 limitConstraints = 100,
                                   keepResults = TRUE)[,6] > 100))

    expect_true(all(comboGeneral(5, 3, FALSE,
                                 constraintFun = "max", comparisonFun = "=<",
                                 limitConstraints = 4,
                                   keepResults = TRUE)[,4] <= 4))

    ## N.B. When there are two comparisons (i.e. comparisonFun = c(">=","<"))
    ## and only one limitConstraint, the first comparison is used. Similarly,
    ## when there are two limitConstraints and one comparison the first
    ## limitConstraint is use and the other is ignored (See next test)
    expect_true(all(comboGeneral(3, 5, TRUE,
                                 constraintFun = "mean", comparisonFun = c(">=","<"),
                                 limitConstraints = 2L,
                                   keepResults = TRUE)[,6] >= 2))

    expect_true(all(comboGeneral(3, 5, TRUE,
                                 constraintFun = "mean", comparisonFun = "==",
                                 limitConstraints = 2,
                                 keepResults = TRUE)[,6] == 2))

    expect_true(all(comboGeneral(3, 5, TRUE,
                                 constraintFun = "mean", comparisonFun = ">=",
                                 limitConstraints = c(2L, 1e10),
                                 keepResults = TRUE)[,6] >= 2))

    expect_true(all(comboGeneral(5, 5, FALSE, constraintFun = "sum", comparisonFun = ">",
                                 limitConstraints = 18,
                                   freqs = c(1,2,1,2,4),
                                   keepResults = TRUE)[,6] > 18))

    expect_equal(comboGeneral(5, 5, TRUE, constraintFun = "sum",
                              comparisonFun = "==",
                              limitConstraints = 25, keepResults = TRUE)[,6], 25)

    all_combs = comboGeneral(10, 5, freqs = rep(1:5, 2), constraintFun = "sum")
    all_combs = all_combs[500:nrow(all_combs), ]
    expect_equal(comboGeneral(10, 5, freqs = rep(1:5, 2), keepResults = TRUE,
                              constraintFun = "sum", comparisonFun = "==",
                              limitConstraints =  25, lower = 500),
                 all_combs[all_combs[,6] == 25, ])


    test = comboGeneral(10, 5, freqs = rep(1:5, 2),
                        constraintFun = "sum", comparisonFun = c(">=", "=<"),
                        limitConstraints = c(27, 23),
                        lower = 1000, tolerance = 1, keepResults = TRUE)
    bench = comboGeneral(10, 5, freqs = rep(1:5, 2),
                         constraintFun = "sum", lower = 1000)
    expect_equal(test, bench[c(which(bench[, 6] > 26 | bench[, 6] < 24)), ])
})

test_that("comboGeneral produces correct results with use of FUN", {

    test <- comboGeneral(10, 5, constraintFun = "sum")
    expect_equal(as.vector(test[,6]), unlist(comboGeneral(10, 5, FUN = sum)))
    expect_equal(dim(comboGeneral(LETTERS, 5, freqs = c(rep(1:4, 6), 1:2),
                                  FUN = function(x) sapply(x, charToRaw),
                                  upper = 8, FUN.VALUE = as.raw(1:5))), c(8, 5))

    expect_equal(class(comboGeneral(
        LETTERS, 5, TRUE, FUN = function(x) {
            sapply(x, function(y) rawToBits(charToRaw(y)))
        }, upper = 8, FUN.VALUE = rawToBits(as.raw(1:5))
    )[1, ]), "raw")
    test <- comboGeneral(10, 4, TRUE)
    testFun <- apply(test, 1, function(x) mean(x) * 2)
    expect_equal(testFun,
                 comboGeneral(10, 4, TRUE, FUN = function(x) {mean(x) * 2},
                              FUN.VALUE = 2.2))

    test <- comboGeneral(8, 4, freqs = rep(1:4, 2))
    testFun <- lapply(1:nrow(test), function(x) cumsum(test[x, ]))
    expect_equal(testFun, comboGeneral(8, 4, freqs = rep(1:4, 2), FUN = cumsum))

    expect_equal(class(comboGeneral(5, 3, FUN = cumsum,
                                    FUN.VALUE = as.numeric(1:3))[1, ]), "numeric")
    expect_equal(class(comboGeneral(5, 3, FUN = cumsum,
                                    FUN.VALUE = as.complex(1:3))[1, ]), "complex")

    expect_equal(testFun[1:100], comboGeneral(8, 4, freqs = rep(1:4, 2), upper = 100, FUN = cumsum))
    expect_equal(testFun[101:length(testFun)],
                 comboGeneral(8, 4, freqs = rep(1:4, 2), lower = 101, FUN = cumsum))
    expect_equal(testFun[121:123],
                 comboGeneral(8, 4, freqs = rep(1:4, 2), lower = 121, upper = 123, FUN = cumsum))

    expect_equal(comboGeneral(as.raw(1:5), 3, FUN = rawToChar),
                 combn(as.raw(1:5), 3, rawToChar, simplify = FALSE))

    expect_equal(unlist(comboGeneral(letters[1:5], 3, FUN = function(x) {
        paste0(x, collapse = "")
    })), apply(comboGeneral(letters[1:5], 3), 1, paste0, collapse = ""))
})

test_that("comboGeneral produces correct results with exotic constraints", {

    a = t(combn(10, 7))
    expect_equal(comboGeneral(10, 7, constraintFun = "sum",
                 comparisonFun = c(">","<"),
                 limitConstraints = c(40, 45)), a[which(rowSums(a) > 40 & rowSums(a) < 45), ])

    set.seed(13)
    rSet = 1:10 + rnorm(10)
    a = comboGeneral(sort(rSet), 7, TRUE)
    b = rowSums(a)
    expect_equal(comboGeneral(rSet, 7, TRUE, constraintFun = "sum",
                              comparisonFun = c(">=","<="),
                              limitConstraints = c(42.50001, 45.76277)),
                 a[which(b >= 42.50001 & b <= 45.76277), ])

    temp1 = comboGeneral(rSet, 7, TRUE, constraintFun = "sum",
                         comparisonFun = c("<=",">="),
                         limitConstraints = c(20.05669, 60.93901),
                         keepResults = TRUE)

    temp2 = cbind(a, b)
    temp2 = temp2[which(b <= 20.05669 | b >= 60.93901), ]

    expect_equal(sort(temp1[,8]), sort(temp2[,8]))
    a = comboGeneral(10, 7, freqs = rep(3, 10))
    b = rowSums(a)
    expect_equal(comboGeneral(10, 7, freqs = rep(3, 10), constraintFun = "sum",
                 comparisonFun = c("<=", ">"),
                 limitConstraints = c(50, 47)), a[which(b > 47 & b <= 50), ])

    b = apply(a, 1, max)
    expect_equal(comboGeneral(10, 7, freqs = rep(3, 10),
                 constraintFun = "max",
                 comparisonFun = c("<=", ">"),
                 limitConstraints = c(9, 7)), a[which(b > 7 & b <= 9), ])

    b = apply(a, 1, min)
    expect_equal(comboGeneral(10, 7, freqs = rep(3, 10),
                              constraintFun = "min",
                              comparisonFun = "==",
                              limitConstraints = 3,
                              lower = 7900, upper = 8500),
                 a[(7900:8500)[b[7900:8500] == 3], ])

    a = comboGeneral(5, 7, TRUE)
    b = apply(a, 1, prod)
    expect_equal(comboGeneral(5, 7, TRUE, constraintFun = "prod",
                 comparisonFun = c(">=","<="),
                 limitConstraints = c(2000, 5000)), a[which(b >= 2000 & b <= 5000), ])

    a = comboGeneral(-5, 7, TRUE)
    b = apply(a, 1, prod)
    expect_equal(nrow(comboGeneral(-5, 7, TRUE, constraintFun = "prod",
                              comparisonFun = c("<=",">="),
                              limitConstraints = c(-2000, 5000),
                              keepResults = TRUE)),
                 nrow(rbind(a[which(b <= -2000),], a[which(b >= 5000), ])))

    set.seed(4321)
    samp = sample(-50:50, 16)
    a = comboGeneral(samp, 6, TRUE)
    b = apply(a, 1, prod)
    expect_equal(nrow(comboGeneral(samp, 6, TRUE, constraintFun = "prod",
                                   comparisonFun = c("<",">"),
                                   limitConstraints = c(-6e9, 6e9),
                                   keepResults = TRUE, nThreads = 2)),
                 nrow(rbind(a[which(b <= -6e9),], a[which(b >= 6e9), ])))

    ## Testing sums in a range
    a = comboGeneral(10, 8, TRUE, lower = 23500, upper = 24000,
                      constraintFun = "sum", keepResults = TRUE)

    expect_equal(comboGeneral(c(NA, 1:10), 8, TRUE, constraintFun = "sum",
                              comparisonFun = c("=>","=<"),
                              limitConstraints = c(72, 78),
                              lower = 23500, upper = 24000),
                 a[a[,9] >= 72 & a[,9] <= 78, 1:8])


    comp1 = c("<", "<=")
    comp2 = c(">", ">=")

    ## Test that unsorted vector is being handled properly
    ## for both numeric and integer type vectors
    # identical(sort(scrambled), 1:10)
    # [1] TRUE
    scrambled = as.integer(c(8, 2, 5, 1, 6, 3, 4, 7))
    scramFreqs = rep(1:4, 2)[scrambled]
    funs <- c("sum", "prod", "mean", "max", "min")
    m <- 7

    allCombs1 = lapply(funs, function(f) {
        comboGeneral(8, m, freqs = rep(1:4, 2), constraintFun = f)
    })

    allCombs2 = lapply(funs, function(f) {
        comboGeneral(8:1, m, freqs = rev(rep(1:4, 2)), constraintFun = f)
    })

    ## ensure the left bound is in the solution space
    tars = lapply(allCombs1, function(x) {
        vals <- sort(x[, m + 1])
        t <- quantile(as.numeric(names(table(vals))),
                      c(0.25, 0.75), names = FALSE)
        t[1] <- vals[findInterval(t[1], vals)]
        t
    })

    theVals1  = lapply(allCombs1, function(x) x[, m + 1])
    theVals2  = lapply(allCombs2, function(x) x[, m + 1])
    allCombs1 = lapply(allCombs1, function(x) x[, 1:m])
    allCombs2 = lapply(allCombs2, function(x) x[, 1:m])

    for (f in seq_along(funs)) {
        for (i in 1:2) {

            if (i == 1) {
                a = comp1
                b = comp2
            } else {
                a = comp2
                b = comp1
            }

            for (j in a) {
                for (k in b) {
                    myComp = c(j, k)
                    myTest = comboGeneral(scrambled, m, freqs = scramFreqs,
                                          constraintFun = funs[f],
                                          comparisonFun = myComp,
                                          limitConstraints = tars[[f]],
                                          tolerance = 0)

                    fun1 = match.fun(j)
                    fun2 = match.fun(k)

                    if (i == 1) {
                        temp1 = allCombs1[[f]][fun1(theVals1[[f]],
                                                    tars[[f]][1]), ]
                        temp2 = allCombs2[[f]][fun2(theVals2[[f]],
                                                    tars[[f]][2]), ]
                        temp = rbind(temp1, temp2)
                    } else {
                        temp = allCombs1[[f]][fun1(theVals1[[f]], tars[[f]][1]) &
                                              fun2(theVals1[[f]], tars[[f]][2]), ]
                    }

                    expect_equal(temp, myTest, info = list(myComp, funs[f], tars[[f]], i))
                }
            }
        }
    }

    allCombs1 = lapply(funs, function(f) {
        comboGeneral(8, m, TRUE, constraintFun = f, keepResults = TRUE)
    })

    allCombs2 = lapply(funs, function(f) {
        comboGeneral(8:1, m, TRUE, constraintFun = f, keepResults = TRUE)
    })

    ## ensure the right bound is in the solution space
    tars = lapply(allCombs1, function(x) {
        vals <- sort(x[, m + 1])
        t <- quantile(as.numeric(names(table(vals))),
                      c(0.25, 0.75), names = FALSE)
        t[2] <- vals[findInterval(t[2], vals)]
        t
    })

    theVals1  = lapply(allCombs1, function(x) x[, m + 1])
    theVals2  = lapply(allCombs2, function(x) x[, m + 1])
    allCombs1 = lapply(allCombs1, function(x) x[, 1:m])
    allCombs2 = lapply(allCombs2, function(x) x[, 1:m])

    for (f in seq_along(funs)) {
        for (i in 1:2) {

            if (i == 1) {
                a = comp1
                b = comp2
            } else {
                a = comp2
                b = comp1
            }

            for (j in a) {
                for (k in b) {
                    myComp = c(j, k)
                    myTest = comboGeneral(8, m, TRUE,
                                          constraintFun = funs[f],
                                          comparisonFun = myComp,
                                          limitConstraints = tars[[f]],
                                          tolerance = 0)

                    fun1 = match.fun(j)
                    fun2 = match.fun(k)

                    if (i == 1) {
                        temp1 = allCombs1[[f]][fun1(theVals1[[f]], tars[[f]][1]), ]
                        temp2 = allCombs2[[f]][fun2(theVals2[[f]], tars[[f]][2]), ]
                        temp = rbind(temp1, temp2)
                    } else {
                        temp = allCombs1[[f]][fun1(theVals1[[f]], tars[[f]][1]) &
                                              fun2(theVals1[[f]], tars[[f]][2]), ]
                    }

                    expect_equal(temp, myTest, info = list(myComp, funs[f], tars[[f]], i))
                }
            }
        }
    }

    allCombs1 = lapply(funs, function(f) {
        comboGeneral(15, m, constraintFun = f, keepResults = TRUE)
    })

    allCombs2 = lapply(funs, function(f) {
        comboGeneral(15:1, m, constraintFun = f, keepResults = TRUE)
    })

    ## ensure the both bounds are in the solution space
    tars = lapply(allCombs1, function(x) {
        vals <- sort(x[, m + 1])
        t <- quantile(as.numeric(names(table(vals))),
                      c(0.25, 0.75), names = FALSE)
        t <- vals[findInterval(t, vals)]
        t
    })

    theVals1  = lapply(allCombs1, function(x) x[, m + 1])
    theVals2  = lapply(allCombs2, function(x) x[, m + 1])
    allCombs1 = lapply(allCombs1, function(x) x[, 1:m])
    allCombs2 = lapply(allCombs2, function(x) x[, 1:m])

    for (f in seq_along(funs)) {
        for (i in 1:2) {

            if (i == 1) {
                a = comp1
                b = comp2
            } else {
                a = comp2
                b = comp1
            }

            for (j in a) {
                for (k in b) {
                    myComp = c(j, k)
                    myTest = comboGeneral(15, m,
                                          constraintFun = funs[f],
                                          comparisonFun = myComp,
                                          limitConstraints = tars[[f]],
                                          tolerance = 0)

                    fun1 = match.fun(j)
                    fun2 = match.fun(k)

                    if (i == 1) {
                        temp1 = allCombs1[[f]][fun1(theVals1[[f]], tars[[f]][1]), ]
                        temp2 = allCombs2[[f]][fun2(theVals2[[f]], tars[[f]][2]), ]
                        temp = rbind(temp1, temp2)
                    } else {
                        temp = allCombs1[[f]][fun1(theVals1[[f]], tars[[f]][1]) &
                                              fun2(theVals1[[f]], tars[[f]][2]), ]
                    }

                    expect_equal(temp, myTest)
                }
            }
        }
    }
})
