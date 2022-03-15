context("testing special subset sum")

test_that("comboGeneral produces correct results for special subset sum", {
    # testing comboGeneral and permuteGeneral results for constrainFun = 'sum',
    # comparisonFun = '==', and with v having the property that if you were to
    # sort v, the difference of each element with it's neighbor is constant.

    testCombFun <- function(v, m, myRep = FALSE, verbose = FALSE,
                            f = "sum", isExact = TRUE) {
        v <- sort(v)
        allSums <- comboGeneral(v, m, myRep, constraintFun = f)
        tbl <- table(allSums[, m + 1])
        possVals <- as.numeric(names(tbl))

        if (verbose) {
            print(possVals)
            print(tbl)
        }

        ## Credit to Johan Larsson: https://stackoverflow.com/a/39175037/4408538
        is_equal_tol <- function(x, y, tol = sqrt(.Machine$double.eps)) {
            abs(x - y) < tol
        }

        t <- sapply(possVals, function(x) {
            a <- comboGeneral(v, m, myRep,
                              constraintFun = f,
                              comparisonFun = "==",
                              limitConstraints = x)
            if (isExact) {
                u <- allSums[allSums[, m + 1] == x, 1:m]
            } else {
                u <- allSums[which(is_equal_tol(allSums[, m + 1], x)), 1:m]
            }

            if (nrow(a) > 1) {
                identical(a, u)
            } else {
                identical(as.vector(a), u)
            }
        })

        if (verbose) {
            print(t)
        }

        all(t)
    }

    expect_true(testCombFun(1:18, 9))
    expect_true(testCombFun(0:17, 9))
    expect_true(testCombFun(-8:8, 9))
    expect_true(testCombFun(1:5, 5))
    expect_true(testCombFun((1e10 + 1):(1e10 + 18), 9))
    expect_true(testCombFun((-1e10 - 1):(-1e10 - 18), 9))
    expect_true(testCombFun(-49:50, 99))
    expect_true(testCombFun(1:100, 99))
    expect_true(testCombFun(-49:50, 2))
    expect_true(testCombFun(1:100, 2))
    expect_true(testCombFun(1:100, 100))
    expect_true(testCombFun((-1e12 - 50):(-1e12 - 1), 49))
    expect_true(testCombFun((-1e12 - 50):(-1e12 - 1), 3))

    expect_true(testCombFun(1:10, 7, myRep = TRUE))
    expect_true(testCombFun(0:9, 7, myRep = TRUE))
    expect_true(testCombFun(-4:5, 7, myRep = TRUE))
    expect_true(testCombFun((1e10 + 1):(1e10 + 10), 7, myRep = TRUE))
    expect_true(testCombFun((-1e10 - 1):(-1e10 - 10), 7, myRep = TRUE))
    expect_true(testCombFun(-49:50, 2, myRep = TRUE))
    expect_true(testCombFun(1:100, 2, myRep = TRUE))
    expect_true(testCombFun((-1e12 - 50):(-1e12 - 1), 3, myRep = TRUE))
    expect_true(testCombFun(1:5, 10, myRep = TRUE))
    expect_true(testCombFun(-1:1, 100, myRep = TRUE))

    expect_true(testCombFun(seq(100, 180, 5), 9))
    expect_true(testCombFun(seq(-80, 80, 10), 9))
    expect_true(testCombFun(seq(1e10, 1e10 + 180, 10), 9))

    expect_true(testCombFun(seq(100, 210, 10), 7, myRep = TRUE))
    expect_true(testCombFun(seq(-140L, -100L, 10L), 10, myRep = TRUE))
    expect_true(testCombFun(seq(-100L, 300L, 100L), 10, myRep = TRUE))
    expect_true(testCombFun(seq(1e10, 1e10 + 500, 100), 10, myRep = TRUE))
    expect_true(testCombFun(seq(-1e10 - 500, -1e10, 100), 10, myRep = TRUE))

    ## Irregular vector input... i.e. the distance between neighbors varies
    ## This will test the BruteNextElem and main constraint functions
    p <- as.integer(c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41))
    expect_true(testCombFun(p, 1))
    expect_true(testCombFun(p, 2))
    expect_true(testCombFun(p, 7))
    expect_true(testCombFun(p, 13))

    pN <- as.numeric(p)
    expect_true(testCombFun(pN, 1))
    expect_true(testCombFun(pN, 2))
    expect_true(testCombFun(pN, 7))
    expect_true(testCombFun(pN, 13))

    expect_true(testCombFun(p, 1, f = "prod"))
    expect_true(testCombFun(p, 2, f = "prod"))
    expect_true(testCombFun(p, 6, f = "prod"))

    expect_true(testCombFun(p, 1, f = "mean", isExact = F))
    expect_true(testCombFun(p, 2, f = "mean", isExact = F))
    expect_true(testCombFun(p, 7, f = "mean", isExact = F))
    expect_true(testCombFun(p, 13, f = "mean", isExact = F))

    pS <- p[1:6]
    expect_true(testCombFun(pS, 1, T))
    expect_true(testCombFun(pS, 2, T))
    expect_true(testCombFun(pS, 6, T))
    expect_true(testCombFun(pS, 8, T))

    pNS <- as.numeric(pS)
    expect_true(testCombFun(pNS, 1, T))
    expect_true(testCombFun(pNS, 2, T))
    expect_true(testCombFun(pNS, 6, T))
    expect_true(testCombFun(pNS, 8, T))

    expect_true(testCombFun(pS, 1, T, f = "prod"))
    expect_true(testCombFun(pS, 2, T, f = "prod"))
    expect_true(testCombFun(pS, 6, T, f = "prod"))
    expect_true(testCombFun(pS, 7, T, f = "prod"))

    expect_true(testCombFun(pS, 1, T, f = "mean", isExact = F))
    expect_true(testCombFun(pS, 2, T, f = "mean", isExact = F))
    expect_true(testCombFun(pS, 6, T, f = "mean", isExact = F))
    expect_true(testCombFun(pS, 8, T, f = "mean", isExact = F))

    ## Standard partitions into k parts
    tempCombs <- comboGeneral(15, 6, TRUE, constraintFun = "sum")
    expect_equal(tempCombs[tempCombs[,ncol(tempCombs)] == 15, 1:6],
                 comboGeneral(15, 6, TRUE,
                              constraintFun = "sum",
                              comparisonFun = "==",
                              limitConstraints = 15))

    ## Standard partitions
    tempCombs <- comboGeneral(0:10, 10, TRUE, constraintFun = "sum")
    expect_equal(tempCombs[tempCombs[,ncol(tempCombs)] == 10, 1:10],
                 comboGeneral(0:10, 10, TRUE,
                              constraintFun = "sum",
                              comparisonFun = "==",
                              limitConstraints = 10))

    expect_true(all(rowSums(comboGeneral(0:20, 20, TRUE,
                                         constraintFun = "sum",
                                         comparisonFun = "==",
                                         limitConstraints = 20,
                                         keepResults = TRUE)) == 40))

    ## Testing cases where no results should be returned
    expect_equal(nrow(comboGeneral(10, 6,
                                   constraintFun = "sum",
                                   comparisonFun = "==",
                                   limitConstraints = 20)), 0)

    expect_equal(nrow(comboGeneral(10, 6,
                                   constraintFun = "sum",
                                   comparisonFun = "==",
                                   limitConstraints = 46)), 0)

    expect_equal(nrow(comboGeneral(6, 10, TRUE,
                                   constraintFun = "sum",
                                   comparisonFun = "==",
                                   limitConstraints = 5)), 0)

    expect_equal(nrow(comboGeneral(6, 10, TRUE,
                                   constraintFun = "sum",
                                   comparisonFun = "==",
                                   limitConstraints = 61)), 0)

    expect_equal(nrow(comboGeneral(seq(10, 100, 5), 8,
                                   constraintFun = "sum",
                                   comparisonFun = "==",
                                   limitConstraints = 402)), 0)

    expect_equal(nrow(comboGeneral(seq(10, 100, 5), 8, TRUE,
                                   constraintFun = "sum",
                                   comparisonFun = "==",
                                   limitConstraints = 402)), 0)


    ## nrow(comboGeneral(4, 5, TRUE,
    ##                   constraintFun = "sum",
    ##                   comparisonFun = "==",
    ##                   limitConstraints = 10))
    ## [1] 5

    expect_equal(nrow(comboGeneral(4, 5, TRUE,
                                   constraintFun = "sum",
                                   comparisonFun = "==",
                                   limitConstraints = 10,
                                   upper = 3)), 3)

    testCombMultiset <- function(v, m, frqs, verbose = FALSE,
                                 f = "sum", isExact = TRUE) {
        v <- sort(v)
        allSums <- comboGeneral(v, m, freqs = frqs, constraintFun = f)
        tbl <- table(allSums[, m + 1])
        possVals <- as.numeric(names(tbl))

        if (verbose) {
            print(possVals)
            print(tbl)
        }

        ## Credit to Johan Larsson: https://stackoverflow.com/a/39175037/4408538
        is_equal_tol <- function(x, y, tol = sqrt(.Machine$double.eps)) {
            abs(x - y) < tol
        }

        t <- sapply(possVals, function(x) {
            t <- comboGeneral(v, m, freqs = frqs,
                              constraintFun = f,
                              comparisonFun = "==",
                              limitConstraints = x)
            if (isExact) {
                u <- allSums[allSums[, m + 1] == x, 1:m]
            } else {
                u <- allSums[which(is_equal_tol(allSums[, m + 1], x)), 1:m]
            }

            if (nrow(t) > 1) {
                identical(t, u)
            } else {
                identical(as.vector(t), u)
            }
        })

        if (verbose) {
            print(t)
        }

        all(t)
    }

    expect_true(testCombMultiset(1:10, 7, rep(1:5, 2)))
    scrambled = as.integer(c(8, 2, 5, 10, 1, 6, 3, 9, 4, 7))
    expect_true(testCombMultiset(scrambled, 7, rep(1:5, 2)[scrambled]))

    expect_true(testCombMultiset(0:9, 7, rep(1:5, 2)))
    expect_true(testCombMultiset(-4:5, 7, rep(1:5, 2)))
    expect_true(testCombMultiset((1e10 + 1):(1e10 + 10), 7, rep(1:5, 2)))
    expect_true(testCombMultiset((-1e10 - 1):(-1e10 - 10), 7, rep(1:5, 2)))
    expect_true(testCombMultiset(-49:50, 2, rep(1:2, 50)))
    expect_true(testCombMultiset(1:100, 2, rep(1:2, 50)))
    expect_true(testCombMultiset((-1e12 - 50):(-1e12 - 1), 3, rep(1:2, 25)))
    expect_true(testCombMultiset(1:5, 10, 1:5))
    expect_true(testCombMultiset(-1:1, 100, c(20, 30, 50)))

    expect_true(testCombMultiset(seq(100, 210, 10), 7, rep(1:4, 3)))
    expect_true(testCombMultiset(seq(-140L, -100L, 10L), 10, c(1, 2, 3, 4, 3)))
    expect_true(testCombMultiset(seq(-100L, 300L, 100L), 9, c(5, 1, 1, 1, 1)))
    expect_true(testCombMultiset(seq(1e10, 1e10 + 500, 100), 10, c(1, 1, 5, 1, 1, 1)))
    expect_true(testCombMultiset(seq(-1e10 - 500, -1e10, 100), 10, c(1, 1, 1, 1, 1, 5)))

    ## Irregular vector input... i.e. the distance between neighbors varies
    ## This will test the BruteNextElem and main constraint functions
    pS <- as.integer(c(2, 3, 5, 7, 11, 13))
    expect_true(testCombMultiset(pS, 1, frqs = 1:6))
    expect_true(testCombMultiset(pS, 2, frqs = 1:6))
    expect_true(testCombMultiset(pS, 6, frqs = 1:6))
    expect_true(testCombMultiset(pS, 8, frqs = 1:6))

    ## This is equivalent to combinations without rep
    expect_true(testCombMultiset(pS, 1, frqs = rep(1, 6)))
    expect_true(testCombMultiset(pS, 2, frqs = rep(1, 6)))
    expect_true(testCombMultiset(pS, 6, frqs = rep(1, 6)))

    pNS <- as.numeric(pS)
    expect_true(testCombMultiset(pNS, 1, frqs = 1:6))
    expect_true(testCombMultiset(pNS, 2, frqs = 1:6))
    expect_true(testCombMultiset(pNS, 6, frqs = 1:6))
    expect_true(testCombMultiset(pNS, 8, frqs = 1:6))

    expect_true(testCombMultiset(pS, 1, frqs = 1:6, f = "prod"))
    expect_true(testCombMultiset(pS, 2, frqs = 1:6, f = "prod"))
    expect_true(testCombMultiset(pS, 7, frqs = 1:6, f = "prod"))

    expect_true(testCombMultiset(pS, 1, frqs = 1:6, f = "mean", isExact = F))
    expect_true(testCombMultiset(pS, 2, frqs = 1:6, f = "mean", isExact = F))
    expect_true(testCombMultiset(pS, 6, frqs = 1:6, f = "mean", isExact = F))
    expect_true(testCombMultiset(pS, 8, frqs = 1:6, f = "mean", isExact = F))


    ## ******************* Testing Distinct Partitions ******************* ##
    expect_equal(nrow(comboGeneral(0:100, 4, freqs = c(20, rep(1, 100)),
                                     constraintFun = "sum",
                                     comparisonFun = "==",
                                     limitConstraints = 100)), 6786)

    expect_equal(nrow(comboGeneral(0:100, freqs = c(12, rep(1, 100)),
                                   constraintFun = "sum",
                                   comparisonFun = "==", limitConstraints = 100)), 444793)

    expect_equal(nrow(comboGeneral(0:100, 5, freqs = c(12, rep(1, 100)),
                                   constraintFun = "sum",
                                   comparisonFun = "==", limitConstraints = 100)), 32123)

    expect_equal(nrow(comboGeneral(0:10, 3,
                                   constraintFun = "sum",
                                   comparisonFun = "==",
                                   limitConstraints = 10)), 8)

    expect_equal(nrow(permuteGeneral(0:10, constraintFun = "sum",
                                     comparisonFun = "==", limitConstraints = 10)), 120)

    expect_equal(nrow(permuteGeneral(0:10, 4, constraintFun = "sum",
                                     comparisonFun = "==", limitConstraints = 10)), 120)

    expect_equal(nrow(permuteGeneral(10, 3, constraintFun = "sum",
                                     comparisonFun = "==", limitConstraints = 10)), 24)

    expect_equal(nrow(permuteGeneral(0:10, repetition = TRUE, constraintFun = "sum",
                                     comparisonFun = "==", limitConstraints = 10)), 92378)

    expect_equal(nrow(permuteGeneral(0:10, 10, repetition = TRUE, constraintFun = "sum",
                                     comparisonFun = "==", limitConstraints = 10)), 92378)
})

test_that("permuteGeneral produces correct results for special subset sum", {
    testPermFun <- function(v, m, myRep = FALSE, verbose = FALSE) {
        v <- sort(v)
        allSums <- permuteGeneral(v, m, myRep, constraintFun = "sum")
        tbl <- table(allSums[, m + 1])
        possVals <- as.numeric(names(tbl))

        if (verbose) {
            print(possVals)
            print(tbl)
        }

        t <- sapply(possVals, function(x) {
            t <- permuteGeneral(v, m, myRep,
                                constraintFun = "sum",
                                comparisonFun = "==",
                                limitConstraints = x)
            u <- allSums[allSums[, m + 1] == x, 1:m]

            if (nrow(t) > 1) {
                identical(t[do.call(order, as.data.frame(t)), ], u)
            } else {
                identical(as.vector(t), u)
            }
        })

        if (verbose)
            print(t)

        all(t)
    }

    expect_true(testPermFun(1:8, 6))
    expect_true(testPermFun(0:7, 6))
    expect_true(testPermFun(-4:3, 6))
    expect_true(testPermFun(1:5, 5))
    expect_true(testPermFun((1e10 + 1):(1e10 + 8), 6))
    expect_true(testPermFun((-1e10 - 1):(-1e10 - 8), 6))
    expect_true(testPermFun(-49:50, 2))
    expect_true(testPermFun(1:100, 2))
    expect_true(testPermFun((-1e12 - 30):(-1e12 - 1), 3))

    expect_true(testPermFun(1:7, 5, myRep = TRUE))
    expect_true(testPermFun(0:6, 5, myRep = TRUE))
    expect_true(testPermFun(-3:3, 5, myRep = TRUE))
    expect_true(testPermFun((1e10 + 1):(1e10 + 7), 5, myRep = TRUE))
    expect_true(testPermFun((-1e10 - 1):(-1e10 - 7), 5, myRep = TRUE))
    expect_true(testPermFun(-49:50, 2, myRep = TRUE))
    expect_true(testPermFun(1:100, 2, myRep = TRUE))
    expect_true(testPermFun((-1e12 - 50):(-1e12 - 1), 3, myRep = TRUE))
    expect_true(testPermFun(1:4, 7, myRep = TRUE))
    expect_true(testPermFun(-1:1, 9, myRep = TRUE))

    expect_true(testPermFun(seq(100, 135, 5), 6))
    expect_true(testPermFun(seq(-40, 30, 10), 6))
    expect_true(testPermFun(seq(1e10, 1e10 + 80, 10), 6))

    expect_true(testPermFun(seq(100, 160, 10), 5, myRep = TRUE))
    expect_true(testPermFun(seq(-130L, -100L, 10L), 7, myRep = TRUE))
    expect_true(testPermFun(seq(-100L, 100L, 100L), 9, myRep = TRUE))
    expect_true(testPermFun(seq(1e10, 1e10 + 300, 100), 7, myRep = TRUE))
    expect_true(testPermFun(seq(-1e10 - 300, -1e10, 100), 7, myRep = TRUE))

    testPermMultiset <- function(v, m, frqs, verbose = FALSE,
                                 f = "sum", isExact = TRUE, my_p = F) {
        v <- sort(v)
        allSums <- permuteGeneral(v, m, freqs = frqs, constraintFun = f)
        tbl <- table(allSums[, m + 1])
        possVals <- as.numeric(names(tbl))

        if (verbose) {
            print(possVals)
            print(tbl)
        }

        ## Credit to Johan Larsson: https://stackoverflow.com/a/39175037/4408538
        is_equal_tol <- function(x, y, tol = sqrt(.Machine$double.eps)) {
            abs(x - y) < tol
        }

        t <- sapply(possVals, function(x) {
            t <- permuteGeneral(v, m, freqs = frqs,
                                constraintFun = f,
                                comparisonFun = "==",
                                limitConstraints = x)

            if (my_p) {
                print(partitionsDesign(v, m, freqs = frqs, target = x))
            }

            if (isExact) {
                u <- allSums[allSums[, m + 1] == x, 1:m]
            } else {
                u <- allSums[which(is_equal_tol(allSums[, m + 1], x)), 1:m]
            }

            if (nrow(t) > 1) {
                identical(t[do.call(order, as.data.frame(t)), ], u)
            } else {
                identical(as.vector(t), u)
            }
        })

        if (verbose) {
            print(t)
        }

        all(t)
    }

    expect_true(testPermMultiset(1:8, 5, rep(1:4, 2)))
    expect_true(testPermMultiset(0:7, 5, rep(1:4, 2)))
    expect_true(testPermMultiset(0:7, 5, rep(4:1, 2)))

    for (i in 2:4) {
        for (m in 1:(i + 5)) {
            expect_true(testPermMultiset(0:5, m, c(i, rep(1, 5))))
            expect_true(testPermMultiset(7L + 3L * 0:5, m, c(i, rep(1, 5)), my_p = F))
        }
    }

    expect_true(testPermMultiset(-3:4, 5, rep(1:4, 2)))
    expect_true(testPermMultiset((1e10 + 1):(1e10 + 8), 5, rep(1:4, 2)))
    expect_true(testPermMultiset((-1e10 - 1):(-1e10 - 8), 5, rep(1:4, 2)))
    expect_true(testPermMultiset(-49:50, 2, rep(1:2, 50)))
    expect_true(testPermMultiset(1:100, 2, rep(1:2, 50)))
    expect_true(testPermMultiset((-1e12 - 30):(-1e12 - 1), 3, rep(1:2, 15)))
    expect_true(testPermMultiset(1:4, 8, 1:4))
    expect_true(testPermMultiset(-1:1, 10, c(2, 3, 5)))

    expect_true(testPermMultiset(seq(100, 180, 10), 5, rep(1:3, 3)))
    expect_true(testPermMultiset(seq(-130L, -100L, 10L), 6, c(1, 2, 3, 2)))
    expect_true(testPermMultiset(seq(-100L, 300L, 100L), 9, c(5, 1, 1, 1, 1)))
    expect_true(testPermMultiset(seq(1e10, 1e10 + 500, 100), 10, c(1, 1, 5, 1, 1, 1)))
    expect_true(testPermMultiset(seq(-1e10 - 500, -1e10, 100), 10, c(1, 1, 1, 1, 1, 5)))

    ## Irregular vector input... i.e. the distance between neighbors varies
    ## This will test the BruteNextElem and main constraint functions
    pS <- as.integer(c(2, 3, 5, 7, 11))
    expect_true(testPermMultiset(pS, 1, frqs = 1:5))
    expect_true(testPermMultiset(pS, 2, frqs = 1:5))
    expect_true(testPermMultiset(pS, 5, frqs = 1:5))
    expect_true(testPermMultiset(pS, 7, frqs = 1:5))

    ## This is equivalent to combinations without rep
    expect_true(testPermMultiset(pS, 1, frqs = rep(1, 5)))
    expect_true(testPermMultiset(pS, 2, frqs = rep(1, 5)))
    expect_true(testPermMultiset(pS, 5, frqs = rep(1, 5)))

    pNS <- as.numeric(pS)
    expect_true(testPermMultiset(pNS, 1, frqs = 1:5))
    expect_true(testPermMultiset(pNS, 2, frqs = 1:5))
    expect_true(testPermMultiset(pNS, 5, frqs = 1:5))
    expect_true(testPermMultiset(pNS, 7, frqs = 1:5))

    expect_true(testPermMultiset(pS, 1, frqs = 1:5, f = "prod"))
    expect_true(testPermMultiset(pS, 2, frqs = 1:5, f = "prod"))
    expect_true(testPermMultiset(pS, 7, frqs = 1:5, f = "prod"))

    expect_true(testPermMultiset(pS, 1, frqs = 1:5, f = "mean", isExact = F))
    expect_true(testPermMultiset(pS, 2, frqs = 1:5, f = "mean", isExact = F))
    expect_true(testPermMultiset(pS, 5, frqs = 1:5, f = "mean", isExact = F))
    expect_true(testPermMultiset(pS, 7, frqs = 1:5, f = "mean", isExact = F))

    ## Standard compositions into k parts
    tempPerms <- permuteGeneral(10, 5, TRUE, constraintFun = "sum")
    expect_equal(nrow(tempPerms[tempPerms[,ncol(tempPerms)] == 10, 1:5]),
                 nrow(permuteGeneral(10, 5, TRUE,
                                     constraintFun = "sum",
                                     comparisonFun = "==",
                                     limitConstraints = 10)))

    ## Standard compositions
    tempPerms <- permuteGeneral(0:6, 6, TRUE, constraintFun = "sum")
    expect_equal(nrow(tempPerms[tempPerms[,ncol(tempPerms)] == 6, 1:6]),
                 nrow(permuteGeneral(0:6, 6, TRUE,
                                     constraintFun = "sum",
                                     comparisonFun = "==",
                                     limitConstraints = 6)))

    ## Testing cases where no results should be returned
    expect_equal(nrow(permuteGeneral(10, 6,
                                     constraintFun = "sum",
                                     comparisonFun = "==",
                                     limitConstraints = 20)), 0)

    expect_equal(nrow(permuteGeneral(10, 6,
                                     constraintFun = "sum",
                                     comparisonFun = "==",
                                     limitConstraints = 46)), 0)

    expect_equal(nrow(permuteGeneral(6, 10, TRUE,
                                     constraintFun = "sum",
                                     comparisonFun = "==",
                                     limitConstraints = 5)), 0)

    expect_equal(nrow(permuteGeneral(6, 10, TRUE,
                                     constraintFun = "sum",
                                     comparisonFun = "==",
                                     limitConstraints = 61)), 0)

    ## nrow(permuteGeneral(10, 6,
    ##                     constraintFun = "sum",
    ##                     comparisonFun = "==",
    ##                     limitConstraints = 30))
    ## [1] 10080
    expect_equal(nrow(permuteGeneral(10, 6,
                                     constraintFun = "sum",
                                     comparisonFun = "==",
                                     limitConstraints = 30,
                                     upper = 10)), 10)

    expect_equal(nrow(permuteGeneral(10, 6,
                                     constraintFun = "sum",
                                     comparisonFun = "==",
                                     limitConstraints = 30,
                                     upper = 10079)), 10079)

    ## nrow(permuteGeneral(4, 5, TRUE,
    ##                     constraintFun = "sum",
    ##                     comparisonFun = "==",
    ##                     limitConstraints = 10))
    ## [1] 101
    expect_equal(nrow(permuteGeneral(4, 5, TRUE,
                                     constraintFun = "sum",
                                     comparisonFun = "==",
                                     limitConstraints = 10,
                                     upper = 11)), 11)

    ## nrow(permuteGeneral(4, 5, TRUE,
    ##                     constraintFun = "sum",
    ##                     comparisonFun = "==",
    ##                     limitConstraints = 8))
    ## [1] 35
    expect_equal(nrow(permuteGeneral(4, 5, TRUE,
                                     constraintFun = "sum",
                                     comparisonFun = "==",
                                     limitConstraints = 8,
                                     upper = 34)), 34)

    ## nrow(permuteGeneral(0:10, 7, TRUE,
    ##                    constraintFun = "sum",
    ##                     comparisonFun = "==",
    ##                     limitConstraints = 10))
    ## [1] 8008
    expect_equal(nrow(permuteGeneral(0:10, 7, TRUE,
                                     constraintFun = "sum",
                                     comparisonFun = "==",
                                     limitConstraints = 10,
                                     upper = 10)), 10)
})
