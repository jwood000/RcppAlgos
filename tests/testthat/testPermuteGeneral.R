context("testing permuteGeneral")

test_that("permuteGeneral produces correct results with no constraints and no repetition", {
    expect_equal(nrow(permuteGeneral(3, 3)), 6)
    expect_equal(as.vector(permuteGeneral(1,1)), 1)

    expect_equal(permuteGeneral(10, 5)[500:600, ], permuteGeneral(10, 5,
                                                              lower = 500,
                                                              upper = 600))

    expect_equal(permuteGeneral(month.abb[1:7], 7)[4000:5000, ],
                 permuteGeneral(month.abb[1:7], 7, lower = 4000, upper = 5000))

    expect_equal(as.vector(permuteGeneral(100,1)), 1:100)

    ## Constraint should not be carried out if no comparisonFun is given
    expect_equal(permuteGeneral(3, 3, constraintFun = "sum",
                                limitConstraints = 100), permuteGeneral(1:3, 3))
    
    expect_equal(nrow(permuteGeneral(8, 8)), factorial(8))

    set.seed(11)
    myNums <- rnorm(5)
    expect_equal(permuteGeneral(myNums, 3),
                 permuteGeneral(myNums, 3, freqs = rep(1, 5)))

    expect_equal(ncol(permuteGeneral(5, 3)), 3)
    expect_equal(ncol(permuteGeneral(5, 3, FALSE, constraintFun = "prod", keepResults = TRUE)), 4)
    expect_equal(nrow(permuteGeneral(5, 3, upper = 20)), 20)
    expect_equal(nrow(permuteGeneral(5, 3, FALSE, constraintFun = "prod",
                                     keepResults = TRUE, upper = 10L)), 10)
})

test_that("permuteGeneral produces correct results with no constraints and has repetition", {
    expect_equal(as.vector(permuteGeneral(1,1,TRUE)), 1)
    expect_equal(as.vector(permuteGeneral(1,5,TRUE)), rep(1, 5))
    expect_equal(permuteGeneral(letters[1:9], 5, TRUE)[59000:permuteCount(9,5,T), ],
                 permuteGeneral(letters[1:9], 5, TRUE, lower = 59000L))
    expect_equal(ncol(permuteGeneral(5, 3, TRUE)), 3)
    expect_equal(nrow(permuteGeneral(2, 2, TRUE)), 4)
    expect_equal(nrow(permuteGeneral(5, 3, TRUE, constraintFun = "prod", upper = 10)), 10)

    set.seed(111)
    myNums <- rnorm(5)
    expect_equal(permuteGeneral(myNums, 3, TRUE),
                 permuteGeneral(myNums, 3, freqs = rep(3, 5)))

    expect_true(all(permuteGeneral(3, 3, TRUE) ==
                        as.matrix(expand.grid(1:3, 1:3, 1:3))[,3:1]))

    expect_equal(nrow(permuteGeneral(5, 3, TRUE, upper = 10)), 10)
    expect_equal(ncol(permuteGeneral(5, 3, TRUE, constraintFun = "prod", keepResults = TRUE)), 4)
    
    ## In older versions the test below would fail b/c it would produce NaNs during prep
    expect_equal(nrow(permuteGeneral(2, 180, freqs = c(180, 2))), 
                 permuteCount(2, 180, freqs = c(180, 2)))
})

test_that("permuteGeneral produces correct results with no constraints for multisets", {
    expect_equal(nrow(permuteGeneral(5, 5, freqs = 1:5, upper = 10)), 10)
    expect_equal(ncol(permuteGeneral(5, 3, FALSE,
                                     constraintFun = "prod", freqs = c(1,2,1,2,4),
                                     keepResults = TRUE)), 4)
    expect_equal(as.vector(permuteGeneral(1, 2, freqs = 2)), c(1, 1))

    expect_equal(permuteGeneral(LETTERS[1:5], 3),
                 permuteGeneral(LETTERS[1:5], 3, freqs = rep(1, 5)))

    expect_equal(permuteGeneral(month.name[1:5], 3, TRUE),
                 permuteGeneral(month.name[1:5], 3, freqs = rep(3, 5)))

    myNums2 <- 1:10 / 3
    expect_equal(permuteGeneral(myNums2, 5, freqs = rep(2, 10))[80000:90000, ],
                 permuteGeneral(myNums2, 5, freqs = rep(2, 10),
                                lower = 80000, upper = 90000))

    expect_equal(sum(permuteGeneral(3, 3, freqs = c(1, 1, 1))),
                 sum(permuteGeneral(3, 3)))

    expect_equal(permuteGeneral(factor(1:5, ordered = TRUE), 5, freqs = rep(3, 5)),
                 permuteSample(factor(1:5, ordered = TRUE), 5, freqs = rep(3, 5),
                               sampleVec = 1:permuteCount(5, 5, freqs = rep(3, 5))))

    expect_equal(permuteGeneral(5, 5),
                 permuteGeneral(5, 5, freqs = rep(1, 5)))
    
    expect_equal(permuteCount(30, freqs = rep(1:2, 15)),
                 permuteCount(30, 45, freqs = rep(1:2, 15)))

    expect_equal(do.call(rbind, lapply(seq(1L, 1680, 168), function(x) {
        permuteGeneral(4, freqs = c(2,1,2,3), lower = x, upper = x+167)
    })), permuteGeneral(4, 20, freqs = c(2,1,2,3)))
})

test_that("permuteGeneral produces correct results with constraints", {
    expect_equal(nrow(permuteGeneral(15, 7,
                                     comparisonFun = "==", constraintFun = "sum",
                                     limitConstraints = 80, upper = 100)), 100)
    
    expect_equal(nrow(permuteGeneral(15, 7, TRUE, 
                                     comparisonFun = "==", constraintFun = "sum",
                                     limitConstraints = 80, upper = 200)), 200)
    
    expect_equal(nrow(permuteGeneral(15, 7, freqs = rep(1:5, 3), 
                                     comparisonFun = "==", constraintFun = "sum",
                                     limitConstraints = 80, upper = 220)), 220)

    expect_equal(nrow(permuteGeneral(3, 3, FALSE, constraintFun = "sum",
                                     comparisonFun = "==", limitConstraints = 6)), 6)

    ## NA should be removed when constraint check is carried out
    expect_equal(unique(permuteGeneral(c(NA,1:5), 5, TRUE,
                                       constraintFun = "sum",
                                       comparisonFun = "==", limitConstraints = 9,
                                       keepResults = TRUE)[,6]), 9)

    expect_true(all(permuteGeneral(5, 5L, TRUE,
                                   constraintFun = "min",
                                   comparisonFun = "<", limitConstraints = 3,
                                       keepResults = TRUE)[,6] < 3))

    expect_true(all(permuteGeneral(5, 5, TRUE,
                                   constraintFun = "prod",
                                   comparisonFun = ">", limitConstraints = 100,
                                   keepResults = TRUE)[,6] > 100))

    expect_true(all(permuteGeneral(5, 3, FALSE,
                                   constraintFun = "max",
                                   comparisonFun = "=<", limitConstraints = 4,
                                keepResults = TRUE)[,4] <= 4))

    expect_true(all(permuteGeneral(3, 5, TRUE,
                                   constraintFun = "mean",
                                   comparisonFun = ">=", limitConstraints = 2,
                                   keepResults = TRUE)[,6] >= 2))

    expect_true(all(permuteGeneral(5, 5, FALSE, constraintFun = "sum",
                                   comparisonFun = ">", limitConstraints = 18,
                                     freqs = c(1,2,1,2,4),
                                     keepResults = TRUE)[,6] > 18))

    expect_equal(sum(permuteGeneral(4, 6, freqs = c(1,3,2,2),
                                   constraintFun = "sum",
                                   keepResults = TRUE)[,7] < 16),
                sum(apply(permuteGeneral(c(3,1,4,2), 6, freqs = c(2,1,2,3),
                                         constraintFun = "sum",
                                         comparisonFun = "<",
                                         limitConstraints = 16), 1, sum) < 16))

    expect_true(min(permuteGeneral(15, 5, constraintFun = "prod",
                                   comparisonFun = ">",
                                   limitConstraints = 3*10^5,
                                   keepResults = TRUE)[,6]) > 3*10^5)

    a <- permuteGeneral(8, 5, freqs = rep(3, 8))
    b <- apply(a, 1, min)
    expect_equal(permuteGeneral(8, 5, freqs = rep(3, 8),
                              constraintFun = "min",
                              comparisonFun = "==",
                              limitConstraints = 3,
                              lower = 17900, upper = 18500), a[(17900:18500)[b[17900:18500] == 3], ])

    a <- permuteGeneral(c(-1L,1:5), 7, T)
    b <- apply(a, 1, prod)
    expect_equal(nrow(permuteGeneral(c(-1L,1:5), 7, TRUE, constraintFun = "prod",
                              comparisonFun = c(">=","<="),
                              limitConstraints = c(2000, 5000))),
                 nrow(a[which(b >= 2000 & b <= 5000), ]))
})

test_that("permuteGeneral produces correct results with exotic constraints", {

    comp1 <- c("<", "<=")
    comp2 <- c(">", ">=")
    allPerms <- permuteGeneral(c(-6:(-1),1:2), 5, freqs = c(rep(1:3, 2), 2:3),
                              constraintFun = "prod", keepResults = TRUE)

    theSum <- allPerms[, 6]
    allPerms <- allPerms[, 1:5]
    q <- quantile(theSum)

    for (i in 1:2) {

        if (i == 1) {
            a <- comp1
            b <- comp2
        } else {
            a <- comp2
            b <- comp1
        }

        for (j in a) {
            for (k in b) {
                myComp <- c(j, k)
                myTest <- permuteGeneral(c(-6:(-1),1:2), 5, freqs = c(rep(1:3, 2), 2:3),
                                       constraintFun = "prod", comparisonFun = myComp,
                                       limitConstraints = c(q[2], q[4]))
                fun1 <- match.fun(j)
                fun2 <- match.fun(k)

                if (i == 1) {
                    temp <- allPerms[fun1(theSum, q[2]) | fun2(theSum, q[4]),]
                } else {
                    temp <- allPerms[fun1(theSum, q[2]) & fun2(theSum, q[4]),]
                }

                expect_equal(temp, myTest)
            }
        }
    }

    allPerms <- permuteGeneral(c(-6:(-1),1:2), 5, TRUE,
                              constraintFun = "prod", keepResults = TRUE)

    theSum <- allPerms[, 6]
    allPerms <- allPerms[, 1:5]
    q <- quantile(theSum)

    for (i in 1:2) {

        if (i == 1) {
            a <- comp1
            b <- comp2
        } else {
            a <- comp2
            b <- comp1
        }

        for (j in a) {
            for (k in b) {
                myComp <- c(j, k)
                myTest <- permuteGeneral(c(-6:(-1),1:2), 5, TRUE,
                                       constraintFun = "prod", comparisonFun = myComp,
                                       limitConstraints = c(q[2], q[4]))
                fun1 <- match.fun(j)
                fun2 <- match.fun(k)

                if (i == 1) {
                    temp <- allPerms[fun1(theSum, q[2]) | fun2(theSum, q[4]),]
                } else {
                    temp <- allPerms[fun1(theSum, q[2]) & fun2(theSum, q[4]),]
                }

                expect_equal(temp, myTest)
            }
        }
    }

    allPerms <- permuteGeneral(c(-6:(-1),1:4), 5,
                              constraintFun = "prod", keepResults = TRUE)

    theSum <- allPerms[, 6]
    allPerms <- allPerms[, 1:5]
    q <- quantile(theSum)

    for (i in 1:2) {

        if (i == 1) {
            a <- comp1
            b <- comp2
        } else {
            a <- comp2
            b <- comp1
        }

        for (j in a) {
            for (k in b) {
                myComp <- c(j, k)
                myTest <- permuteGeneral(c(-6:(-1),1:4), 5,
                                       constraintFun = "prod", comparisonFun = myComp,
                                       limitConstraints = c(q[2], q[4]))
                fun1 <- match.fun(j)
                fun2 <- match.fun(k)

                if (i == 1) {
                    temp <- allPerms[fun1(theSum, q[2]) | fun2(theSum, q[4]),]
                } else {
                    temp <- allPerms[fun1(theSum, q[2]) & fun2(theSum, q[4]),]
                }
                expect_equal(temp, myTest)
            }
        }
    }
})

test_that("permuteGeneral produces correct results with use of FUN", {
    test <- permuteGeneral(6, 6, constraintFun = "mean")[, 7]
    expect_equal(as.vector(test), unlist(permuteGeneral(6, 6, FUN = mean)))
    
    expect_equal(sum(unlist(permuteGeneral(as.complex(c(1, -1, -1i, 1i)), 3,
                                           FUN = function(x) sum(Re(x))))), 0)
    
    test <- permuteGeneral(6, 6, lower = 100, constraintFun = "prod")[, 7]
    expect_equal(as.vector(test), unlist(permuteGeneral(6, 6, lower = 100, FUN = prod)))
    
    test <- permuteGeneral(10, 5, constraintFun = "sum", keepResults = TRUE)
    expect_equal(as.vector(test[,6]), unlist(permuteGeneral(10, 5, FUN = sum)))

    test <- permuteGeneral(10, 4, TRUE)
    testFun <- apply(test, 1, function(x) mean(x) * 2)
    expect_equal(testFun, unlist(permuteGeneral(10, 4, T, FUN = function(x) {mean(x) * 2})))

    test <- permuteGeneral(8, 4, freqs = rep(1:4, 2))
    testFun <- lapply(1:nrow(test), function(x) cumsum(test[x, ]))
    expect_equal(testFun, permuteGeneral(8, 4, freqs = rep(1:4, 2), FUN = cumsum))
    
    test <- apply(permuteGeneral(4, 8, freqs = c(1,3,1,3)), 1, {
        function(x) paste0(cumprod(x), collapse = "")
    })
    
    testFun <- unlist(permuteGeneral(4, 8, freqs = c(1,3,1,3), FUN = function(x) {
        paste0(cumprod(x), collapse = "")
    }))
    
    expect_equal(test, testFun)
    
    expect_equal(permuteGeneral(4, 8, freqs = c(1,3,1,3), 
                                constraintFun = "sum", nThreads = 2)[, 9], 
                 rowSums(permuteGeneral(4, 8, freqs = c(1,3,1,3), nThreads = 2)))
})

test_that("permuteGeneral produces correct results with very large results", {
    ##******** BIG TESTS *********##

    ## NO REPETITION
    numR <- permuteCount(1000, 10)
    n1 <- gmp::sub.bigz(numR, 99)

    ## accepts raw values
    expect_equal(nrow(permuteGeneral(1000, 10, lower = n1)), 100)
    ## accepts characters
    expect_equal(nrow(permuteGeneral(1000, 10, lower = as.character(n1))), 100)
    expect_equal(as.vector(permuteGeneral(1000, 10, lower = numR)),
                 1000:991)

    ## WITH REPETITION
    numR <- permuteCount(1000, 10, TRUE)
    n1 <- gmp::sub.bigz(numR, 99)

    expect_equal(nrow(permuteGeneral(1000, 10, TRUE, lower = n1)), 100)
    expect_equal(nrow(permuteGeneral(1000, 10, TRUE, lower = as.character(n1))), 100)
    expect_equal(as.vector(permuteGeneral(1000, 10, TRUE, lower = numR)),
                 rep(1000, 10))

    ## MULTISETS
    numR <- permuteCount(100, 10, freqs = rep(1:4, 25))
    n1 <- gmp::sub.bigz(numR, 99)

    expect_equal(nrow(permuteGeneral(100, 10, freqs = rep(1:4, 25), lower = n1)), 100)
    expect_equal(nrow(permuteGeneral(100, 10, freqs = rep(1:4, 25), lower = as.character(n1))), 100)
    expect_equal(as.vector(permuteGeneral(100, 10, freqs = rep(1:4, 25), lower = numR)),
                 rep(100:97, times = 4:1))
})
