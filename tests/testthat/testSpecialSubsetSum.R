context("testing special subset sum")

test_that("comboGeneral produces correct results for special subset sum", {
    # testing comboGeneral and permuteGeneral results for constrainFun = 'sum', 
    # comparisonFun = '==', and with v having the property that if you were to
    # sort v, the difference of each element with it's neighbor is constant.
    
    testFun <- function(v, m, myRep = FALSE, verbose = FALSE, f = "sum", isExact = TRUE) {
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
            t <- comboGeneral(v, m, myRep,
                              constraintFun = f, 
                              comparisonFun = "==",
                              limitConstraints = x)
            if (isExact)
                u <- allSums[allSums[, m + 1] == x, 1:m]
            else
                u <- allSums[which(is_equal_tol(allSums[, m + 1], x)), 1:m]
            
            if (nrow(t) > 1) {
                identical(t, u)
            } else {
                identical(as.vector(t), u)
            }
        })
        
        if (verbose)
            print(t)
        
        all(t)
    }
    
    expect_true(testFun(1:18, 9))
    expect_true(testFun(0:17, 9))
    expect_true(testFun(-8:8, 9))
    expect_true(testFun(1:5, 5))
    expect_true(testFun((1e10 + 1):(1e10 + 18), 9))
    expect_true(testFun((-1e10 - 1):(-1e10 - 18), 9))
    expect_true(testFun(-49:50, 99))
    expect_true(testFun(1:100, 99))
    expect_true(testFun(-49:50, 2))
    expect_true(testFun(1:100, 2))
    expect_true(testFun(1:100, 100))
    expect_true(testFun((-1e12 - 50):(-1e12 - 1), 49))
    expect_true(testFun((-1e12 - 50):(-1e12 - 1), 3))
    
    expect_true(testFun(1:10, 7, myRep = TRUE))
    expect_true(testFun(0:9, 7, myRep = TRUE))
    expect_true(testFun(-4:5, 7, myRep = TRUE))
    expect_true(testFun((1e10 + 1):(1e10 + 10), 7, myRep = TRUE))
    expect_true(testFun((-1e10 - 1):(-1e10 - 10), 7, myRep = TRUE))
    expect_true(testFun(-49:50, 2, myRep = TRUE))
    expect_true(testFun(1:100, 2, myRep = TRUE))
    expect_true(testFun((-1e12 - 50):(-1e12 - 1), 3, myRep = TRUE))
    expect_true(testFun(1:5, 10, myRep = TRUE))
    expect_true(testFun(-1:1, 100, myRep = TRUE))
    
    expect_true(testFun(seq(100, 180, 5), 9))
    expect_true(testFun(seq(-80, 80, 10), 9))
    expect_true(testFun(seq(1e10, 1e10 + 180, 10), 9))
    
    expect_true(testFun(seq(100, 210, 10), 7, myRep = TRUE))
    expect_true(testFun(seq(-140L, -100L, 10L), 10, myRep = TRUE))
    expect_true(testFun(seq(-100L, 300L, 100L), 10, myRep = TRUE))
    expect_true(testFun(seq(1e10, 1e10 + 500, 100), 10, myRep = TRUE))
    expect_true(testFun(seq(-1e10 - 500, -1e10, 100), 10, myRep = TRUE))
    
    ## Irregular vector input... i.e. the distance between neighbors varies
    ## This will test the BruteNextElem and main constraint functions
    p <- as.integer(c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41))
    expect_true(testFun(p, 1))
    expect_true(testFun(p, 2))
    expect_true(testFun(p, 7))
    expect_true(testFun(p, 13))
    
    pN <- as.numeric(p)
    expect_true(testFun(pN, 1))
    expect_true(testFun(pN, 2))
    expect_true(testFun(pN, 7))
    expect_true(testFun(pN, 13))
    
    expect_true(testFun(p, 1, f = "prod"))
    expect_true(testFun(p, 2, f = "prod"))
    expect_true(testFun(p, 6, f = "prod"))
    
    expect_true(testFun(p, 1, f = "mean", isExact = F))
    expect_true(testFun(p, 2, f = "mean", isExact = F))
    expect_true(testFun(p, 7, f = "mean", isExact = F))
    expect_true(testFun(p, 13, f = "mean", isExact = F))
    
    pS <- p[1:6]
    expect_true(testFun(pS, 1, T))
    expect_true(testFun(pS, 2, T))
    expect_true(testFun(pS, 6, T))
    expect_true(testFun(pS, 8, T))
    
    pNS <- as.numeric(pS)
    expect_true(testFun(pNS, 1, T))
    expect_true(testFun(pNS, 2, T))
    expect_true(testFun(pNS, 6, T))
    expect_true(testFun(pNS, 8, T))
    
    expect_true(testFun(pS, 1, T, f = "prod"))
    expect_true(testFun(pS, 2, T, f = "prod"))
    expect_true(testFun(pS, 6, T, f = "prod"))
    expect_true(testFun(pS, 7, T, f = "prod"))
    
    expect_true(testFun(pS, 1, T, f = "mean", isExact = F))
    expect_true(testFun(pS, 2, T, f = "mean", isExact = F))
    expect_true(testFun(pS, 6, T, f = "mean", isExact = F))
    expect_true(testFun(pS, 8, T, f = "mean", isExact = F))
    
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
    
    testMultiset <- function(v, m, frqs, verbose = FALSE, f = "sum", isExact = TRUE) {
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
            if (isExact)
                u <- allSums[allSums[, m + 1] == x, 1:m]
            else
                u <- allSums[which(is_equal_tol(allSums[, m + 1], x)), 1:m]
            
            if (nrow(t) > 1) {
                identical(t, u)
            } else {
                identical(as.vector(t), u)
            }
        })
        
        if (verbose)
            print(t)
        
        all(t)
    }
    
    expect_true(testMultiset(1:10, 7, rep(1:5, 2)))
    expect_true(testMultiset(0:9, 7, rep(1:5, 2)))
    expect_true(testMultiset(-4:5, 7, rep(1:5, 2)))
    expect_true(testMultiset((1e10 + 1):(1e10 + 10), 7, rep(1:5, 2)))
    expect_true(testMultiset((-1e10 - 1):(-1e10 - 10), 7, rep(1:5, 2)))
    expect_true(testMultiset(-49:50, 2, rep(1:2, 50)))
    expect_true(testMultiset(1:100, 2, rep(1:2, 50)))
    expect_true(testMultiset((-1e12 - 50):(-1e12 - 1), 3, rep(1:2, 25)))
    expect_true(testMultiset(1:5, 10, 1:5))
    expect_true(testMultiset(-1:1, 100, c(20, 30, 50)))
    
    expect_true(testMultiset(seq(100, 210, 10), 7, rep(1:4, 3)))
    expect_true(testMultiset(seq(-140L, -100L, 10L), 10, c(1, 2, 3, 4, 3)))
    expect_true(testMultiset(seq(-100L, 300L, 100L), 9, c(5, 1, 1, 1, 1)))
    expect_true(testMultiset(seq(1e10, 1e10 + 500, 100), 10, c(1, 1, 5, 1, 1, 1)))
    expect_true(testMultiset(seq(-1e10 - 500, -1e10, 100), 10, c(1, 1, 1, 1, 1, 5)))
    
    ## Irregular vector input... i.e. the distance between neighbors varies
    ## This will test the BruteNextElem and main constraint functions
    pS <- as.integer(c(2, 3, 5, 7, 11, 13))
    expect_true(testMultiset(pS, 1, frqs = 1:6))
    expect_true(testMultiset(pS, 2, frqs = 1:6))
    expect_true(testMultiset(pS, 6, frqs = 1:6))
    expect_true(testMultiset(pS, 8, frqs = 1:6))
    
    ## This is equivalent to combinations without rep
    expect_true(testMultiset(pS, 1, frqs = rep(1, 6)))
    expect_true(testMultiset(pS, 2, frqs = rep(1, 6)))
    expect_true(testMultiset(pS, 6, frqs = rep(1, 6)))
    
    pNS <- as.numeric(pS)
    expect_true(testMultiset(pNS, 1, frqs = 1:6))
    expect_true(testMultiset(pNS, 2, frqs = 1:6))
    expect_true(testMultiset(pNS, 6, frqs = 1:6))
    expect_true(testMultiset(pNS, 8, frqs = 1:6))
    
    expect_true(testMultiset(pS, 1, frqs = 1:6, f = "prod"))
    expect_true(testMultiset(pS, 2, frqs = 1:6, f = "prod"))
    expect_true(testMultiset(pS, 7, frqs = 1:6, f = "prod"))
    
    expect_true(testMultiset(pS, 1, frqs = 1:6, f = "mean", isExact = F))
    expect_true(testMultiset(pS, 2, frqs = 1:6, f = "mean", isExact = F))
    expect_true(testMultiset(pS, 6, frqs = 1:6, f = "mean", isExact = F))
    expect_true(testMultiset(pS, 8, frqs = 1:6, f = "mean", isExact = F))
    
    
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
                                     comparisonFun = "==", limitConstraints = 10)), 48)
    
    expect_equal(nrow(permuteGeneral(0:10, 4, constraintFun = "sum",
                                     comparisonFun = "==", limitConstraints = 10)), 120)
    
    expect_equal(nrow(permuteGeneral(10, 3, constraintFun = "sum",
                                     comparisonFun = "==", limitConstraints = 10)), 24)
    
    expect_equal(nrow(permuteGeneral(0:10, repetition = TRUE, constraintFun = "sum",
                                     comparisonFun = "==", limitConstraints = 10)), 512)
    
    expect_equal(nrow(permuteGeneral(0:10, 10, repetition = TRUE, constraintFun = "sum",
                                     comparisonFun = "==", limitConstraints = 10)), 92378)
})

test_that("permuteGeneral produces correct results for special subset sum", {
    testFun <- function(v, m, myRep = FALSE, verbose = FALSE) {
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
    
    expect_true(testFun(1:8, 6))
    expect_true(testFun(0:7, 6))
    expect_true(testFun(-4:3, 6))
    expect_true(testFun(1:5, 5))
    expect_true(testFun((1e10 + 1):(1e10 + 8), 6))
    expect_true(testFun((-1e10 - 1):(-1e10 - 8), 6))
    expect_true(testFun(-49:50, 2))
    expect_true(testFun(1:100, 2))
    expect_true(testFun((-1e12 - 30):(-1e12 - 1), 3))
    
    expect_true(testFun(1:7, 5, myRep = TRUE))
    expect_true(testFun(0:6, 5, myRep = TRUE))
    expect_true(testFun(-3:3, 5, myRep = TRUE))
    expect_true(testFun((1e10 + 1):(1e10 + 7), 5, myRep = TRUE))
    expect_true(testFun((-1e10 - 1):(-1e10 - 7), 5, myRep = TRUE))
    expect_true(testFun(-49:50, 2, myRep = TRUE))
    expect_true(testFun(1:100, 2, myRep = TRUE))
    expect_true(testFun((-1e12 - 50):(-1e12 - 1), 3, myRep = TRUE))
    expect_true(testFun(1:4, 7, myRep = TRUE))
    expect_true(testFun(-1:1, 9, myRep = TRUE))
    
    expect_true(testFun(seq(100, 135, 5), 6))
    expect_true(testFun(seq(-40, 30, 10), 6))
    expect_true(testFun(seq(1e10, 1e10 + 80, 10), 6))
    
    expect_true(testFun(seq(100, 160, 10), 5, myRep = TRUE))
    expect_true(testFun(seq(-130L, -100L, 10L), 7, myRep = TRUE))
    expect_true(testFun(seq(-100L, 100L, 100L), 9, myRep = TRUE))
    expect_true(testFun(seq(1e10, 1e10 + 300, 100), 7, myRep = TRUE))
    expect_true(testFun(seq(-1e10 - 300, -1e10, 100), 7, myRep = TRUE))
    
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
