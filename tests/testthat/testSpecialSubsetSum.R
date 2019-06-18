context("testing special subset sum")

test_that("comboGeneral produces correct results for special subset sum", {
    # testing comboGeneral and permuteGeneral results for constrainFun = 'sum', 
    # comparisonFun = '==', and with v having the property that if you were to
    # sort v, the difference of each element with it's neighbor is constant.
    
    testFun <- function(v, m, myRep = FALSE, verbose = FALSE) {
        v <- sort(v)
        allSums <- comboGeneral(v, m, myRep, constraintFun = "sum")
        n <- length(v)
        tbl <- table(allSums[,(m+1)])
        possVals <- as.numeric(names(tbl))
        
        if (verbose) {
            print(possVals)
            print(tbl)
        }
        
        t <- sapply(possVals, function(x) {
            t <- comboGeneral(v, m, myRep,
                              constraintFun = "sum", 
                              comparisonFun = "==",
                              limitConstraints = x)
            u <- allSums[allSums[, m + 1] == x, 1:m]
            
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
    
    ## Testing cases where no results should be returned
    expect_equal(nrow(comboGeneral(10, 6, constraintFun = "sum", comparisonFun = "==", limitConstraints = 20)), 0)
    expect_equal(nrow(comboGeneral(10, 6, constraintFun = "sum", comparisonFun = "==", limitConstraints = 46)), 0)
    expect_equal(nrow(comboGeneral(6, 10, TRUE, constraintFun = "sum", comparisonFun = "==", limitConstraints = 5)), 0)
    expect_equal(nrow(comboGeneral(6, 10, TRUE, constraintFun = "sum", comparisonFun = "==", limitConstraints = 61)), 0)
})

test_that("permuteGeneral produces correct results for special subset sum", {
    testFun <- function(v, m, myRep = FALSE, verbose = FALSE) {
        v <- sort(v)
        allSums <- permuteGeneral(v, m, myRep, constraintFun = "sum")
        n <- length(v)
        tbl <- table(allSums[,(m+1)])
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
    
    ## Testing cases where no results should be returned
    expect_equal(nrow(permuteGeneral(10, 6, constraintFun = "sum", comparisonFun = "==", limitConstraints = 20)), 0)
    expect_equal(nrow(permuteGeneral(10, 6, constraintFun = "sum", comparisonFun = "==", limitConstraints = 46)), 0)
    expect_equal(nrow(permuteGeneral(6, 10, TRUE, constraintFun = "sum", comparisonFun = "==", limitConstraints = 5)), 0)
    expect_equal(nrow(permuteGeneral(6, 10, TRUE, constraintFun = "sum", comparisonFun = "==", limitConstraints = 61)), 0)
})
