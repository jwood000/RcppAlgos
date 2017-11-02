library(RcppAlgos)
context("testing comboGeneral with constraints")

test_that("comboGeneral produces correct combinations with constraints and no repetition", {
    WithConstNoRepNum <- function() {
        funs <- c("sum", "prod", "mean", "min", "max")
        comps <- c("<", "<=", ">", ">=", "==")
        allCombs1 <- t(combn(5,3))
        allCombs2 <- t(combn(5:1,3))
        
        results <- sapply(funs, function(f) {
            myfun <- match.fun(FUN = f)
            t <- myfun(2:4)
            sapply(comps, function(e) {
                operator <- match.fun(FUN = e)
                myCombo <- comboGeneral(5,3,FALSE,f,e,t)
                if (nrow(myCombo) == 1) {myCombo <- as.vector(myCombo)}
                if (e %in% c(">", ">=")) {
                    testCombs <- allCombs2[apply(allCombs2, 1, function(x) operator(myfun(x), t)), ]
                } else {
                    testCombs <- allCombs1[apply(allCombs1, 1, function(x) operator(myfun(x), t)), ]
                }
                all.equal(myCombo, testCombs)
            })
        })
        
        all(results)
    }
    expect_true(WithConstNoRepNum())
    expect_equal(nrow(comboGeneral(100,3,FALSE,"sum","<",100,10)), 10)
})

test_that("comboGeneral produces correct combinations with constraints and with repetition", {
    WithConstWithRepNum <- function() {
        funs <- c("sum", "prod", "mean", "min", "max")
        comps <- c("<", "<=", ">", ">=", "==")
        allCombs1 <- gtools::combinations(5,3,1:5,FALSE,TRUE)
        allCombs2 <- gtools::combinations(5,3,5:1,FALSE,TRUE)
        
        results <- sapply(funs, function(f) {
            myfun <- match.fun(FUN = f)
            t <- myfun(2:4)
            sapply(comps, function(e) {
                operator <- match.fun(FUN = e)
                myCombo <- comboGeneral(5,3,TRUE,f,e,t)
                if (nrow(myCombo) == 1) {myCombo <- as.vector(myCombo)}
                if (e %in% c(">", ">=")) {
                    testCombs <- allCombs2[apply(allCombs2, 1, function(x) operator(myfun(x), t)), ]
                } else {
                    testCombs <- allCombs1[apply(allCombs1, 1, function(x) operator(myfun(x), t)), ]
                }
                all.equal(myCombo, testCombs)
            })
        })
        
        all(results)
    }
    expect_true(WithConstWithRepNum(), TRUE)
    expect_equal(nrow(comboGeneral(100,3,TRUE,"sum","<",100,10)), 10)
})

test_that("comboGeneral produces appropriate error messages", {
    expect_error(comboGeneral(9,4,TRUE,"summ","<",10), "prod, sum, mean, max, or min")
    expect_error(comboGeneral(9,4,TRUE,"sum","=<",10), ">, >=, <, <=, or ==")
    expect_error(comboGeneral(9,4,TRUE,"sum",60,10), "must be passed as a character")
    expect_error(comboGeneral(9,4,FALSE,sum,"<",10), "must be passed as a character")
    expect_error(comboGeneral(9,4,TRUE,"sum","<",10,-1), "must be positive")
    expect_error(comboGeneral(170,7,FALSE,"sum","<",100), "The number of rows cannot exceed")
    expect_error(comboGeneral(170,7,FALSE,"sum","<",100,10^10), "The number of rows cannot exceed")
})