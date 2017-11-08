library(RcppAlgos)
library(gtools)
context("testing permuteGeneral with constraints")

test_that("permuteGeneral produces correct permutations with constraints and no repetition", {
    WithConstNoRepNum <- function() {
        funs <- c("sum", "prod", "mean", "min", "max")
        comps <- c("<", "<=", ">", ">=", "==")
        allCombs <- permutations(5,3)
        
        results <- sapply(funs, function(f) {
            myfun <- match.fun(FUN = f)
            t <- myfun(2:4)
            sapply(comps, function(e) {
                operator <- match.fun(FUN = e)
                myPermute <- permuteGeneral(5,3,FALSE,f,e,t)
                myPermute <- myPermute[order(myPermute[,1], myPermute[,2], myPermute[,3]), ]
                if (nrow(myPermute) == 1) {myPermute <- as.vector(myPermute)}
                testPermute <- allCombs[apply(allCombs, 1, function(x) operator(myfun(x), t)), ]
                testPermute <- testPermute[order(testPermute[,1], testPermute[,2], testPermute[,3]), ]
                all.equal(myPermute, testPermute)
            })
        })
        
        all(results)
    }
    expect_true(WithConstNoRepNum())
    expect_equal(nrow(permuteGeneral(100,3,FALSE,"sum","<",100,10)), 10)
})

test_that("permuteGeneral produces correct permutations with constraints and with repetition", {
    WithConstWithRepNum <- function() {
        funs <- c("sum", "prod", "mean", "min", "max")
        comps <- c("<", "<=", ">", ">=", "==")
        allCombs <- gtools::permutations(5,3,1:5,FALSE,TRUE)
        
        results <- sapply(funs, function(f) {
            myfun <- match.fun(FUN = f)
            t <- myfun(2:4)
            sapply(comps, function(e) {
                operator <- match.fun(FUN = e)
                myPermute <- permuteGeneral(5,3,TRUE,f,e,t)
                if (nrow(myPermute) == 1) {myPermute <- as.vector(myPermute)}
                myPermute <- myPermute[order(myPermute[,1], myPermute[,2], myPermute[,3]), ]
                if (nrow(myPermute) == 1) {myPermute <- as.vector(myPermute)}
                testPermute <- allCombs[apply(allCombs, 1, function(x) operator(myfun(x), t)), ]
                testPermute <- testPermute[order(testPermute[,1], testPermute[,2], testPermute[,3]), ]
                all.equal(myPermute, testPermute)
            })
        })
        
        all(results)
    }
    expect_true(WithConstWithRepNum(), TRUE)
    expect_equal(nrow(permuteGeneral(100,3,TRUE,"sum","<",100,10)), 10)
})

test_that("permuteGeneral produces appropriate error messages", {
    expect_error(permuteGeneral(9,4,TRUE,"summ","<",10), "prod, sum, mean, max, or min")
    expect_error(permuteGeneral(9,4,TRUE,"sum","=<>",10), ">, >=, <, <=, or ==")
    expect_error(permuteGeneral(9,4,TRUE,"sum",60,10), "must be passed as a character")
    expect_error(permuteGeneral(9,4,FALSE,sum,"<",10), "must be passed as a character")
    expect_error(permuteGeneral(9,4,TRUE,"sum","<",10,-1), "must be positive")
    expect_error(permuteGeneral(170,7,FALSE,"sum","<",100), "The number of rows cannot exceed")
    expect_error(permuteGeneral(170,7,FALSE,"sum","<",100,10^10), "The number of rows cannot exceed")
})