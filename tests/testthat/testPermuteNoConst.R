library(testthat)
library(RcppAlgos)
library(gtools)
context("testing permuteGeneral no constraints")

test_that("permuteGeneral produces correct combinations with no constraints and no repetition", {
    NoConstNoRepNum1 <- function(n) {
        results <- do.call(c, sapply(1:n, function(x) {
            sapply(1:x, function(y) {
                    t1 <- permuteGeneral(x, y)
                    t1 <- t1[do.call(order,as.data.frame(t1)),]
                    t2 <- permutations(x, y)
                    t2 <- t2[do.call(order,as.data.frame(t2)),]
                    all.equal(t1,t2)
                })
        }))
        all(results)
    }
    NoConstNoRepNum2 <- function(n) {
        set.seed(101)
        results <- do.call(c, sapply(2:n, function(x) {
            myV <- sort(runif(x, sample(-100:0, 1), sample(1:100, 1)))
            sapply(1:x, function(y) {
                t1 <- permuteGeneral(myV, y)
                t1 <- t1[do.call(order,as.data.frame(t1)),]
                t2 <- permutations(x, y, myV)
                t2 <- t2[do.call(order,as.data.frame(t2)),]
                all.equal(t1,t2)
            })
        }))
        all(results)
    }
    NoConstNoRepChar <- function() {
        myStr1 <- c("VxFVwte6K2","d6nLA724c4","N64jydNVPa","0Co4liTeCF","tIx6aDPNbD","wlGjg63lxj")
        results <- sapply(1:length(myStr1), function(y) {
                            t1 <- permuteGeneral(myStr1, y)
                            t1 <- t1[do.call(order,as.data.frame(t1)),]
                            t2 <- permutations(length(myStr1), y, myStr1)
                            t2 <- t2[do.call(order,as.data.frame(t2)),]
                            identical(t1, t2)
                    })
        all(results)
    }
    expect_true(NoConstNoRepNum1(7))
    expect_true(NoConstNoRepNum2(7))
    expect_true(NoConstNoRepChar())
})

test_that("permuteGeneral produces correct combinations with no constraints and with repetition", {
    NoConstWithRepNum1 <- function(n) {
        results <- do.call(c, sapply(1:n, function(x) {
            sapply(1:x, function(y) {
                t1 <- permuteGeneral(x, y, TRUE)
                t1 <- t1[do.call(order,as.data.frame(t1)),]
                t2 <- permutations(x, y, 1:x, FALSE, TRUE)
                t2 <- t2[do.call(order,as.data.frame(t2)),]
                all.equal(t1,t2)
            })
        }))
        all(results)
    }
    NoConstWithRepNum2 <- function(n) {
        set.seed(103)
        results <- do.call(c, sapply(2:n, function(x) {
            myV <- sort(runif(x, sample(-100:0, 1), sample(1:100, 1)))
            sapply(1:x, function(y) {
                t1 <- permuteGeneral(myV, y, TRUE)
                t1 <- t1[do.call(order,as.data.frame(t1)),]
                t2 <- permutations(x, y, myV, FALSE, TRUE)
                t2 <- t2[do.call(order,as.data.frame(t2)),]
                all.equal(t1,t2)
            })
        }))
        all(results)
    }
    NoConstWithRepChar <- function() {
        myStr2 <- c("sDXjdHpzjY","PJJ4c4opEi","JJUvO8NKA7","jVfiEMeqrm","xrPLmD7kvM","FpbTuOrpXJ")
        results <- sapply(1:length(myStr2), function(y) {
                            t1 <- permuteGeneral(myStr2, y, TRUE)
                            t1 <- t1[do.call(order,as.data.frame(t1)),]
                            t2 <- permutations(length(myStr2), y, myStr2, FALSE, TRUE)
                            t2 <- t2[do.call(order,as.data.frame(t2)),]
                            identical(t1, t2)
                        })
        all(results)
    }
    expect_true(NoConstWithRepNum1(6))
    expect_true(NoConstWithRepNum2(6))
    expect_true(NoConstWithRepChar())
})

test_that("permuteGeneral produces appropriate error messages", {
    expect_error(permuteGeneral(7,-4), "m must be positive")
    expect_error(permuteGeneral(170, 7), "The number of rows cannot exceed")
    expect_error(permuteGeneral(100, 10, TRUE), "The number of rows cannot exceed")
    expect_error(permuteGeneral(c(2,3,5,7),-4,TRUE), "m must be positive")
})