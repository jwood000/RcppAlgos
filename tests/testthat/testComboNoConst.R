library(testthat)
library(RcppAlgos)
context("testing comboGeneral no constraints")

## Must convert strings to upper case as C++ uses the ASCII value
## and R doesn't. From ?Comparison: "lexicographic within the 
## strings using the collating sequence of the locale in use"

test_that("comboGeneral produces correct combinations with no constraints and no repetition", {
    NoConstNoRepNum1 <- function(n) {
        results <- do.call(c, sapply(1:n, function(x) {
            sapply(1:x, function(y) all.equal(comboGeneral(x, y), 
                                              t(combn(x, y))))
        }))
        all(results)
    }
    NoConstNoRepNum2 <- function(n) {
        set.seed(101)
        results <- do.call(c, sapply(2:n, function(x) {
            myV <- sort(runif(x, sample(-100:0, 1), sample(1:100, 1)))
            sapply(1:x, function(y) all.equal(comboGeneral(myV, y), 
                                              t(combn(myV, y))))
        }))
        all(results)
    }
    NoConstNoRepChar <- function() {
        myStr1 <- sort(toupper(c("VxFVwte6K2","d6nLA724c4","N64jydNVPa","0Co4liTeCF","tIx6aDPNbD",
                                 "wlGjg63lxj","E77pznzQRE","KXDbiyEGgp","2dnP37k8nr","zQvvOLanqn")))
        results <- sapply(1:10, function(y) identical(comboGeneral(myStr1, y), t(combn(myStr1, y))))
        all(results)
    }
    expect_true(NoConstNoRepNum1(10))
    expect_true(NoConstNoRepNum2(10))
    expect_true(NoConstNoRepChar())
})

test_that("comboGeneral produces correct combinations with no constraints and with repetition", {
    NoConstWithRepNum1 <- function(n) {
        results <- do.call(c, sapply(1:n, function(x) {
            sapply(1:x, function(y) all.equal(comboGeneral(x,y,TRUE), 
                                              gtools::combinations(x,y,1:x,FALSE,TRUE)))
        }))
        all(results)
    }
    NoConstWithRepNum2 <- function(n) {
        set.seed(103)
        results <- do.call(c, sapply(2:n, function(x) {
            myV <- sort(runif(x, sample(-100:0, 1), sample(1:100, 1)))
            sapply(1:x, function(y) all.equal(comboGeneral(myV,y,TRUE), 
                                              gtools::combinations(x,y,myV,FALSE,TRUE)))
        }))
        all(results)
    }
    NoConstWithRepChar <- function() {
        myStr2 <- sort(toupper(c("sDXjdHpzjY","PJJ4c4opEi","JJUvO8NKA7",
                                 "jVfiEMeqrm","xrPLmD7kvM","FpbTuOrpXJ")))
        results <- sapply(1:6, function(y) identical(comboGeneral(myStr2, y, TRUE), 
                                                     gtools::combinations(6, y, myStr2,FALSE,TRUE)))
        all(results)
    }
    expect_true(NoConstWithRepNum1(5))
    expect_true(NoConstWithRepNum2(5))
    expect_true(NoConstWithRepChar())
})

test_that("comboGeneral produces appropriate error messages", {
    expect_error(comboGeneral(-9,4), "v cannot be less than m")
    expect_error(comboGeneral(-9,4,TRUE), "v cannot be less than m")
    expect_error(comboGeneral(c(1,4,5),4,TRUE), "m must be less than or equal to the length of v")
    expect_error(comboGeneral(LETTERS[1:3],4,TRUE), "m must be less than or equal to the length of v")
    expect_error(comboGeneral(7,-4), "m must be positive")
    expect_error(comboGeneral(170, 7), "The number of rows cannot exceed")
    expect_error(comboGeneral(100, 10, TRUE), "The number of rows cannot exceed")
    expect_error(comboGeneral(c(2,3,5,7),-4,TRUE), "m must be positive")
})