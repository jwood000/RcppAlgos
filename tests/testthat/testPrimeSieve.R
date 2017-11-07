library(RcppAlgos)
context("testing primeSieve")

test_that("primeSieve generates correct numbers with 1 argument", {
    EqualityTest1 <- function() {
        set.seed(17)
        samp <- sample(10^6, 10^2)
        results <- sapply(samp, function(x) all.equal(numbers::Primes(x), primeSieve(x)))
        all(results)
    }
    expect_true(EqualityTest1())
})

test_that("generates correct numbers with 2 argument", {
    EqualityTest2 <- function() {
        set.seed(19)
        b1 <- c(sample(10^3, 50), sample(10^7, 50))
        b2 <- b1 + c(sample(10^5, 50), sample(10^7, 50))
        results <- sapply(1:100, function(x) all.equal(numbers::Primes(b1[x], b2[x]),
                                                      primeSieve(b1[x], b2[x])))
        b3 <- sample(10^13, 10)
        b4 <- b3 + sample(10^5, 10)
        results <- c(results, sapply(1:10, function(x) all.equal(numbers::Primes(b3[x], b4[x]),
                                                                 primeSieve(b3[x], b4[x]))))
        all(results)
    }
    expect_true(EqualityTest2())
})

test_that("primeSieve produces appropriate error messages", {
    expect_error(primeSieve(-1), "must be a positive number less")
    expect_error(primeSieve(1,-1), "must be a positive number less")
    expect_error(primeSieve(1,2^53), "must be a positive number less")
    expect_error(primeSieve(2^31), "must be less than")
    expect_error(primeSieve(2^4, "1"), "must be of type numeric or integer")
    expect_error(primeSieve("500"), "must be of type numeric or integer")
})