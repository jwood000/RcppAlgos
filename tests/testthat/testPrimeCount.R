context("testing primeCount")

test_that("primeCount generates correct numbers", {
    expect_equal(primeCount(10), 4)
    expect_equal(primeCount(10^6), 78498)
    expect_equal(primeCount(1), 0)
    expect_equal(primeCount(10.1), primeCount(10))
    expect_equal(primeCount(10.9), primeCount(10))
})

test_that("primeCount produces appropriate error messages", {
    expect_error(primeCount(0), "must be a positive")
    expect_error(primeCount(-1), "must be a positive")
    expect_error(primeCount(2^53), "must be a positive number less")
    expect_error(primeCount("100000"), "must be of type numeric or integer")
})
