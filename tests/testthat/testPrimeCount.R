context("testing primeCount")

test_that("primeCount generates correct numbers", {
    expect_equal(primeCount(10), 4)
    
    ## Pi(x) obtained from https://en.wikipedia.org/wiki/Prime-counting_function
    expect_equal(primeCount(10^6), 78498L)
    expect_equal(primeCount(as.integer(1e8)), 5761455L)
    expect_equal(primeCount(10^10), 455052511)
    expect_equal(primeCount(10^10, nThreads = 2), 455052511)
    expect_equal(primeCount(1e11, nThreads = 2), 4118054813)
    expect_equal(primeCount(1e12, nThreads = 4), 37607912018)
    expect_equal(primeCount(1), 0)
    expect_equal(sapply(1:9, primeCount), c(0,1,2,2,3,3,4,4,4))
})

test_that("primeCount produces appropriate error messages", {
    expect_error(primeCount(0), "must be a positive")
    expect_error(primeCount(-1), "must be a positive")
    expect_error(primeCount(2^53), "must be a positive number less")
    expect_error(primeCount("100000"), "must be of type numeric or integer")
})
