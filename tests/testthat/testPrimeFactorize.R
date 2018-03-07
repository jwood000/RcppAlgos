context("testing primeFactorize")

test_that("primeFactorize generates correct numbers", {
    options(scipen = 999)
    expect_equal(primeFactorize(100), c(2, 2, 5, 5))
    expect_equal(length(primeFactorize(1:100)), 100)
    expect_equal(primeFactorize(2), 2)
    expect_equal(primeFactorize(1), numeric(0))
    expect_equal(primeFactorize(1000), c(2,2,2,5,5,5))
    expect_equal(primeFactorize(1000, TRUE), c(2,2,2,5,5,5))
    
    ## Test Names
    expect_equal(as.integer(names(primeFactorize(100, namedList = TRUE))), integer(0))
    expect_equal(as.numeric(names(primeFactorize((10^12):(10^12 + 100),
                                                      namedList = TRUE))), (10^12):(10^12 + 100))
})

test_that("primeFactorize produces appropriate error messages", {
    expect_error(primeFactorize(2^53), "each element must be less than")
    expect_error(primeFactorize(-2^53), "each element must be less than")
    expect_error(primeFactorize("10"), "must be of type numeric or integer")
    expect_error(primeFactorize(100, namedList = "TRUE"), "Not compatible with requested type")
})