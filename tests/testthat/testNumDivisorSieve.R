context("testing numDivisorSieve")

test_that("numDivisorSieve generates correct numbers", {
    options(scipen = 999)
    expect_equal(numDivisorSieve(10)[10], 4)
    expect_equal(length(numDivisorSieve(1000)), 1000)
    expect_equal(numDivisorSieve(2), 1:2)
    expect_equal(numDivisorSieve(1), 1)
    expect_equal(numDivisorSieve(99, 100), c(6, 9))
    
    ## Test Names
    expect_equal(as.integer(names(numDivisorSieve(100, namedVector = TRUE))), 1:100)
    expect_equal(as.numeric(names(numDivisorSieve(10^12, 10^12 + 100,
                                                  namedVector = TRUE))), (10^12):(10^12 + 100))
})

test_that("numDivisorSieve produces appropriate error messages", {
    expect_error(numDivisorSieve(-1), "bound1 must be a positive number less than")
    expect_error(numDivisorSieve(0), "bound1 must be a positive number less than")
    expect_error(numDivisorSieve(2^53), "bound1 must be a positive number less than")
    expect_error(numDivisorSieve(2^53, 1), "must be a positive number less")
    expect_error(numDivisorSieve(1, 2^53), "must be a positive number less")
    expect_error(numDivisorSieve("10"), "must be of type numeric or integer")
    expect_error(numDivisorSieve(2, "10"), "must be of type numeric or integer")
    expect_error(numDivisorSieve(100, namedVector = "TRUE"), "Not compatible with requested type")
})
