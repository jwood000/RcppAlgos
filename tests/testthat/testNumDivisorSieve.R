context("testing numDivisorSieve")

test_that("numDivisorSieve generates correct numbers", {
    expect_equal(numDivisorSieve(10)[10], 4)
    expect_equal(length(numDivisorSieve(1000)), 1000)
    expect_equal(numDivisorSieve(2), 1:2)
    expect_equal(numDivisorSieve(1), 1)
})

test_that("numDivisorSieve produces appropriate error messages", {
    expect_error(numDivisorSieve(-1), "must be positive")
    expect_error(numDivisorSieve(0), "must be positive")
    expect_error(numDivisorSieve(1,10), "unused argument")
    expect_error(numDivisorSieve(2^31), "must be less than")
    expect_error(numDivisorSieve("10"), "must be of type numeric or integer")
})