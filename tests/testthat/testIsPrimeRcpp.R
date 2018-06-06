context("testing isPrimeRcpp")

test_that("isPrimeRcpp generates correct numbers", {
    options(scipen = 999)
    expect_equal(isPrimeRcpp(10), FALSE)
    expect_equal(isPrimeRcpp(c(999983, 10^6)), c(TRUE, FALSE))
    expect_equal(isPrimeRcpp(1), FALSE)
    expect_equal(all(isPrimeRcpp(c(2,3,5,7,11))), TRUE)
    expect_equal(isPrimeRcpp(11.1), FALSE)
    expect_equal(isPrimeRcpp(c(11, 11.1, 12.9999, 13)), c(TRUE, FALSE, FALSE, TRUE))
    
    expect_true(isPrimeRcpp(.Machine$integer.max))
    
    ## Test Names
    expect_equal(as.integer(names(isPrimeRcpp(100, namedVector = TRUE))), 100)
    expect_equal(as.numeric(names(isPrimeRcpp((10^12):(10^12 + 100),
                                                 namedVector = TRUE))), (10^12):(10^12 + 100))
})

test_that("isPrimeRcpp produces appropriate error messages", {
    expect_error(isPrimeRcpp(0), "each element must be positive")
    expect_error(isPrimeRcpp(-1), "each element must be positive")
    expect_error(isPrimeRcpp(2^53), "each element must be less than")
    expect_error(isPrimeRcpp("100000"), "must be of type numeric or integer")
})
