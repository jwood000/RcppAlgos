context("testing primeFactorizeSieve")

test_that("primeFactorizeSieve generates correct numbers", {
    expect_equal(primeFactorizeSieve(100)[[100]], c(2, 2, 5, 5))
    expect_equal(length(primeFactorizeSieve(1000)), 1000)
    expect_equal(primeFactorizeSieve(2)[[2]], 2)
    expect_equal(primeFactorizeSieve(2, 2)[[1]], 2)
    expect_equal(primeFactorizeSieve(1000, 1000)[[1]], c(2,2,2,5,5,5))
})

test_that("primeFactorizeSieve produces appropriate error messages", {
    expect_error(primeFactorizeSieve(-1), "must be a positive number")
    expect_error(primeFactorizeSieve(2^53), "must be a positive number less")
    expect_error(primeFactorizeSieve(2^53, 1), "must be a positive number less")
    expect_error(primeFactorizeSieve(1, 2^53), "must be a positive number less")
    expect_error(primeFactorizeSieve("10"), "must be of type numeric or integer")
    expect_error(primeFactorizeSieve(2, "10"), "must be of type numeric or integer")
    expect_error(primeFactorizeSieve(100, namedList = "TRUE"), "Not compatible with requested type")
})
