context("testing primeFactorizationList")

test_that("primeFactorizationList generates correct numbers", {
    expect_equal(primeFactorizationList(100)[[100]], c(2, 2, 5, 5))
    expect_equal(length(primeFactorizationList(1000)), 1000)
    expect_equal(primeFactorizationList(2)[[2]], 2)
})

test_that("primeFactorizationList produces appropriate error messages", {
    expect_error(primeFactorizationList(-1), "must be positive")
    expect_error(primeFactorizationList(1,10), "unused argument")
    expect_error(primeFactorizationList(2^31), "must be less than")
    expect_error(primeFactorizationList("10"), "must be of type numeric or integer")
})