context("testing divisorsList")

test_that("divisorsList generates correct numbers", {
    expect_equal(divisorsList(10)[[10]], c(1, 2, 5, 10))
    expect_equal(length(divisorsList(1000)), 1000)
    expect_equal(divisorsList(1)[[1]], 1)
})

test_that("divisorsList produces appropriate error messages", {
    expect_error(divisorsList(-1), "must be positive")
    expect_error(divisorsList(0), "must be positive")
    expect_error(divisorsList(1,10), "unused argument")
    expect_error(divisorsList(2^31), "must be less than")
    expect_error(divisorsList("10"), "must be of type numeric or integer")
})