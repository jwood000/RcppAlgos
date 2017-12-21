context("testing eulerPhiSieve")

test_that("eulerPhiSieve generates correct numbers", {
    expect_equal(eulerPhiSieve(10)[10], 4)
    expect_equal(length(eulerPhiSieve(1000)), 1000)
    expect_equal(eulerPhiSieve(2), c(1, 1))
    expect_equal(eulerPhiSieve(1), 1)
})

test_that("eulerPhiSieve produces appropriate error messages", {
    expect_error(eulerPhiSieve(-1), "must be positive")
    expect_error(eulerPhiSieve(0), "must be positive")
    expect_error(eulerPhiSieve(1,10), "unused argument")
    expect_error(eulerPhiSieve(2^31), "must be less than")
    expect_error(eulerPhiSieve("10"), "must be of type numeric or integer")
})