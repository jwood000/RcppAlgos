context("testing eulerPhiSieve")

test_that("eulerPhiSieve generates correct numbers", {
    options(scipen = 999)
    expect_equal(eulerPhiSieve(10)[10], 4)
    expect_equal(length(eulerPhiSieve(1000)), 1000)
    expect_equal(eulerPhiSieve(2), c(1, 1))
    expect_equal(eulerPhiSieve(2, 3), c(1, 2))
    expect_true(eulerPhiSieve(1, 1, TRUE)==1)
    expect_true(eulerPhiSieve(1L, namedVector = TRUE)==1)
    expect_equal(eulerPhiSieve(11, 20), c(10, 4, 12, 6, 8, 8, 16, 6, 18, 8))
    expect_equal(eulerPhiSieve(1e13, 1e13), 4000000000000)
    ## Test Names
    expect_equal(as.integer(names(eulerPhiSieve(100, namedVector = TRUE))), 1:100)
    expect_equal(as.numeric(names(eulerPhiSieve(10^12, 10^12 + 100,
                                                namedVector = TRUE))), (10^12):(10^12 + 100))
})

test_that("eulerPhiSieve produces appropriate error messages", {
    expect_error(eulerPhiSieve(-1), "must be a positive number less")
    expect_error(eulerPhiSieve(0), "must be a positive number less")
    expect_error(eulerPhiSieve(2^53), "bound1 must be a positive number less than")
    expect_error(eulerPhiSieve(2^53, 1), "must be a positive number less")
    expect_error(eulerPhiSieve(1, 2^53), "must be a positive number less")
    expect_error(eulerPhiSieve("10"), "must be of type numeric or integer")
    expect_error(eulerPhiSieve(2, "10"), "must be of type numeric or integer")
    expect_error(eulerPhiSieve(100, namedVector = "TRUE"), "Not compatible with requested type")
})
