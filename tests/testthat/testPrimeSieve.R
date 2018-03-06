context("testing primeSieve")

test_that("primeSieve generates correct numbers", {
    expect_equal(primeSieve(10), c(2L, 3L, 5L, 7L))
    expect_equal(max(primeSieve(1000)), 997L)
    expect_equal(min(primeSieve(1000)), 2L)
    expect_equal(primeSieve(1), integer(0))
    expect_equal(primeSieve(2), 2L)
    
    expect_equal(primeSieve(6,8), 7)
    expect_equal(primeSieve(999982,10^6), 999983)
    expect_equal(primeSieve(2, 7), c(2, 3, 5, 7))
    
    expect_equal(primeSieve(10.1), primeSieve(10))
    expect_equal(primeSieve(2.1, 10.9), primeSieve(3, 7))
})

test_that("primeSieve produces appropriate error messages", {
    expect_error(primeSieve(-1), "must be a positive")
    expect_error(primeSieve(1,-1), "must be a positive number less")
    expect_error(primeSieve(1,2^53), "must be a positive number less")
    expect_error(primeSieve(2^53), "must be a positive number less")
    expect_error(primeSieve(2^53, 1), "must be a positive number less")
    expect_error(primeSieve(2^4, "1"), "must be of type numeric or integer")
    expect_error(primeSieve("500"), "must be of type numeric or integer")
})
