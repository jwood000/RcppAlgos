context("testing isPrimeRcpp")

test_that("isPrimeRcpp preserves empty result types and names", {
    for (v in list(integer(0), numeric(0), logical(0))) {
        expect_identical(isPrimeRcpp(v), logical(0))
        expect_identical(isPrimeRcpp(v, namedVector = TRUE),
                         structure(logical(0), names = character(0)))
    }

    expect_identical(isPrimeRcpp(numeric(0), nThreads = 2L), logical(0))
})

test_that("isPrimeRcpp generates correct numbers", {
    options(scipen = 999)
    expect_equal(isPrimeRcpp(10), FALSE)
    expect_equal(isPrimeRcpp(c(999983, 10^6)), c(TRUE, FALSE))
    expect_equal(isPrimeRcpp(1), FALSE)
    expect_equal(all(isPrimeRcpp(c(2,3,5,7,11))), TRUE)
    expect_equal(1e10 + which(isPrimeRcpp((1e10 + 1):(1e10 + 1e4))), primeSieve(1e10, 1e10 + 1e4))
    expect_equal(isPrimeRcpp(11:13), c(TRUE, FALSE, TRUE))

    expect_true(isPrimeRcpp(.Machine$integer.max))

    ## Test Names
    expect_equal(as.integer(names(isPrimeRcpp(100, namedVector = TRUE))), 100)
    expect_equal(as.numeric(names(isPrimeRcpp((10^12):(10^12 + 100),
                                                 namedVector = TRUE))), (10^12):(10^12 + 100))

    expect_equal(10000L + which(isPrimeRcpp(10001:15000)), primeSieve(10001, 15000))
    expect_equal(1e10 + which(isPrimeRcpp((1e10 + 1):(1e10 + 5000))), primeSieve(1e10, 1e10 + 5000))

    ## Test Parallel
    set.seed(150)
    samp <- sample(1e12, 10000)
    expect_equal(isPrimeRcpp(samp), isPrimeRcpp(samp, nThreads = 2))

    ## This hardware query launches no workers; do not use it as a test budget.
    expect_true(stdThreadMax() > 0)
})

test_that("primality agrees across the native and GMP dispatch boundary", {
    ## 31607 is prime, so its square reaches the native primality test after
    ## trial division. The primes straddle the 1e9 dispatch boundary.
    values <- c(31607^2, 999999937, 1e9 - 1, 1e9, 1e9 + 1, 1000000007)
    expected <- c(FALSE, TRUE, FALSE, FALSE, FALSE, TRUE)

    expect_identical(isPrimeRcpp(values, nThreads = 1), expected)
    expect_identical(isPrimeRcpp(values, nThreads = 2), expected)
})
