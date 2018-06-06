context("testing primeSieve")

test_that("primeSieve generates correct numbers", {
    options(scipen = 999)
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
    expect_equal(primeSieve(1, 1), integer(0))
    
    expect_equal(2L + length(primeSieve(5, 1e5)), length(primeSieve(1e5)))
    
    ## Primes obtained from https://primes.utm.edu
    expect_equal(primeSieve(1e8, 1e8 + 100), c(100000007, 100000037, 100000039, 100000049, 100000073, 100000081))
    
    ## Primes obtained from http://compoasso.free.fr/primelistweb/page/prime/liste_online_en.php
    expect_equal(primeSieve(1e10, 1e10 + 100), c(10000000019, 10000000033, 10000000061, 10000000069, 10000000097))
    
    ## Primes obtained from http://compoasso.free.fr/primelistweb/page/prime/liste_online_en.php
    expect_equal(primeSieve(1e12 - 100, 1e12), c(999999999937, 999999999959,	999999999961, 999999999989))
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
