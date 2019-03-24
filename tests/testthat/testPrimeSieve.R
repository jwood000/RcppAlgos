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
    
    expect_equal(primeSieve(31625, 31630), 31627)
    
    expect_equal(2L + length(primeSieve(5, 1e5)), length(primeSieve(1e5)))
    
    ## Primes obtained from https://primes.utm.edu
    expect_equal(primeSieve(1e8, 1e8 + 100), c(100000007, 100000037, 100000039, 
                                               100000049, 100000073, 100000081))
    
    ## Primes obtained from http://compoasso.free.fr/primelistweb/page/prime/liste_online_en.php
    expect_equal(primeSieve(1e10, 1e10 + 100), c(10000000019, 10000000033, 
                                                 10000000061, 10000000069, 10000000097))
    
    ## Primes obtained from http://compoasso.free.fr/primelistweb/page/prime/liste_online_en.php
    expect_equal(primeSieve(1e12 - 100, 1e12), c(999999999937, 999999999959, 
                                                 999999999961, 999999999989))
    
    ## Large prime gap greater than 1e15 http://www.trnicely.net/gaps/gaps3.html
    expect_equal(primeSieve(1693182318746371, 1693182318746371 + 1132), 
                                    c(1693182318746371, 1693182318747503))
    
    ## Test Parallel
    expect_equal(primeSieve(1e7), primeSieve(1e7, nThreads = 2))
    expect_equal(primeSieve(1e12, 1e12 + 1e5), primeSieve(1e12, 1e12 + 1e5, nThreads = 2))
    expect_equal(primeSieve(1e15, 1e15 + 1e5), primeSieve(1e15, 1e15 + 1e5, nThreads = 2))
    
    expect_equal(primeSieve(1e12, 1e12 + 1e7), primeSieve(1e12, 1e12 + 1e7, nThreads = 3))
    expect_equal(primeSieve(1e15, 1e15 + 1e8), primeSieve(1e15, 1e15 + 1e8, nThreads = 3))
    
    ## The number 2916073, was obtained from Kim Walish's primesieve library
    expect_equal(length(primeSieve(769166929090560, 769167029090560, nThreads = 2)), 2916073)
    
    ## This test is for full code coverage... It takes about 1.5 seconds
    ## system.time(a <- primeSieve(6, 1.1e9, nThreads = 1))
})

test_that("primeSieve produces appropriate error messages", {
    expect_error(primeSieve(-1), "must be a positive")
    expect_error(primeSieve(1,-1), "must be a positive number")
    expect_error(primeSieve(1,2^53), "must be less than")
    expect_error(primeSieve(2^53), "must be less than")
    expect_error(primeSieve(2^53, 1), "must be less than")
    expect_error(primeSieve(2^4, "1"), "must be of type numeric or integer")
    expect_error(primeSieve("500"), "must be of type numeric or integer")
})
