context("testing primeSieve")

test_that("primeSieve generates correct numbers", {
    options(scipen = 999)
    expect_equal(primeSieve(10), c(2L, 3L, 5L, 7L))
    expect_equal(max(primeSieve(1000)), 997L)
    expect_equal(min(primeSieve(1000)), 2L)
    expect_equal(primeSieve(1), integer(0))
    expect_equal(primeSieve(2), 2L)
    
    expect_equal(primeSieve(6,8), 7)
    expect_equal(primeSieve(999982, 10^6), 999983)
    expect_equal(primeSieve(1e10, 1e10 + 20), 10000000019)
    expect_equal(primeSieve(2, 7), c(2, 3, 5, 7))
    
    expect_equal(primeSieve(10.1), primeSieve(10))
    expect_equal(primeSieve(2.1, 10.9), primeSieve(3, 7))
    expect_equal(primeSieve(1, 1), integer(0))
    
    expect_equal(primeSieve(31625, 31630), 31627)
    expect_equal(2L + length(primeSieve(5, 1e5)), length(primeSieve(1e5)))
    
    ## Primes obtained from https://primes.utm.edu
    expect_equal(primeSieve(1e8, 1e8 + 100), c(100000007, 100000037, 100000039, 100000049, 100000073, 100000081))
    
    ## Primes obtained from http://compoasso.free.fr/primelistweb/page/prime/liste_online_en.php
    expect_equal(primeSieve(1e10, 1e10 + 100), c(10000000019, 10000000033, 10000000061, 10000000069, 10000000097))
    
    ## Primes obtained from http://compoasso.free.fr/primelistweb/page/prime/liste_online_en.php
    expect_equal(primeSieve(1e12 - 100, 1e12), c(999999999937, 999999999959, 999999999961, 999999999989))
    
    ## Large prime gap greater than 1e15 http://www.trnicely.net/gaps/gaps3.html
    expect_equal(primeSieve(1693182318746371, 1693182318746371 + 1132), 
                                    c(1693182318746371, 1693182318747503))
    
    ## Test Parallel
    funTestPar <- function(b1, b2, nT = 2) {
        ser <- primeSieve(b1, b2)
        par <- primeSieve(b1, b2, nThreads = nT)
        all.equal(par, ser)
    }

    expect_true(funTestPar(1, 1e7))
    expect_true(funTestPar(1e12, 1e12 + 1e7))
    expect_true(funTestPar(1e15, 1e15 + 1e7))
    
    # The following test can be used for fuller coverage but are commented
    # because they take a bit of time. All results were confirmed with
    # primesieve by Kim Walisch
    # gc()
    # expect_equal(length(primeSieve(633318687598976, 633318707598976)), 587185)
    # gc()
    # expect_equal(length(primeSieve(12e8)), 60454705)
    # gc()
    # expect_equal(length(primeSieve(39582118599936, 39582718599936, nThreads = 4)), 19161558)
})
