context("testing nthPrime")

test_that("nthPrime agrees with an independent sieve", {
    ## A base-R Eratosthenes sieve is independent of the native estimate,
    ## prime-counting kernels, and segmented selection routine.
    sieve <- rep(TRUE, 200000L)
    sieve[1L] <- FALSE
    for (p in 2:floor(sqrt(length(sieve)))) {
        if (sieve[p]) sieve[seq.int(p * p, length(sieve), by = p)] <- FALSE
    }
    primes <- which(sieve)
    ranks <- sort(unique(c(1:256, 999:1001, 1999:2001, 9999:10001,
                           seq.int(300L, length(primes), by = 137L))))
    actual <- vapply(ranks, nthPrime, integer(1))
    expect_identical(actual, primes[ranks])
})

test_that("nthPrime preserves known prime values", {
    ## Reference values were checked with the independent primecount CLI.
    ranks <- 10^(0:9)
    expected <- c(2, 29, 541, 7919, 104729, 1299709, 15485863,
                  179424673, 2038074743, 22801763489)
    expect_equal(vapply(ranks, nthPrime, numeric(1)), expected)
})

test_that("nthPrime preserves integer-double result boundaries", {
    ## primecount reports pi(2^31 - 1) = 105097565. Independent trial
    ## division confirms 2147483659 is the next prime after 2147483647.
    expect_identical(nthPrime(105097565), .Machine$integer.max)
    expect_identical(nthPrime(105097566), 2147483659)
    expect_identical(nthPrime(1L), 2L)
    expect_identical(nthPrime(25), 97L)
})

test_that("nthPrime returns the first integer with the requested prime count", {
    ## These are integration checks in addition to the independent sieve and
    ## reference values above. Counting both sides catches a composite on a
    ## pi plateau being mistaken for the requested prime.
    ranks <- c(1L, 2L, 25L, 65L, 1001L, 12345L, 1000000L)
    for (n in ranks) {
        p <- nthPrime(n)
        expect_equal(primeCount(p), n)
        expect_equal(primeCount(p - 1), n - 1L)
    }
})

test_that("nthPrime preserves adjacent ranks with opposite search directions", {
    ## Independent sieve references. The initial search anchors for these
    ## adjacent ranks straddle 1e7: 9999998 and 10000014. The former is a
    ## composite plateau searched backward; the latter is searched forward.
    ranks <- c(664579L, 664580L)
    expected <- c(9999991L, 10000019L)
    threads <- unique(c(1L, pmin(c(2L, 8L), RcppAlgos:::pkgEnv$maxThreads)))
    for (nThreads in threads) {
        expect_identical(vapply(ranks, nthPrime, integer(1),
                                nThreads = nThreads), expected)
    }
})

test_that("nthPrime agrees across serial and parallel requests", {
    ranks <- c(1e6, 1e9)
    expected <- list(15485863L, 22801763489)
    threads <- unique(c(1L, pmin(c(2L, 4L), RcppAlgos:::pkgEnv$maxThreads)))
    for (nThreads in threads) {
        expect_identical(lapply(ranks, nthPrime, nThreads = nThreads), expected)
    }
    ## An oversized request is capped before any workers are launched.
    expect_identical(nthPrime(1e9, .Machine$integer.max), expected[[2L]])
})

test_that("nthPrime preserves defaults and shared logical conversion", {
    expect_identical(nthPrime(25L, NULL), 97L)
    expect_identical(nthPrime(25L, 2L), 97L)
    expect_identical(nthPrime(TRUE), 2L)
    expect_identical(nthPrime(25L, TRUE), 97L)
})

test_that("nthPrime preserves larger independently checked results", {
    ## These ranks exercise bounded final selection after exact counting.
    ## References were obtained independently with the primecount CLI.
    ranks <- c(1e10, 1e12)
    expected <- c(252097800623, 29996224275833)
    threads <- unique(c(1L, pmin(8L, RcppAlgos:::pkgEnv$maxThreads)))
    for (nThreads in threads) {
        expect_identical(vapply(ranks, nthPrime, numeric(1),
                                nThreads = nThreads), expected)
    }
})

test_that("nthPrime rejects invalid ranks before preparing large tables", {
    invalid <- list(0, -1, 1.5, FALSE, NA, NA_integer_, NA_real_, NaN,
                    Inf, -Inf, numeric(), logical(), c(1, 2), c(TRUE, FALSE),
                    "25", NULL, list(1), 1i, as.raw(1),
                    252252704148405, 2^53 - 1, 2^53)
    for (n in invalid) {
        expect_error(nthPrime(n))
    }
})

test_that("nthPrime validates thread requests even for lookup results", {
    invalid <- list(0L, -1L, 1.5, FALSE, NA, NA_integer_, NaN, Inf, -Inf,
                    2^31, numeric(), logical(), c(1, 2), c(TRUE, FALSE),
                    "2", list(1), 1i)
    for (nThreads in invalid) {
        expect_error(nthPrime(1L, nThreads))
        expect_error(nthPrime(1000000L, nThreads))
    }
    expect_identical(nthPrime(1L, .Machine$integer.max), 2L)
})
