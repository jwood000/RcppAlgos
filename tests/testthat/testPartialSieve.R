context("testing partial prime sieves")

test_that("partial prime sieves agree with independent trial factorization", {
    powers <- as.vector(outer(c(2, 3, 7), 2:4, "^"))
    products <- c(powers, 17^2, 17^3, 17^2 * 19, 17 * 19^2)
    boundaries <- sort(unique(as.vector(outer(products, -1:1, "+"))))
    maximum <- max(boundaries)

    ## Trial division by consecutive integers is independent of the native
    ## prime tables, Phi recurrence, and P2/P3 summation identities.
    omega <- integer(maximum)
    least <- rep(Inf, maximum)
    for (n in 2:maximum) {
        remaining <- n
        divisor <- 2L
        while (divisor <= remaining %/% divisor) {
            while (remaining %% divisor == 0) {
                if (is.infinite(least[n])) least[n] <- divisor
                omega[n] <- omega[n] + 1L
                remaining <- remaining %/% divisor
            }
            divisor <- divisor + 1L
        }
        if (remaining > 1) {
            if (is.infinite(least[n])) least[n] <- remaining
            omega[n] <- omega[n] + 1L
        }
    }
    primes <- which(omega == 1L)
    reference <- function(a) {
        cutoff <- if (a == 0) 1 else if (a <= length(primes)) primes[a] else Inf
        survives <- least > cutoff
        survives[1] <- TRUE
        list(phi = c(0L, cumsum(survives)),
             pairs = c(0L, cumsum(survives & omega == 2L)),
             triples = c(0L, cumsum(survives & omega == 3L)))
    }

    inputs <- 0:120
    for (a in c(0:(sum(primes <= 120) + 1L), 2^53 - 1)) {
        expected <- reference(a)
        expect_identical(vapply(inputs, function(n) primePhi(n, a), integer(1)),
                         expected$phi[inputs + 1L])
        expect_identical(vapply(inputs, function(n) semiPrimeCount(n, a), integer(1)),
                         expected$pairs[inputs + 1L])
        expect_identical(vapply(inputs, function(n) triPrimeCount(n, a), integer(1)),
                         expected$triples[inputs + 1L])
    }
    for (a in c(0, 1, 3, 6, 7, 8, 25)) {
        expected <- reference(a)
        expect_identical(vapply(boundaries, function(n) primePhi(n, a), integer(1)),
                         expected$phi[boundaries + 1L])
        expect_identical(vapply(boundaries, function(n) semiPrimeCount(n, a), integer(1)),
                         expected$pairs[boundaries + 1L])
        expect_identical(vapply(boundaries, function(n) triPrimeCount(n, a), integer(1)),
                         expected$triples[boundaries + 1L])
    }
})

test_that("partial prime sieves preserve defaults and positional arguments", {
    expect_identical(semiPrimeCount(100), 34L)
    expect_identical(triPrimeCount(100), 22L)
    expect_identical(primePhi(100L, 1L), 50L)
    expect_identical(primePhi(100, 1, NULL), 50L)
    expect_identical(semiPrimeCount(100, 0L, NULL), 34L)
    expect_identical(triPrimeCount(100, 0L, NULL), 22L)
    expect_identical(primePhi(100, 1L, 2L), 50L)
    expect_identical(semiPrimeCount(100, 0L, 2L), 34L)
    expect_identical(triPrimeCount(100, 0L, 2L), 22L)
    expect_error(primePhi(100), "a")
})

test_that("primePhi preserves exact values and integer-double result boundaries", {
    maxExact <- 2^53 - 1
    inputs <- list(0L, 1L, 2L, .Machine$integer.max, 2^31, maxExact)
    expected <- list(0L, 1L, 2L, .Machine$integer.max, 2147483648, maxExact)
    expect_identical(lapply(inputs, primePhi, a = 0L), expected)

    ## Removing only multiples of two leaves exactly ceiling(n/2).
    expect_identical(primePhi(4294967294, 1), .Machine$integer.max)
    expect_identical(primePhi(4294967295, 1), 2147483648)
    expect_identical(primePhi(maxExact - 1, 1), 4503599627370495)
    expect_identical(primePhi(maxExact, 1), 4503599627370496)
})

test_that("cutoffs beyond pi(n) and large empty cases are supported", {
    maxExact <- 2^53 - 1
    for (n in c(1, 100, maxExact)) {
        expect_identical(primePhi(n, maxExact), 1L)
        expect_identical(semiPrimeCount(n, maxExact), 0L)
        expect_identical(triPrimeCount(n, maxExact), 0L)
    }
    for (fun in list(primePhi, semiPrimeCount, triPrimeCount)) {
        expect_identical(fun(0, maxExact), 0L)
    }

    ## Here a < n, so Phi must subtract beyond pi(sqrt(n)) and clamp at one.
    expect_identical(primePhi(100, 10), 16L)
    expect_identical(primePhi(100, 25), 1L)
    expect_identical(primePhi(100, 26), 1L)
    expect_identical(semiPrimeCount(100, 5), 0L)
    expect_identical(triPrimeCount(100, 3), 0L)
})

test_that("partial prime sieves satisfy known prime-count identities", {
    cases <- list(c(n = 1e6, cubeA = 25, fourthA = 11, pi = 78498),
                  c(n = 1e9, cubeA = 168, fourthA = 40, pi = 50847534))
    for (case in cases) {
        n <- unname(case["n"])
        a <- unname(case["cubeA"])
        expect_equal(primePhi(n, a) + a - 1 - semiPrimeCount(n, a),
                     unname(case["pi"]))
        a <- unname(case["fourthA"])
        expect_equal(primePhi(n, a) + a - 1 - semiPrimeCount(n, a) - triPrimeCount(n, a),
                     unname(case["pi"]))
    }
})

test_that("partial prime sieves agree across serial and parallel calls", {
    ## Inclusion-exclusion independently supplies the large Phi expectation.
    firstPrimes <- c(2, 3, 5, 7, 11, 13, 17, 19)
    expectedPhi <- 0
    for (mask in 0:255) {
        selected <- as.logical(intToBits(mask)[1:8])
        expectedPhi <- expectedPhi + (-1)^sum(selected) * floor(1e9 / prod(firstPrimes[selected]))
    }
    expected <- list(phi = as.integer(expectedPhi),
                     pairs = semiPrimeCount(1e10, 325, nThreads = 1L),
                     triples = triPrimeCount(1e11, 102, nThreads = 1L))
    ## These P2/P3 cases pass their real-worker thresholds with modest tables.
    ## Their serial references test scheduling; independent small-factor tests
    ## and known-count identities above establish separate mathematical checks.
    threads <- unique(c(1L, pmin(c(2L, 4L), as.integer(RcppAlgos:::pkgEnv$maxThreads))))
    for (nThreads in threads) {
        actual <- list(phi = primePhi(1e9, 8, nThreads),
                       pairs = semiPrimeCount(1e10, 325, nThreads),
                       triples = triPrimeCount(1e11, 102, nThreads))
        expect_identical(actual, expected)
    }
    expect_identical(primePhi(1e9, 8, nThreads = 2L), expected$phi)
})

test_that("semiPrimeCount handles short tails and prime-square boundaries", {
    ## This independent pi table covers every quotient used below. The large
    ## cutoff leaves a short tail, including pairs with repeated prime factors.
    sieve <- rep(TRUE, 100000L)
    sieve[1L] <- FALSE
    for (p in 2:floor(sqrt(length(sieve)))) {
        if (sieve[p]) sieve[seq.int(p * p, length(sieve), by = p)] <- FALSE
    }
    primes <- which(sieve)
    primeCounts <- cumsum(sieve)
    reference <- function(n, a) {
        b <- sum(primes <= sqrt(n))
        if (a >= b) return(0L)
        indices <- seq.int(a + 1L, b)
        quotients <- floor(n / primes[indices])
        stopifnot(max(quotients) <= length(primeCounts))
        as.integer(sum(primeCounts[quotients] - (indices - 1L)))
    }

    p <- max(primes[primes <= sqrt(3e9)])
    a <- sum(primes <= p) - 2048L
    inputs <- c(as.double(p)^2 + c(-1, 0, 1), 3e9)
    expected <- vapply(inputs, reference, integer(1), a = a)
    threads <- unique(c(1L, pmin(2L, as.integer(RcppAlgos:::pkgEnv$maxThreads))))
    for (nThreads in threads) {
        actual <- vapply(inputs, semiPrimeCount, integer(1),
                         a = a, nThreads = nThreads)
        expect_identical(actual, expected)
    }

    ## Keep the last two eligible smaller primes, then just the last, then
    ## none. This checks the same cutoff convention as the longer tail above.
    b <- sum(primes <= sqrt(3e9))
    cutoffs <- b + c(-2L, -1L, 0L, 1L)
    expected <- vapply(cutoffs, function(a) reference(3e9, a), integer(1))
    actual <- vapply(cutoffs, function(a) semiPrimeCount(3e9, a), integer(1))
    expect_identical(actual, expected)
})

test_that("semiPrimeCount preserves large counts across cutoffs and threads", {
    ## Independent per-prime quotient counts are retained in
    ## .github/development/p2-hybrid/results-matched-profile/reference-terms.csv.gz.
    ## They use a separate sieve and primecount, not this package's Phi or P2.
    cutoffs <- c(0L, 590L)
    expected <- list(11238880589, 2069687562L)
    threads <- unique(c(pmin(c(2L, 4L),
                            as.integer(RcppAlgos:::pkgEnv$maxThreads)), 1L))
    for (nThreads in threads) {
        actual <- lapply(cutoffs, function(a) semiPrimeCount(8e10, a, nThreads))
        expect_identical(actual, expected)
    }
})

test_that("partial prime sieves validate scalar inputs before shortcuts", {
    invalid <- list(-1, 1.5, NA_real_, NA_integer_, NA, NaN, Inf, -Inf, 2^53,
                    numeric(), logical(), c(1, 2), c(TRUE, FALSE), "10", NULL, list(1))
    for (fun in list(primePhi, semiPrimeCount, triPrimeCount)) {
        for (value in invalid) {
            expect_error(fun(value, 0))
            expect_error(fun(0, value))
        }
        expect_identical(fun(0L, 0L), 0L)
    }
})

test_that("partial prime sieves use the shared logical input conversion", {
    for (fun in list(primePhi, semiPrimeCount, triPrimeCount)) {
        expect_identical(fun(FALSE, FALSE), 0L)
        expect_identical(fun(TRUE, FALSE), fun(1L, 0L))
        expect_identical(fun(100, FALSE), fun(100, 0L))
        expect_identical(fun(100, TRUE, nThreads = TRUE), fun(100, 1L, 1L))
        expect_identical(fun(0, 0, nThreads = TRUE), 0L)
    }
})

test_that("partial prime sieves always validate nThreads", {
    invalid <- list(0L, -1L, 1.5, NA_integer_, NA, Inf, 2^31,
                    numeric(), logical(), c(1, 2), c(TRUE, FALSE), "2", FALSE)
    for (fun in list(primePhi, semiPrimeCount, triPrimeCount)) {
        for (value in invalid) {
            expect_error(fun(0, 0, nThreads = value), "nThreads")
            expect_error(fun(10, 2^53 - 1, nThreads = value), "nThreads")
        }
        ## Unsupported types use the shared converter's generic error message.
        expect_error(fun(0, 0, nThreads = list(1)))
        expect_error(fun(10, 2^53 - 1, nThreads = list(1)))
        expect_identical(fun(0, 0, nThreads = .Machine$integer.max), 0L)
    }
})

test_that("triPrimeCount rejects unsupported table coverage before allocation", {
    expect_error(triPrimeCount(379625064, 0), "94906265")
    expect_error(triPrimeCount(379625065, 0), "94906265")
    ## Check valid neighboring n without allocating the huge a=0 table.
    expect_identical(triPrimeCount(379625063, a = 100000000), 0L)
    expect_identical(triPrimeCount(379625064, a = 100000000), 0L)
})
