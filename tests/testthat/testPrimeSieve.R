context("testing primeSieve")

test_that("primeSieve handles the precomputed prime table endpoint", {
    # 31627 is the final entry in smallPrimeBase. Exercise its exact endpoint
    # so sanitizer builds also check that traversal stops within the table.
    primes <- primeSieve(31627)
    expect_identical(length(primes), 3402L)
    expect_identical(tail(primes, 1L), 31627L)

    expect_identical(primeSieve(31625, 31626), integer(0))
    expect_identical(primeSieve(31625, 31627), 31627L)
    expect_identical(primeSieve(31625, 31628), 31627L)
    expect_identical(primeSieve(31627, 31627), 31627L)

    # Rounding an interval with no integers puts the lower bound past the table.
    expect_identical(primeSieve(31627.1, 31627.2), integer(0))

    # Prime upper bounds must remain inclusive on the table path.
    expect_identical(primeSieve(13), c(2L, 3L, 5L, 7L, 11L, 13L))
})

# An independent reference for modest values and short intervals.
primeSieveTrialPrimes <- function(lower, upper) {
    values <- seq.int(lower, upper)
    keep <- vapply(values, function(n) {
        if (n < 2) return(FALSE)
        if (n == 2) return(TRUE)
        if (n %% 2 == 0) return(FALSE)
        limit <- floor(sqrt(n))
        limit < 3 || all(n %% seq.int(3, limit, by = 2) != 0)
    }, logical(1))
    primes <- values[keep]
    if (upper <= .Machine$integer.max) as.integer(primes) else as.double(primes)
}

test_that("primeSieve agrees with trial division across segment boundaries", {
    # Almost2310L1Cache in NumbersUtils/Wheel.h.
    segmentSize <- 32340
    for (boundary in c(segmentSize, 2 * segmentSize)) {
        lower <- boundary - 50
        upper <- boundary + 50
        expected <- primeSieveTrialPrimes(lower, upper)
        for (offset in -1:1) {
            endpoint <- boundary + offset
            expect_identical(primeSieve(lower, endpoint),
                             expected[expected <= endpoint])
            expect_identical(primeSieve(endpoint, upper),
                             expected[expected >= endpoint])
        }
    }

    # Include full intervening segments as well as partial first/last segments.
    lower <- segmentSize - 17
    upper <- 4 * segmentSize + 23
    expect_identical(primeSieve(lower, upper),
                     primeSieveTrialPrimes(lower, upper))
})

test_that("primeSieve carries offsets across the small-sieve marking threshold", {
    # PrimeSieveSmall switches marking strategies at segSize^2 / 4.
    # Span both branches with full segments and partial first/last segments.
    segmentSize <- 32340
    threshold <- segmentSize^2 / 4
    lower <- threshold - 2 * segmentSize - 17
    upper <- threshold + 2 * segmentSize + 23

    values <- seq.int(lower, upper)
    expected <- as.integer(values[isPrimeRcpp(values)])
    expect_identical(primeSieve(lower, upper, nThreads = 1), expected)
})

test_that("primeSieve preserves results when cropping small-sieve segments", {
    segmentSize <- 32340
    # Prime factors with only the smaller factor within the sieving bound.
    # 101 marks densely; 20011 has a stride larger than one Small segment.
    for (composite in c(101 * 1009, 20011 * 40009)) {
        boundary <- floor(composite / segmentSize) * segmentSize
        intervals <- list(
            odd_start = c(composite, composite + 100),
            even_start = c(composite - 1, composite + 100),
            composite_at_end = c(composite - 100, composite),
            aligned_start = c(boundary, boundary + 100),
            full_segment_carry = c(composite, composite + 4 * segmentSize + 101)
        )

        # Short interior ranges crop both ends and must still strike the
        # composite at either endpoint. The longer range includes full
        # segments, including one skipped by the sparse prime's carry.
        for (caseName in names(intervals)) {
            bounds <- intervals[[caseName]]
            values <- seq.int(bounds[1], bounds[2])
            expected <- as.integer(values[isPrimeRcpp(values)])
            expect_identical(primeSieve(bounds[1], bounds[2], nThreads = 1),
                             expected,
                             info = paste("composite =", composite, caseName))
        }
    }
})

test_that("primeSieve excludes prime squares at interval boundaries", {
    # Exercise the small and medium kernels; these are all prime bases.
    for (p in c(181, 257, 65537)) {
        square <- p * p
        expected <- primeSieveTrialPrimes(square - 20, square + 20)
        for (offset in -1:1) {
            endpoint <- square + offset
            expect_identical(primeSieve(square - 20, endpoint),
                             expected[expected <= endpoint])
            expect_identical(primeSieve(endpoint, square + 20),
                             expected[expected >= endpoint])
        }
    }
})

test_that("primeSieve preserves results across the small/medium cutoff", {
    # PrimeSieveMain defines smallCut from Almost2310L1Cache.
    cutoff <- (32340 - 100)^2
    expected <- primeSieveTrialPrimes(cutoff - 100, cutoff + 100)
    for (offset in c(-1, 0, 1, 100)) {
        upper <- cutoff + offset
        expect_identical(primeSieve(cutoff - 100, upper),
                         expected[expected <= upper])
    }
    expect_identical(primeSieve(cutoff + 1, cutoff + 100),
                     expected[expected > cutoff])
})

test_that("primeSieve preserves values and types at the integer limit", {
    limit <- .Machine$integer.max
    expected <- primeSieveTrialPrimes(limit - 50, limit + 50)
    for (offset in c(-1, 0, 1, 50)) {
        upper <- limit + offset
        reference <- expected[expected <= upper]
        if (upper <= limit) reference <- as.integer(reference)
        expect_identical(primeSieve(limit - 50, upper), reference)
    }
})

test_that("primeSieve preserves square starts and carry when cropping medium segments", {
    # Fix dispatch inputs so these intervals use PrimeSieveMedium on any host.
    nativeSieve <- function(lower, upper) {
        .Call(RcppAlgos:::`_RcppAlgos_PrimeSieveCpp`,
              lower, upper, 1L, 1L, 1L)
    }

    # Prime squares in the medium range, with integer and double outputs.
    for (p in c(32771, 1000003)) {
        square <- p * p
        segmentSize <- ceiling(sqrt(square + 101) / 30030) * 30030
        boundary <- floor(square / segmentSize) * segmentSize
        intervals <- list(
            square_at_end = c(square - 100, square),
            square_at_start = c(square, square + 100),
            after_square = c(square + 1, square + 101),
            aligned_start = c(boundary, boundary + 100),
            crossing_boundary = c(boundary - 50, boundary + 50)
        )

        # Exercise minNum below, equal to, and above p^2. The square must
        # still be marked at either inclusive endpoint. The short ranges
        # also cover a single cropped segment and carry into a second one.
        for (caseName in names(intervals)) {
            bounds <- intervals[[caseName]]
            values <- seq.int(bounds[1], bounds[2])
            expected <- values[isPrimeRcpp(values)]
            if (bounds[2] <= .Machine$integer.max) {
                expected <- as.integer(expected)
            }
            expect_identical(nativeSieve(bounds[1], bounds[2]), expected,
                             info = paste("p =", p, caseName))
        }
    }
})

test_that("primeSieve carries cropped medium offsets through several full segments", {
    # Fix dispatch inputs to exercise Medium with integer and double outputs.
    nativeSieve <- function(lower, upper) {
        .Call(RcppAlgos:::`_RcppAlgos_PrimeSieveCpp`,
              lower, upper, 1L, 1L, 1L)
    }

    for (base in c(1.5e9, 4e9)) {
        # Medium's segment size stays constant over each selected interval.
        segmentSize <- ceiling(sqrt(base) / 30030) * 30030
        boundary <- floor(base / segmentSize) * segmentSize
        lower <- boundary + segmentSize / 2 + 17
        upper <- boundary + 4 * segmentSize + 101

        # A partial first segment, three full segments, and a partial last
        # segment. Use independent primality checks, not another sieve run.
        values <- seq.int(lower, upper)
        expected <- values[isPrimeRcpp(values)]
        if (upper <= .Machine$integer.max) {
            expected <- as.integer(expected)
        }
        expect_identical(nativeSieve(lower, upper), expected,
                         info = paste("base =", base))
    }
})

test_that("primeSieve preserves results across the medium/large cutoff", {
    # Fix the cache heuristic inputs rather than depending on host core counts.
    # This invokes the same registered entry point as the public R wrapper.
    nativeSieve <- function(lower, upper) {
        .Call(RcppAlgos:::`_RcppAlgos_PrimeSieveCpp`,
              lower, upper, 1L, 1L, 1L)
    }
    cutoff <- (6 * 1024^2)^2
    values <- seq.int(cutoff - 100, cutoff + 100)
    expected <- values[isPrimeRcpp(values)]
    for (offset in c(-1, 0, 1, 100)) {
        upper <- cutoff + offset
        expect_identical(nativeSieve(cutoff - 100, upper),
                         expected[expected <= upper])
    }
    expect_identical(nativeSieve(cutoff + 1, cutoff + 100),
                     expected[expected > cutoff])
})

test_that("primeSieve bounds bucket scheduling at the final segment", {
    # Fix nThreads = 1, maxCores = 1, maxThreads = 10 so PrimeSieveBig
    # uses 640 wheels of 30030 integers, regardless of the host core count.
    nativeSieve <- function(lower, upper) {
        .Call(RcppAlgos:::`_RcppAlgos_PrimeSieveCpp`,
              lower, upper, 1L, 1L, 10L)
    }
    segmentSize <- 640 * 30030
    composite <- 10000019 * 103751579
    boundary <- composite - 1
    expect_equal(boundary %% segmentSize, 0)

    # The first factor belongs to the bucketed primes. Its product lies at
    # offset segSize + 1 from the preceding segment's start. Previously an
    # upper bound of boundary - 1 allocated one bucket but scheduled this
    # product in bucket 1 (past the end). Sanitizers detect the invalid access.
    lower <- boundary - 1001
    values <- seq.int(lower, boundary + 100)
    expected <- values[isPrimeRcpp(values)]
    expect_false(isPrimeRcpp(composite))

    # Include an exact composite endpoint to ensure the final valid strike
    # is retained, and a prime endpoint to check inclusive result extraction.
    endpoints <- c(boundary + c(-2, -1, 0, 1, 2), tail(expected, 1L))
    for (upper in endpoints) {
        expect_identical(nativeSieve(lower, upper),
                         expected[expected <= upper])
    }
})

test_that("primeSieve preserves results when cropping large-sieve segments", {
    # Force PrimeSieveBig with 640 wheels per segment, independent of the host.
    nativeSieve <- function(lower, upper) {
        .Call(RcppAlgos:::`_RcppAlgos_PrimeSieveCpp`,
              lower, upper, 1L, 1L, 10L)
    }
    segmentSize <- 640 * 30030
    boundary <- floor(1e14 / segmentSize) * segmentSize

    # Both factors are prime. Only the smaller factor is within the sieving
    # bound: 1000003 is dense (2*p < segmentSize), while 10000019 is bucketed.
    # Neither product is removed by the 30030 wheel itself.
    denseComposite <- 1000003 * 100000007
    bucketComposite <- 10000019 * 100000007
    intervals <- list(
        interior = c(boundary + segmentSize / 2 + 17,
                     boundary + segmentSize / 2 + 1017),
        aligned_start = c(boundary, boundary + 1000),
        crossing_boundary = c(boundary - 50, boundary + 50),
        dense_odd_start = c(denseComposite, denseComposite + 1000),
        dense_even_start = c(denseComposite - 1, denseComposite + 1000),
        bucket_odd_start = c(bucketComposite, bucketComposite + 1000),
        bucket_even_start = c(bucketComposite - 1, bucketComposite + 1000)
    )

    # Odd composite starts must be marked inclusively. The preceding even
    # starts exercise finding that multiple from a non-divisible lower bound.
    # Short interior ranges crop both ends; crossing a boundary tests carry.
    for (caseName in names(intervals)) {
        bounds <- intervals[[caseName]]
        values <- seq.int(bounds[1], bounds[2])
        expected <- values[isPrimeRcpp(values)]
        expect_identical(nativeSieve(bounds[1], bounds[2]), expected,
                         info = caseName)
    }
})

test_that("primeSieve handles uneven parallel workloads", {
    maxThreads <- get("pkgEnv", asNamespace("RcppAlgos"))$maxThreads
    skip_if(maxThreads < 2L, "At least two supported threads are required")

    # Seventeen segments supply enough work for both two and three workers
    # under PrimeSieveMain's four-segments-per-worker rule. Neither worker
    # count divides the segment count; the bounds also leave partial segments.
    segmentSize <- 32340
    lower <- segmentSize + 17
    upper <- lower + 17 * segmentSize + 101
    expected <- primeSieve(lower, upper, nThreads = 1)
    for (threads in c(2L, 3L)[c(2L, 3L) <= maxThreads]) {
        expect_identical(primeSieve(lower, upper, nThreads = threads), expected)
    }
})

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
})

test_that("primeSieve agrees between serial and parallel generation", {
    funTestPar <- function(b1, b2) {
        ser <- primeSieve(b1, b2)
        par <- primeSieve(b1, b2, nThreads = 2)
        all.equal(par, ser)
    }

    expect_true(funTestPar(1, 1e7))
    expect_true(funTestPar(1e12, 1e12 + 1e7))
    expect_true(funTestPar(1e15, 1e15 + 1e7))
})

test_that("primeSieve matches expensive reference counts", {
    # These active tests take more time and memory than the boundary checks.
    # Reference counts were confirmed with primesieve by Kim Walisch.
    gc()
    expect_equal(length(primeSieve(633318687598976, 633318707598976)), 587185)
    gc()
    expect_equal(length(primeSieve(12e8)), 60454705)
    gc()
    expect_equal(length(primeSieve(39582118599936, 39582718599936,
                                   nThreads = 2)), 19161558)
})
