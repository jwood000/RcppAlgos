context("testing primeCount")

test_that("parallel Legendre handles cache growth and reuse", {
    skip_if(RcppAlgos:::pkgEnv$maxThreads < 2L,
            "At least two available threads are needed")

    ## Run before the serial cases so a standalone run starts with a cold cache.
    ## Known counts avoid warming the cache with a serial reference calculation.
    inputs <- c(1e7, 1e8, 1e9, 1e10, 1e11, 1e12)
    counts <- c(664579, 5761455, 50847534, 455052511, 4118054813, 37607912018)
    order <- c(3L, 1L, 5L, 2L, 6L, 4L, 3L, 1L)

    for (nThreads in c(2L, 3L, 4L)) {
        actual <- vapply(inputs[order], primeCount, numeric(1),
                         nThreads = nThreads, method = "legendre")
        expect_equal(actual, counts[order])
    }
})

test_that("parallel Legendre handles the batch boundary and repeated calls", {
    maxThreads <- RcppAlgos:::pkgEnv$maxThreads
    skip_if(maxThreads < 2L, "At least two available threads are needed")

    ## Cross the single/multiple-batch boundary, then repeat the multiple-batch call.
    inputs <- c(1e10 + 1, 1e10 - 1, 1e10, 1e10 + 1)
    threadCounts <- c(2L, 3L, 4L)
    threadCounts <- threadCounts[threadCounts <= maxThreads]

    for (nThreads in threadCounts) {
        actual <- vapply(inputs, primeCount, numeric(1),
                         nThreads = nThreads, method = "legendre")
        expect_equal(actual, rep(455052511, length(inputs)))
    }
})

test_that("primeCount generates correct numbers", {
    expect_equal(primeCount(10), 4)

    ## Pi(x) obtained from https://en.wikipedia.org/wiki/Prime-counting_function
    expect_equal(primeCount(10^6), 78498L)
    expect_equal(primeCount(as.integer(1e8)), 5761455L)
    expect_equal(primeCount(10^10), 455052511)
    expect_equal(primeCount(10^10, nThreads = 2), 455052511)
    expect_equal(primeCount(1e11, nThreads = 2), 4118054813)
    expect_equal(primeCount(1e12, nThreads = 2), 37607912018)
    expect_equal(primeCount(1), 0)
    expect_equal(sapply(1:9, primeCount), c(0,1,2,2,3,3,4,4,4))
})

test_that("primeCount methods agree with an independent sieve at boundaries", {
    methods <- c("auto", "legendre", "meissel", "lehmer")
    powers <- as.vector(outer(c(2, 3, 17, 37), 2:4, "^"))
    inputs <- sort(unique(c(1:12, 99999:100001,
                            as.vector(outer(powers, -1:1, "+")))))

    ## A base-R Eratosthenes sieve supplies the expected counts independently
    ## of primeCount, primeSieve, and the shared native Phi infrastructure.
    sieve <- rep(TRUE, max(inputs))
    sieve[1] <- FALSE
    for (p in 2:floor(sqrt(length(sieve)))) {
        if (sieve[p]) sieve[seq.int(p * p, length(sieve), by = p)] <- FALSE
    }
    expected <- as.integer(cumsum(sieve)[inputs])

    for (method in methods) {
        actual <- vapply(inputs, primeCount, integer(1), method = method)
        expect_identical(actual, expected)
    }
})

test_that("primeCount methods preserve known counts and result types", {
    inputs <- list(1000000L, 1e10, 1e11)
    expected <- list(78498L, 455052511L, 4118054813)

    for (method in c("auto", "legendre", "meissel", "lehmer")) {
        expect_identical(lapply(inputs, primeCount, method = method), expected)
    }
})

test_that("auto primeCount preserves counts at method boundaries", {
    inputs <- c(1e10 - 1, 1e10, 1e10 + 1, 1e12 - 1, 1e12, 1e12 + 1)
    expected <- c(rep(list(455052511L), 3), rep(list(37607912018), 3))
    threads <- c(1L, 2L, 6L)
    threads <- threads[threads <= RcppAlgos:::pkgEnv$maxThreads]

    ## The lower neighbors are divisible by 9; the upper neighbors are
    ## divisible by 101 and 73, respectively. Thus the known counts at the
    ## powers of ten also apply on either side of each dispatch boundary.
    for (nThreads in threads) {
        actual <- lapply(inputs, primeCount, nThreads = nThreads)
        expect_identical(actual, expected)
    }
})

test_that("primeCount appends method without changing positional calls", {
    inputs <- list(1, 10, 99999, 100000, 1e10, 1e11)
    expect_identical(lapply(inputs, primeCount),
                     lapply(inputs, primeCount, method = "auto"))
    expect_identical(primeCount(1e10, 2L), 455052511L)

    for (method in c("auto", "legendre", "meissel", "lehmer")) {
        expect_identical(primeCount(1e10, 2L, method), 455052511L)
    }
})

test_that("primeCount uses standard method matching", {
    abbreviations <- c("a", "leg", "m", "leh")
    actual <- lapply(abbreviations, function(method) {
        primeCount(1e6, method = method)
    })
    expect_identical(actual, rep(list(78498L), length(abbreviations)))
    expect_identical(primeCount(1e6, method = NULL),
                     primeCount(1e6, method = "auto"))

    invalidMethods <- list("le", "l", "unknown", "Legendre", "", NA_character_,
                           character(), 1, TRUE, c("auto", "lehmer"))
    for (method in invalidMethods) {
        expect_error(primeCount(1e6, method = method), "arg")
    }
})

test_that("explicit primeCount methods support parallel calls and reuse", {
    maxThreads <- RcppAlgos:::pkgEnv$maxThreads
    skip_if(maxThreads < 2L, "At least two available threads are needed")
    threadCounts <- unique(pmin(c(2L, 4L), as.integer(maxThreads)))
    inputs <- c(1e10 + 1, 1e7, 1e11, 1e10 + 1)
    expected <- c(455052511, 664579, 4118054813, 455052511)

    for (method in c("legendre", "meissel", "lehmer")) {
        for (nThreads in threadCounts) {
            actual <- vapply(inputs, primeCount, numeric(1),
                             nThreads = nThreads, method = method)
            expect_identical(actual, expected)
        }
    }
})

test_that("Meissel and Lehmer preserve larger known counts across threads", {
    threads <- unique(c(pmin(4L, as.integer(RcppAlgos:::pkgEnv$maxThreads)), 1L))
    for (method in c("meissel", "lehmer")) {
        for (nThreads in threads) {
            expect_identical(primeCount(1e12, nThreads = nThreads, method = method),
                             37607912018)
        }
    }
})

test_that("primeCount methods retain input validation", {
    invalidInputs <- list(0, -1, 2^53, 1.5, NA_real_, Inf, "100000")
    messages <- c("positive", "positive", "less than", "whole number",
                  "NA or NaN", "less than", "numeric or integer")

    for (method in c("auto", "legendre", "meissel", "lehmer")) {
        for (i in seq_along(invalidInputs)) {
            expect_error(primeCount(invalidInputs[[i]], method = method), messages[i])
        }
        expect_identical(primeCount(100, method = method), 25L)
    }
})

test_that("primeCount methods retain the thread-validation threshold", {
    invalidThreads <- list(0L, -1L, 1.5, NA_integer_, Inf, "two", numeric())

    for (method in c("auto", "legendre", "meissel", "lehmer")) {
        ## Historically these inputs ignore nThreads and execute serially.
        actual <- vapply(invalidThreads, function(nThreads) {
            primeCount(99999, nThreads = nThreads, method = method)
        }, integer(1))
        expect_identical(actual, rep(9592L, length(invalidThreads)))

        for (nThreads in invalidThreads) {
            expect_error(primeCount(100000, nThreads = nThreads, method = method),
                         "nThreads")
        }
        expect_identical(primeCount(100000, nThreads = .Machine$integer.max,
                                   method = method), 9592L)
    }
})
