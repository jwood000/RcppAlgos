context("testing primeCount")

test_that("parallel primeCount handles cache growth and reuse", {
    skip_if(RcppAlgos:::pkgEnv$maxThreads < 2L,
            "At least two available threads are needed")

    ## Run before the serial cases so a standalone run starts with a cold cache.
    ## Known counts avoid warming the cache with a serial reference calculation.
    inputs <- c(1e7, 1e8, 1e9, 1e10, 1e11, 1e12)
    counts <- c(664579, 5761455, 50847534, 455052511, 4118054813, 37607912018)
    order <- c(3L, 1L, 5L, 2L, 6L, 4L, 3L, 1L)

    for (nThreads in c(2L, 3L, 4L)) {
        actual <- vapply(inputs[order], primeCount, numeric(1),
                         nThreads = nThreads)
        expect_equal(actual, counts[order])
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
