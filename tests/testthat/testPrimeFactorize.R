context("testing primeFactorize")

test_that("primeFactorize generates correct numbers", {
    options(scipen = 50)
    expect_equal(primeFactorize(100), c(2, 2, 5, 5))
    expect_equal(length(primeFactorize(1:100)), 100)
    expect_equal(primeFactorize(2), 2)
    expect_equal(primeFactorize(1), numeric(0))
    expect_equal(primeFactorize(1000), c(2,2,2,5,5,5))
    expect_equal(primeFactorize(1000, TRUE), c(2,2,2,5,5,5))
    
    expect_equal(primeFactorize(0), integer(0))
    
    a <- primeFactorize(-10:10)
    b <- primeFactorizeSieve(10)
    b1 <- lapply(rev(b), function(x) c(-1L, x))
    expect_equal(a, c(b1, list(integer(0)), b))
    
    expect_equal(rle(primeFactorize(-1e10))$lengths, c(1, 10, 10))
    expect_equal(rle(primeFactorize(1e15))$lengths, c(15, 15))
    
    expect_equal(primeFactorize((1e6):(1e6 + 1e4)), primeFactorizeSieve(1e6, 1e6 + 1e4))
    expect_equal(primeFactorize((1e12):(1e12 + 1e2)), primeFactorizeSieve(1e12, 1e12 + 1e2))
    
    ## Test Names
    expect_equal(as.integer(names(primeFactorize(100, namedList = TRUE))), integer(0))
    expect_equal(as.numeric(names(primeFactorize((10^12):(10^12 + 100),
                                                      namedList = TRUE))), (10^12):(10^12 + 100))
    ## Test Parallel
    set.seed(567)
    samp <- sample(1e12, 100)
    expect_equal(primeFactorize(samp), primeFactorize(samp, nThreads = 2))
})
