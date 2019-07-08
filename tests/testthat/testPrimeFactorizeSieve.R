context("testing primeFactorizeSieve")

test_that("primeFactorizeSieve generates correct numbers", {
    options(scipen = 50)
    expect_equal(primeFactorizeSieve(100)[[100]], c(2, 2, 5, 5))
    expect_equal(length(primeFactorizeSieve(1000)), 1000)
    expect_equal(primeFactorizeSieve(2)[[2]], 2)
    expect_equal(primeFactorizeSieve(2, 2)[[1]], 2)
    expect_equal(primeFactorizeSieve(1000, 1000)[[1]], c(2,2,2,5,5,5))
    expect_equal(length(primeFactorizeSieve(1L, namedList = TRUE)[[1]]), 0)
    expect_equal(length(primeFactorizeSieve(1L, 1L, namedList = TRUE)[[1]]), 0)
    
    a <- primeFactorizeSieve(1e12+10, 1e12, namedList = TRUE)
    expect_true(all(sapply(a, prod) == as.numeric(names(a))))
    
    ## Test Names
    expect_equal(as.integer(names(primeFactorizeSieve(100, namedList = TRUE))), 1:100)
    expect_equal(as.numeric(names(primeFactorizeSieve(10^12, 10^12 + 100,
                                          namedList = TRUE))), (10^12):(10^12 + 100))
    
    ## Test Parallel
    expect_equal(primeFactorizeSieve(117, 5e4), 
                    primeFactorizeSieve(117, 5e4, nThreads = 2))
    
    ## two threads will only be used
    expect_equal(primeFactorizeSieve(1e7, 1e7 + 25000), 
                 primeFactorizeSieve(1e7, 1e7 + 25000, nThreads = 3))
    
    expect_equal(primeFactorizeSieve(1e12, 1e12 + 2e4), 
                    primeFactorizeSieve(1e12, 1e12 + 2e4, nThreads = 2))
})
