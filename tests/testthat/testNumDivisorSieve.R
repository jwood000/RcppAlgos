context("testing numDivisorSieve")

test_that("numDivisorSieve generates correct numbers", {
    options(scipen = 50)
    expect_equal(numDivisorSieve(10)[10], 4)
    expect_equal(length(numDivisorSieve(1000)), 1000)
    expect_equal(numDivisorSieve(2), 1:2)
    expect_equal(numDivisorSieve(1), 1)
    expect_equal(numDivisorSieve(99, 100), c(6, 9))
    
    expect_equal(numDivisorSieve(100L), 
                 sapply(1:100, function(x) length((1:x)[x %% (1:x) == 0])))
    
    ## lower bound less than sqrt(100) and greater than 1
    expect_equal(numDivisorSieve(100L, 5), 
                 sapply(5:100, function(x) length((1:x)[x %% (1:x) == 0])))
    
    expect_true(all.equal(numDivisorSieve(1000000L, 1000005L, TRUE), 
                 sapply(1000000:1000005, function(x) length((1:x)[x %% (1:x) == 0])), 
                 check.attributes = FALSE))
    
    expect_true(numDivisorSieve(1, namedVector = TRUE)==1)
    expect_true(numDivisorSieve(1, 1, TRUE)==1)
    
    ## Test Names
    expect_equal(as.integer(names(numDivisorSieve(100, namedVector = TRUE))), 1:100)
    expect_equal(as.numeric(names(numDivisorSieve(10^12, 10^12 + 100,
                                                  namedVector = TRUE))), (10^12):(10^12 + 100))
    
    expect_equal(numDivisorSieve(1e6), numDivisorSieve(1e6, nThreads = 2))
    expect_equal(numDivisorSieve(1e12, 1e12 + 1e5), 
                 numDivisorSieve(1e12, 1e12 + 1e5, nThreads = 2))
})
