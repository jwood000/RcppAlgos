context("testing divisorsSieve")

test_that("divisorsSieve generates correct numbers", {
    options(scipen = 50)
    expect_equal(divisorsSieve(10)[[10]], c(1, 2, 5, 10))
    expect_equal(length(divisorsSieve(1000)), 1000)
    expect_equal(divisorsSieve(1000, 1009)[[10]], c(1, 1009))
    expect_equal(divisorsSieve(1)[[1]], 1)
    expect_equal(divisorsSieve(2, 2)[[1]], c(1, 2))
    expect_equal(divisorsSieve(997, 997)[[1]], c(1, 997))
    expect_equal(divisorsSieve(1000, 1000)[[1]], c(1,2,4,5,8,10,20,
                                                   25,40,50,100,125,
                                                   200,250,500,1000))
    expect_equal(divisorsSieve(100L), 
                 lapply(1:100, function(x) (1:x)[x %% (1:x) == 0]))
    
    ## lower bound less than sqrt(100) and greater than 1
    expect_equal(divisorsSieve(100L, 5), 
                 lapply(5:100, function(x) (1:x)[x %% (1:x) == 0]))
    
    expect_true(divisorsSieve(1, namedList = TRUE) == 1)
    expect_true(divisorsSieve(1, 1, TRUE) == 1)
    
    ## Test against brute force
    expect_equal(divisorsSieve(1000000L, 1000005L), 
                 lapply(1000000:1000005, function(x) (1:x)[x %% (1:x) == 0]))
    
    expect_equal(divisorsSieve(1e12, 1e12 + 1e2), divisorsRcpp((1e12):(1e12 + 1e2)))
    
    ## Test Names
    expect_equal(as.integer(names(divisorsSieve(100, namedList = TRUE))), 1:100)
    expect_equal(as.numeric(names(divisorsSieve(10^12, 10^12 + 100,
                                                namedList = TRUE))), (10^12):(10^12 + 100))
    
    ## Test Parallel
    expect_equal(divisorsSieve(1e9, 1e9 + 25000),
                 divisorsSieve(1e9, 1e9 + 25000, nThreads = 3))
    expect_equal(divisorsSieve(3e4), divisorsSieve(3e4, nThreads = 2))
})
