context("testing divisorsRcpp")

test_that("divisorsRcpp generates correct numbers", {
    options(scipen = 50)
    expect_equal(divisorsRcpp(0), numeric(0))
    expect_equal(length(divisorsRcpp(1:100)), 100)
    expect_equal(divisorsRcpp(2), c(1, 2))
    expect_equal(divisorsRcpp(1), 1)
    expect_equal(divisorsRcpp(-1), c(-1, 1))
    expect_equal(divisorsRcpp(-10), c(-10, -5, -2, -1, 1, 2, 5, 10))
    
    expect_equal(divisorsRcpp(-10:10), 
                 lapply(-10:10, function(x) {
                        if (x == 0) {
                            div <- vector("integer")
                        } else {
                            y <- abs(x)
                            div <- (1:y)[y %% (1:y) == 0]
                            if (x < 0)
                                div <- c(-1L * rev(div), div)
                        }
                        div
                     }))
    
    expect_equal(divisorsRcpp((1e6):(1e6 + 1e3)), divisorsSieve(1e6, 1e6 + 1e3))
    
    expect_equal(divisorsRcpp(-1*1e10), c(-1* rev(divisorsRcpp(1e10)), divisorsRcpp(1e10)))
    
    expect_equal(divisorsRcpp(1000), c(1,2,4,5,8,10,20,
                                       25,40,50,100,125,
                                       200,250,500,1000))
    expect_equal(divisorsRcpp(1000, TRUE), c(1,2,4,5,8,10,20,
                                             25,40,50,100,125,
                                             200,250,500,1000))
    
    ## Test Names
    expect_equal(as.integer(names(divisorsRcpp(100, namedList = TRUE))), integer(0))
    expect_equal(as.numeric(names(divisorsRcpp((10^12):(10^12 + 100),
                                                 namedList = TRUE))), (10^12):(10^12 + 100))
    
    ## Test Parallel
    mySample = sample(10^9, 20)
    expect_equal(divisorsRcpp(mySample), divisorsRcpp(mySample, nThreads = 2))
    options(scipen = 0)
})
