context("testing comboGeneral")

test_that("comboGeneral produces correct results with no constraints", {
    
    expect_equal(comboGeneral(5, 3), t(combn(5, 3)))
    
    expect_equal(nrow(comboGeneral(6, 3)), choose(6, 3))
    expect_equal(as.vector(comboGeneral(1,1)), 1)
    expect_equal(as.vector(comboGeneral(1,1,TRUE)), 1)
    
    expect_equal(sum(comboGeneral(3, 3, freqs = c(1, 1, 1))),
                 sum(comboGeneral(3, 3)))
    
    expect_equal(as.vector(comboGeneral(1, 2, freqs = 2)), c(1, 1))
    
    expect_equal(ncol(comboGeneral(5, 3)), 3)
    expect_equal(ncol(comboGeneral(5, 3, TRUE)), 3)
    expect_equal(ncol(comboGeneral(5, 3, FALSE, "prod", keepResults = TRUE)), 4)
    expect_equal(ncol(comboGeneral(5, 3, TRUE, "prod", keepResults = TRUE)), 4)
    
    expect_equal(ncol(comboGeneral(5, 3, FALSE,
                                     "prod", freqs = c(1,2,1,2,4), 
                                     keepResults = TRUE)), 4)
    
    expect_equal(nrow(comboGeneral(10, 3, TRUE, rowCap = 10)), 10)
    expect_equal(nrow(comboGeneral(10, 3, rowCap = 10)), 10)
    
    expect_equal(nrow(comboGeneral(10, 3, FALSE, "prod", 
                                     keepResults = TRUE, rowCap = 10)), 10)
    
    expect_equal(nrow(comboGeneral(5, 5, freqs = 1:5, rowCap = 10)), 10)
    
    expect_equal(nrow(comboGeneral(10, 3, TRUE, "prod", 
                                     keepResults = TRUE, rowCap = 10)), 10)
})

test_that("comboGeneral produces correct results with constraints", {
    
    expect_equal(nrow(comboGeneral(-5:5, 4, FALSE, "sum", "==", 6)), 
                 length(which(apply(combn(-5:5, 4), 2, sum)==6)))
    
    expect_equal(unique(comboGeneral(5, 5, TRUE,
                                       "sum", "==", 9, 
                                       keepResults = TRUE)[,6]), 9)
    
    expect_true(all(comboGeneral(5, 5, TRUE,
                                   "min", "<", 3, 
                                   keepResults = TRUE)[,6] < 3))
    
    expect_true(all(comboGeneral(5, 5, TRUE,
                                   "prod", ">", 100, 
                                   keepResults = TRUE)[,6] > 100))
    
    expect_true(all(comboGeneral(5, 3, FALSE,
                                   "max", "=<", 4, 
                                   keepResults = TRUE)[,4] <= 4))
    
    expect_true(all(comboGeneral(3, 5, TRUE,
                                   "mean", ">=", 2, 
                                   keepResults = TRUE)[,6] >= 2))
    
    expect_true(all(comboGeneral(5, 5, FALSE, "sum", ">", 18,
                                   freqs = c(1,2,1,2,4), 
                                   keepResults = TRUE)[,6] > 18))
})

test_that("comboGeneral produces appropriate error messages", {
    expect_error(comboGeneral(9,4,TRUE,"summ","<",10), "prod, sum, mean, max, or min")
    expect_error(comboGeneral(9,4,TRUE,"sum","=<>",10), ">, >=, <, <=, or ==")
    expect_error(comboGeneral(9,4,TRUE,"sum",60,10), "must be passed as a character")
    expect_error(comboGeneral(9,4,FALSE,sum,"<",10), "must be passed as a character")
    expect_error(comboGeneral(9,4,TRUE,"sum","<",10,-1), "must be positive")
    expect_error(comboGeneral(170,7,FALSE,"sum","<",100), "The number of rows cannot exceed")
    expect_error(comboGeneral(170,7,FALSE,"sum","<",100,10^10), "number of rows cannot exceed")
    
    expect_error(comboGeneral(5,3,freqs = c(1,2,3,"2",1)), "freqs must be of type numeric")
    expect_error(comboGeneral(5,3,freqs = c(1,2,3,-2,1)), "in freqs must be a positive")
    expect_error(comboGeneral(5,1000,freqs = rep(5000, 5)), "number of rows cannot exceed")
    expect_error(comboGeneral(5,freqs = rep(1,6)), "the length of freqs must equal the")
})