context("testing comboSample and permuteSample")

test_that("comboSample produces correct results", {
    
    expect_equal(comboSample(5, 3, sampleVec = 1:10), comboGeneral(5, 3))
    set.seed(2)
    expect_equal(nrow(comboSample(5, 3, n = 5)), 5)
    expect_equal(as.vector(comboSample(1, 1, n = 1)), 1)
    
    totalCombs <- comboCount(5, 8, freqs = rep(3, 5))
    expect_equal(comboSample(5, 8, freqs = rep(3, 5), sampleVec = 1:totalCombs), 
                 comboGeneral(5, 8, freqs = rep(3, 5)))
    
    totalCombs <- comboCount(5, 8, TRUE)
    expect_equal(comboSample(5, 8, TRUE, sampleVec = 1:totalCombs), 
                 comboGeneral(5, 8, TRUE))
    
    expect_equal(comboSample(5, 3, TRUE, sampleVec = 1:35), comboGeneral(5, 3, TRUE))
    
    myNums <- rnorm(5)
    expect_equal(comboSample(myNums, 3, sampleVec = 1:10), comboGeneral(myNums, 3))
    expect_equal(nrow(comboSample(myNums, 3, n = 5)), 5)
    
    totalCombs <- comboCount(5, 8, freqs = rep(3, 5))
    expect_equal(comboSample(myNums, 8, freqs = rep(3, 5), sampleVec = 1:totalCombs), 
                 comboGeneral(myNums, 8, freqs = rep(3, 5)))
    
    totalCombs <- comboCount(5, 8, TRUE)
    expect_equal(comboSample(myNums, 8, TRUE, sampleVec = 1:totalCombs), 
                 comboGeneral(myNums, 8, TRUE))
    
    totalCombs <- comboCount(LETTERS, 4)
    samp <- sample(totalCombs, 10)
    expect_equal(comboSample(LETTERS, 4, sampleVec = samp), 
             comboGeneral(LETTERS, 4)[samp, ])
    
    expect_equal(comboSample(myNums, 3, TRUE, sampleVec = 1:35), comboGeneral(myNums, 3, TRUE))
    
    mySamp <- seq.int(1, comboCount(5, 5, freqs = c(2,1,2,1,2)), 3)
    expect_equal(comboSample(5, 5, freqs = c(2,1,2,1,2), sampleVec = mySamp),
                 comboGeneral(5, 5, freqs = c(2,1,2,1,2))[mySamp, ])
})

test_that("permuteSample produces correct results", {
    expect_equal(permuteSample(5, 3, sampleVec = 1:(permuteCount(5, 3))), permuteGeneral(5, 3))
    set.seed(2)
    expect_equal(nrow(permuteSample(5, 3, n = 5)), 5)
    expect_equal(as.vector(permuteSample(1, 1, n = 1)), 1)
    
    expect_equal(permuteSample(5, 6, freqs = rep(2, 5), 
                             sampleVec = 1:(permuteCount(5, 6, freqs = rep(2, 5)))), 
                 permuteGeneral(5, 6, freqs = rep(2, 5)))
    
    expect_equal(permuteSample(5, 6, TRUE, freqs = rep(2, 5), 
                             sampleVec = 1:(permuteCount(5, 6, freqs = rep(2, 5)))), 
                 permuteGeneral(5, 6, freqs = rep(2, 5)))
    
    expect_equal(permuteSample(5, 3, TRUE, 
                               sampleVec = 1:(permuteCount(5, 3, T))), 
                 permuteGeneral(5, 3, TRUE))
    
    totalCombs <- permuteCount(LETTERS, 3)
    samp <- as.numeric(sample(totalCombs, 10))
    expect_equal(permuteSample(LETTERS, 3, sampleVec = samp), 
                 permuteGeneral(LETTERS, 3)[samp, ])
    
    expect_equal(permuteSample(c(1,2,3,NA), 3, sampleVec = 1:24),
                 permuteGeneral(c(1,2,3,NA), 3))
    
    expect_equal(permuteSample(5, freqs = c(2,1,2,1,2), 
                               sampleVec = 1:(permuteCount(5, freqs = c(2,1,2,1,2)))), 
                 permuteGeneral(5, freqs = c(2,1,2,1,2)))
})

test_that("comboSample produces appropriate error messages", {
    expect_error(comboSample(5, 3), "n and sampleVec cannot both be NULL")
    expect_error(comboSample(5,3,freqs = c(1,2,3,-2,1)), "in freqs must be a positive")
    expect_error(comboSample(5,3, n = 100), "n exceeds the maximum number of possible results")
    expect_error(comboSample(5,3, comboSample(5,3, sampleVec = 1:200)), "exceeds the maximum number of possible results")
    expect_error(comboSample(5,freqs = rep(1,6)), "the length of freqs must equal the")
    
    expect_error(comboSample(5,3, n = "5"), 
                 "n must be a number")
    expect_error(comboSample(5,3, n = 1:5), 
                 "length of n must be 1")
})

test_that("permuteSample produces appropriate error messages", {
    expect_error(permuteSample(5, 3), 
                 "n and sampleVec cannot both be NULL")
    expect_error(permuteSample(5,3,freqs = c(1,2,3,-2,1)), 
                 "in freqs must be a positive")
    expect_error(permuteSample(5,3, n = 100), 
                 "n exceeds the maximum number of possible results")
    expect_error(permuteSample(5,3, n = "5"), 
                 "n must be a number")
    expect_error(permuteSample(5,3, n = 1:5), 
                 "length of n must be 1")
    expect_error(permuteSample(5,3, permuteSample(5,3, sampleVec = 1:200)), 
                 "exceeds the maximum number of possible results")
    expect_error(permuteSample(5,freqs = rep(1,6)), 
                 "the length of freqs must equal the")
    expect_error(permuteSample(5, 4, sampleVec = "adfs"), 
                 "sampleVec must be of type numeric or integer")
})