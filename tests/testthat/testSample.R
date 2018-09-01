context("testing comboSample and permuteSample")

test_that("comboSample produces correct results", {
    
    expect_equal(comboSample(5, 3, sampleVec = 1:10), comboGeneral(5, 3))
    expect_equal(nrow(comboSample(5, 3, n = 5, seed = 2)), 5)
    expect_equal(as.vector(comboSample(1, 1, n = 1)), 1)
    
    totalCombs = comboCount(5, 8, freqs = rep(3, 5))
    expect_equal(comboSample(5, 8, freqs = rep(3, 5), sampleVec = 1:totalCombs), 
                 comboGeneral(5, 8, freqs = rep(3, 5)))
    
    totalCombs = comboCount(5, 8, TRUE)
    expect_equal(comboSample(5, 8, TRUE, sampleVec = 1:totalCombs), 
                 comboGeneral(5, 8, TRUE))
    
    expect_equal(comboCount(100, 25), gmp::chooseZ(100, 25))
    expect_equal(comboCount(50, 25, TRUE), gmp::chooseZ(50 + 25 - 1, 25))
    expect_equal(comboCount(100, 20), comboCount(100, 20, freqs = rep(1, 100)))
    
    expect_equal(comboSample(5, 3, TRUE, sampleVec = 1:35), comboGeneral(5, 3, TRUE))
    
    myNums = rnorm(5)
    expect_equal(comboSample(myNums, 3, sampleVec = 1:10), comboGeneral(myNums, 3))
    expect_equal(nrow(comboSample(myNums, 3, n = 5)), 5)
    
    totalCombs = comboCount(5, 8, freqs = rep(3, 5))
    expect_equal(comboSample(myNums, 8, freqs = rep(3, 5), sampleVec = 1:totalCombs), 
                 comboGeneral(myNums, 8, freqs = rep(3, 5)))
    
    totalCombs = comboCount(5, 8, TRUE)
    expect_equal(comboSample(myNums, 8, TRUE, sampleVec = 1:totalCombs), 
                 comboGeneral(myNums, 8, TRUE))
    
    totalCombs = comboCount(LETTERS, 4)
    samp = sample(totalCombs, 10)
    expect_equal(comboSample(LETTERS, 4, sampleVec = samp), 
             comboGeneral(LETTERS, 4)[samp, ])
    
    expect_equal(comboSample(myNums, 3, TRUE, sampleVec = 1:35), comboGeneral(myNums, 3, TRUE))
    
    mySamp = seq.int(1, comboCount(5, 5, freqs = c(2,1,2,1,2)), 3)
    expect_equal(comboSample(5, 5, freqs = c(2,1,2,1,2), sampleVec = mySamp),
                 comboGeneral(5, 5, freqs = c(2,1,2,1,2))[mySamp, ])
    
    expect_equal(sum(comboGeneral(c(T,F), 10, TRUE)), sum(1:10))
    
    expect_equal(comboSample(500, 200, n = 100, seed = 42),
                 comboSample(500, 200, n = 100, seed = 42))
    
    expect_equal(comboSample(500, 200, TRUE, n = 100, seed = 42),
                 comboSample(500, 200, TRUE, n = 100, seed = 42))
    
    ## **** MULTISETS ******
    lastRow = comboCount(100, 20, freqs = rep(1:4, 25))
    myChars = as.character(1:100)
    
    expect_equal(comboSample(as.character(1:100), 20, freqs = rep(1:4, 25), 
                        sampleVec = c(gmp::as.bigz(1), lastRow)), 
                 rbind(rep(as.character(1:8), times = c(1:4,1:4)),
                       rep(as.character(93:100), times = c(1:4,1:4))))
    
    expect_equal(comboSample(100, 20, freqs = rep(1:4, 25), n = 50, seed = 42),
                 comboSample(100, 20, freqs = rep(1:4, 25), n = 50, seed = 42))
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

    ## Count Test
    expect_equal(gmp::as.bigz(permuteCount(100, 25)),
                 gmp::as.bigz(gmp::div.bigz(gmp::factorialZ(100),
                                            gmp::factorialZ(75))))
    expect_equal(permuteCount(50, 25, TRUE), gmp::pow.bigz(50, 25))
    expect_equal(permuteCount(100, 20),
                 permuteCount(100, 20, freqs = rep(1, 100)))

    totalPerms = permuteCount(LETTERS, 3)
    samp = as.numeric(sample(totalPerms, 10))
    expect_equal(permuteSample(LETTERS, 3, sampleVec = samp),
                 permuteGeneral(LETTERS, 3)[samp, ])

    expect_equal(permuteSample(c(1,2,3,NA), 3, sampleVec = 1:24),
                 permuteGeneral(c(1,2,3,NA), 3))

    expect_equal(permuteSample(5, freqs = c(2,1,2,1,2),
                               sampleVec = 1:(permuteCount(5, freqs = c(2,1,2,1,2)))),
                 permuteGeneral(5, freqs = c(2,1,2,1,2)))

    expect_equal(matrix(as.integer(permuteSample(c(T,F),
                                                 freqs = c(50, 50),
                                                 n = 3, seed = 22)), nrow = 3),
                 permuteSample(1:0,
                               freqs = c(50, 50),
                               n = 3, seed = 22))

    expect_equal(nrow(permuteSample(c(T,F), 10, TRUE, n = 5, seed = 5)), 5)
    expect_equal(permuteSample(factor(state.name), 20, sampleVec = 1e12),
                 permuteGeneral(factor(state.name), 20, lower = 1e12, upper = 1e12))

    expect_equal(permuteSample(500, 100, n = 100, seed = 42),
                 permuteSample(500, 100, n = 100, seed = 42))

    expect_equal(permuteSample(100, 20, n = 100, seed = 42),
                 permuteSample(100, 20, n = 100, seed = 42))

    expect_equal(permuteSample(100, 20, TRUE, n = 100, seed = 42),
                 permuteSample(100, 20, TRUE, n = 100, seed = 42))

    expect_equal(permuteSample(75, 10, freqs = rep(1:3, 25), n = 20, seed = 42),
                 permuteSample(75, 10, freqs = rep(1:3, 25), n = 20, seed = 42))
})

test_that("comboSample produces correct results when FUN is applied", {
    
    num <- comboCount(80, 40) - 10
    vec <- do.call(c, lapply(0:10, function(x)  gmp::add.bigz(x, num)))
    expect_equal(comboSample(80, 40, FUN = sd, sampleVec = vec),
                 comboGeneral(80, 40, lower = num, FUN = sd))
    
    num <- comboCount(80, 40) - 10
    vec <- do.call(c, lapply(0:10, function(x)  gmp::add.bigz(x, num)))
    expect_equal(comboSample(80, 40, FUN = sd, sampleVec = vec),
                 comboGeneral(80, 40, lower = num, FUN = sd))
    
    num <- comboCount(80, 40) - 10
    vec <- do.call(c, lapply(0:10, function(x)  gmp::add.bigz(x, num)))
    expect_equal(comboSample(80, 40, FUN = sd, sampleVec = vec),
                 comboGeneral(80, 40, lower = num, FUN = sd))
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
                 "sampleVec must be convertible to a real number")
    expect_error(permuteSample(100, 10, n = 2^50),
                 "The number of rows cannot exceed")
    expect_error(permuteSample(100, 10, n = -200),
                 "n must be a positive number")
    expect_error(permuteSample(100, 10, sampleVec = c("62815650955529472001")),
                 "One or more of the requested values in sampleVec")
})