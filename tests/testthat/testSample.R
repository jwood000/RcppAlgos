context("testing comboSample and permuteSample")

test_that("comboSample produces correct results", {
    
    expect_equal(comboSample(5, 3, sampleVec = 1:10), comboGeneral(5, 3))
    expect_equal(nrow(comboSample(5, 3, n = 5, seed = 2)), 5)
    expect_equal(as.vector(comboSample(1, 1, n = 1)), 1)
    expect_equal(class(comboSample(c(1:5, NA), 5, n = 3)[1, ]), "numeric")
    
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
    
    expect_equal(comboSample(500, 200, n = 10, seed = 42),
                 comboSample(500, 200, n = 10, seed = 42, nThreads = 2))
    
    expect_equal(comboSample(500, 200, TRUE, n = 10, seed = 42),
                 comboSample(500, 200, TRUE, n = 10, seed = 42, nThreads = 2))
    
    ## **** MULTISETS ******
    lastRow = comboCount(100, 20, freqs = rep(1:4, 25))
    
    expect_equal(comboSample(as.character(1:100), 20, freqs = rep(1:4, 25), 
                        sampleVec = c(gmp::as.bigz(1), lastRow)), 
                 rbind(rep(as.character(1:8), times = c(1:4,1:4)),
                       rep(as.character(93:100), times = c(1:4,1:4))))
    
    expect_equal(comboSample(100, 20, freqs = rep(1:4, 25), sampleVec = seq(1L, 100L, 10L)),
                 comboGeneral(100, 20, freqs = rep(1:4, 25), upper = 100)[seq(1L, 100L, 10L), ])
    
    set.seed(123)
    v <- rnorm(100)
    expect_equal(comboSample(v, 20, freqs = rep(1:4, 25), n = 10, seed = 42),
                 comboSample(v, 20, freqs = rep(1:4, 25), n = 10, seed = 42, nThreads = 2))
    
    ## only used 2 threads
    expect_equal(comboSample(as.factor(LETTERS), 12, TRUE, n = 2, seed = 17),
                 comboSample(as.factor(LETTERS), 12, TRUE, n = 2, seed = 17, Parallel = TRUE))
})

test_that("permuteSample produces correct results", {
    expect_equal(permuteSample(5, 3, sampleVec = 1:(permuteCount(5, 3))), permuteGeneral(5, 3))
    set.seed(2)
    expect_equal(nrow(permuteSample(5, 3, n = 5)), 5)
    expect_equal(as.vector(permuteSample(1, 1, n = 1)), 1)
    expect_equal(ncol(permuteSample(5, 20, freqs = 1:5, n = 3)), sum(1:5))

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

    expect_equal(permuteSample(500, 100, n = 10, seed = 42),
                 permuteSample(500, 100, n = 10, seed = 42, nThreads = 2))

    expect_equal(permuteSample(100, 20, n = 10, seed = 42),
                 permuteSample(100, 20, n = 10, seed = 42, nThreads = 2))

    expect_equal(permuteSample(100, 20, TRUE, n = 10, seed = 42),
                 permuteSample(100, 20, TRUE, n = 10, seed = 42, nThreads = 2))

    expect_equal(permuteSample(75, 10, freqs = rep(1:3, 25), n = 10, seed = 42),
                 permuteSample(75, 10, freqs = rep(1:3, 25), n = 10, seed = 42, nThreads = 2))
    
    expect_equal(permuteSample(c(TRUE, FALSE), 20, freqs = c(10, 15), seed = 97, n = 20), 
                 permuteSample(c(TRUE, FALSE), 20, freqs = c(10, 15), seed = 97, n = 20, nThreads = 2))
})

test_that("comboSample produces correct results when FUN is applied", {
    
    num <- comboCount(80, 40) - 10
    vec <- do.call(c, lapply(0:10, function(x)  gmp::add.bigz(x, num)))
    expect_equal(comboSample(80, 40, FUN = sd, sampleVec = vec),
                 comboGeneral(80, 40, lower = num, FUN = sd))
    
    expect_equal(apply(comboGeneral(letters[1:8], 5), 1, paste0, collapse = "")[30:40],
                 unlist(comboSample(letters[1:8], 5, sampleVec = 30:40,
                             FUN = function(x) paste0(x, collapse = ""))))
    
    expect_equal(comboSample(c(NA, 1, 2, 5, 10), 3, FUN = sum, sampleVec = 1:10),
                 comboGeneral(c(NA, 1, 2, 5, 10), 3, FUN = sum))
    
    set.seed(20)
    vec <- runif(10, -1e4, 1e4)
    num <- comboCount(10, 7, TRUE)
    samp <- sample(num, 100)
    expect_equal(comboSample(vec, 7, TRUE, sampleVec = samp, FUN = function(x) {
        mean(dcauchy(x))
    }), comboGeneral(vec, 7, TRUE, FUN = function(x) mean(dcauchy(x)))[samp])

    set.seed(19)
    vec <- runif(8, -1e4, 1e4)
    num <- comboCount(8, 8, freqs = rep(1:4, 2))
    samp <- sample(num, 100)
    expect_equal(comboSample(vec, 8, freqs = rep(1:4, 2), FUN = cospi, sampleVec = samp),
                 comboGeneral(vec, 8, freqs = rep(1:4, 2), FUN = cospi)[samp])
})

test_that("permuteSample produces correct results when FUN is applied", {
    
    num <- permuteCount(80, 40) - 10
    vec <- do.call(c, lapply(0:10, function(x)  gmp::add.bigz(x, num)))
    expect_equal(permuteSample(80, 40, FUN = sd, sampleVec = vec),
                 permuteGeneral(80, 40, lower = num, FUN = sd))
    
    set.seed(18)
    vec <- runif(6, -1e4, 1e4)
    num <- permuteCount(6, 5, TRUE)
    samp <- sample(num, 100)
    expect_equal(permuteSample(vec, 7, TRUE, sampleVec = samp, FUN = function(x) {
        mean(dcauchy(x))
    }), permuteGeneral(vec, 7, TRUE, FUN = function(x) mean(dcauchy(x)))[samp])
    
    set.seed(17)
    vec <- runif(8, -1e4, 1e4)
    num <- permuteCount(8, 5, freqs = rep(1:4, 2))
    samp <- sample(num, 100)
    expect_equal(permuteSample(vec, 5, freqs = rep(1:4, 2), FUN = cospi, sampleVec = samp),
                 permuteGeneral(vec, 5, freqs = rep(1:4, 2), FUN = cospi)[samp])
})

test_that("comboSample produces appropriate error messages", {
    expect_error(comboSample(5, 3), "n and sampleVec cannot both be NULL")
    expect_error(comboSample(5,3,freqs = c(1,2,3,-2,1)), "in freqs must be a positive")
    expect_error(comboSample(5,3, n = 100), "n exceeds the maximum number of possible results")
    expect_error(comboSample(5,3, comboSample(5,3, sampleVec = 1:200)), "exceeds the maximum number of possible results")
    expect_error(comboSample(5,freqs = rep(1,6)), "the length of freqs must equal the")
    expect_error(comboSample(5,3, n = 5, FUN = "sum"), "FUN must be a function!")
    expect_error(comboSample(5,3, n = "5"), "n must be of type numeric or integer")
    expect_error(comboSample(5,3, n = 1:5), "length of n must be 1")
    expect_error(comboSample(100000, 10, sampleVec = -1L),
                 "sampleVec must be a positive number")
    expect_error(comboSample(100000, 10, sampleVec = NaN),
                 "sampleVec cannot be NA or NaN")
    expect_error(comboSample(100000, 10, sampleVec = 2^1000),
                 "Number is too large for double precision. Consider using gmp::as.bigz or as.character for sampleVec")
    expect_error(comboSample(100000, 10, sampleVec = as.complex(1)),
                 "This type is not supported! No conversion possible for sampleVec")
    
    expect_error(comboSample(1000, 100, sampleVec = gmp::as.bigz(c(1, -1))),
                 "Each element in sampleVec must be a positive number")
    expect_error(comboSample(1000, 100, sampleVec = NA_character_), "sampleVec cannot be NA or NaN")
    expect_error(comboSample(1000, 100, sampleVec = 2^53), 
                 "Number is too large for double precision. Consider using gmp::as.bigz or as.character for sampleVec")
    expect_error(comboSample(1000, 15, sampleVec = gmp::as.bigz(c(NA, 1))),
                 "Each element in sampleVec cannot be NA or NaN")
})

test_that("permuteSample produces appropriate error messages", {
    expect_error(permuteSample(5, 6, n = 3), "m must be less than or equal to the length of v")
    expect_error(permuteSample(5, repetition = TRUE), "m and freqs cannot both be NULL")
    expect_error(permuteSample(5, 1:5), "length of m must be 1")
    expect_error(permuteSample(5, -4), "m must be a positive number")
    expect_error(permuteSample(5, 3), "n and sampleVec cannot both be NULL")
    expect_error(permuteSample(5,3,freqs = c(1,2,3,-2,1)), "in freqs must be a positive")
    expect_error(permuteSample(5,3, n = 100), 
                 "n exceeds the maximum number of possible results")
    expect_error(permuteSample(5,3, n = "5"), "n must be of type numeric or integer")
    expect_error(permuteSample(5,3, n = 1:5), "length of n must be 1")
    expect_error(permuteSample(5,3, permuteSample(5,3, sampleVec = 1:200)), 
                 "exceeds the maximum number of possible results")
    expect_error(permuteSample(5,freqs = rep(1,6)), "the length of freqs must equal the")
    expect_error(permuteSample(5, 4, sampleVec = "adfs"), 
                 "sampleVec must be a positive whole number")
    expect_error(permuteSample(5, 4, sampleVec = 1.1), 
                 "Each element in sampleVec must be a whole number")
    expect_error(permuteSample(5000, 40, sampleVec = 1.1),
                 "sampleVec must be a whole number")
    expect_error(permuteSample(100, 10, n = 2^50), "n must be less than or equal to 2147483647")
    expect_error(permuteSample(100, 10, n = -200), "n must be a positive number")
    expect_error(permuteSample(100, 10, sampleVec = c("62815650955529472001")),
                 "One or more of the requested values in sampleVec")
    expect_error(permuteSample(3.3, 2, n = 1), 
                 "If v is not a character and of length 1, it must be a whole number")
    expect_error(permuteSample(3, 2.2, n = 1), "m must be a whole number")
    expect_error(permuteSample(3, 2, n = 1.1), "n must be a whole number")
    expect_error(permuteSample(3, 2, sampleVec = as.complex(1)), 
                 "This type is not supported! No conversion possible for sampleVec")
    expect_error(permuteSample(3, 2, Parallel = "TRUE"), 
                 "Only logical values are supported for Parallel")
    expect_error(permuteSample(5, 3, n = 5, nThreads = "5"), 
                 "nThreads must be of type numeric or integer")
    expect_error(permuteSample(5, 3, n = 5, nThreads = 3.2), "nThreads must be a whole number")
    expect_error(permuteSample(5000, 10, n = 5.5), "n must be a whole number")
    expect_error(permuteSample(5000, 10, n = NA_integer_), "n cannot be NA or NaN")
    expect_error(permuteSample(5000, 10, n = NA), 
                 "This type is not supported! No conversion possible for n")
    expect_error(permuteSample(5000, 10, sampleVec = -1), 
                 "sampleVec must be a positive number")
    expect_error(permuteSample(5000, 10, sampleVec = c(1e9, -1)), 
                 "Each element in sampleVec must be a positive number")
    expect_error(permuteSample(NA_integer_, 3, n = 5), 
                 "If v is not a character and of length 1, it cannot be NA or NaN")
    expect_error(permuteSample(1000, 20, sampleVec = c(NA, "1234567890")),
                 "Each element in sampleVec cannot be NA or NaN")
})