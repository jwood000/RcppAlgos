context("testing comboRank and permuteRank")

test_that("comboRank produces correct results", {

    expect_equal(comboRank(), integer(0))
    expect_equal(comboRank(comboGeneral(5, 3), v = 5), 1:10)
    expect_equal(comboRank(1, v = 1), 1)

    totalCombs = comboCount(5, 8, freqs = rep(3, 5))
    expect_equal(comboRank(comboGeneral(5, 8, freqs = rep(3, 5)), v = 5,
                           freqs = rep(3, 5)), 1:totalCombs)

    totalCombs = comboCount(5, 8, TRUE)
    expect_equal(comboRank(comboGeneral(5, 8, TRUE),
                           v = 5, repetition = TRUE), 1:totalCombs)

    expect_equal(comboRank(comboGeneral(letters[1:5], 3, TRUE),
                           v = letters[1:5], repetition = TRUE), 1:35)

    myNums = rnorm(5)
    expect_equal(comboRank(comboGeneral(myNums, 3), v = myNums), 1:10)

    totalCombs = comboCount(LETTERS, 4)
    samp = sample(totalCombs, 10)
    expect_equal(comboRank(comboSample(LETTERS, 4, sampleVec = samp),
                           v = LETTERS), samp)

    mySamp = seq.int(1, comboCount(5, 5, freqs = c(2,1,2,1,2)), 3)
    expect_equal(comboRank(comboSample(5, 5, freqs = c(2,1,2,1,2),
                                       sampleVec = mySamp),
                           v = 5, freqs = c(2,1,2,1,2)), mySamp)

    mySamp = comboSample(500, 200, n = 10, seed = 42, namedSample = TRUE)
    expect_equal(as.character(comboRank(mySamp, v = 500)),
                 rownames(mySamp))

    mySamp = comboSample(500, 200, TRUE, n = 10, seed = 42, namedSample = TRUE)
    expect_equal(as.character(comboRank(mySamp, v = 500, repetition = TRUE)),
                 rownames(mySamp))

    ## **** MULTISETS ******
    lastRow = comboCount(100, 20, freqs = rep(1:4, 25))
    mySamp = comboSample(as.character(1:100), 20, freqs = rep(1:4, 25),
                         sampleVec = c(gmp::as.bigz(1), lastRow))
    expect_equal(comboRank(mySamp, v = 100, freqs = rep(1:4, 25)),
                 c(gmp::as.bigz(1), lastRow))

    mySamp = comboSample(100, 20, freqs = rep(1:4, 25),
                         sampleVec = seq(1L, 100L, 10L))
    expect_equal(
        as.integer(comboRank(mySamp, v = 100, freqs = rep(1:4, 25))),
        seq(1L, 100L, 10L)
    )

    set.seed(123)
    v <- rnorm(100)
    mySamp = comboSample(v, 20, freqs = rep(1:4, 25),
                         n = 10, seed = 42, namedSample = TRUE)
    expect_equal(as.character(comboRank(mySamp, v = v, freqs = rep(1:4, 25))),
                 rownames(mySamp))

    mySamp = comboSample(as.factor(LETTERS), 12, TRUE, n = 2,
                         seed = 17, namedSample = TRUE)
    expect_equal(comboRank(mySamp, v = as.factor(LETTERS), repetition = TRUE),
                 as.integer(rownames(mySamp)))

    ## Multiple inputs
    mat1 = comboSample(100, 5, n = 5, seed = 123, namedSample = TRUE)
    mat2 = comboSample(100, 7, n = 5, seed = 987, namedSample = TRUE)
    mat3 = comboSample(100, 15, n = 5, seed = 97, namedSample = TRUE)
    myRank = comboRank(int = mat1, dbl = mat2, bigz = mat3, v = 100)
    expect_named(myRank, c("int", "dbl", "bigz"))
    expect_equal(unname(sapply(myRank, class)), c("integer", "numeric", "bigz"))
    expect_true(all(mapply(function(x, y) {
        identical(as.character(x), rownames(y))
    }, myRank, list(mat1, mat2, mat3))))
})

test_that("permuteRank produces correct results", {
    expect_equal(permuteRank(permuteGeneral(5, 3), v = 5),
                 1:(permuteCount(5, 3)))
    set.seed(2)
    mySamp = permuteSample(5, 3, n = 5, namedSample = TRUE)
    expect_equal(permuteRank(mySamp, v = 5), as.integer(rownames(mySamp)))
    expect_equal(permuteRank(1, v = 1), 1)
    mySamp = permuteSample(5, 20, freqs = 1:5, n = 3, namedSample = TRUE)
    expect_equal(permuteRank(mySamp, v = 5, freqs = 1:5),
                 as.integer(rownames(mySamp)))

    expect_equal(permuteRank(permuteGeneral(5, 6, freqs = rep(2, 5)),
                             v = 5, freqs = rep(2, 5)),
                 1:(permuteCount(5, 6, freqs = rep(2, 5))))

    expect_equal(permuteRank(permuteGeneral(5, 6, TRUE),
                             v = 5, repetition = TRUE),
                 1:(permuteCount(5, 6, TRUE)))

    expect_equal(permuteRank(permuteGeneral(5, 3, TRUE),
                             v = 5, repetition = TRUE),
                 1:(permuteCount(5, 3, TRUE)))

    mySamp1 = permuteSample(as.raw(1:4), 8, T,
                            sampleVec = c(31788, 59688, 3780))
    mySamp2 = permuteSample(as.complex(c(1, -1, 1i, -1i)),
                            8, T, sampleVec = c(31788, 59688, 3780))
    expect_equal(permuteRank(mySamp1, v = as.raw(1:4), repetition = TRUE),
                 permuteRank(mySamp2, v = as.complex(c(1, -1, 1i, -1i)),
                             repetition = TRUE))

    totalPerms = permuteCount(LETTERS, 3)
    samp = as.numeric(sample(totalPerms, 10))
    mySamp = permuteSample(LETTERS, 3, sampleVec = samp)
    expect_equal(permuteRank(mySamp, v = LETTERS), samp)

    expect_equal(permuteRank(permuteGeneral(5, freqs = c(2,1,2,1,2)),
                             v = 5, freqs = c(2,1,2,1,2)),
                 1:(permuteCount(5, freqs = c(2,1,2,1,2))))

    mySamp1 = permuteSample(c(T,F), freqs = c(50, 50), n = 3, seed = 22)
    mySamp2 = permuteSample(1:0, freqs = c(50, 50), n = 3, seed = 22)
    expect_equal(permuteRank(mySamp1, v = c(T, F), freqs = c(50, 50)),
                 permuteRank(mySamp2, v = 1:0, freqs = c(50, 50)))

    mySamp = permuteSample(factor(state.name), 20, sampleVec = 1e12)
    expect_equal(permuteRank(mySamp, v = factor(state.name)), 1e12)

    mySamp = permuteSample(500, 100, n = 5, seed = 42, namedSample = TRUE)
    expect_equal(as.character(permuteRank(mySamp, v = 500)),
                 rownames(mySamp))

    mySamp = permuteSample(100, 20, TRUE, n = 5, seed = 42, namedSample = TRUE)
    expect_equal(as.character(permuteRank(mySamp, v = 100, repetition = TRUE)),
                 rownames(mySamp))

    mySamp = permuteSample(75, 10, freqs = rep(1:3, 25),
                           n = 5, seed = 42, namedSample = TRUE)
    expect_equal(as.character(permuteRank(mySamp, v = 75,
                                          freqs = rep(1:3, 25))),
                 rownames(mySamp))

    mySamp = permuteSample(c(TRUE, FALSE), 20, freqs = c(10, 15),
                           seed = 97, n = 5, namedSample = TRUE)
    expect_equal(as.character(permuteRank(mySamp, v = c(TRUE, FALSE),
                                          freqs = c(10, 15))),
                 rownames(mySamp))
})
