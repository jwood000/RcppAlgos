context("testing partitionsIter")

test_that("partitionsIter produces correct results", {

    partitionClassTest <- function(v_pass, m_pass = NULL, rep = FALSE,
                                   fr = NULL, tar = NULL, testRand = TRUE) {

        myResults <- vector(mode = "logical")
        myRows <- partitionsCount(v_pass, m_pass, rep, fr, tar)

        a <- partitionsIter(v_pass, m_pass, rep, fr, tar)
        b <- partitionsGeneral(v_pass, m_pass, rep, fr, tar)
        myResults <- c(myResults, isTRUE(all.equal(
            a@summary()$totalResults, myRows)
        ))

        if (length(v_pass) == 1) {
            myResults <- c(myResults, isTRUE(
                all.equal(v_pass, length(a@sourceVector()))
            ))
        } else {
            myResults <- c(myResults, isTRUE(
                all.equal(sort(v_pass), a@sourceVector())
            ))
        }

        if (testRand) {
            myResults <- c(myResults, isTRUE(
                all.equal(a@front(), b[1 ,])
            ))
            myResults <- c(myResults, isTRUE(all.equal(a@currIter(),
                                                       b[1 ,])))
            myResults <- c(myResults, isTRUE(all.equal(a@back(),
                                                       b[myRows, ])))
            myResults <- c(myResults, isTRUE(all.equal(a@currIter(),
                                                       b[myRows, ])))
        }

        a@startOver()
        msg <- capture.output(noMore <- a@currIter())
        myResults <- c(myResults, is.null(noMore))
        myResults <- c(myResults, grepl("Iterator Initialized. To see the first", msg[1]))
        a1 <- b

        for (i in 1:myRows) {
            a1[i, ] <- a@nextIter()
        }

        myResults <- c(myResults, isTRUE(all.equal(a1, b)))
        a@startOver()
        numTest <- as.integer(myRows / 3);

        s <- 1L
        e <- numTest

        for (i in 1:3) {
            myResults <- c(myResults, isTRUE(all.equal(a@nextNIter(numTest),
                                                       b[s:e, ])))
            s <- e + 1L
            e <- e + numTest
        }

        a@startOver()
        myResults <- c(myResults, isTRUE(all.equal(a@nextRemaining(), b)))
        msg <- capture.output(noMore <- a@nextIter())
        myResults <- c(myResults, is.null(noMore))

        if (testRand) {
            a@back()
            msg <- capture.output(noMore <- a@nextNIter(1))
            myResults <- c(myResults, is.null(noMore))
            myResults <- c(myResults, "No more results." == msg[1])
            msg <- capture.output(noMore <- a@currIter())
            myResults <- c(myResults, "No more results." == msg[1])

            a@startOver()
            a@back()
            msg <- capture.output(noMore <- a@nextRemaining())
            myResults <- c(myResults, is.null(noMore))
            myResults <- c(myResults, "No more results." == msg[1])

            a@startOver()
            a@back()
            msg <- capture.output(noMore <- a@nextIter())
            myResults <- c(myResults, is.null(noMore))
            myResults <- c(myResults, "No more results." == msg[1])

            samp <- sample(myRows, numTest)
            myResults <- c(myResults, isTRUE(all.equal(a[[samp]], b[samp, ])))
            one_samp <- sample(myRows, 1)
            myResults <- c(myResults, isTRUE(all.equal(a[[one_samp]], b[one_samp, ])))
        }

        rm(a, a1, b)
        gc()
        all(myResults)
    }

    #### Distinct; Length determined internally; No zero;
    expect_true(partitionClassTest(189))

    #### Distinct; Length determined internally; One zero;
    expect_true(partitionClassTest(0:50))

    #### Distinct; Specific Length; No zero
    expect_true(partitionClassTest(50, 5))
    #### Mapped version
    ## 50 * 3 + 6 * 5 = 180
    expect_true(partitionClassTest(6 + (1:50) * 3, 5, tar = 180))

    #### Distinct; Specific Length; One zero
    expect_true(partitionClassTest(0:50, 5))
    #### Mapped version
    expect_true(partitionClassTest(6 + (0:50) * 3, 5, tar = 180))

    #### Distinct; Specific Length; Multiple Zeros; Not enough to maximize
    expect_true(partitionClassTest(0:50, 9, fr = c(4, rep(1, 50))))
    #### Mapped version
    ## 50 * 13 + 7 * 9 = 713
    expect_true(partitionClassTest(7 + (0:50) * 13, 9,
                                   fr = c(4, rep(1, 50)), tar = 713))

    #### Distinct; Specific Length; Multiple Zeros; Enough to maximize;
    #### Length is restrictive
    expect_true(partitionClassTest(0:50, 5, fr = c(8, rep(1, 50))))
    #### Mapped version
    ## 50 * 13 + 7 * 5 = 713
    expect_true(partitionClassTest(7 + (0:50) * 13, 5,
                                   fr = c(4, rep(1, 50)), tar = 685))

    #### Distinct; Length determined internally; Multiple Zeros;
    #### Enough to maximize; N.B. There is no mapped version of this case
    expect_true(partitionClassTest(0:50, fr = c(50, rep(1, 50))))

    #### Distinct; Specific Length; No Zeros; Specific Target;
    expect_true(partitionClassTest(30, 8, tar = 75))

    #### Distinct; Specific Length; Multi Zeros; Specific Target;
    expect_true(partitionClassTest(0:30, 6, tar = 75, fr = c(3, rep(1, 30))))

    #### Repetition; Length determined internally; Multiple Zero;
    expect_true(partitionClassTest(0:30, rep = TRUE))
    #### Mapped version
    ## 19 * 30 + 30 * 3e9 = 90000000570
    expect_true(partitionClassTest(19 + (0:30) * 3e9, 30,
                                   rep = TRUE, tar = 90000000570))

    #### Repetition; Specific Length; No zero
    expect_true(partitionClassTest(50, 5, TRUE))
    #### Mapped version
    ## 19 * 5 + 50 * 3 = 245
    expect_true(partitionClassTest(19 + (1:50) * 3, 5, TRUE, tar = 245))

    #### Repetition; Specific Length; Zero included
    expect_true(partitionClassTest(0:30, 10, rep = TRUE))
    #### Mapped version
    ## 19 * 10 + 30 * 3 = 280
    expect_true(partitionClassTest(19 + (0:30) * 3, 10,
                                   rep = TRUE, tar = 280))

    #### Repetition; Specific Length; No Zeros; Specific Target;
    expect_true(partitionClassTest(20, 10, rep = TRUE, tar = 45))

    #### Multiset; Specific Length;
    expect_true(partitionClassTest(50, 6, fr = rep(4, 50),
                                   testRand = FALSE))

    #### Multiset; Mapped;
    # $num_partitions
    # [1] 15080
    #
    # $mapped_vector
    # [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
    # [27] 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52
    # [53] 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70
    #
    # $mapped_target
    # [1] 853
    #
    # $first_index_vector
    # [1] 41 66 66 67 67 67 67 68 68 68 69 69 70
    #
    # $eqn_check
    # [1] TRUE
    #
    # $partition_type
    # [1] "Multiset"
    expect_true(partitionClassTest(79L + -2L * (1:70), 13, fr = rep(1:10, 7),
                                   tar = 887L, testRand = FALSE))

    ## N.B. In the above we see the mapped target is 853. We must remember to
    ## also reverse freqs as the mapped vector is 1:70
    # RcppAlgos:::partitionsDesign(70, 13, freqs = rep(10:1, 7), target = 853L)
    # $num_partitions
    # [1] 15080
    #
    # $mapped_vector
    # [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
    # [28] 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54
    # [55] 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70
    #
    # $mapped_target
    # [1] 853
    #
    # $first_index_vector
    # [1] 41 66 66 67 67 67 67 68 68 68 69 69 70
    #
    # $eqn_check
    # [1] TRUE
    #
    # $partition_type
    # [1] "Multiset"
    expect_true(partitionClassTest(70, 13, fr = rep(10:1, 7),
                                   tar = 853L, testRand = FALSE))

    #### Multiset; Mapped; Double Precision
    expect_true(partitionClassTest((1:50) * 1e10, 13, fr = rep(1:10, 5),
                                   tar = 6.17e12, testRand = FALSE))

    #### Multiset; zero included; random freqs; non-standard target
    set.seed(123)
    expect_true(partitionClassTest(0:50, 6, fr = sample(1:8, 51, TRUE),
                                   tar = 60, testRand = FALSE))

    ##******** BIG TESTS *********##
    partitionClassBigZTest <- function(v_pass, m_pass = NULL, rep = FALSE,
                                       fr = NULL, tar = NULL, lenCheck = 1000) {

        myResults <- vector(mode = "logical")
        myRows <- partitionsCount(v_pass, m_pass, rep, fr, tar)

        a <- partitionsIter(v_pass, m_pass, rep, fr, tar)
        b1 <- partitionsGeneral(v_pass, m_pass, rep, fr, tar, upper = lenCheck)
        b2 <- partitionsGeneral(v_pass, m_pass, rep, fr, tar,
                                lower = gmp::sub.bigz(myRows, lenCheck - 1))

        myResults <- c(myResults, isTRUE(all.equal(
            a@summary()$totalResults, myRows)
        ))

        if (length(v_pass) == 1) {
            myResults <- c(myResults, isTRUE(
                all.equal(v_pass, length(a@sourceVector()))
            ))
        } else {
            myResults <- c(myResults, isTRUE(
                all.equal(sort(v_pass), a@sourceVector())
            ))
        }

        myResults <- c(myResults, isTRUE(
            all.equal(a@front(), b1[1 ,])
        ))
        myResults <- c(myResults, isTRUE(all.equal(a@currIter(),
                                                   b1[1 ,])))
        myResults <- c(myResults, isTRUE(all.equal(a@back(),
                                                   b2[lenCheck, ])))
        myResults <- c(myResults, isTRUE(all.equal(a@currIter(),
                                                   b2[lenCheck, ])))

        a@startOver()
        a1 <- b1

        for (i in 1:lenCheck) {
            a1[i, ] <- a@nextIter()
        }

        myResults <- c(myResults, isTRUE(all.equal(a1, b1)))
        a@startOver()
        numTest <- as.integer(lenCheck / 3);
        s <- 1L
        e <- numTest

        for (i in 1:3) {
            myResults <- c(myResults, isTRUE(all.equal(a@nextNIter(numTest),
                                                       b1[s:e, ])))
            s <- e + 1L
            e <- e + numTest
        }

        a@startOver()
        a[[gmp::sub.bigz(myRows, lenCheck)]]
        myResults <- c(myResults, isTRUE(all.equal(a@nextRemaining(), b2)))

        t <- capture.output(a@nextIter())
        myResults <- c(myResults, is.null(a@nextIter()))
        myResults <- c(myResults, is.null(a@nextNIter(1)))
        myResults <- c(myResults, is.null(a@nextRemaining()))

        samp1 <- sample(lenCheck, 5)
        samp2 <- gmp::sub.bigz(myRows, lenCheck) + gmp::as.bigz(samp1)
        myResults <- c(myResults, isTRUE(all.equal(a[[samp1]], b1[samp1, ])))
        myResults <- c(myResults, isTRUE(all.equal(a[[samp2]], b2[samp1, ])))
        rm(a, a1, b1, b2)
        gc()
        all(myResults)
    }

    expect_true(partitionClassBigZTest(2000, 10, TRUE))
    #### Mapped version
    ## 17 * 10 + 2000 * 123456789 = 280
    expect_true(partitionClassBigZTest(17 + (1:2000) * 123456789,
                                       10, TRUE, tar = 246913578170))
    expect_true(partitionClassBigZTest(2000, 10))
})
