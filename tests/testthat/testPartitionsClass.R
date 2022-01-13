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

        if (testRand) {
            samp <- sample(myRows, numTest)
            myResults <- c(myResults, isTRUE(all.equal(a[[samp]], b[samp, ])))
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
    ## 19 * 30 + 30 * 3 = 660
    expect_true(partitionClassTest(19 + (0:30) * 3, 30,
                                   rep = TRUE, tar = 660))
    
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
    expect_true(partitionClassTest(79L + -2L * (1:70), 13, fr = rep(1:10, 7),
                                   tar = 887L, testRand = FALSE))
    
    #### Multiset; zero included; random freqs; non-standard target
    set.seed(123)
    expect_true(partitionClassTest(0:50, 6, fr = sample(1:8, 51, TRUE),
                                   tar = 60, testRand = FALSE))
})