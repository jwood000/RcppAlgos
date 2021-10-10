context("testing partitionsIter")

test_that("partitionsIter produces correct results", {
    
    partitionClassTest <- function(v1, m1 = NULL, rep1 = FALSE,
                                   freqs1 = NULL, limit1 = NULL) {
        
        myResults <- vector(mode = "logical")
        myRows <- partitionsCount(v1, m1, rep1, freqs1, limit1)

        a <- partitionsIter(v1, m1, rep1, freqs1, limit1)
        b <- partitionsGeneral(v1, m1, rep1, freqs1, limit1)
        myResults <- c(myResults, isTRUE(all.equal(a@summary()$totalResults, myRows)))

        if (length(v1) == 1) {
            myResults <- c(myResults, isTRUE(all.equal(v1, length(a@sourceVector()))))
        } else {
            myResults <- c(myResults, isTRUE(all.equal(v1, a@sourceVector())))
        }

        myResults <- c(myResults, isTRUE(all.equal(a@front(), b[1 ,])))
        myResults <- c(myResults, isTRUE(all.equal(a@currIter(), b[1 ,])))
        myResults <- c(myResults, isTRUE(all.equal(a@back(), b[myRows, ])))
        myResults <- c(myResults, isTRUE(all.equal(a@currIter(), b[myRows, ])))

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
            myResults <- c(myResults, isTRUE(all.equal(a@nextNIter(numTest), b[s:e, ])))
            s <- e + 1L
            e <- e + numTest
        }

        a@startOver()
        myResults <- c(myResults, isTRUE(all.equal(a@nextRemaining(), b)))

        samp <- sample(myRows, numTest)
        myResults <- c(myResults, isTRUE(all.equal(a[[samp]], b[samp, ])))
        rm(a, a1, b)
        gc()
        all(myResults)
    }
    
    expect_true(partitionClassTest(50, 5))
    expect_true(partitionClassTest(50, 5, TRUE))
    
    expect_true(partitionClassTest(0:50))
    expect_true(partitionClassTest(0:50, freqs1 = c(50, rep(1, 50))))
    expect_true(partitionClassTest(0:30, rep1 = TRUE))
    expect_true(partitionClassTest(0:30, 10, rep1 = TRUE))
})
