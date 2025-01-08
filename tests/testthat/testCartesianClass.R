context("testing expandGridIter")

test_that("expandGridIter produces correct results", {

    CartesianClassTest <- function(lst) {

        myResults <- vector(mode = "logical")
        myRows <- expandGridCount(lst)

        a <- expandGridIter(lst)
        b <- expandGrid(lst)

        # .method("summary", &CartesianClass::summary)
        myResults <- c(myResults, isTRUE(
            all.equal(a@summary()$totalResults, myRows)
        ))

        # .method("sourceVector", &CartesianClass::sourceVector)
        myResults <- c(
            myResults,
            isTRUE(all.equal(unname(lst), unname(a@sourceVector())))
        )

        # .method("front", &CartesianClass::front)
        # .method("currIter", &CartesianClass::currIter)
        # .method("back", &CartesianClass::back)
        # .method("currIter", &CartesianClass::currIter)
        myResults <- c(myResults, isTRUE(all.equal(a@front(), b[1, ])))
        myResults <- c(myResults, isTRUE(all.equal(a@currIter(), b[1, ])))

        temp <- b[myRows, ]
        rownames(temp) <- NULL

        myResults <- c(myResults, isTRUE(all.equal(a@back(), temp)))
        myResults <- c(myResults, isTRUE(all.equal(a@currIter(), temp)))

        # .method("startOver", &CartesianClass::startOver)
        a@startOver()
        a@back()
        msg <- capture.output(noMore <- a@nextNIter(1))
        myResults <- c(myResults, is.null(noMore))
        myResults <- c(myResults, grepl("No more results. To see the last ", msg[1]))
        msg <- capture.output(noMore <- a@currIter())
        myResults <- c(myResults, grepl("No more results. To see the last ", msg[1]))

        a@startOver()
        a@back()
        msg <- capture.output(noMore <- a@nextRemaining())
        myResults <- c(myResults, is.null(noMore))
        myResults <- c(myResults, grepl("No more results. To see the last ", msg[1]))

        a@startOver()
        a1 <- b

        for (i in 1:myRows) {
            a1[i, ] <- a@nextIter()
        }

        # .method("nextIter", &CartesianClass::nextIter)
        myResults <- c(myResults, isTRUE(all.equal(a1, b)))

        # .method("startOver", &CartesianClass::startOver)
        a@startOver()
        numTest <- as.integer(myRows / 3);

        s <- 1L
        e <- numTest

        # .method("nextNIter", &CartesianClass::nextNumIters)
        for (i in 1:3) {
            temp <- b[s:e, ]
            rownames(temp) <- NULL
            iter_df <- a@nextNIter(numTest)

            myResults <- c(
                myResults, isTRUE(all.equal(iter_df, temp))
            )

            s <- e + 1L
            e <- e + numTest
        }

        # .method("nextRemaining", &CartesianClass::nextGather)
        a@startOver()
        myResults <- c(myResults, isTRUE(all.equal(a@nextRemaining(), b)))

        # .method("[[", &CartesianClass::randomAccess)
        samp <- sample(myRows, numTest)
        temp <- b[samp, ]
        rownames(temp) <- NULL
        myResults <- c(myResults, isTRUE(all.equal(a[[samp]], temp)))

        rm(a, a1,  b, temp, iter_df)
        gc()
        all(myResults)
    }

    ## INTEGERS
    myList <- list(1:5, 2:6, 3:7)
    expect_true(CartesianClassTest(myList))
    expect_true(CartesianClassTest(lapply(myList, sample)))

    ## RAW
    myList <- list(as.raw(1:5), as.raw(2:6), as.raw(3:7))
    expect_true(CartesianClassTest(myList))
    expect_true(CartesianClassTest(lapply(myList, sample)))

    ## CHARACTER
    myList <- list(letters[1:5], letters[2:6], letters[3:7])
    expect_true(CartesianClassTest(myList))
    expect_true(CartesianClassTest(lapply(myList, sample)))

    ## COMPLEX
    myList <- list(as.complex(1:5 + 1i), as.complex(2:6 + 1i),
                   as.complex(3:7 + 1i))
    expect_true(CartesianClassTest(myList))
    expect_true(CartesianClassTest(lapply(myList, sample)))

    ## DOUBLE
    myList <- list(1:5 + 0.1, 2:6 + 0.1)
    expect_true(CartesianClassTest(myList))
    expect_true(CartesianClassTest(lapply(myList, sample)))

    ## BOOLEAN
    myList <- rep(list(c(TRUE, FALSE)), 10)
    expect_true(CartesianClassTest(myList))
    expect_true(CartesianClassTest(lapply(myList, sample)))

    ## DATA.FRAME
    myList <- list(1:5, as.complex(2:6 + 1i), letters[1:4], letters[2:4])
    expect_true(CartesianClassTest(myList))
    expect_true(CartesianClassTest(lapply(myList, sample)))

    ## DATA.FRAME w/ doubles and factors
    myList <- list(1:5 + 0.1, factor(2:6 + 0.1))
    expect_true(CartesianClassTest(myList))
    expect_true(CartesianClassTest(lapply(myList, sample)))

    ## DATA.FRAME w/ integers and factors
    myList <- list(1:5, factor(2:6), as.raw(sample(10, 4)), c(TRUE, FALSE))
    expect_true(CartesianClassTest(myList))
    expect_true(CartesianClassTest(lapply(myList, sample)))

    ## FACTORS with same levels. Should return a matrix.
    myList <- list("v1" = factor(1:5, levels = 1:10),
                   "v2" = factor(2:6, levels = 1:10))
    expect_true(CartesianClassTest(myList))
    expect_true(CartesianClassTest(lapply(myList, sample)))

    ## FACTORS with different levels. Should return a data.frame.
    myList <- list("v1" = factor(1:5), "v2" = factor(2:6))
    expect_true(CartesianClassTest(myList))
    expect_true(CartesianClassTest(lapply(myList, sample)))

    ## INTEGERS w/ custom names
    myList <- list("name1" = 1:5, 2:6)
    expect_true(CartesianClassTest(myList))
    expect_true(CartesianClassTest(lapply(myList, sample)))

    ##******** BIG TESTS *********##
    CartesianClassBigZTest <- function(lst, lenCheck = 500) {

        myResults <- vector(mode = "logical")
        myRows <- expandGridCount(lst)

        a <- expandGridIter(lst)
        b1 <- expandGrid(lst, upper = lenCheck)
        b2 <- expandGrid(lst, lower = gmp::sub.bigz(myRows, lenCheck - 1))

        # .method("summary", &CartesianClass::summary)
        myResults <- c(myResults, isTRUE(
            all.equal(a@summary()$totalResults, myRows)
        ))

        # .method("sourceVector", &CartesianClass::sourceVector)
        myResults <- c(
            myResults,
            isTRUE(all.equal(unname(lst), unname(a@sourceVector())))
        )

        # .method("front", &CartesianClass::front)
        # .method("currIter", &CartesianClass::currIter)
        # .method("back", &CartesianClass::back)
        # .method("currIter", &CartesianClass::currIter)
        myResults <- c(myResults, isTRUE(all.equal(a@front(), b1[1, ])))
        myResults <- c(myResults, isTRUE(all.equal(a@currIter(), b1[1, ])))

        temp <- b2[lenCheck, ]
        rownames(temp) <- NULL

        myResults <- c(myResults, isTRUE(all.equal(a@back(), temp)))
        myResults <- c(myResults, isTRUE(all.equal(a@currIter(), temp)))

        # .method("startOver", &CartesianClass::startOver)
        a@startOver()
        a@back()
        msg <- capture.output(noMore <- a@nextNIter(1))
        myResults <- c(myResults, is.null(noMore))
        myResults <- c(myResults, grepl("No more results. To see the last ", msg[1]))
        msg <- capture.output(noMore <- a@currIter())
        myResults <- c(myResults, grepl("No more results. To see the last ", msg[1]))

        a@startOver()
        a@back()
        msg <- capture.output(noMore <- a@nextRemaining())
        myResults <- c(myResults, is.null(noMore))
        myResults <- c(myResults, grepl("No more results. To see the last ", msg[1]))

        a@startOver()
        a1 <- b1

        for (i in 1:lenCheck) {
            a1[i, ] <- a@nextIter()
        }

        # .method("nextIter", &CartesianClass::nextIter)
        myResults <- c(myResults, isTRUE(all.equal(a1, b1)))

        # .method("startOver", &CartesianClass::startOver)
        a@startOver()
        numTest <- as.integer(lenCheck / 3);

        s <- 1L
        e <- numTest

        # .method("nextNIter", &CartesianClass::nextNumIters)
        for (i in 1:3) {
            temp <- b1[s:e, ]
            rownames(temp) <- NULL
            iter_df <- a@nextNIter(numTest)

            myResults <- c(
                myResults, isTRUE(all.equal(iter_df, temp))
            )

            s <- e + 1L
            e <- e + numTest
        }

        # .method("nextRemaining", &CartesianClass::nextGather)
        a@startOver()
        a[[gmp::sub.bigz(myRows, lenCheck)]]
        myResults <- c(myResults, isTRUE(all.equal(a@nextRemaining(), b2)))

        # .method("[[", &CartesianClass::randomAccess)
        samp1 <- sample(lenCheck, 5)
        samp2 <- gmp::sub.bigz(myRows, lenCheck) + gmp::as.bigz(samp1)

        temp <- b1[samp1, ]
        rownames(temp) <- NULL
        myResults <- c(myResults, isTRUE(all.equal(a[[samp1]], temp)))

        temp <- b2[samp1, ]
        rownames(temp) <- NULL
        myResults <- c(myResults, isTRUE(all.equal(a[[samp2]], temp)))

        rm(a, a1, b1, b2, temp, iter_df)
        gc()
        all(myResults)
    }

    set.seed(8675309)

    ## INTEGERS
    myList <- Map(\(x, y) x:y, 10:30, 30:50)
    expect_true(CartesianClassBigZTest(myList))
    expect_true(CartesianClassBigZTest(lapply(myList, sample)))

    ## DOUBLES
    myList <- Map(\(x, y) x:y + rnorm(1), 10:30, 30:50)
    expect_true(CartesianClassBigZTest(myList))
    expect_true(CartesianClassBigZTest(lapply(myList, sample)))

    ## Non-GMP
    myList <- Map(\(x, y) x:y + rnorm(1), 10:15, 16:21)
    expect_true(CartesianClassBigZTest(myList))
    expect_true(CartesianClassBigZTest(lapply(myList, sample)))

    ## BOOLEANS
    myList <- rep(list(c(TRUE, FALSE)), 60)
    expect_true(CartesianClassBigZTest(myList))
    expect_true(CartesianClassBigZTest(lapply(myList, sample)))

    ## CHARACTERS
    myList <- rep(list(letters, LETTERS, state.abb), 4)
    expect_true(CartesianClassBigZTest(myList))
    expect_true(CartesianClassBigZTest(lapply(myList, sample)))

    ## DATA.FRAME
    myList <- rep(
        list(LETTERS, unique(state.region), c(TRUE, FALSE),
             as.raw(1:10), as.complex(10:1 + 1i), 1:100), 4
    )
    expect_true(CartesianClassBigZTest(myList))
    expect_true(CartesianClassBigZTest(lapply(myList, sample)))
})
