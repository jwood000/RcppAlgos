test_that(paste("partitionsGeneral and partitionsIter produces empty",
                "matrix when there are no partitions"), {

    ### *************************** Partitions **************************** ###
    ## Distinct case
    expect_identical(partitionsGeneral(10, 5),
                     matrix(integer(0), nrow = 0, ncol = 5))
    iter <- partitionsIter(10, 5)
    msg <- capture.output(noMore <- iter@nextIter())
    expect_null(noMore)

    ## Mapped version
    expect_identical(partitionsGeneral((1:10) * 3L + 13L, 5, target = 95L),
                     matrix(integer(0), nrow = 0, ncol = 5))
    iter <- partitionsIter((1:10) * 3L + 13L, 5, target = 95L)
    msg <- capture.output(noMore <- iter@nextIter())
    expect_null(noMore)

    ## Distinct case with comboGeneral
    expect_identical(comboGeneral(10, 5, constraintFun = "prod",
                                  comparisonFun = "==",
                                  limitConstraints = 7),
                     matrix(integer(0), nrow = 0, ncol = 5))
    iter <- comboIter(10, 5, constraintFun = "prod",
                      comparisonFun = "==", limitConstraints = 7)
    msg <- capture.output(noMore <- iter@nextIter())
    expect_null(noMore)

    ## Mapped version
    expect_identical(comboGeneral((1:10) * 3L + 13L, 5, constraintFun = "prod",
                                  comparisonFun = "==",
                                  limitConstraints = 7),
                     matrix(integer(0), nrow = 0, ncol = 5))
    iter <- comboIter((1:10) * 3L + 13L, 5, constraintFun = "prod",
                      comparisonFun = "==", limitConstraints = 7)
    msg <- capture.output(noMore <- iter@nextIter())
    expect_null(noMore)

    ## Repetition case
    expect_identical(partitionsGeneral(10, 11, TRUE),
                     matrix(integer(0), nrow = 0, ncol = 11))
    iter <- partitionsIter(10, 11, TRUE)
    msg <- capture.output(noMore <- iter@nextIter())
    expect_null(noMore)

    ## Mapped version
    expect_identical(partitionsGeneral((1:10) * 3L + 13L, 11, TRUE),
                     matrix(integer(0), nrow = 0, ncol = 11))
    iter <- partitionsIter((1:10) * 3L + 13L, 11, TRUE)
    msg <- capture.output(noMore <- iter@nextIter())
    expect_null(noMore)

    ## Repetition case with comboGeneral
    expect_identical(comboGeneral(10, 5, TRUE,
                                  constraintFun = "prod",
                                  comparisonFun = "==",
                                  limitConstraints = 17),
                     matrix(integer(0), nrow = 0, ncol = 5))
    iter <- comboIter(10, 5, TRUE, constraintFun = "prod",
                      comparisonFun = "==", limitConstraints = 17)
    msg <- capture.output(noMore <- iter@nextIter())
    expect_null(noMore)

    ## Mapped version
    expect_identical(comboGeneral((1:10) * 3L + 13L, 5, TRUE,
                                  constraintFun = "prod",
                                  comparisonFun = "==",
                                  limitConstraints = 17),
                     matrix(integer(0), nrow = 0, ncol = 5))
    iter <- comboIter((1:10) * 3L + 13L, 5, TRUE, constraintFun = "prod",
                      comparisonFun = "==", limitConstraints = 17)
    msg <- capture.output(noMore <- iter@nextIter())
    expect_null(noMore)

    ## Multiset case
    expect_identical(partitionsGeneral(10, 5, freqs = rep(1:2, 5)),
                     matrix(integer(0), nrow = 0, ncol = 5))
    iter <- partitionsIter(10, 5, freqs = rep(1:2, 5))
    msg <- capture.output(noMore <- iter@nextIter())
    expect_null(noMore)

    ## Mapped version
    expect_identical(partitionsGeneral((1:10) * 3L + 13L, 5,
                                       freqs = rep(1:2, 5)),
                     matrix(integer(0), nrow = 0, ncol = 5))
    iter <- partitionsIter((1:10) * 3L + 13L, 5, freqs = rep(1:2, 5))
    msg <- capture.output(noMore <- iter@nextIter())
    expect_null(noMore)

    ## Multiset case with comboGeneral
    expect_identical(comboGeneral(10, 5, freqs = rep(1:2, 5),
                                  constraintFun = "prod",
                                  comparisonFun = "==",
                                  limitConstraints = 7),
                     matrix(integer(0), nrow = 0, ncol = 5))
    iter <- comboIter(10, 5, freqs = rep(1:2, 5), constraintFun = "prod",
                      comparisonFun = "==", limitConstraints = 7)
    msg <- capture.output(noMore <- iter@nextIter())
    expect_null(noMore)

    ## Mapped version
    expect_identical(comboGeneral((1:10) * 3L + 13L, 5, freqs = rep(1:2, 5),
                                  constraintFun = "prod",
                                  comparisonFun = "==",
                                  limitConstraints = 7),
                     matrix(integer(0), nrow = 0, ncol = 5))
    iter <- comboIter((1:10) * 3L + 13L, 5, freqs = rep(1:2, 5),
                      constraintFun = "prod",
                      comparisonFun = "==", limitConstraints = 7)
    msg <- capture.output(noMore <- iter@nextIter())
    expect_null(noMore)

    ### ************************** Compositions *************************** ###
    ## Distinct case
    expect_identical(compositionsGeneral(10, 5),
                     matrix(integer(0), nrow = 0, ncol = 5))
    iter <- compositionsIter(10, 5)
    msg <- capture.output(noMore <- iter@nextIter())
    expect_null(noMore)

    ## Mapped version
    expect_identical(compositionsGeneral((1:10) * 3L + 13L, 5, target = 95L),
                     matrix(integer(0), nrow = 0, ncol = 5))
    iter <- compositionsIter((1:10) * 3L + 13L, 5, target = 95L)
    msg <- capture.output(noMore <- iter@nextIter())
    expect_null(noMore)

    ## Distinct case with permuteGeneral
    expect_identical(permuteGeneral(10, 5, constraintFun = "prod",
                                  comparisonFun = "==",
                                  limitConstraints = 7),
                     matrix(integer(0), nrow = 0, ncol = 5))
    iter <- permuteIter(10, 5, constraintFun = "prod",
                      comparisonFun = "==", limitConstraints = 7)
    msg <- capture.output(noMore <- iter@nextIter())
    expect_null(noMore)

    ## Mapped version
    expect_identical(permuteGeneral((1:10) * 3L + 13L, 5, constraintFun = "prod",
                                  comparisonFun = "==",
                                  limitConstraints = 7),
                     matrix(integer(0), nrow = 0, ncol = 5))
    iter <- permuteIter((1:10) * 3L + 13L, 5, constraintFun = "prod",
                      comparisonFun = "==", limitConstraints = 7)
    msg <- capture.output(noMore <- iter@nextIter())
    expect_null(noMore)

    ## Repetition case
    expect_identical(compositionsGeneral(10, 11, TRUE),
                     matrix(integer(0), nrow = 0, ncol = 11))
    iter <- compositionsIter(10, 11, TRUE)
    msg <- capture.output(noMore <- iter@nextIter())
    expect_null(noMore)

    ## Mapped version
    expect_identical(compositionsGeneral((1:10) * 3L + 13L, 11, TRUE),
                     matrix(integer(0), nrow = 0, ncol = 11))
    iter <- compositionsIter((1:10) * 3L + 13L, 11, TRUE)
    msg <- capture.output(noMore <- iter@nextIter())
    expect_null(noMore)

    ## Repetition case with permuteGeneral
    expect_identical(permuteGeneral(10, 5, TRUE,
                                  constraintFun = "prod",
                                  comparisonFun = "==",
                                  limitConstraints = 17),
                     matrix(integer(0), nrow = 0, ncol = 5))
    iter <- permuteIter(10, 5, TRUE, constraintFun = "prod",
                      comparisonFun = "==", limitConstraints = 17)
    msg <- capture.output(noMore <- iter@nextIter())
    expect_null(noMore)

    ## Mapped version
    expect_identical(permuteGeneral((1:10) * 3L + 13L, 5, TRUE,
                                  constraintFun = "prod",
                                  comparisonFun = "==",
                                  limitConstraints = 17),
                     matrix(integer(0), nrow = 0, ncol = 5))
    iter <- permuteIter((1:10) * 3L + 13L, 5, TRUE, constraintFun = "prod",
                      comparisonFun = "==", limitConstraints = 17)
    msg <- capture.output(noMore <- iter@nextIter())
    expect_null(noMore)

    ## Multiset case with permuteGeneral
    expect_identical(permuteGeneral(10, 5, freqs = rep(1:2, 5),
                                  constraintFun = "prod",
                                  comparisonFun = "==",
                                  limitConstraints = 7),
                     matrix(integer(0), nrow = 0, ncol = 5))
    iter <- permuteIter(10, 5, freqs = rep(1:2, 5), constraintFun = "prod",
                      comparisonFun = "==", limitConstraints = 7)
    msg <- capture.output(noMore <- iter@nextIter())
    expect_null(noMore)

    ## Mapped version
    expect_identical(permuteGeneral((1:10) * 3L + 13L, 5, freqs = rep(1:2, 5),
                                  constraintFun = "prod",
                                  comparisonFun = "==",
                                  limitConstraints = 7),
                     matrix(integer(0), nrow = 0, ncol = 5))
    iter <- permuteIter((1:10) * 3L + 13L, 5, freqs = rep(1:2, 5),
                      constraintFun = "prod",
                      comparisonFun = "==", limitConstraints = 7)
    msg <- capture.output(noMore <- iter@nextIter())
    expect_null(noMore)

    ## ********** Multiset case doesn't exist yet for compositions **********
    skip(paste("Multiset compositions not implemented yet;",
               "kept as a spec for future work."))
    expect_identical(compositionsGeneral(10, 5, freqs = rep(1:2, 5)),
                     matrix(integer(0), nrow = 0, ncol = 5))
    iter <- compositionsIter(10, 5, freqs = rep(1:2, 5))
    msg <- capture.output(noMore <- iter@nextIter())
    expect_null(noMore)

    expect_identical(compositionsGeneral((1:10) * 3L + 13L, 5, freqs = rep(1:2, 5)),
                     matrix(integer(0), nrow = 0, ncol = 5))
    iter <- compositionsIter((1:10) * 3L + 13L, 5, freqs = rep(1:2, 5))
    msg <- capture.output(noMore <- iter@nextIter())
    expect_null(noMore)
})

test_that("partitionsIter produces correct results", {

    partitionClassTest <- function(
        v_pass, m_pass = NULL, rep = FALSE, fr = NULL, tar = NULL,
        testRand = TRUE, IsComposition = FALSE, IsWeak = FALSE,
        sanity = TRUE, requiresWidthRebuild = FALSE
    ) {

        myResults <- vector(mode = "logical")

        if (IsComposition) {
            myRows <- compositionsCount(v_pass, m_pass, rep, fr, tar, IsWeak)
            a <- compositionsIter(v_pass, m_pass, rep, fr, tar, IsWeak)
            b <- compositionsGeneral(v_pass, m_pass, rep, fr, tar, IsWeak)

            if (sanity) {
                perm_tar <- sum(b[1, ])

                b1 <- if (!requiresWidthRebuild) {
                    permuteGeneral(
                        v = v_pass, m = ncol(b), repetition = rep, freqs = fr,
                        constraintFun = "sum", comparisonFun = "==",
                        limitConstraints = perm_tar
                    )
                } else {
                    min_m <- max(1, ncol(b) - sum(0 == b[1, ]))
                    max_m <- ncol(b)
                    new_v <- if (is.null(fr)) {
                        v_pass[-which(v_pass == 0)]
                    } else {
                        v_pass[-which(fr == max(fr))]
                    }

                    do.call(
                        rbind,
                        lapply(min_m:max_m, function(width) {
                            perms <- permuteGeneral(
                                v = new_v, m = width, repetition = rep,
                                constraintFun = "sum", comparisonFun = "==",
                                limitConstraints = perm_tar
                            )

                            if (width < max_m) {
                                if (class(perms[1, ]) == "integer") {
                                    myMat <- cbind(0L, perms)
                                    count <- 1L

                                    while (ncol(myMat) < max_m) {
                                        myMat <- cbind(0L, myMat)

                                    }

                                    return(myMat)
                                } else {
                                    myMat <- cbind(0, perms)

                                    while (ncol(myMat) < max_m) {
                                        myMat <- cbind(0, myMat)
                                    }

                                    return(myMat)
                                }
                            } else {
                                return(perms)
                            }
                        })
                    )
                }

                b1 <- b1[do.call(order, as.data.frame(b1)), , drop = FALSE]
                myResults <- c(myResults, identical(b1, b))
            }
        } else {
            if (class(v_pass) != "table") {
                a <- partitionsIter(v_pass, m_pass, rep, fr, tar)
                b <- partitionsGeneral(v_pass, m_pass, rep, fr, tar)
                myRows <- partitionsCount(v_pass, m_pass, rep, fr, tar)
            } else {
                a <- partitionsIter(v_pass, m_pass, tar)
                b <- partitionsGeneral(v_pass, m_pass, tar)
                myRows <- partitionsCount(v_pass, m_pass, tar)
            }
        }

        myResults <- c(myResults, isTRUE(all.equal(
            a@summary()$totalResults, myRows)
        ))

        if (length(v_pass) == 1 && v_pass == 0) {
            myResults <- c(myResults, v_pass == a@sourceVector())
        } else if (length(v_pass) == 1) {
            myResults <- c(myResults, isTRUE(
                all.equal(abs(v_pass), length(a@sourceVector()))
            ))
        } else if (class(v_pass) != "table") {
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

        if (myRows) {
            for (i in 1:myRows) {
                a1[i, ] <- a@nextIter()
            }

            myResults <- c(myResults, isTRUE(all.equal(a1, b)))
            a@startOver()
            num_iters <- if (myRows > 10) 3L else 1L
            numTest   <- as.integer(myRows / num_iters);

            s <- 1L
            e <- numTest

            for (i in 1:num_iters) {
                myResults <- c(myResults, isTRUE(all.equal(a@nextNIter(numTest),
                                                           b[s:e, , drop = FALSE])))
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
        } else {
            a@startOver()
            msg <- capture.output(noMore <- a@nextNIter(1))
            myResults <- c(myResults, is.null(noMore))
            myResults <- c(myResults, "No more results." == msg[1])
            msg <- capture.output(noMore <- a@currIter())
            myResults <- c(myResults, "No more results." == msg[1])

            a@startOver()
            msg <- capture.output(noMore <- a@nextIter())
            myResults <- c(myResults, is.null(noMore))
            myResults <- c(myResults, "No more results." == msg[1])

            a@startOver()
            msg <- capture.output(noMore <- a@nextRemaining())
            myResults <- c(myResults, is.null(noMore))
            myResults <- c(myResults, "No more results." == msg[1])
        }

        rm(a, a1, b)
        gc()
        all(myResults)
    }

    #### Trivial Cases
    expect_true(partitionClassTest(0, testRand = FALSE))
    expect_true(partitionClassTest(1, testRand = FALSE))
    expect_true(partitionClassTest(2, testRand = FALSE))
    expect_true(partitionClassTest(0, IsComposition = TRUE, testRand = FALSE))
    expect_true(partitionClassTest(1, IsComposition = TRUE, testRand = FALSE))
    expect_true(partitionClassTest(2, IsComposition = TRUE, testRand = FALSE))
    expect_true(partitionClassTest(0, rep = TRUE, testRand = FALSE))
    expect_true(partitionClassTest(1, rep = TRUE, testRand = FALSE))
    expect_true(partitionClassTest(2, rep = TRUE, testRand = FALSE))
    expect_true(partitionClassTest(2, rep = TRUE, testRand = FALSE,
                                   IsComposition = TRUE))
    expect_true(partitionClassTest(0:1, rep = TRUE, testRand = FALSE))
    expect_true(partitionClassTest(0:1, rep = TRUE, testRand = FALSE,
                                   IsComposition = TRUE))
    expect_true(partitionClassTest(0:1, rep = TRUE, testRand = FALSE,
                                   IsComposition = TRUE, IsWeak = TRUE))
    expect_true(partitionClassTest(0:2, rep = TRUE, testRand = FALSE))
    expect_true(partitionClassTest(
        0:2, rep = TRUE, testRand = FALSE, IsComposition = TRUE,
        requiresWidthRebuild = TRUE
    ))
    expect_true(partitionClassTest(0:2, rep = TRUE, testRand = FALSE,
                                   IsComposition = TRUE, IsWeak = TRUE))
    expect_true(partitionClassTest(-1, testRand = FALSE))
    expect_true(partitionClassTest(-2, testRand = FALSE))
    expect_true(partitionClassTest(-1, IsComposition = TRUE, testRand = FALSE))
    expect_true(partitionClassTest(-2, IsComposition = TRUE, testRand = FALSE))
    expect_true(partitionClassTest(-1, rep = TRUE, testRand = FALSE))
    expect_true(partitionClassTest(-2, rep = TRUE, testRand = FALSE))
    expect_true(partitionClassTest(-1:0, rep = TRUE, testRand = FALSE))
    expect_true(partitionClassTest(-2:0, rep = TRUE, testRand = FALSE))
    expect_true(partitionClassTest(-1:0, rep = TRUE, tar = -1,
                                   testRand = FALSE))
    expect_true(partitionClassTest(-2:0, rep = TRUE, tar = -2,
                                   testRand = FALSE))
    expect_true(partitionClassTest(-2:0, 2, rep = TRUE,
                                   tar = -2, testRand = FALSE,
                                   IsComposition = TRUE))
    expect_true(partitionClassTest(-2:0, 2, rep = TRUE,
                                   tar = -2, testRand = FALSE,
                                   IsComposition = TRUE, IsWeak = TRUE))
    expect_true(partitionClassTest(0:3, fr = c(2, rep(1, 3)),
                                   IsComposition = TRUE,
                                   requiresWidthRebuild = TRUE))

    #### Distinct; Length determined internally; No zero;
    expect_true(partitionClassTest(189))
    expect_true(partitionClassTest(35, IsComposition = TRUE))
    expect_true(partitionClassTest(0:10, fr = c(2, rep(1, 10)),
                                   IsComposition = TRUE,
                                   requiresWidthRebuild = TRUE))
    expect_true(partitionClassTest(0:16, fr = c(4, rep(1, 16)),
                                   IsComposition = TRUE,
                                   requiresWidthRebuild = TRUE))

    #### Distinct; Length determined internally; One zero;
    expect_true(partitionClassTest(0:50))
    expect_true(partitionClassTest(0:30, IsComposition = TRUE,
                                   requiresWidthRebuild = TRUE))
    expect_true(partitionClassTest(0:25, IsComposition = TRUE, IsWeak = TRUE))

    #### Distinct; Specific Length; No zero
    expect_true(partitionClassTest(50, 5))
    expect_true(partitionClassTest(40, 5, IsComposition = TRUE))

    #### Mapped version
    ## 50 * 3 + 6 * 5 = 180
    expect_true(partitionClassTest(6 + (1:50) * 3, 5, tar = 180))

    #### Compositions Mapped version
    ## 40 * 3 + 6 * 5 = 150
    expect_true(partitionClassTest(6 + (1:40) * 3, 5, tar = 150,
                                   IsComposition = TRUE))

    #### Distinct; Specific Length; One zero
    expect_true(partitionClassTest(0:50, 5))
    expect_true(partitionClassTest(0:40, 5, IsComposition = TRUE,
                                   requiresWidthRebuild = TRUE))
    expect_true(partitionClassTest(0:40, 5, IsComposition = TRUE,
                                   IsWeak = TRUE))
    expect_true(partitionClassTest(0:15, 5, tar = 30,
                                   IsComposition = TRUE, IsWeak = TRUE))

    #### Mapped version
    expect_true(partitionClassTest(6 + (0:50) * 3, 5, tar = 180))
    expect_true(partitionClassTest(6 + (0:40) * 3, 5, tar = 150,
                                   IsComposition = TRUE))
    expect_true(partitionClassTest(6 + (0:40) * 3, 5, tar = 150,
                                   IsComposition = TRUE, IsWeak = TRUE))

    #### Distinct; Specific Length; Multiple Zeros; Not enough to maximize
    expect_true(partitionClassTest(0:50, 9, fr = c(4, rep(1, 50))))
    expect_true(partitionClassTest(0:30, 7, fr = c(4, rep(1, 30)),
                                   IsComposition = TRUE,
                                   requiresWidthRebuild = TRUE))
    expect_true(partitionClassTest(0:25, 6, fr = c(4, rep(1, 25)),
                                   IsComposition = TRUE, IsWeak = TRUE,
                                   testRand = FALSE))

    #### Mapped version
    ## 50 * 13 + 7 * 9 = 713
    expect_true(partitionClassTest(7 + (0:50) * 13, 9,
                                   fr = c(4, rep(1, 50)), tar = 713))

    ## Currently, we don't have an algorithm for this case. We have to think
    ## carefully here... Since there technically isn't an actual zero, the
    ## idea of weakness doesn't come into play even though we can map a value
    ## isomorphically to zero. Because of this coupled with the fact we are
    ## dealing with compositions where order matters, we technically have
    ## a strange compositions multiset case. Our current algorithm may in
    ## fact work, but testing is still needed.
    ##
    ## This will throw the error:
    ##
    ## Error: Currently, there is no composition algorithm for this case.
    ##  Use permuteCount, permuteIter, permuteGeneral, permuteSample, or
    ##  permuteRank instead.
    expect_error(partitionClassTest(7 + (0:30) * 13, 7, fr = c(4, rep(1, 30)),
                                    tar = 439, IsComposition = TRUE),
                 "Currently, there is no composition algorithm for this case")

    ## If we don't shift the results, we should be able to yield results.
    ## - target will be determined internally to 13 * 30 = 390
    expect_true(partitionClassTest((0:30) * 13, 7, fr = c(4, rep(1, 30)),
                                   IsComposition = TRUE,
                                   requiresWidthRebuild = TRUE))
    expect_true(partitionClassTest((0:25) * 13, 6, fr = c(4, rep(1, 25)),
                                   IsComposition = TRUE, IsWeak = TRUE,
                                   testRand = FALSE))

    #### Distinct; Specific Length; Multiple Zeros; Enough to maximize;
    #### Length is restrictive
    expect_true(partitionClassTest(0:50, 5, fr = c(8, rep(1, 50))))
    expect_true(partitionClassTest(0:30, 5, fr = c(7, rep(1, 30)),
                                   IsComposition = TRUE,
                                   requiresWidthRebuild = TRUE))
    expect_true(partitionClassTest(0:30, 5, fr = c(7, rep(1, 30)),
                                   IsComposition = TRUE, IsWeak = TRUE,
                                   testRand = FALSE))

    #### Mapped version
    ## 50 * 13 + 7 * 5 = 713
    ##
    ## 30 * 13 + 7 * 5 = 425
    expect_true(partitionClassTest(7 + (0:50) * 13, 5,
                                   fr = c(8, rep(1, 50)), tar = 685))

    ## Same as above... this will throw an error
    expect_error(partitionClassTest(7 + (0:30) * 13, 5,
                                    fr = c(7, rep(1, 30)), tar = 425,
                                    IsComposition = TRUE),
                 "Currently, there is no composition algorithm for this case")

    ## Again not shifting the source vector yields results.
    ## 30 * 13 = 390
    expect_true(partitionClassTest((0:30) * 13, 5, fr = c(7, rep(1, 30)),
                                   IsComposition = TRUE,
                                   requiresWidthRebuild = TRUE))
    expect_true(partitionClassTest((0:30) * 13, 5, fr = c(7, rep(1, 30)),
                                   IsComposition = TRUE, IsWeak = TRUE,
                                   testRand = FALSE))

    #### Distinct; Length determined internally; Multiple Zeros;
    #### Enough to maximize;
    expect_true(partitionClassTest(0:50, fr = c(50, rep(1, 50))))
    expect_true(partitionClassTest(0:30, fr = c(30, rep(1, 30)),
                                   IsComposition = TRUE,
                                   requiresWidthRebuild = TRUE))
    expect_true(partitionClassTest(0:25, fr = c(25, rep(1, 25)),
                                   IsComposition = TRUE, IsWeak = TRUE,
                                   testRand = FALSE))

    #### Mapped Versions
    ## N.B. We don't shift
    ## We have to explicitly set the width as the internal code will try
    ## maximize given the inputs. When we have 25, the max width is 6
    ## (i.e. sum(1:7) = 28 > 25 > sum(1:76) = 21). When we multiply by 17,
    ## we will get a huge width (40 to be exact). This causes major issues
    ## when we are dealing with the weak case.
    expect_true(partitionClassTest((0:50) * 17, 9, fr = c(50, rep(1, 50))))
    expect_true(partitionClassTest((0:30) * 17, 7, fr = c(30, rep(1, 30)),
                                   IsComposition = TRUE,
                                   requiresWidthRebuild = TRUE))
    expect_true(partitionClassTest((0:25) * 17, 6, fr = c(25, rep(1, 25)),
                                   IsComposition = TRUE, IsWeak = TRUE,
                                   testRand = FALSE))

    #### Distinct; Specific Length; No Zeros; Specific Target;
    expect_true(partitionClassTest(30, 8, tar = 75))
    expect_true(
        partitionClassTest(
            20, 5, tar = 55, IsComposition = TRUE, testRand = FALSE
        )
    )

    #### Mapped Versions
    ## 8 * 97 + 75 * 3 = 1001
    ##
    ## 5 * 97 + 55 * 3 = 650
    expect_true(partitionClassTest(97L + (1:30) * 3L, 8, tar = 1001L))
    expect_true(
        partitionClassTest(
            97L + (1:20) * 3L, 5, tar = 650L,
            IsComposition = TRUE, testRand = FALSE
        )
    )

    #### Distinct; Specific Length; Multi-Zeros; Specific Target;
    expect_true(partitionClassTest(0:30, 6, tar = 75, fr = c(3, rep(1, 30)),
                                   requiresWidthRebuild = TRUE))
    expect_true(partitionClassTest(0:10, 6, tar = 33, fr = c(3, rep(1, 10)),
                                   IsComposition = TRUE,
                                   requiresWidthRebuild = TRUE))
    expect_true(partitionClassTest(0:10, 6, tar = 33, fr = c(3, rep(1, 10)),
                                   IsComposition = TRUE, IsWeak = TRUE,
                                   testRand = FALSE))
    expect_true(partitionClassTest(0:15, 7, fr = c(4, rep(1, 15)), tar = 32,
                                   IsComposition = TRUE, testRand = FALSE,
                                   requiresWidthRebuild = TRUE))
    expect_true(partitionClassTest(0:15, 7, fr = c(4, rep(1, 15)), tar = 32,
                                   IsComposition = TRUE, IsWeak = TRUE,
                                   testRand = FALSE))

    #### Mapped Versions
    ## 13 * 75 = 975
    ##
    ## 13 * 33 = 429
    expect_true(
        partitionClassTest((0:30) * 13L, 6, tar = 975L, fr = c(3, rep(1, 30)),
                           requiresWidthRebuild = TRUE)
    )
    expect_true(
        partitionClassTest((0:10) * 13L, 6, tar = 429L, fr = c(3, rep(1, 10)),
                           IsComposition = TRUE, requiresWidthRebuild = TRUE)
    )
    expect_true(
        partitionClassTest((0:10) * 13L, 6, tar = 429L, fr = c(3, rep(1, 10)),
                           IsComposition = TRUE, IsWeak = TRUE,
                           testRand = FALSE)
    )

    #### Repetition; Length determined internally; Multiple Zero;
    expect_true(partitionClassTest(0:30, rep = TRUE))
    #### Mapped version
    ## 19 * 30 + 30 * 3e9 = 90000000570
    expect_true(partitionClassTest(19 + (0:30) * 3e9, 30,
                                   rep = TRUE, tar = 90000000570))

    expect_true(partitionClassTest(0:15, rep = TRUE,
                                   IsComposition = TRUE,
                                   requiresWidthRebuild = TRUE))
    expect_true(partitionClassTest(0:7, rep = TRUE, IsWeak = TRUE,
                                   IsComposition = TRUE))
    #### Mapped version
    ## 15 * 3e9 = 45000000000
    comp <- compositionsGeneral((0:15) * 3e9, 15, repetition = TRUE,
                                target = 45000000000)
    expect_equal(nrow(comp), compositionsCount(0:15, repetition = TRUE))
    expect_equal(comp[nrow(comp), ], rep(3e9, 15))
    expect_equal(comp[1, ], c(rep(0, 14), 45000000000))
    expect_true(partitionClassTest((0:15) * 3e9, 15,
                                   rep = TRUE, tar = 45000000000,
                                   IsComposition = TRUE,
                                   requiresWidthRebuild = TRUE))

    #### Mapped version
    ## 7 * 3e9 = 21000000000
    comp <- compositionsGeneral((0:7) * 3e9, 7, repetition = TRUE,
                                weak = TRUE, target = 21000000000)
    expect_equal(nrow(comp), compositionsCount(0:7, repetition = TRUE,
                                               weak = TRUE))
    expect_equal(comp[nrow(comp), ], c(21000000000, rep(0, 6)))
    expect_equal(comp[1, ], c(rep(0, 6), 21000000000))
    expect_true(partitionClassTest((0:7) * 3e9, 7, rep = TRUE,
                                   IsWeak = TRUE, tar = 21000000000,
                                   IsComposition = TRUE))

    #### Repetition; Specific Length; No zero
    expect_true(partitionClassTest(50, 5, TRUE))
    #### Mapped version
    ## 19 * 5 + 50 * 3 = 245
    expect_true(partitionClassTest(19 + (1:50) * 3, 5, TRUE, tar = 245))

    expect_true(partitionClassTest(20, 5, TRUE, IsComposition = TRUE))
    #### Mapped version
    ## 20 * 3 = 60
    expect_true(partitionClassTest((1:20) * 3, 5, TRUE, tar = 60,
                                   IsComposition = TRUE))

    #### Repetition; Specific Length; Zero included
    expect_true(partitionClassTest(0:30, 10, rep = TRUE))
    #### Mapped version
    ## 19 * 10 + 30 * 3 = 280
    expect_true(partitionClassTest(19 + (0:30) * 3, 10,
                                   rep = TRUE, tar = 280))

    expect_true(partitionClassTest(0:20, 5, rep = TRUE,
                                   IsComposition = TRUE,
                                   requiresWidthRebuild = TRUE))
    expect_true(partitionClassTest(0:20, 5, rep = TRUE,
                                   IsComposition = TRUE,
                                   IsWeak = TRUE))

    #### Mapped version
    ## 20 * 3 = 60
    expect_true(partitionClassTest((0:20) * 3, 5,
                                   rep = TRUE, tar = 60,
                                   IsComposition = TRUE,
                                   requiresWidthRebuild = TRUE))
    expect_true(partitionClassTest((0:20) * 3, 5, IsWeak = TRUE,
                                   rep = TRUE, tar = 60,
                                   IsComposition = TRUE))

    #### Repetition; Specific Length; No Zeros; Specific Target;
    expect_true(partitionClassTest(20, 10, rep = TRUE, tar = 45))

    #### Multiset; class table;
    expect_true(partitionClassTest(table(sample(10, 100, TRUE)),
                                   15, tar = 55, testRand = FALSE))

    #### Multiset: Specific Length;
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
                                  tar = 60, testRand = FALSE,
                                  requiresWidthRebuild = TRUE))

    ##******** BIG TESTS *********##
    partitionClassBigZTest <- function(v_pass, m_pass = NULL, rep = FALSE,
                                       fr = NULL, tar = NULL, lenCheck = 1000,
                                       IsComposition = FALSE, IsWeak = FALSE) {

        myResults <- vector(mode = "logical")

        if (IsComposition) {
            myRows <- compositionsCount(v_pass, m_pass, rep, fr, tar, IsWeak)
            a  <- compositionsIter(
                v_pass, m_pass, rep, fr, tar, IsWeak, nThreads = 2
            )
            b1 <- compositionsGeneral(v_pass, m_pass, rep, fr, tar, IsWeak,
                                      upper = lenCheck)
            b2 <- compositionsGeneral(
                v_pass, m_pass, rep, fr, tar, IsWeak,
                lower = gmp::sub.bigz(myRows, lenCheck - 1)
            )
        } else {
            myRows <- partitionsCount(v_pass, m_pass, rep, fr, tar)
            a  <- partitionsIter(v_pass, m_pass, rep, fr, tar, nThreads = 2)
            b1 <- partitionsGeneral(v_pass, m_pass, rep, fr, tar, upper = lenCheck)
            b2 <- partitionsGeneral(v_pass, m_pass, rep, fr, tar,
                                    lower = gmp::sub.bigz(myRows, lenCheck - 1))
        }

        myResults <- c(myResults, isTRUE(all.equal(
            a@summary()$totalResults, myRows)
        ))

        myResults <- c(myResults, class(myRows) == "bigz")

        if (!is.null(tar)) {
            myResults <- c(myResults, all(rowSums(b1) == tar))
            myResults <- c(myResults, all(rowSums(b2) == tar))
        } else {
            myResults <- c(myResults, all(rowSums(b1) == sum(b1[1, ])))
            myResults <- c(myResults, all(rowSums(b2) == sum(b1[1, ])))
        }

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

        samp1 <- sample(lenCheck, 2)
        samp2 <- gmp::sub.bigz(myRows, lenCheck) + gmp::as.bigz(samp1)
        myResults <- c(myResults, isTRUE(all.equal(a[[samp1]], b1[samp1, ])))
        myResults <- c(myResults, isTRUE(all.equal(a[[samp2]], b2[samp1, ])))
        rm(a, a1, b1, b2)
        gc()
        all(myResults)
    }

    expect_true(partitionClassBigZTest(2000, 10, TRUE))
    #### Mapped version
    ## 17 * 10 + 2000 * 123456789 = 246913578170
    expect_true(partitionClassBigZTest(17 + (1:2000) * 123456789,
                                       10, TRUE, tar = 246913578170))

    expect_true(partitionClassBigZTest(2000, 10, TRUE, IsComposition = TRUE))
    #### Mapped version
    ## 2000 * 123456789 = 246913578000
    expect_true(partitionClassBigZTest((1:2000) * 123456789, 10, TRUE,
                                       IsComposition = TRUE,
                                       tar = 246913578000))

    expect_true(partitionClassBigZTest(0:150, rep = TRUE, IsComposition = TRUE))
    #### Mapped version
    ## 150 * 123456789 = 18518518350
    expect_true(partitionClassBigZTest((0:150) * 123456789, rep = TRUE,
                                       IsComposition = TRUE,
                                       tar = 18518518350))
    expect_true(partitionClassBigZTest(2000, 10))
    expect_true(partitionClassBigZTest(300, 10, IsComposition = TRUE))
    expect_true(partitionClassBigZTest(0:300, 10, IsComposition = TRUE))
    expect_true(partitionClassBigZTest(0:300, 10, fr = c(7, rep(1, 300)),
                                       IsComposition = TRUE))
    expect_true(partitionClassBigZTest(0:200, IsComposition = TRUE))

    # Lots of results: weak compositions with repetition
    expect_true(partitionClassBigZTest(
        0:250, 20, rep = TRUE, IsWeak = TRUE,
        IsComposition = TRUE, lenCheck = 1000)
    )

    # Same, mapped (slope big); target = 40 * 98765431
    expect_true(partitionClassBigZTest((0:250) * 98765431, 20, rep = TRUE,
                                       IsComposition = TRUE, IsWeak = TRUE,
                                       lenCheck = 1000))

    # Distinct weak comps with 0 included (lots of results)
    expect_true(partitionClassBigZTest(0:400, 15, rep = FALSE, IsWeak = TRUE,
                                       IsComposition = TRUE, lenCheck = 1000))

    # Mapped distinct weak comps (0 included)
    # slope=13, shift=7; choose tar divisible-ish
    expect_true(partitionClassBigZTest((0:400) * 13, 15, IsWeak = TRUE,
                                       IsComposition = TRUE,
                                       lenCheck = 1000))

    # Distinct partitions, no zero, fixed width, nontrivial target
    expect_true(partitionClassBigZTest(500, 20, tar = 1000, lenCheck = 1000))

    # Mapped distinct partitions, no zero
    expect_true(partitionClassBigZTest(19 + (1:500) * 3, 20,
                                       tar = 19 * 20 + 3 * 1000,
                                       lenCheck = 1000))

    expect_true(partitionClassBigZTest(
        220, 25, rep = TRUE, tar = 440, lenCheck = 1000
    ))
    expect_true(partitionClassBigZTest(
        300, 26, rep = TRUE, tar = 600, lenCheck = 1000
    ))

    expect_true(partitionClassBigZTest(
        0:280, 15, rep = TRUE, tar = 600, lenCheck = 1000
    ))

    # tar = shift*m + slope*baseTar
    # baseTar = 240, m=12, slope=7, shift=5  => tar = 5*25 + 7*440 = 1740
    expect_true(partitionClassBigZTest(
        5 + (0:220) * 7, 25, rep = TRUE, tar = 3205, lenCheck = 1000
    ))

    # No-zero mapped (start at 1)
    # baseTar = 240, m=12, slope=5, shift=9 => tar = 9*26 + 5*600 = 3234
    expect_true(partitionClassBigZTest(
        9 + (1:300) * 5, 26, rep = TRUE, tar = 3234, lenCheck = 1000)
    )

    expect_true(partitionClassBigZTest(
        220, 25, rep = FALSE, tar = 840, lenCheck = 1000)
    )
    expect_true(partitionClassBigZTest(
        300, 26, rep = FALSE, tar = 999, lenCheck = 1000)
    )

    expect_true(partitionClassBigZTest(0:220, 25, tar = 840, lenCheck = 1000))
    expect_true(partitionClassBigZTest(0:300, 26, tar = 999, lenCheck = 1000))

    expect_true(
        partitionClassBigZTest(0:25, 15, IsComposition = TRUE, IsWeak = FALSE,
                               tar = 200, lenCheck = 1000)
    )

    # No zero compositions (start at 1)
    expect_true(partitionClassBigZTest(
        25, 15, IsComposition = TRUE, IsWeak = FALSE, tar = 200
    ))

    expect_true(partitionClassBigZTest(
        0:20, 17, IsComposition = TRUE, IsWeak = TRUE, tar = 180
    ))
    expect_true(partitionClassBigZTest(
        0:100, 10, IsComposition = TRUE, IsWeak = TRUE, tar = 300
    ))

    # baseTar=360, m=10, slope=7, shift=6 => tar = 6*10 + 7*360 = 2580
    expect_true(partitionClassBigZTest(6 + (0:80) * 7, 10,
                                       IsComposition = TRUE, IsWeak = FALSE,
                                       tar = 2580, lenCheck = 1000))

    # Weak mapped
    # baseTar=300, m=10, slope=5 => tar = 5*300 = 1500
    expect_true(partitionClassBigZTest((0:100) * 5, 10,
                                       IsComposition = TRUE, IsWeak = TRUE,
                                       tar = 1500, lenCheck = 1000))

    expect_true(partitionClassBigZTest(0:200, 15, fr = c(9, rep(1, 200)),
                                       tar = 700, lenCheck = 1000))

    expect_true(partitionClassBigZTest(0:800, 15, lenCheck = 1000))
    expect_true(partitionClassBigZTest(0:800, 15, fr = c(8, rep(1, 800))))
    expect_true(partitionClassBigZTest(0:600, fr = c(600, rep(1, 600))))
    expect_true(partitionClassBigZTest(0:600, 15, rep = TRUE))
    expect_true(partitionClassBigZTest(0:300, rep = TRUE))
})
