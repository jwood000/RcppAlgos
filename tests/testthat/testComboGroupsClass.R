test_that("comboGroupsIter produces correct results", {

    comboGroupsClassTest <- function(
        v_pass, n_grps = NULL, grp_sizes = NULL,
        ret = "matrix", testRand = TRUE
    ) {

        myResults <- vector(mode = "logical")

        myRows <- comboGroupsCount(v_pass, n_grps, grp_sizes)
        a <- comboGroupsIter(v_pass, n_grps, grp_sizes, ret)
        b <- comboGroups(v_pass, n_grps, grp_sizes, ret)

        myResults <- c(myResults, isTRUE(all.equal(
            a@summary()$totalResults, myRows)
        ))

        if (length(v_pass) == 1 && v_pass == 0) {
            myResults <- c(myResults, v_pass == a@sourceVector())
        } else if (length(v_pass) == 1) {
            myResults <- c(myResults, isTRUE(
                all.equal(abs(v_pass), length(a@sourceVector()))
            ))
        } else if (is.integer(v_pass) || is.numeric(v_pass)) {
            myResults <- c(myResults, isTRUE(
                all.equal(sort(v_pass), sort(a@sourceVector()))
            ))
        } else {
            myResults <- c(myResults, isTRUE(
                all.equal(v_pass, a@sourceVector())
            ))
        }

        if (testRand && ret == "matrix") {
            myResults <- c(myResults, isTRUE(
                all.equal(a@front(), b[1, ])
            ))
            myResults <- c(myResults, isTRUE(all.equal(a@currIter(),
                                                       b[1, ])))
            myResults <- c(myResults, isTRUE(all.equal(a@back(),
                                                       b[myRows, ])))
            myResults <- c(myResults, isTRUE(all.equal(a@currIter(),
                                                       b[myRows, ])))
        } else if (testRand) {
            myResults <- c(myResults, isTRUE(
                all.equal(a@front(), b[1, , ])
            ))
            myResults <- c(myResults, isTRUE(all.equal(a@currIter(),
                                                       b[1, , ])))
            myResults <- c(myResults, isTRUE(all.equal(a@back(),
                                                       b[myRows, , ])))
            myResults <- c(myResults, isTRUE(all.equal(a@currIter(),
                                                       b[myRows, , ])))
        }

        a@startOver()
        msg <- capture.output(noMore <- a@currIter())
        myResults <- c(myResults, is.null(noMore))
        myResults <- c(myResults, grepl("Iterator Initialized. To see the first", msg[1]))
        a1 <- b

        if (myRows) {
            if (ret == "matrix") {
                for (i in 1:myRows) {
                    a1[i, ] <- a@nextIter()
                }
            } else {
                for (i in 1:myRows) {
                    a1[i, , ] <- a@nextIter()
                }
            }

            myResults <- c(myResults, isTRUE(all.equal(a1, b)))
            a@startOver()
            num_iters <- if (myRows > 10) 3L else 1L
            numTest   <- as.integer(myRows / num_iters);

            s <- 1L
            e <- numTest

            if (ret == "matrix") {
                for (i in 1:num_iters) {
                    myResults <- c(
                        myResults, isTRUE(
                            all.equal(a@nextNIter(numTest),
                                      b[s:e, , drop = FALSE])
                            )
                    )
                    s <- e + 1L
                    e <- e + numTest
                }
            } else {
                for (i in 1:num_iters) {
                    myResults <- c(
                        myResults, isTRUE(
                            all.equal(a@nextNIter(numTest),
                                      b[s:e, , ])
                        )
                    )
                    s <- e + 1L
                    e <- e + numTest
                }
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
                one_samp <- sample(myRows, 1)

                if (ret == "matrix") {
                    myResults <- c(
                        myResults, isTRUE(all.equal(a[[samp]], b[samp, ]))
                    )

                    myResults <- c(
                        myResults, isTRUE(all.equal(a[[one_samp]],
                                                    b[one_samp, ]))
                    )
                } else {
                    myResults <- c(
                        myResults, isTRUE(all.equal(a[[samp]], b[samp, , ]))
                    )

                    myResults <- c(
                        myResults, isTRUE(all.equal(a[[one_samp]],
                                                    b[one_samp, , ]))
                    )
                }
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

    expect_true(comboGroupsClassTest(12, 3))
    expect_true(comboGroupsClassTest(runif(10) + rnorm(10) * 1i, 2))
    expect_true(comboGroupsClassTest(as.raw(1:10), grp_sizes = 1:4))
    expect_true(comboGroupsClassTest(sample(state.name, 10), grp_sizes = 4:1))
    expect_true(comboGroupsClassTest(12, 3, ret = "3Darray"))
    expect_true(comboGroupsClassTest(LETTERS[1:8], 2, ret = "3Darray"))
    expect_true(comboGroupsClassTest(as.factor(LETTERS[1:8]),
                                     4, ret = "3Darray"))
    expect_true(comboGroupsClassTest(8, 1))
    expect_true(comboGroupsClassTest(8, 8))

    ## Not tested right now. Sent email to r-devel
    # expect_true(comboGroupsClassTest(8, 1, ret = "3Darray"))
    # expect_true(comboGroupsClassTest(8, 8, ret = "3Darray"))

    expect_true(comboGroupsClassTest(12, grp_sizes = c(3, 3, 6)))
    expect_true(comboGroupsClassTest(12, grp_sizes = c(3, 4, 5)))
    expect_true(comboGroupsClassTest(11, grp_sizes = c(1, 2, 2, 3, 3)))
    expect_true(comboGroupsClassTest(1:4 + 0.1, grp_sizes = c(1, 1, 1, 1)))
    expect_true(comboGroupsClassTest(LETTERS[1:4], grp_sizes = c(1, 1, 2)))
    expect_true(comboGroupsClassTest(as.factor(state.abb[1:3]),
                                     grp_sizes = c(1, 2)))
    expect_true(comboGroupsClassTest(3, 1))
    expect_true(comboGroupsClassTest(1, 1))


    ###****************************** BIG TESTS *****************************###
    comboGroupsClassBigZTest <- function(
        v_pass, n_grps = NULL, grp_sizes = NULL,
        ret = "matrix", lenCheck = 1000, testRand = TRUE
    ) {

        myResults <- vector(mode = "logical")

        myRows <- comboGroupsCount(v_pass, n_grps, grp_sizes)
        a  <- comboGroupsIter(v_pass, n_grps, grp_sizes, ret)
        b1 <- comboGroups(v_pass, n_grps, grp_sizes, ret, upper = lenCheck)
        b2 <- comboGroups(v_pass, n_grps, grp_sizes, ret,
                          lower = gmp::sub.bigz(myRows, lenCheck - 1))

        myResults <- c(myResults, isTRUE(all.equal(
            a@summary()$totalResults, myRows)
        ))

        if (length(v_pass) == 1) {
            myResults <- c(myResults, isTRUE(
                all.equal(v_pass, length(a@sourceVector()))
            ))
        } else if (is.integer(v_pass) || is.numeric(v_pass)) {
            myResults <- c(myResults, isTRUE(
                all.equal(sort(v_pass), sort(a@sourceVector()))
            ))
        } else {
            myResults <- c(myResults, isTRUE(
                all.equal(v_pass, a@sourceVector())
            ))
        }

        if (ret == "matrix") {
            myResults <- c(myResults, isTRUE(
                all.equal(a@front(), b1[1 ,])
            ))
            myResults <- c(myResults, isTRUE(all.equal(a@currIter(),
                                                       b1[1 ,])))
            myResults <- c(myResults, isTRUE(all.equal(a@back(),
                                                       b2[lenCheck, ])))
            myResults <- c(myResults, isTRUE(all.equal(a@currIter(),
                                                       b2[lenCheck, ])))
        } else {
            myResults <- c(myResults, isTRUE(
                all.equal(a@front(), b1[1, ,])
            ))
            myResults <- c(myResults, isTRUE(all.equal(a@currIter(),
                                                       b1[1, ,])))
            myResults <- c(myResults, isTRUE(all.equal(a@back(),
                                                       b2[lenCheck, ,])))
            myResults <- c(myResults, isTRUE(all.equal(a@currIter(),
                                                       b2[lenCheck, ,])))
        }

        a@startOver()
        a1 <- b1

        if (ret == "matrix") {
            for (i in 1:lenCheck) {
                a1[i, ] <- a@nextIter()
            }
        } else {
            for (i in 1:lenCheck) {
                a1[i, ,] <- a@nextIter()
            }
        }

        myResults <- c(myResults, isTRUE(all.equal(a1, b1)))
        a@startOver()
        numTest <- as.integer(lenCheck / 3);
        s <- 1L
        e <- numTest

        if (ret == "matrix") {
            for (i in 1:3) {
                myResults <- c(myResults, isTRUE(all.equal(a@nextNIter(numTest),
                                                           b1[s:e, ])))
                s <- e + 1L
                e <- e + numTest
            }
        } else {
            for (i in 1:3) {
                myResults <- c(myResults, isTRUE(all.equal(a@nextNIter(numTest),
                                                           b1[s:e, ,])))
                s <- e + 1L
                e <- e + numTest
            }
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

        if (ret == "matrix") {
            myResults <- c(myResults, isTRUE(all.equal(a[[samp1]], b1[samp1, ])))
            myResults <- c(myResults, isTRUE(all.equal(a[[samp2]], b2[samp1, ])))
        } else {
            myResults <- c(myResults, isTRUE(all.equal(a[[samp1]], b1[samp1, ,])))
            myResults <- c(myResults, isTRUE(all.equal(a[[samp2]], b2[samp1, ,])))
        }

        rm(a, a1, b1, b2)
        gc()
        all(myResults)
    }

    expect_true(comboGroupsClassBigZTest(50, 10))
    expect_true(comboGroupsClassBigZTest(
        runif(50) + rnorm(50) * 1i, 5)
    )
    expect_true(comboGroupsClassBigZTest(
        as.raw(1:45), grp_sizes = 1:9)
    )
    expect_true(comboGroupsClassBigZTest(
        sample(state.name, 45), grp_sizes = 9:1)
    )
    expect_true(
        comboGroupsClassBigZTest(50, 10, ret = "3Darray")
    )
    expect_true(comboGroupsClassBigZTest(state.abb, 10, ret = "3Darray"))
    expect_true(comboGroupsClassBigZTest(
        as.factor(state.abb), 5, ret = "3Darray"
    ))
    expect_true(
        comboGroupsClassBigZTest(40, grp_sizes = c(9, 9, 22))
    )
    expect_true(
        comboGroupsClassBigZTest(42, grp_sizes = c(10, 10, 22))
    )
    expect_true(
        comboGroupsClassBigZTest(48, grp_sizes = 15:17)
    )
    expect_true(
        comboGroupsClassBigZTest(55, grp_sizes = rep(1:5, time = 1:5))
    )
    expect_true(
        comboGroupsClassBigZTest(state.name[1:40], grp_sizes = rep(1:4, 4))
    )
    expect_true(
        comboGroupsClassBigZTest(
            as.factor(state.abb), grp_sizes =
                rep(1:7, times = c(1, 4, 1, 2, 1, 3, 1))
            )
    )
})
