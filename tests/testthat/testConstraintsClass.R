context("testing ConstraintsClass")

test_that("ConstraintsClass produces correct results", {

    constraintsClassTest <- function(v_pass, m_pass = NULL, rep = FALSE,
                                     fr = NULL, tar, fun, comp,
                                     keep = FALSE, tol = NULL,
                                     expect_zero_res = FALSE, IsComb = TRUE) {

        myResults <- vector(mode = "logical")

        if (IsComb) {
            a <- comboIter(v_pass, m_pass, rep, fr, constraintFun = fun,
                           comparisonFun = comp, limitConstraints = tar,
                           keepResults = keep, tolerance = tol)
            b <- comboGeneral(v_pass, m_pass, rep, fr, constraintFun = fun,
                              comparisonFun = comp, limitConstraints = tar,
                              keepResults = keep, tolerance = tol)
        } else {
            a <- permuteIter(v_pass, m_pass, rep, fr, constraintFun = fun,
                             comparisonFun = comp, limitConstraints = tar,
                             keepResults = keep, tolerance = tol)
            b <- permuteGeneral(v_pass, m_pass, rep, fr, constraintFun = fun,
                                comparisonFun = comp, limitConstraints = tar,
                                keepResults = keep, tolerance = tol)

            if (fun == "sum" && length(comp) == 1 && comp == "==") {
                parts <- partitionsGeneral(
                    v_pass, m_pass, rep, freqs = fr, target = tar
                )

                res <- do.call(
                    rbind,
                    lapply(
                        seq_len(nrow(parts)), function(i) {
                            permuteGeneral(table(parts[i, ]))
                        }
                    )
                )

                myResults <- c(myResults, identical(res, b))
            }
        }

        myRows <- nrow(b)
        myResults <- c(myResults,
                       if (expect_zero_res) myRows == 0 else myRows > 0)
        myResults <- c(myResults, is.na(a@summary()$totalResults))

        if (length(v_pass) == 1) {
            myResults <- c(myResults, isTRUE(
                all.equal(abs(v_pass), length(a@sourceVector()))
            ))
        } else {
            myResults <- c(myResults, isTRUE(
                all.equal(sort(v_pass), a@sourceVector())
            ))
        }

        a1 <- b

        if (myRows) {
            a@nextIter()
            myResults <- c(myResults, all.equal(a@currIter(), b[1, ]))
            a@startOver()

            for (i in 1:myRows) {
                a1[i, ] <- a@nextIter()
            }

            myResults <- c(myResults, isTRUE(all.equal(a1, b)))
        }

        msg <- capture.output(noMore <- a@nextIter())
        myResults <- c(myResults, is.null(noMore))
        myResults <- c(myResults, msg[1] == "No more results.")
        myResults <- c(myResults, is.null(a@nextNIter(1)))
        myResults <- c(myResults, is.null(a@nextRemaining()))
        myResults <- c(myResults, is.null(a@nextIter()))
        a@startOver()

        if (myRows) {
            numTest <- as.integer(myRows / 3);

            s <- 1L
            e <- numTest

            for (i in 1:3) {
                myResults <- c(myResults, isTRUE(all.equal(a@nextNIter(numTest),
                                                           b[s:e, ])))
                s <- e + 1L
                e <- e + numTest
            }
        }

        idx <- a@summary()$currentIndex

        while (idx < myRows) {
            a@nextNIter(1)
            idx <- idx + 1L
        }

        capture.output(noMore <- a@nextNIter(1))
        myResults <- c(myResults, is.null(noMore))
        a@startOver()
        if (myRows) myResults <- c(myResults, isTRUE(all.equal(a@nextRemaining(), b)))

        a@startOver()
        if (myRows) tmp <- a@nextNIter(myRows)
        msg <- capture.output(noMore <- a@nextNIter(1))
        myResults <- c(myResults, is.null(noMore))
        myResults <- c(myResults, msg[1] == "No more results.")

        a@startOver()
        if (myRows) tmp <- a@nextNIter(myRows)
        msg <- capture.output(noMore <- a@nextRemaining())
        myResults <- c(myResults, is.null(noMore))
        myResults <- c(myResults, msg[1] == "No more results.")

        rm(a, a1, b)
        gc()
        all(myResults)
    }

    ## no viable results 1st
    expect_true(constraintsClassTest(10, 7, fun = "sum",
                                     comp = c(">","<"), tar = c(49, 51),
                                     expect_zero_res = TRUE))
    expect_true(constraintsClassTest(10, 7, fun = "sum",
                                     comp = c(">","<"), tar = c(40, 45)))

    set.seed(13)
    rSet = 1:10 + rnorm(10)
    ## no viable results 1st
    expect_true(constraintsClassTest(rSet, 7, TRUE, fun = "sum",
                                     comp = ">=", tar = 100,
                                     expect_zero_res = TRUE))
    expect_true(constraintsClassTest(rSet, 7, TRUE, fun = "sum",
                                     comp = c(">=","<="),
                                     tar = c(42.50001, 45.76277)))
    expect_true(constraintsClassTest(rSet, 7, TRUE, fun = "sum",
                                     comp = c("<=",">="),
                                     tar = c(20.05669, 60.93901), keep = TRUE))

    ## This test fails for version prior to 2.8.1
    expect_true(constraintsClassTest(rSet, 7, TRUE, fun = "sum",
                                     comp = c("<=",">="),
                                     tar = c(20.05669, 60.93901), keep = TRUE,
                                     IsComb = FALSE))

    ## no viable results 1st
    expect_true(constraintsClassTest(10, 7, fr = rep(2:3, 5), fun = "sum",
                                     comp = ">", tar = 100,
                                     expect_zero_res = TRUE))

    expect_true(constraintsClassTest(10, 7, fr = rep(3, 10), fun = "sum",
                                     comp = c("<=",">"),
                                     tar = c(50, 47), keep = TRUE))

    ## This test fails for version prior to 2.8.1
    expect_true(constraintsClassTest(10, 7, fr = rep(3, 10), fun = "sum",
                                     comp = c("<=",">"),
                                     tar = c(50, 47), keep = TRUE,
                                     IsComb = FALSE))

    expect_true(constraintsClassTest(10, 7, fr = rep(3, 10), fun = "max",
                                     comp = c("<=",">"), tar = c(9, 7)))

    expect_true(constraintsClassTest(10, 7, fr = rep(3, 10), fun = "min",
                                     comp = "==", tar = 3))

    ## no viable results 1st
    expect_true(constraintsClassTest(5, 7, TRUE, fun = "prod",
                                     comp = "<", tar = -1,
                                     expect_zero_res = TRUE))
    expect_true(constraintsClassTest(5, 7, TRUE, fun = "prod",
                                     comp = c(">=","<="),
                                     tar = c(2000, 5000)))
    ## Need special
    expect_true(constraintsClassTest(-5, 7, TRUE, fun = "prod",
                                     comp = c(">=","<="),
                                     tar = c(2000, 5000),
                                     keep = TRUE, expect_zero_res = TRUE))
    expect_true(constraintsClassTest(-5, 7, TRUE, fun = "prod",
                                     comp = c("<=",">="),
                                     tar = c(-2000, -5000),
                                     keep = TRUE))

    set.seed(42)
    s <- runif(10, -5, 5)
    expect_true(constraintsClassTest(
        s, 5, fr = rep(2:3, 5), fun = "prod",
        comp = ">", tar = 1000, keep = TRUE
    ))

    expect_true(constraintsClassTest(
        s, 5, fun = "prod", comp = "==", tar = 100, tol = 10, keep = TRUE
    ))

    ## Testing sums in a range
    expect_true(constraintsClassTest(
        c(NA, 1:10), 8, TRUE, fun = "sum", comp = c("=>","=<"), tar = c(72, 78)
    ))

    ## This test fails for version prior to 2.8.1
    expect_true(constraintsClassTest(
        c(NA, 1:10), 8, TRUE, fun = "sum",
        comp = c("=>","=<"), tar = c(72, 78), IsComb = FALSE
    ))

    ## This test fails for version prior to 2.8.1
    expect_true(constraintsClassTest(0:10, 8, TRUE, fun = "sum",
                                     comp = c(">","<"), tar = c(9, 11),
                                     IsComb = FALSE))

    ## Permutation Partition tests
    ## "PrmRepPartNoZ"
    expect_true(
        constraintsClassTest(
            20, 5, TRUE, fun = "sum", comp = "==", tar = 20, IsComb = FALSE
        )
    )

    ## "PrmRepPart"
    expect_true(
        constraintsClassTest(
            0:20, 5, TRUE, fun = "sum", comp = "==", tar = 20, IsComb = FALSE
        )
    )

    ## "PrmRepPartNoZ"
    expect_true(
        constraintsClassTest(
            12, 5, TRUE, fun = "sum", comp = "==", tar = 30, IsComb = FALSE
        )
    )

    ## "PrmDstPrtOneZ"
    expect_true(
        constraintsClassTest(
            0:25, 5, fun = "sum", comp = "==", tar = 25, IsComb = FALSE
        )
    )

    ## "PrmDstPartMZ"
    expect_true(
        constraintsClassTest(
            0:25, 5, fr = c(3, rep(1, 25)),
            fun = "sum", comp = "==", tar = 25, IsComb = FALSE
        )
    )

    ## Test table method
    expect_identical(
        permuteGeneral(0:25, 5, freqs = c(3, rep(1, 25)),
                       constraintFun = "sum", comparisonFun = "==",
                       limitConstraints = 25),
        permuteGeneral(table(c(0L, 0L, 0L, 1:25)), 5,
                       constraintFun = "sum", comparisonFun = "==",
                       limitConstraints = 25)
    )

    ## "PrmDstPrtCap"
    expect_true(
        constraintsClassTest(
            0:15, 5, fun = "sum", comp = "==", tar = 35, IsComb = FALSE
        )
    )

    ## "PrmDstPrtCap"
    expect_true(
        constraintsClassTest(
            15, 5, fun = "sum", comp = "==", tar = 35, IsComb = FALSE
        )
    )

    ## "PrmDstPrtCapMZ"
    expect_true(
        constraintsClassTest(
            0:15, 5, fr = c(3, rep(1, 15)),
            fun = "sum", comp = "==", tar = 35, IsComb = FALSE
        )
    )

    ## Test table method
    expect_identical(
        permuteGeneral(0:15, 5, freqs = c(3, rep(1, 15)),
                       constraintFun = "sum", comparisonFun = "==",
                       limitConstraints = 35),
        permuteGeneral(table(c(0L, 0L, 0L, 1:15)), 5,
                       constraintFun = "sum", comparisonFun = "==",
                       limitConstraints = 35)
    )

    ## "PrmMultiset"
    expect_true(
        constraintsClassTest(
            0:20, 6, fr = c(5, rep(1:4, 5)),
            fun = "sum", comp = "==", tar = 20, IsComb = FALSE
        )
    )

    ## "PrmMultiset"
    expect_true(
        constraintsClassTest(
            20, 6, fr = rep(1:4, 5),
            fun = "sum", comp = "==", tar = 20, IsComb = FALSE
        )
    )

    ## "PrmMultiset"
    expect_true(
        constraintsClassTest(
            12, 5, fr = rep(1:4, 3),
            fun = "sum", comp = "==", tar = 30, IsComb = FALSE
        )
    )

    ## "PrmMultiset"
    expect_true(
        constraintsClassTest(
            0:12, 5, fr = c(4, rep(1:4, 3)),
            fun = "sum", comp = "==", tar = 30, IsComb = FALSE
        )
    )

    expect_true(constraintsClassTest(0:10, 8, TRUE, fun = "sum",
                                     comp = c(">","<"), tar = c(9, 11),
                                     IsComb = FALSE))


    comp1 = c("<", "<=")
    comp2 = c(">", ">=")

    ## Test that unsorted vector is being handled properly
    ## for both numeric and integer type vectors
    # identical(sort(scrambled), 1:10)
    # [1] TRUE
    scrambled = as.integer(c(8, 2, 5, 1, 6, 3, 4, 7))
    scramFreqs = rep(1:5, 2)[scrambled]
    funs <- c("sum", "prod", "mean", "max", "min")
    m <- 7

    allCombs1 = lapply(funs, function(f) {
        comboGeneral(8, m, freqs = rep(1:4, 2), constraintFun = f)
    })

    ## ensure the left bound is in the solution space
    tars = lapply(allCombs1, function(x) {
        vals <- sort(x[, m + 1])
        t <- quantile(as.numeric(names(table(vals))),
                      c(0.25, 0.75), names = FALSE)
        t[1] <- vals[findInterval(t[1], vals)]
        t
    })

    for (f in seq_along(funs)) {
        for (i in 1:2) {

            if (i == 1) {
                a = comp1
                b = comp2
            } else {
                a = comp2
                b = comp1
            }

            for (j in a) {
                for (k in b) {
                    expect_true(constraintsClassTest(scrambled, m,
                                                     fr = scramFreqs,
                                                     fun = funs[f],
                                                     comp = c(j, k),
                                                     tar = tars[[f]],
                                                     tol = 0))
                }
            }
        }
    }

    allCombs1 = lapply(funs, function(f) {
        comboGeneral(8, m, TRUE, constraintFun = f, keepResults = TRUE)
    })

    ## ensure the right bound is in the solution space
    tars = lapply(allCombs1, function(x) {
        vals <- sort(x[, m + 1])
        t <- quantile(as.numeric(names(table(vals))),
                      c(0.25, 0.75), names = FALSE)
        t[2] <- vals[findInterval(t[2], vals)]
        t
    })

    for (f in seq_along(funs)) {
        for (i in 1:2) {

            if (i == 1) {
                a = comp1
                b = comp2
            } else {
                a = comp2
                b = comp1
            }

            for (j in a) {
                for (k in b) {
                    expect_true(constraintsClassTest(8, m, TRUE,
                                                     fun = funs[f],
                                                     comp = c(j, k),
                                                     tar = tars[[f]],
                                                     tol = 0))
                }
            }
        }
    }

    allCombs1 = lapply(funs, function(f) {
        comboGeneral(15, m, constraintFun = f, keepResults = TRUE)
    })

    ## ensure the both bounds are in the solution space
    tars = lapply(allCombs1, function(x) {
        vals <- sort(x[, m + 1])
        t <- quantile(as.numeric(names(table(vals))),
                      c(0.25, 0.75), names = FALSE)
        t <- vals[findInterval(t, vals)]
        t
    })

    for (f in seq_along(funs)) {
        for (i in 1:2) {

            if (i == 1) {
                a = comp1
                b = comp2
            } else {
                a = comp2
                b = comp1
            }

            for (j in a) {
                for (k in b) {
                    expect_true(constraintsClassTest(15, m,
                                                     fun = funs[f],
                                                     comp = c(j, k),
                                                     tar = tars[[f]],
                                                     tol = 0))
                }
            }
        }
    }
})
