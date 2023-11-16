context("testing comboIter & permuteIter")

test_that("comboIter & permuteIter produces correct results", {

    comboClassTest <- function(v1, m1 = NULL, rep1 = FALSE, freqs1 = NULL,
                               constr1 = NULL, compar1 = NULL, limit1 = NULL,
                               FUN1 = NULL, FUN.VALUE1 = NULL, tol1 = NULL,
                               IsComb = TRUE) {

        myResults <- vector(mode = "logical")

        myRows <- if (IsComb && class(v1) == "table") {
            comboCount(v1, m1)
        } else if (IsComb) {
            comboCount(v1, m1, rep1, freqs1)
        } else if (class(v1) == "table") {
            permuteCount(v1, m1)
        } else {
            permuteCount(v1, m1, rep1, freqs1)
        }

        a <- if (IsComb && class(v1) == "table") {
            comboIter(v1, m1, constr1, compar1, limit1, NULL,
                      FUN1, FALSE, NULL, tol1, FUN.VALUE1)
        } else if (IsComb && is.numeric(v1)) {
            comboIter(v1, m1, rep1, freqs1, constr1, compar1, limit1,
                      NULL, FUN1, FALSE, NULL, tol1, FUN.VALUE1)
        } else if (IsComb && is.list(v1)) {
            comboIter(v1, m1, rep1, freqs1)
        } else if (IsComb) {
            comboIter(v1, m1, rep1, freqs1, FUN1, FUN.VALUE1)
        } else if (class(v1) == "table") {
            permuteIter(v1, m1, constr1, compar1, limit1,
                        NULL, FUN1, FALSE, NULL, tol1, FUN.VALUE1)
        } else if (is.numeric(v1)) {
            permuteIter(v1, m1, rep1, freqs1, constr1, compar1, limit1,
                        NULL, FUN1, FALSE, NULL, tol1, FUN.VALUE1)
        } else if (is.list(v1)) {
            permuteIter(v1, m1, rep1, freqs1)
        } else {
            permuteIter(v1, m1, rep1, freqs1, FUN1, FUN.VALUE1)
        }

        b <- if (IsComb && class(v1) == "table") {
            comboGeneral(v1, m1, NULL, NULL, constr1, compar1, limit1,
                         NULL, FUN1, FALSE, NULL, tol1, FUN.VALUE1)
        } else if (IsComb && is.numeric(v1)) {
            comboGeneral(v1, m1, rep1, freqs1, NULL, NULL, constr1, compar1,
                         limit1, NULL, FUN1, FALSE, NULL, tol1, FUN.VALUE1)
        } else if (IsComb && is.list(v1)) {
            comboGeneral(v1, m1, rep1, freqs1)
        } else if (IsComb) {
            comboGeneral(v1, m1, rep1, freqs1, NULL, NULL, FUN1, FUN.VALUE1)
        } else if (class(v1) == "table") {
            permuteGeneral(v1, m1, NULL, NULL, constr1, compar1, limit1,
                           NULL, FUN1, FALSE, NULL, tol1, FUN.VALUE1)
        } else if (is.numeric(v1)) {
            permuteGeneral(v1, m1, rep1, freqs1, NULL, NULL, constr1, compar1,
                           limit1, NULL, FUN1, FALSE, NULL, tol1, FUN.VALUE1)
        } else if (is.list(v1)) {
            permuteGeneral(v1, m1, rep1, freqs1)
        } else {
            permuteGeneral(v1, m1, rep1, freqs1, NULL, NULL, FUN1, FUN.VALUE1)
        }

        # .method("summary", &Combo::summary)
        myResults <- c(myResults, isTRUE(all.equal(a@summary()$totalResults, myRows)))

        # .method("sourceVector", &Combo::sourceVector)
        if (length(v1) == 1 && is.atomic(v1) && class(v1) != "table") {
            myResults <- c(myResults, isTRUE(all.equal(v1, length(a@sourceVector()))))
        } else if (is.atomic(v1) && class(v1) != "table") {
            myResults <- c(myResults, isTRUE(all.equal(v1, a@sourceVector())))
        }

        # .method("front", &Combo::front)
        # .method("currIter", &Combo::currComb)
        # .method("back", &Combo::back)
        # .method("currIter", &Combo::currComb)
        if (is.atomic(b) && !is.matrix(b)) {
            myResults <- c(myResults, isTRUE(all.equal(a@front(), b[1])))
            myResults <- c(myResults, isTRUE(all.equal(a@currIter(), b[1])))
            myResults <- c(myResults, isTRUE(all.equal(a@back(), b[myRows])))
            myResults <- c(myResults, isTRUE(all.equal(a@currIter(), b[myRows])))
        } else if (is.list(b)) {
            myResults <- c(myResults, isTRUE(all.equal(a@front(), b[[1]])))
            myResults <- c(myResults, isTRUE(all.equal(a@currIter(), b[[1]])))
            myResults <- c(myResults, isTRUE(all.equal(a@back(), b[[myRows]])))
            myResults <- c(myResults, isTRUE(all.equal(a@currIter(), b[[myRows]])))
        } else {
            myResults <- c(myResults, isTRUE(all.equal(a@front(), b[1, ])))
            myResults <- c(myResults, isTRUE(all.equal(a@currIter(), b[1, ])))
            myResults <- c(myResults, isTRUE(all.equal(a@back(), b[myRows, ])))
            myResults <- c(myResults, isTRUE(all.equal(a@currIter(), b[myRows, ])))
        }

        # .method("startOver", &Combo::startOver)
        a@startOver()
        a@nextIter()
        msg <- capture.output(noMore <- a@prevIter())
        myResults <- c(myResults, is.null(noMore))
        myResults <- c(myResults, grepl("Iterator Initialized. To see the first", msg[1]))
        msg <- capture.output(noMore <- a@currIter())
        myResults <- c(myResults, grepl("Iterator Initialized. To see the first", msg[1]))
        myResults <- c(myResults, is.null(a@prevIter()))

        a@nextIter()
        msg <- capture.output(noMore <- a@prevNIter(1))
        myResults <- c(myResults, is.null(noMore))
        myResults <- c(myResults, grepl("Iterator Initialized. To see the first", msg[1]))
        myResults <- c(myResults, is.null(a@prevNIter(1)))

        a@nextIter()
        msg <- capture.output(noMore <- a@prevRemaining())
        myResults <- c(myResults, is.null(noMore))
        myResults <- c(myResults, grepl("Iterator Initialized. To see the first", msg[1]))
        myResults <- c(myResults, is.null(a@prevRemaining()))

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

        if (is.atomic(b) && !is.matrix(b)) {
            for (i in 1:myRows) {
                a1[i] <- a@nextIter()
            }
        } else if (is.list(b)) {
            for (i in 1:myRows) {
                a1[[i]] <- a@nextIter()
            }
        } else{
            for (i in 1:myRows) {
                a1[i, ] <- a@nextIter()
            }
        }

        # .method("nextIter", &Combo::nextComb)
        myResults <- c(myResults, isTRUE(all.equal(a1, b)))

        # .method("startOver", &Combo::startOver)
        a@startOver()
        numTest <- as.integer(myRows / 3);

        s <- 1L
        e <- numTest

        # .method("nextNIter", &Combo::nextNumCombs)
        for (i in 1:3) {
            if (is.atomic(b) && !is.matrix(b)) {
                myResults <- c(myResults, isTRUE(all.equal(a@nextNIter(numTest), b[s:e])))
            } else if (is.list(b)) {
                myResults <- c(myResults, isTRUE(all.equal(a@nextNIter(numTest), b[s:e])))
            } else {
                myResults <- c(myResults, isTRUE(all.equal(a@nextNIter(numTest), b[s:e, ])))
            }

            s <- e + 1L
            e <- e + numTest
        }

        # .method("nextRemaining", &Combo::nextGather)
        a@startOver()
        myResults <- c(myResults, isTRUE(all.equal(a@nextRemaining(), b)))

        ## Prepare a for reverse iteration
        temp <- a@back()
        t <- capture.output(a@nextIter())
        myResults <- c(myResults, is.null(a@nextIter()))
        a2 <- b

        if (is.atomic(b) && !is.matrix(b)) {
            for (i in myRows:1) {
                a2[i] <- a@prevIter()
            }
        } else if (is.list(b)) {
            for (i in myRows:1) {
                a2[[i]] <- a@prevIter()
            }
        } else {
            for (i in myRows:1) {
                a2[i, ] <- a@prevIter()
            }
        }

        # .method("prevIter", &Combo::prevComb)
        myResults <- c(myResults, isTRUE(all.equal(a2, b)))
        a@startOver()

        s <- myRows
        e <- myRows - numTest + 1L

        ## Prepare a for reverse iteration
        temp <- a@back()
        t <- capture.output(a@nextIter())
        myResults <- c(myResults, is.null(a@nextNIter()))

        # .method("prevNIter", &Combo::prevNumCombs)
        for (i in 1:3) {
            if (is.atomic(b) && !is.matrix(b)) {
                myResults <- c(myResults, isTRUE(all.equal(a@prevNIter(numTest), b[s:e])))
            } else if (is.list(b)) {
                myResults <- c(myResults, isTRUE(all.equal(a@prevNIter(numTest), b[s:e])))
            } else {
                myResults <- c(myResults, isTRUE(all.equal(a@prevNIter(numTest), b[s:e, ])))
            }

            s <- e - 1L
            e <- e - numTest
        }

        ## Prepare a for reverse iteration
        temp <- a@back()
        t <- capture.output(a@nextIter())
        myResults <- c(myResults, is.null(a@nextRemaining()))
        a3 <- a@prevRemaining()

        # .method("prevRemaining", &Combo::prevGather)
        if (is.atomic(b) && !is.matrix(b)) {
            myResults <- c(myResults, all(sapply(1:myRows, function(x) {
                isTRUE(all.equal(a3[[myRows + 1 - x]], b[x]))
            })))
        } else if (is.list(b)) {
            myResults <- c(myResults, all(sapply(1:myRows, function(x) {
                                                isTRUE(all.equal(a3[[myRows + 1 - x]], b[[x]]))
                                            })))
        } else {
            myResults <- c(myResults, all(sapply(1:myRows, function(x) {
                                                isTRUE(all.equal(a3[myRows + 1 - x, ], b[x, ]))
                                            })))
        }

        # .method("[[", &Combo::combIndex)
        samp <- sample(myRows, numTest)

        if (is.atomic(b) && !is.matrix(b)) {
            myResults <- c(myResults, isTRUE(all.equal(a[[samp]], b[samp])))
        } else if (is.list(b)) {
            myResults <- c(myResults, isTRUE(all.equal(a[[samp]], b[samp])))
        } else {
            myResults <- c(myResults, isTRUE(all.equal(a[[samp]], b[samp, ])))
        }
        rm(a, a1, a2, a3, b)
        gc()
        all(myResults)
    }

    expect_true(comboClassTest(5, 3))
    expect_true(comboClassTest(as.raw(1:5), 3))
    expect_true(comboClassTest(factor(1:5), 3))

    ## S3 table mehtod
    set.seed(32)
    s <- sample(letters[1:5], 10, TRUE)
    expect_true(comboClassTest(table(s), 3))

    ## S3 list method
    expect_true(comboClassTest(as.list(1:5), 3))

    expect_true(comboClassTest(as.complex(1:5), 3, TRUE))
    expect_true(comboClassTest(c(TRUE, FALSE), 20, TRUE))
    expect_true(comboClassTest(letters[1:5], 6, freqs1 = 1:5))

    set.seed(103)
    myNums = rnorm(5)
    expect_true(comboClassTest(myNums, 4, freqs1 = c(1, 2, 3, 2, 1)))

    ## Permutations
    expect_true(comboClassTest(5, 3, IsComb = FALSE))
    expect_true(comboClassTest(as.raw(1:5), 3, IsComb = FALSE))
    expect_true(comboClassTest(factor(1:5), 3, IsComb = FALSE))

    ## S3 table mehtod
    expect_true(comboClassTest(table(s), 3, IsComb = FALSE))

    ## S3 list method
    expect_true(comboClassTest(as.list(1:5), 3, IsComb = FALSE))

    expect_true(comboClassTest(as.complex(1:5), 3, TRUE, IsComb = FALSE))
    expect_true(comboClassTest(c(TRUE, FALSE), 3, TRUE, IsComb = FALSE))
    expect_true(comboClassTest(letters[1:4], 5, freqs1 = 1:4, IsComb = FALSE))

    set.seed(103)
    myNums = rnorm(5)
    expect_true(comboClassTest(myNums, 4, freqs1 = c(1, 2, 3, 2, 1), IsComb = FALSE))

    ## With FUN
    expect_true(comboClassTest(5, 3, FUN1 = sum))
    expect_true(comboClassTest(as.complex(1:5), 3, TRUE, FUN1 = mean))
    expect_true(comboClassTest(as.raw(1:5), 3, FUN1 = rawToChar))

    set.seed(103)
    myNums = rnorm(5)
    expect_true(comboClassTest(myNums, 4, freqs1 = c(1, 2, 3, 2, 1),
                               FUN1 = cumprod))
    expect_true(comboClassTest(5, 3, IsComb = FALSE, FUN1 = sd))
    expect_true(comboClassTest(as.complex(1:5), 3, TRUE,
                               IsComb = FALSE, FUN1 = var))
    expect_true(comboClassTest(letters[1:4], 5, freqs1 = 1:4, IsComb = FALSE,
                               FUN1 = function(x) {paste0(x, collapse = "")}))
    expect_true(comboClassTest(letters[1:4], 5, freqs1 = 1:4, IsComb = TRUE,
                               FUN1 = function(x) {paste0(x, collapse = "")}))
    expect_true(comboClassTest(letters[1:3], 5, TRUE, IsComb = FALSE,
                               FUN1 = function(x) {charToRaw(paste0(x, collapse = ""))},
                               FUN.VALUE1 = charToRaw("aaaaa")))
    expect_true(comboClassTest(as.character(1:3), 5, TRUE,
                               FUN1 = function(x) {sum(as.integer(x))},
                               FUN.VALUE1 = 2L))
    expect_true(comboClassTest(letters[1:5], 5, IsComb = FALSE,
                               FUN1 = function(x) {which(x == "b")}))

    ## With FUN and FUN.VALUE
    expect_true(comboClassTest(5, 3, FUN1 = sum, FUN.VALUE1 = 1L))
    expect_true(comboClassTest(as.complex(1:5), 3, TRUE, FUN1 = mean,
                               FUN.VALUE1 = as.complex(1.1)))
    expect_true(comboClassTest(as.complex(1:3 + 1i), 5, TRUE,
                               FUN1 = function(x) {sum(x / (as.complex(1:5 + 1i)))},
                               FUN.VALUE1 = as.complex(1), IsComb = FALSE))

    set.seed(103)
    myNums = rnorm(5)
    expect_true(comboClassTest(myNums, 4, freqs1 = c(1, 2, 3, 2, 1),
                               FUN1 = cumprod, FUN.VALUE1 = rnorm(4)))

    expect_true(comboClassTest(5, 3, IsComb = FALSE, FUN1 = sd,
                               FUN.VALUE1 = 1.1))
    expect_true(comboClassTest(as.complex(1:5), 3, TRUE, IsComb = FALSE,
                               FUN1 = var, FUN.VALUE1 = 1.1))
    expect_true(comboClassTest(letters[1:3], freqs1 = 3:1, IsComb = FALSE,
                               FUN1 = function(x) {paste0(x, collapse = "")},
                               FUN.VALUE1 = "a"))
    expect_true(comboClassTest(letters[1:4], 5, freqs1 = 1:4, IsComb = FALSE,
                               FUN1 = function(x) {paste0(x, collapse = "")},
                               FUN.VALUE1 = "a"))
    expect_true(comboClassTest(letters[1:3], 5, TRUE, IsComb = FALSE,
                               FUN1 = function(x) {charToRaw(paste0(x, collapse = ""))},
                               FUN.VALUE1 = charToRaw("aaaaa")))
    expect_true(comboClassTest(as.character(1:3), 5, TRUE,
                               FUN1 = function(x) {sum(as.integer(x))},
                               FUN.VALUE1 = 2L))
    expect_true(comboClassTest(as.character(1:3), 5, TRUE,
                               FUN1 = function(x) {mean(as.numeric(x) + 0.1234)},
                               FUN.VALUE1 = 2.0))
    expect_true(comboClassTest(as.character(1:5), 5, freqs1 = c(2:4, 2:3),
                               FUN1 = function(x) {sum(as.numeric(x) + 1i)},
                               FUN.VALUE1 = as.complex(2.0)))
    expect_true(comboClassTest(letters[1:5], 5, IsComb = FALSE,
                               FUN1 = function(x) {list(x)},
                               FUN.VALUE1 = list(letters[1:5])))
    expect_true(comboClassTest(letters[1:5], 3, IsComb = TRUE,
                               FUN1 = function(x) {x == "a"},
                               FUN.VALUE1 = rep(TRUE, 3)))

    expect_true(comboClassTest(4, 5, freqs1 = 1:4, IsComb = FALSE,
                               FUN1 = function(x) {paste0(letters[x], collapse = "")},
                               FUN.VALUE1 = "a"))
    expect_true(comboClassTest(3, 5, TRUE, IsComb = FALSE,
                               FUN1 = function(x) {
                                   charToRaw(paste0(letters[x], collapse = ""))
                               }, FUN.VALUE1 = charToRaw("aaaaa")))
    expect_true(comboClassTest(5, 3, IsComb = TRUE,
                               FUN1 = function(x) {letters[x] == "a"},
                               FUN.VALUE1 = rep(TRUE, 3)))
    expect_true(comboClassTest(5, 5, IsComb = FALSE,
                               FUN1 = function(x) {list(letters[x])},
                               FUN.VALUE1 = list(letters[1:5])))

    ## With constraintFun
    expect_true(comboClassTest(myNums, 4, freqs1 = c(1, 2, 3, 2, 1),
                               constr1 = "mean"))
    expect_true(comboClassTest(myNums, 4, freqs1 = c(1, 2, 3, 2, 1),
                               constr1 = "mean", IsComb = FALSE))
    expect_true(comboClassTest(6, 4, rep1 = TRUE,
                               constr1 = "sum"))
    expect_true(comboClassTest(6, 4, rep1 = TRUE,
                               constr1 = "sum", IsComb = FALSE))
    expect_true(comboClassTest(10, 4, constr1 = "prod"))
    expect_true(comboClassTest(10, 4, constr1 = "prod", IsComb = FALSE))
    expect_true(comboClassTest(myNums, 3, constr1 = "max"))
    expect_true(comboClassTest(myNums, constr1 = "max", IsComb = FALSE))
    expect_true(comboClassTest(sample(5), 3, constr1 = "min"))
    expect_true(comboClassTest(sample(5), constr1 = "min", IsComb = FALSE))

    ##******** BIG TESTS *********##
    comboClassBigZTest <- function(
        v1, m1 = NULL, rep1 = FALSE, freqs1 = NULL, constr1 = NULL,
        compar1 = NULL, limit1 = NULL, FUN1 = NULL, tol1 = NULL,
        IsComb = TRUE, lenCheck = 500
    ) {

        myResults <- vector(mode = "logical")

        myRows <- if (IsComb && class(v1) == "table") {
            comboCount(v1, m1)
        } else if (IsComb) {
            comboCount(v1, m1, rep1, freqs1)
        } else if (class(v1) == "table") {
            permuteCount(v1, m1)
        } else {
            permuteCount(v1, m1, rep1, freqs1)
        }

        a <- if (IsComb && class(v1) == "table") {
            comboIter(v1, m1, constr1, compar1, limit1, NULL,
                      FUN1, FALSE, NULL, tol1, NULL)
        } else if (IsComb && is.numeric(v1)) {
            comboIter(v1, m1, rep1, freqs1, constr1, compar1, limit1,
                      NULL, FUN1, FALSE, NULL, tol1, NULL)
        } else if (IsComb && is.list(v1)) {
            comboIter(v1, m1, rep1, freqs1)
        } else if (IsComb) {
            comboIter(v1, m1, rep1, freqs1, FUN1, NULL)
        } else if (class(v1) == "table") {
            permuteIter(v1, m1, constr1, compar1, limit1, NULL,
                        FUN1, FALSE, NULL, tol1, NULL)
        } else if (is.numeric(v1)) {
            permuteIter(v1, m1, rep1, freqs1, constr1, compar1, limit1,
                        NULL, FUN1, FALSE, NULL, tol1, NULL)
        } else if (is.list(v1)) {
            permuteIter(v1, m1, rep1, freqs1)
        } else {
            permuteIter(v1, m1, rep1, freqs1, FUN1, NULL)
        }

        b1 <- if (IsComb && class(v1) == "table") {
            comboGeneral(v1, m1, NULL, lenCheck, constr1, compar1, limit1,
                         NULL, FUN1, FALSE, NULL, tol1, NULL)
        } else if (IsComb && is.numeric(v1)) {
            comboGeneral(v1, m1, rep1, freqs1, NULL, lenCheck, constr1, compar1,
                         limit1, NULL, FUN1, FALSE, NULL, tol1, NULL)
        } else if (IsComb && is.list(v1)) {
            comboGeneral(v1, m1, rep1, freqs1, NULL, lenCheck)
        } else if (IsComb) {
            comboGeneral(v1, m1, rep1, freqs1, NULL, lenCheck, FUN1, NULL)
        } else if (class(v1) == "table") {
            permuteGeneral(v1, m1, NULL, lenCheck, constr1, compar1,
                           limit1, NULL, FUN1, FALSE, NULL, tol1, NULL)
        } else if (is.numeric(v1)) {
            permuteGeneral(v1, m1, rep1, freqs1, NULL, lenCheck,
                           constr1, compar1, limit1, NULL, FUN1,
                           FALSE, NULL, tol1, NULL)
        } else if (is.list(v1)) {
            permuteGeneral(v1, m1, rep1, freqs1, NULL, lenCheck)
        } else {
            permuteGeneral(v1, m1, rep1, freqs1, NULL,
                           lenCheck, FUN1, NULL)
        }

        b2 <- if (IsComb && class(v1) == "table") {
            comboGeneral(v1, m1, gmp::sub.bigz(myRows, lenCheck - 1), NULL,
                         constr1, compar1, limit1, NULL, FUN1, FALSE, NULL,
                         tol1, NULL)
        } else if (IsComb && is.numeric(v1)) {
            comboGeneral(v1, m1, rep1, freqs1,
                         gmp::sub.bigz(myRows, lenCheck - 1), NULL, constr1,
                         compar1, limit1, NULL, FUN1, FALSE, NULL,
                         tol1, NULL)
        } else if (IsComb && is.list(v1)) {
            comboGeneral(v1, m1, rep1, freqs1,
                         gmp::sub.bigz(myRows, lenCheck - 1), NULL)
        } else if (IsComb) {
            comboGeneral(v1, m1, rep1, freqs1,
                         gmp::sub.bigz(myRows, lenCheck - 1),
                         NULL, FUN1, NULL)
        } else if (class(v1) == "table") {
            permuteGeneral(v1, m1,
                           gmp::sub.bigz(myRows, lenCheck - 1), NULL, constr1,
                           compar1, limit1, NULL, FUN1, FALSE, NULL,
                           tol1, NULL)
        } else if (is.numeric(v1)) {
            permuteGeneral(v1, m1, rep1, freqs1,
                           gmp::sub.bigz(myRows, lenCheck - 1), NULL, constr1,
                           compar1, limit1, NULL, FUN1, FALSE, NULL,
                           tol1, NULL)
        } else if (is.list(v1)) {
            permuteGeneral(v1, m1, rep1, freqs1,
                           gmp::sub.bigz(myRows, lenCheck - 1), NULL)
        } else {
            permuteGeneral(v1, m1, rep1, freqs1,
                           gmp::sub.bigz(myRows, lenCheck - 1),
                           NULL, FUN1, NULL)
        }

        # .method("summary", &Combo::summary)
        myResults <- c(myResults, isTRUE(all.equal(a@summary()$totalResults, myRows)))

        # .method("sourceVector", &Combo::sourceVector)
        if (length(v1) == 1) {
            myResults <- c(myResults, isTRUE(all.equal(v1, length(a@sourceVector()))))
        } else {
            myResults <- c(myResults, isTRUE(all.equal(v1, a@sourceVector())))
        }

        # .method("front", &Combo::front)
        # .method("currIter", &Combo::currComb)
        # .method("back", &Combo::back)
        # .method("currIter", &Combo::currComb)
        if (is.list(b1)) {
            myResults <- c(myResults, isTRUE(all.equal(a@front(), b1[[1]])))
            myResults <- c(myResults, isTRUE(all.equal(a@currIter(), b1[[1]])))
            myResults <- c(myResults, isTRUE(all.equal(a@back(), b2[[lenCheck]])))
            myResults <- c(myResults, isTRUE(all.equal(a@currIter(), b2[[lenCheck]])))
        } else {
            myResults <- c(myResults, isTRUE(all.equal(a@front(), b1[1, ])))
            myResults <- c(myResults, isTRUE(all.equal(a@currIter(), b1[1, ])))
            myResults <- c(myResults, isTRUE(all.equal(a@back(), b2[lenCheck, ])))
            myResults <- c(myResults, isTRUE(all.equal(a@currIter(), b2[lenCheck, ])))
        }

        # .method("startOver", &Combo::startOver)
        a@startOver()
        a1 <- b1

        if (is.list(b1)) {
            for (i in 1:length(b1)) {
                a1[[i]] <- a@nextIter()
            }
        } else{
            for (i in 1:lenCheck) {
                a1[i, ] <- a@nextIter()
            }
        }

        # .method("nextIter", &Combo::nextComb)
        myResults <- c(myResults, isTRUE(all.equal(a1, b1)))

        # .method("startOver", &Combo::startOver)
        a@startOver()
        numTest <- as.integer(lenCheck / 3);

        s <- 1L
        e <- numTest

        # .method("nextNIter", &Combo::nextNumCombs)
        for (i in 1:3) {
            if (is.list(b1)) {
                myResults <- c(myResults, isTRUE(all.equal(a@nextNIter(numTest), b1[s:e])))
            } else {
                myResults <- c(myResults, isTRUE(all.equal(a@nextNIter(numTest), b1[s:e, ])))
            }

            s <- e + 1L
            e <- e + numTest
        }

        # .method("nextRemaining", &Combo::nextGather)
        a@startOver()
        a[[gmp::sub.bigz(myRows, lenCheck)]]
        myResults <- c(myResults, isTRUE(all.equal(a@nextRemaining(), b2)))

        ## Prepare a for reverse iteration
        temp <- a@back()
        t <- capture.output(a@nextIter())

        a2 <- b2

        if (is.list(b1)) {
            for (i in lenCheck:1) {
                a2[[i]] <- a@prevIter()
            }
        } else {
            for (i in lenCheck:1) {
                a2[i, ] <- a@prevIter()
            }
        }

        # .method("prevIter", &Combo::prevComb)
        myResults <- c(myResults, isTRUE(all.equal(a2, b2)))
        a@startOver()

        s <- lenCheck
        e <- lenCheck - numTest + 1L

        ## Prepare a for reverse iteration
        temp <- a@back()
        t <- capture.output(a@nextIter())

        # .method("prevNIter", &Combo::prevNumCombs)
        for (i in 1:3) {
            if (is.list(b1)) {
                myResults <- c(myResults, isTRUE(all.equal(a@prevNIter(numTest), b2[s:e])))
            } else {
                myResults <- c(myResults, isTRUE(all.equal(a@prevNIter(numTest), b2[s:e, ])))
            }

            s <- e - 1L
            e <- e - numTest
        }

        ## Prepare a for reverse iteration
        temp <- a[[lenCheck]]
        t <- capture.output(a@nextIter())
        a3 <- a@prevRemaining()

        # .method("prevRemaining", &Combo::prevGather)
        if (is.list(b1)) {
            myResults <- c(myResults, all(sapply(1:lenCheck, function(x) {
                isTRUE(all.equal(a3[[lenCheck + 1 - x]], b1[[x]]))
            })))
        } else {
            myResults <- c(myResults, all(sapply(1:lenCheck, function(x) {
                isTRUE(all.equal(a3[lenCheck + 1 - x, ], b1[x, ]))
            })))
        }

        # .method("[[", &Combo::combIndex)
        samp1 <- sample(lenCheck, 5)
        samp2 <- gmp::sub.bigz(myRows, lenCheck) + gmp::as.bigz(samp1)

        if (is.list(b1)) {
            myResults <- c(myResults, isTRUE(all.equal(a[[samp1]], b1[samp1])))
            myResults <- c(myResults, isTRUE(all.equal(a[[samp2]], b2[samp1])))
        } else {
            myResults <- c(myResults, isTRUE(all.equal(a[[samp1]], b1[samp1, ])))
            myResults <- c(myResults, isTRUE(all.equal(a[[samp2]], b2[samp1, ])))
        }

        rm(a, a1, a2, a3, b1, b2)
        gc()

        all(myResults)
    }

    expect_true(comboClassBigZTest(100, 20, lenCheck = 100))
    expect_true(comboClassBigZTest(as.raw(1:100), 10, TRUE, lenCheck = 100))
    expect_true(comboClassBigZTest(factor(1:50), 20, freqs1 = rep(1:10, 5),
                                   lenCheck = 100))

    expect_true(comboClassBigZTest(as.complex(1:50), 10, IsComb = FALSE,
                                   lenCheck = 100))
    expect_true(comboClassBigZTest(state.abb[1:30], 15, freqs1 = rep(1:10, 3),
                                   IsComb = FALSE, lenCheck = 100))
    set.seed(103)
    myNums = rnorm(20)
    expect_true(comboClassBigZTest(myNums, 20, T, IsComb = FALSE,
                                   lenCheck = 100))

    ## With FUN
    expect_true(comboClassBigZTest(100, 20, lenCheck = 100, FUN1 = sum))
    expect_true(comboClassBigZTest((1:100) / 3, 10, TRUE, lenCheck = 100, FUN1 = mean))
    expect_true(comboClassBigZTest(50, 20, freqs1 = rep(1:10, 5),
                                   lenCheck = 100, FUN1 = cumprod))

    expect_true(comboClassBigZTest(as.complex(1:50), 10, IsComb = FALSE,
                                   lenCheck = 100, FUN1 = sd))
    expect_true(comboClassBigZTest(state.abb[1:30], 15,
                                   freqs1 = rep(1:10, 3),
                                   IsComb = FALSE, lenCheck = 100,
                                   FUN1 = function(x) {paste0(x, collapse = "")}))
    set.seed(103)
    myNums = rnorm(20)
    expect_true(comboClassBigZTest(myNums, 20, T, IsComb = FALSE,
                                   lenCheck = 100, FUN1 = cumsum))

    ## With constraintFun
    expect_true(comboClassBigZTest(myNums, 20, T, IsComb = FALSE,
                                   lenCheck = 100, constr1 = "mean"))
    expect_true(comboClassBigZTest(myNums, 20, T, lenCheck = 100,
                                   constr1 = "mean"))
    expect_true(comboClassBigZTest(50, 20, freqs1 = rep(1:10, 5), lenCheck = 100,
                                   constr1 = "sum"))
    expect_true(comboClassBigZTest(50, freqs1 = rep(1:10, 5), IsComb = FALSE,
                                   lenCheck = 100, constr1 = "sum"))

    set.seed(42)
    myNums = rnorm(100)
    expect_true(comboClassBigZTest(myNums, 20, lenCheck = 100, constr1 = "prod"))
    expect_true(comboClassBigZTest(myNums, 20, lenCheck = 100,
                                   constr1 = "prod", IsComb = FALSE))
    expect_true(comboClassBigZTest(sample(100), 30, lenCheck = 100,
                                   constr1 = "max"))
    expect_true(comboClassBigZTest(sample(100), 30, lenCheck = 100,
                                   constr1 = "max", IsComb = FALSE))
    expect_true(comboClassBigZTest(myNums, 30, TRUE,
                                   lenCheck = 100, constr1 = "min"))
    expect_true(comboClassBigZTest(myNums, 30, TRUE,
                                   lenCheck = 100, constr1 = "min", IsComb = FALSE))
})
