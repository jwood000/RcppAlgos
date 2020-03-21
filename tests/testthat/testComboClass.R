context("testing comboIter & permuteIter")

test_that("comboIter & permuteIter produces correct results with no constraints", {
    
    comboClassTest <- function(v1, m1 = NULL, rep1 = FALSE, freqs1 = NULL,
                               constr1 = NULL, compar1 = NULL, limit1 = NULL,
                               FUN1 = NULL, tol1 = NULL, IsComb = TRUE) {
        
        myResults <- vector(mode = "logical")
        
        if (IsComb) {
            myRows <- comboCount(v1, m1, rep1, freqs1)
        } else {
            myRows <- permuteCount(v1, m1, rep1, freqs1)
        }
        
        if (IsComb) {
            a <- comboIter(v1, m1, rep1, freqs1,
                           constr1, compar1, limit1, FUN1, tol1)
            
            b <- comboGeneral(v = v1, m = m1, repetition = rep1, freqs = freqs1,
                              constraintFun = constr1, comparisonFun = compar1, 
                              limitConstraints = limit1, FUN = FUN1, tolerance = tol1)
        } else {
            a <- permuteIter(v1, m1, rep1, freqs1,
                             constr1, compar1, limit1, FUN1, tol1)
            
            b <- permuteGeneral(v = v1, m = m1, repetition = rep1, freqs = freqs1,
                                constraintFun = constr1, comparisonFun = compar1, 
                                limitConstraints = limit1, FUN = FUN1, tolerance = tol1)
        }
        
        # .method("summary", &Combo::summary)
        myResults <- c(myResults, isTRUE(all.equal(a$summary()$totalResults, myRows)))
        
        # .method("sourceVector", &Combo::sourceVector)
        if (length(v1) == 1) {
            myResults <- c(myResults, isTRUE(all.equal(v1, length(a$sourceVector()))))
        } else {
            myResults <- c(myResults, isTRUE(all.equal(v1, a$sourceVector())))
        }
        
        # .method("front", &Combo::front)
        # .method("currIter", &Combo::currComb)
        # .method("back", &Combo::back)
        # .method("currIter", &Combo::currComb)
        if (is.list(b)) {
            myResults <- c(myResults, isTRUE(all.equal(a$front(), b[[1]])))
            myResults <- c(myResults, isTRUE(all.equal(a$currIter(), b[[1]])))
            myResults <- c(myResults, isTRUE(all.equal(a$back(), b[[myRows]])))
            myResults <- c(myResults, isTRUE(all.equal(a$currIter(), b[[myRows]])))
        } else {
            myResults <- c(myResults, isTRUE(all.equal(a$front(), b[1, ])))
            myResults <- c(myResults, isTRUE(all.equal(a$currIter(), b[1, ])))
            myResults <- c(myResults, isTRUE(all.equal(a$back(), b[myRows, ])))
            myResults <- c(myResults, isTRUE(all.equal(a$currIter(), b[myRows, ])))
        }
        
        # .method("startOver", &Combo::startOver)
        a$startOver()
        a1 <- b
        
        if (is.list(b)) {
            for (i in 1:myRows) {
                a1[[i]] <- a$nextIter()
            }
        } else{
            for (i in 1:myRows) {
                a1[i, ] <- a$nextIter()
            }
        }

        # .method("nextIter", &Combo::nextComb)
        myResults <- c(myResults, isTRUE(all.equal(a1, b)))
        
        # .method("startOver", &Combo::startOver)
        a$startOver()
        numTest <- as.integer(myRows / 3);
        
        s <- 1L
        e <- numTest
        
        # .method("nextNIter", &Combo::nextNumCombs)
        for (i in 1:3) {
            if (is.list(b)) {
                myResults <- c(myResults, isTRUE(all.equal(a$nextNIter(numTest), b[s:e])))
            } else {
                myResults <- c(myResults, isTRUE(all.equal(a$nextNIter(numTest), b[s:e, ])))
            }
            
            s <- e + 1L
            e <- e + numTest
        }
        
        # .method("nextRemaining", &Combo::nextGather)
        a$startOver()
        myResults <- c(myResults, isTRUE(all.equal(a$nextRemaining(), b)))
        
        ## Prepare a for reverse iteration
        temp <- a$back()
        t <- capture.output(a$nextIter())
        a2 <- b
        
        if (is.list(b)) {
            for (i in myRows:1) {
                a2[[i]] <- a$prevIter()
            }
        } else {
            for (i in myRows:1) {
                a2[i, ] <- a$prevIter()
            }
        }
        
        # .method("prevIter", &Combo::prevComb)
        myResults <- c(myResults, isTRUE(all.equal(a2, b)))
        a$startOver()
        
        s <- myRows
        e <- myRows - numTest + 1L
        
        ## Prepare a for reverse iteration
        temp <- a$back()
        t <- capture.output(a$nextIter())
        
        # .method("prevNIter", &Combo::prevNumCombs)
        for (i in 1:3) {
            if (is.list(b)) {
                myResults <- c(myResults, isTRUE(all.equal(a$prevNIter(numTest), b[s:e])))
            } else {
                myResults <- c(myResults, isTRUE(all.equal(a$prevNIter(numTest), b[s:e, ])))
            }
            
            s <- e - 1L
            e <- e - numTest
        }
        
        ## Prepare a for reverse iteration
        temp <- a$back()
        t <- capture.output(a$nextIter())
        a3 <- a$prevRemaining()
        
        # .method("prevRemaining", &Combo::prevGather)
        if (is.list(b)) {
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
        
        if (is.list(b)) {
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
    
    expect_true(comboClassTest(as.complex(1:5), 3, TRUE))
    expect_true(comboClassTest(c(TRUE, FALSE), 20, TRUE))
    expect_true(comboClassTest(letters[1:5], 6, freqs1 = 1:5))
    
    set.seed(103)
    myNums = rnorm(5)
    expect_true(comboClassTest(myNums, 4, freqs1 = c(1, 2, 3, 2, 1)))
    
    expect_true(comboClassTest(5, 3, IsComb = FALSE))
    expect_true(comboClassTest(as.raw(1:5), 3, IsComb = FALSE))
    expect_true(comboClassTest(factor(1:5), 3, IsComb = FALSE))
    
    expect_true(comboClassTest(as.complex(1:5), 3, TRUE, IsComb = FALSE))
    expect_true(comboClassTest(c(TRUE, FALSE), 3, TRUE, IsComb = FALSE))
    expect_true(comboClassTest(letters[1:4], 5, freqs1 = 1:4, IsComb = FALSE))
    
    set.seed(103)
    myNums = rnorm(5)
    expect_true(comboClassTest(myNums, 4, freqs1 = c(1, 2, 3, 2, 1), IsComb = FALSE))
    
    ## With FUN
    expect_true(comboClassTest(5, 3, FUN1 = sum))
    expect_true(comboClassTest(as.complex(1:5), 3, TRUE, FUN1 = mean))

    set.seed(103)
    myNums = rnorm(5)
    expect_true(comboClassTest(myNums, 4, freqs1 = c(1, 2, 3, 2, 1), FUN1 = cumprod))
    
    expect_true(comboClassTest(5, 3, IsComb = FALSE, FUN1 = sd))
    expect_true(comboClassTest(as.complex(1:5), 3, TRUE, IsComb = FALSE, FUN1 = var))
    expect_true(comboClassTest(letters[1:4], 5, freqs1 = 1:4, IsComb = FALSE, 
                               FUN1 = function(x) {paste0(x, collapse = "")}))
    
    
    ##******** BIG TESTS *********##
    comboClassBigZTest <- function(v1, m1 = NULL, rep1 = FALSE, freqs1 = NULL,
                                   constr1 = NULL, compar1 = NULL, limit1 = NULL,
                                   FUN1 = NULL, tol1 = NULL, IsComb = TRUE, lenCheck = 500) {
        
        myResults <- vector(mode = "logical")
        
        if (IsComb) {
            myRows <- comboCount(v1, m1, rep1, freqs1)
        } else {
            myRows <- permuteCount(v1, m1, rep1, freqs1)
        }
        
        if (IsComb) {
            a <- comboIter(v1, m1, rep1, freqs1,
                           constr1, compar1, limit1, FUN1, tol1)
            
            b1 <- comboGeneral(v = v1, m = m1, repetition = rep1, freqs = freqs1,
                               constraintFun = constr1, comparisonFun = compar1, 
                               limitConstraints = limit1, FUN = FUN1, tolerance = tol1,
                               upper = lenCheck)
            
            b2 <- comboGeneral(v = v1, m = m1, repetition = rep1, freqs = freqs1,
                               constraintFun = constr1, comparisonFun = compar1, 
                               limitConstraints = limit1, FUN = FUN1, tolerance = tol1,
                               lower = gmp::sub.bigz(myRows, lenCheck - 1))
        } else {
            a <- permuteIter(v1, m1, rep1, freqs1,
                             constr1, compar1, limit1, FUN1, tol1)
            
            b1 <- permuteGeneral(v = v1, m = m1, repetition = rep1, freqs = freqs1,
                                 constraintFun = constr1, comparisonFun = compar1, 
                                 limitConstraints = limit1, FUN = FUN1, tolerance = tol1,
                                 upper = lenCheck)
            
            b2 <- permuteGeneral(v = v1, m = m1, repetition = rep1, freqs = freqs1,
                                 constraintFun = constr1, comparisonFun = compar1, 
                                 limitConstraints = limit1, FUN = FUN1, tolerance = tol1,
                                 lower = gmp::sub.bigz(myRows, lenCheck - 1))
        }
        
        # .method("summary", &Combo::summary)
        myResults <- c(myResults, isTRUE(all.equal(a$summary()$totalResults, myRows)))
        
        # .method("sourceVector", &Combo::sourceVector)
        if (length(v1) == 1) {
            myResults <- c(myResults, isTRUE(all.equal(v1, length(a$sourceVector()))))
        } else {
            myResults <- c(myResults, isTRUE(all.equal(v1, a$sourceVector())))
        }
        
        # .method("front", &Combo::front)
        # .method("currIter", &Combo::currComb)
        # .method("back", &Combo::back)
        # .method("currIter", &Combo::currComb)
        if (is.list(b1)) {
            myResults <- c(myResults, isTRUE(all.equal(a$front(), b1[[1]])))
            myResults <- c(myResults, isTRUE(all.equal(a$currIter(), b1[[1]])))
            myResults <- c(myResults, isTRUE(all.equal(a$back(), b2[[lenCheck]])))
            myResults <- c(myResults, isTRUE(all.equal(a$currIter(), b2[[lenCheck]])))
        } else {
            myResults <- c(myResults, isTRUE(all.equal(a$front(), b1[1, ])))
            myResults <- c(myResults, isTRUE(all.equal(a$currIter(), b1[1, ])))
            myResults <- c(myResults, isTRUE(all.equal(a$back(), b2[lenCheck, ])))
            myResults <- c(myResults, isTRUE(all.equal(a$currIter(), b2[lenCheck, ])))
        }
        
        # .method("startOver", &Combo::startOver)
        a$startOver()
        a1 <- b1
        
        if (is.list(b1)) {
            for (i in 1:length(b1)) {
                a1[[i]] <- a$nextIter()
            }
        } else{
            for (i in 1:lenCheck) {
                a1[i, ] <- a$nextIter()
            }
        }
        
        # .method("nextIter", &Combo::nextComb)
        myResults <- c(myResults, isTRUE(all.equal(a1, b1)))
        
        # .method("startOver", &Combo::startOver)
        a$startOver()
        numTest <- as.integer(lenCheck / 3);
        
        s <- 1L
        e <- numTest
        
        # .method("nextNIter", &Combo::nextNumCombs)
        for (i in 1:3) {
            if (is.list(b1)) {
                myResults <- c(myResults, isTRUE(all.equal(a$nextNIter(numTest), b1[s:e])))
            } else {
                myResults <- c(myResults, isTRUE(all.equal(a$nextNIter(numTest), b1[s:e, ])))
            }
            
            s <- e + 1L
            e <- e + numTest
        }
        
        # .method("nextRemaining", &Combo::nextGather)
        a$startOver()
        a[[gmp::sub.bigz(myRows, lenCheck)]]
        myResults <- c(myResults, isTRUE(all.equal(a$nextRemaining(), b2)))
        
        ## Prepare a for reverse iteration
        temp <- a$back()
        t <- capture.output(a$nextIter())
        
        a2 <- b2
        
        if (is.list(b1)) {
            for (i in lenCheck:1) {
                a2[[i]] <- a$prevIter()
            }
        } else {
            for (i in lenCheck:1) {
                a2[i, ] <- a$prevIter()
            }
        }
        
        # .method("prevIter", &Combo::prevComb)
        myResults <- c(myResults, isTRUE(all.equal(a2, b2)))
        a$startOver()
        
        s <- lenCheck
        e <- lenCheck - numTest + 1L
        
        ## Prepare a for reverse iteration
        temp <- a$back()
        t <- capture.output(a$nextIter())
        
        # .method("prevNIter", &Combo::prevNumCombs)
        for (i in 1:3) {
            if (is.list(b1)) {
                myResults <- c(myResults, isTRUE(all.equal(a$prevNIter(numTest), b2[s:e])))
            } else {
                myResults <- c(myResults, isTRUE(all.equal(a$prevNIter(numTest), b2[s:e, ])))
            }
            
            s <- e - 1L
            e <- e - numTest
        }
        
        ## Prepare a for reverse iteration
        temp <- a[[lenCheck]]
        t <- capture.output(a$nextIter())
        a3 <- a$prevRemaining()
        
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
    expect_true(comboClassBigZTest(factor(1:50), 20, freqs1 = rep(1:10, 5), lenCheck = 100))
    
    expect_true(comboClassBigZTest(as.complex(1:50), 10, IsComb = FALSE, lenCheck = 100))
    expect_true(comboClassBigZTest(state.abb[1:30], 15, 
                                   freqs1 = rep(1:10, 3), IsComb = FALSE, lenCheck = 100))
    set.seed(103)
    myNums = rnorm(20)
    expect_true(comboClassBigZTest(myNums, 20, T, IsComb = FALSE, lenCheck = 100))
    
    ## With FUN
    expect_true(comboClassBigZTest(100, 20, lenCheck = 100, FUN1 = sum))
    expect_true(comboClassBigZTest((1:100) / 3, 10, TRUE, lenCheck = 100, FUN1 = mean))
    expect_true(comboClassBigZTest(50, 20, freqs1 = rep(1:10, 5),
                                   lenCheck = 100, FUN1 = cumprod))
    
    expect_true(comboClassBigZTest(as.complex(1:50), 10, IsComb = FALSE,
                                   lenCheck = 100, FUN1 = sd))
    expect_true(comboClassBigZTest(state.abb[1:30], 15, 
                                   freqs1 = rep(1:10, 3), IsComb = FALSE, lenCheck = 100, 
                                   FUN1 = function(x) {paste0(x, collapse = "")}))
    set.seed(103)
    myNums = rnorm(20)
    expect_true(comboClassBigZTest(myNums, 20, T, IsComb = FALSE, lenCheck = 100, FUN1 = cumsum))
})
