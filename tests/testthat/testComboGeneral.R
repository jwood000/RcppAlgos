context("testing comboGeneral")

test_that("comboGeneral produces correct results with no constraints", {
    
    expect_equal(comboGeneral(5, 3), t(combn(5, 3)))
    expect_equal(comboGeneral(factor(1:5, ordered = TRUE), 3), 
                 t(combn(factor(1:5, ordered = TRUE), 3)))
    
    expect_equal(comboGeneral(factor(1:5, ordered = TRUE), 5, freqs = rep(3, 5)), 
                 comboSample(factor(1:5, ordered = TRUE), 5, freqs = rep(3, 5), 
                             sampleVec = 1:comboCount(5, 5, freqs = rep(3, 5))))
    
    expect_equal(nrow(comboGeneral(6, 3)), choose(6, 3))
    expect_equal(as.vector(comboGeneral(1,1)), 1)
    expect_equal(as.vector(comboGeneral(1,1,TRUE)), 1)
    
    expect_equal(comboGeneral(15, 8)[500:600, ], comboGeneral(15, 8,
                                                              lower = 500,
                                                              upper = 600))
    
    expect_equal(comboGeneral(5, 5), 
                 comboGeneral(5, 5, freqs = rep(1, 5)))
    
    set.seed(103)
    myNums <- rnorm(5)
    expect_equal(comboGeneral(myNums, 3), t(combn(myNums, 3)))
    
    myNums2 <- 1:15 / 3
    expect_equal(comboGeneral(myNums2, 8, freqs = rep(2, 15))[150000:157000, ],
                 comboGeneral(myNums2, 8, freqs = rep(2, 15), lower = 150000, upper = 157000))
    
    expect_equal(comboGeneral(letters[1:15], 8, TRUE)[319000:comboCount(15,8,T), ], 
                 comboGeneral(letters[1:15], 8, TRUE, lower = 319000))
    
    expect_equal(nrow(comboGeneral(myNums, 3, TRUE)), 35)
    expect_equal(comboGeneral(myNums, 3, freqs = rep(3, 5)), 
                 comboGeneral(myNums, 3, TRUE))
    
    expect_equal(comboGeneral(LETTERS[1:5], 3, freqs = rep(3, 5)), 
                 comboGeneral(LETTERS[1:5], 3, TRUE))
    
    expect_equal(comboGeneral(LETTERS[1:5], 3), t(combn(LETTERS[1:5], 3)))
    
    expect_equal(sum(comboGeneral(3, 3, freqs = c(1, 1, 1))),
                 sum(comboGeneral(3, 3)))
    
    expect_equal(as.vector(comboGeneral(1, 2, freqs = 2)), c(1, 1))
    
    expect_equal(ncol(comboGeneral(5, 3)), 3)
    expect_equal(ncol(comboGeneral(5, 3, TRUE)), 3)
    expect_equal(ncol(comboGeneral(5, 3, FALSE, constraintFun = "prod", keepResults = TRUE)), 4)
    expect_equal(ncol(comboGeneral(5, 3, TRUE, constraintFun = "prod", keepResults = TRUE)), 4)
    
    expect_equal(ncol(comboGeneral(5, 3, FALSE,
                                    constraintFun = "prod", freqs = c(1,2,1,2,4), 
                                     keepResults = TRUE)), 4)
    
    expect_equal(nrow(comboGeneral(10, 3, TRUE, upper = 20)), 20)
    expect_equal(nrow(comboGeneral(10, 3, upper = 10)), 10)
    
    expect_equal(nrow(comboGeneral(1:10 + .01, 3, FALSE, constraintFun = "prod", 
                                     keepResults = TRUE, upper = 10)), 10)
    
    expect_equal(nrow(comboGeneral(5, 5, freqs = 1:5, upper = 10)), 10)
    
    expect_equal(nrow(comboGeneral(10, 3, TRUE, constraintFun = "prod", 
                                     keepResults = TRUE, upper = 10)), 10)
})

test_that("comboGeneral produces correct results with constraints", {
    
    expect_equal(nrow(comboGeneral(-5:5, 4, FALSE, constraintFun = "sum", 
                                   comparisonFun = "==", limitConstraints = 6)), 
                 length(which(apply(combn(-5:5, 4), 2, sum)==6)))
    
    expect_equal(unique(comboGeneral(5, 5, TRUE,
                                       constraintFun = "sum", comparisonFun = "==", 
                                     limitConstraints = 9, 
                                       keepResults = TRUE)[,6]), 9)
    
    expect_true(all(comboGeneral(5, 5, TRUE,
                                 constraintFun = "min", comparisonFun = "<", 
                                 limitConstraints = 3, 
                                   keepResults = TRUE)[,6] < 3))
    expect_equal(as.vector(comboGeneral(5, 5, freqs = 1:5, 
                                        constraintFun = "sum", 
                                        comparisonFun = "==", 
                                        limitConstraints = 25)), rep(5, 5))
    
    expect_equal(as.vector(comboGeneral(5, 5, constraintFun = "sum", 
                                        comparisonFun = "==", 
                                        limitConstraints = 15)), 1:5)
    
    expect_true(all(comboGeneral(5, 5, TRUE,
                                 constraintFun = "prod", comparisonFun = ">", 
                                 limitConstraints = 100, 
                                   keepResults = TRUE)[,6] > 100))
    
    expect_true(all(comboGeneral(5, 3, FALSE,
                                 constraintFun = "max", comparisonFun = "=<",
                                 limitConstraints = 4, 
                                   keepResults = TRUE)[,4] <= 4))
    
    ## N.B. When there are two comparisons (i.e. comparisonFun = c(">=","<"))
    ## and only one limitConstraint, the first comparison is used. Similarly,
    ## when there are two limitConstraints and one comparison the first
    ## limitConstraint is use and the other is ignored (See next test)
    expect_true(all(comboGeneral(3, 5, TRUE,
                                 constraintFun = "mean", comparisonFun = c(">=","<"),
                                 limitConstraints = 2L, 
                                   keepResults = TRUE)[,6] >= 2))
    
    expect_true(all(comboGeneral(3, 5, TRUE,
                                 constraintFun = "mean", comparisonFun = ">=",
                                 limitConstraints = c(2L, 1e10), 
                                 keepResults = TRUE)[,6] >= 2))
    
    expect_true(all(comboGeneral(5, 5, FALSE, constraintFun = "sum", comparisonFun = ">",
                                 limitConstraints = 18,
                                   freqs = c(1,2,1,2,4), 
                                   keepResults = TRUE)[,6] > 18))
    
    expect_equal(comboGeneral(5, 5, TRUE, constraintFun = "sum", 
                              comparisonFun = "==", 
                              limitConstraints = 25, keepResults = TRUE)[,6], 25)
})

test_that("constraint functions with no comparison produce correct reults", {
    set.seed(12345)
    myNums <- rnorm(15)
    expect_equal(comboGeneral(myNums, 8, constraintFun = "sum", keepResults = TRUE)[,9], 
                 rowSums(comboGeneral(myNums, 8)))
    
    expect_equal(comboGeneral(myNums, 8, TRUE, constraintFun = "sum", keepResults = TRUE)[,9], 
                 rowSums(comboGeneral(myNums, 8, TRUE)))
    
    expect_equal(comboGeneral(myNums, 8, freqs = rep(2, 15), 
                              constraintFun = "sum", keepResults = TRUE)[,9], 
                 rowSums(comboGeneral(myNums, 8, freqs = rep(2, 15))))
    
    ## test product
    expect_equal(comboGeneral(myNums, 8, constraintFun = "prod", keepResults = TRUE)[,9], 
                 apply(comboGeneral(myNums, 8), 1, prod))
    
    expect_equal(comboGeneral(myNums, 8, TRUE, constraintFun = "prod", keepResults = TRUE)[,9], 
                 apply(comboGeneral(myNums, 8, TRUE), 1, prod))
    
    expect_equal(comboGeneral(myNums, 8, freqs = rep(2, 15), 
                              constraintFun = "prod", keepResults = TRUE)[,9], 
                 apply(comboGeneral(myNums, 8, freqs = rep(2, 15)), 1, prod))
    
    ## test mean
    expect_equal(comboGeneral(myNums, 8, constraintFun = "mean", keepResults = TRUE)[,9], 
                 rowMeans(comboGeneral(myNums, 8)))
    
    expect_equal(comboGeneral(myNums, 8, TRUE, constraintFun = "mean", keepResults = TRUE)[,9], 
                 rowMeans(comboGeneral(myNums, 8, TRUE)))
    
    expect_equal(comboGeneral(myNums, 8, freqs = rep(2, 15), 
                              constraintFun = "mean", keepResults = TRUE)[,9], 
                 rowMeans(comboGeneral(myNums, 8, freqs = rep(2, 15))))
    
    ## test max
    expect_equal(comboGeneral(myNums, 8, constraintFun = "max", keepResults = TRUE)[,9], 
                 apply(comboGeneral(myNums, 8), 1, max))
    
    expect_equal(comboGeneral(myNums, 8, TRUE, constraintFun = "max", keepResults = TRUE)[,9], 
                 apply(comboGeneral(myNums, 8, TRUE), 1, max))
    
    expect_equal(comboGeneral(myNums, 8, freqs = rep(2, 15), 
                              constraintFun = "max", keepResults = TRUE)[,9], 
                 apply(comboGeneral(myNums, 8, freqs = rep(2, 15)), 1, max))
    
    ## test min
    expect_equal(comboGeneral(myNums, 8, constraintFun = "min", keepResults = TRUE)[,9], 
                 apply(comboGeneral(myNums, 8), 1, min))
    
    expect_equal(comboGeneral(myNums, 8, TRUE, constraintFun = "min", keepResults = TRUE)[,9], 
                 apply(comboGeneral(myNums, 8, TRUE), 1, min))
    
    expect_equal(comboGeneral(myNums, 8, freqs = rep(2, 15), 
                              constraintFun = "min", keepResults = TRUE)[,9], 
                 apply(comboGeneral(myNums, 8, freqs = rep(2, 15)), 1, min))
})

test_that("comboGeneral produces correct results with exotic constraints", {
    
    a <- t(combn(10, 7))
    expect_equal(comboGeneral(10, 7, constraintFun = "sum",
                 comparisonFun = c(">","<"),
                 limitConstraints = c(40, 45)), a[which(rowSums(a) > 40 & rowSums(a) < 45), ])
    
    a <- comboGeneral(10, 7, TRUE)
    b <- rowSums(a)
    expect_equal(comboGeneral(10, 7, TRUE, constraintFun = "sum",
                              comparisonFun = c(">=","<"),
                              limitConstraints = c(43, 45)), a[which(b >= 43 & b < 45), ])
    
    a <- comboGeneral(10, 7, freqs = rep(3, 10))
    b <- rowSums(a)
    expect_equal(comboGeneral(10, 7, freqs = rep(3, 10), constraintFun = "sum",
                 comparisonFun = c("<=", ">"),
                 limitConstraints = c(50, 47)), a[which(b > 47 & b <= 50), ])
    
    b <- apply(a, 1, max)
    expect_equal(comboGeneral(10, 7, freqs = rep(3, 10), 
                 constraintFun = "max",
                 comparisonFun = c("<=", ">"),
                 limitConstraints = c(9, 7)), a[which(b > 7 & b <= 9), ])
    
    b <- apply(a, 1, min)
    expect_equal(comboGeneral(10, 7, freqs = rep(3, 10), 
                              constraintFun = "min",
                              comparisonFun = "==",
                              limitConstraints = 3, 
                              lower = 7900, upper = 8500),
                 a[(7900:8500)[b[7900:8500] == 3], ])
    
    a <- comboGeneral(5, 7, T)
    b <- apply(a, 1, prod)
    expect_equal(comboGeneral(5, 7, TRUE, constraintFun = "prod",
                 comparisonFun = c(">=","<="),
                 limitConstraints = c(2000, 5000)), a[which(b >= 2000 & b <= 5000), ])
    
    a <- comboGeneral(-5, 7, T)
    b <- apply(a, 1, prod)
    expect_equal(nrow(comboGeneral(-5, 7, TRUE, constraintFun = "prod",
                              comparisonFun = c("<=",">="),
                              limitConstraints = c(-2000, 5000), 
                              keepResults = TRUE)), 
                 nrow(rbind(a[which(b <= -2000),], a[which(b >= 5000), ])))
    
    ## Testing sums in a range
    a <- comboGeneral(10, 8, TRUE, lower = 23500, upper = 24000, 
                      constraintFun = "sum", keepResults = TRUE)
    
    expect_equal(comboGeneral(c(NA, 1:10), 8, TRUE, constraintFun = "sum", 
                              comparisonFun = c("=>","=<"), 
                              limitConstraints = c(72, 78), 
                              lower = 23500, upper = 24000), 
                 a[a[,9] >= 72 & a[,9] <= 78, 1:8])
    
    
    comp1 <- c("<", "<=")
    comp2 <- c(">", ">=")
    allCombs1 <- comboGeneral(10, 8, freqs = c(rep(1:3, 3), 3), 
                             constraintFun = "sum", keepResults = TRUE)
    allCombs2 <- comboGeneral(10:1, 8, freqs = rev(c(rep(1:3, 3), 3)), 
                              constraintFun = "sum", keepResults = TRUE)
    theSums1 <- allCombs1[, 9]
    theSums2 <- allCombs2[, 9]
    allCombs1 <- allCombs1[, 1:8]
    allCombs2 <- allCombs2[, 1:8]
    
    for (i in 1:2) {
        
        if (i == 1) {
            a <- comp1
            b <- comp2
        } else {
            a <- comp2
            b <- comp1
        }
        
        for (j in a) {
            for (k in b) {
                myComp <- c(j, k)
                myTest <- comboGeneral(10, 8, freqs = c(rep(1:3, 3), 3),
                                       constraintFun = "sum", comparisonFun = myComp,
                                       limitConstraints = c(42, 53))
                fun1 <- match.fun(j)
                fun2 <- match.fun(k)
                
                if (i == 1) {
                    temp1 <- allCombs1[fun1(theSums1, 42), ]
                    temp2 <- allCombs2[fun2(theSums2, 53), ]
                    temp <- rbind(temp1, temp2)
                } else {
                    temp <- allCombs1[fun1(theSums1, 42) & fun2(theSums1, 53),]
                }
               
                expect_equal(temp, myTest)
            }
        }
    }
    
    allCombs1 <- comboGeneral(10, 8, TRUE, 
                              constraintFun = "sum", keepResults = TRUE)
    allCombs2 <- comboGeneral(10:1, 8, TRUE, 
                              constraintFun = "sum", keepResults = TRUE)
    theSums1 <- allCombs1[, 9]
    theSums2 <- allCombs2[, 9]
    allCombs1 <- allCombs1[, 1:8]
    allCombs2 <- allCombs2[, 1:8]
    
    for (i in 1:2) {
        
        if (i == 1) {
            a <- comp1
            b <- comp2
        } else {
            a <- comp2
            b <- comp1
        }
        
        for (j in a) {
            for (k in b) {
                myComp <- c(j, k)
                myTest <- comboGeneral(10, 8, TRUE,
                                       constraintFun = "sum", comparisonFun = myComp,
                                       limitConstraints = c(42, 53))
                fun1 <- match.fun(j)
                fun2 <- match.fun(k)
                
                if (i == 1) {
                    temp1 <- allCombs1[fun1(theSums1, 42), ]
                    temp2 <- allCombs2[fun2(theSums2, 53), ]
                    temp <- rbind(temp1, temp2)
                } else {
                    temp <- allCombs1[fun1(theSums1, 42) & fun2(theSums1, 53),]
                }
                
                expect_equal(temp, myTest)
            }
        }
    }
    
    allCombs1 <- comboGeneral(18, 8, 
                              constraintFun = "sum", keepResults = TRUE)
    allCombs2 <- comboGeneral(18:1, 8, 
                              constraintFun = "sum", keepResults = TRUE)
    theSums1 <- allCombs1[, 9]
    theSums2 <- allCombs2[, 9]
    allCombs1 <- allCombs1[, 1:8]
    allCombs2 <- allCombs2[, 1:8]
    
    for (i in 1:2) {
        
        if (i == 1) {
            a <- comp1
            b <- comp2
        } else {
            a <- comp2
            b <- comp1
        }
        
        for (j in a) {
            for (k in b) {
                myComp <- c(j, k)
                myTest <- comboGeneral(18, 8,
                                       constraintFun = "sum", comparisonFun = myComp,
                                       limitConstraints = c(68, 84))
                fun1 <- match.fun(j)
                fun2 <- match.fun(k)
                
                if (i == 1) {
                    temp1 <- allCombs1[fun1(theSums1, 68), ]
                    temp2 <- allCombs2[fun2(theSums2, 84), ]
                    temp <- rbind(temp1, temp2)
                } else {
                    temp <- allCombs1[fun1(theSums1, 68) & fun2(theSums1, 84),]
                }
                expect_equal(temp, myTest)
            }
        }
    }
})

test_that("comboGeneral produces appropriate error messages", {
    expect_error(comboGeneral(9,4,TRUE,NULL,NULL,NULL,"summ","<",10), "prod, sum, mean, max, or min")
    expect_error(comboGeneral(9,4,TRUE,NULL,NULL,NULL,"sum","=<>",10), ">, >=, <, <=, or ==")
    expect_error(comboGeneral(9,4,TRUE,NULL,NULL,NULL,"sum",60,10), "must be passed as a character")
    expect_error(comboGeneral(9,4,FALSE,NULL,NULL,NULL,sum,"<",10), "must be passed as a character")
    expect_error(comboGeneral(9,4,TRUE,NULL,NULL,-1,"sum","<",10), "must be positive")
    expect_error(comboGeneral(170,7,FALSE,NULL,NULL,NULL,"sum","<",100), "The number of rows cannot exceed")
    expect_error(comboGeneral(170,7,FALSE,NULL,NULL,10^10,"sum","<",100), "number of rows cannot exceed")
    
    expect_error(comboGeneral(50, 5, lower = "1e6"), "bounds must be of type numeric or integer")
    expect_error(comboGeneral(50, 5, upper = "1e6"), "bounds must be of type numeric or integer")
    expect_error(comboGeneral(50, 5, lower = -100), "bounds must be positive")
    expect_error(comboGeneral(50, 5, upper = -100), "bounds must be positive")
    
    expect_error(comboGeneral(5, 50), "m must be less than or equal to the length of v")
    
    expect_error(comboGeneral(5, 3, upper = 100), "bound exceeds the maximum number of possible results")
    expect_error(comboGeneral(5, 3, lower = 100), "bound exceeds the maximum number of possible results")
    
    expect_error(comboGeneral(50, 3, lower = 100, upper = 10), "The number of rows must be positive")
    
    expect_error(comboGeneral(10,7,FALSE,NULL,NULL,NULL,"sum","<","100"), 
                 "limitConstraints must be of type numeric or integer")
    
    expect_error(comboGeneral(10,7,FALSE,NULL,NULL,NULL,"sum","<",c(20,30,40)), 
                 "there cannot be more than 2 limitConstraints")
    expect_error(comboGeneral(10,7,FALSE,NULL,NULL,NULL,"sum","<",c(20,NA)), 
                 "limitConstraints cannot be NA")
    
    expect_error(comboGeneral(10,7,FALSE,NULL,NULL,NULL,"sum",c("<",">","=="),c(20,30)), 
                 "there cannot be more than 2 comparison operator")
    
    expect_error(comboGeneral(10,7,FALSE,NULL,NULL,NULL,"sum",c("<","=="),c(20,30)), 
                 "If comparing against two limitConstraints, the")
    
    expect_error(comboGeneral(10,7,FALSE,NULL,NULL,NULL,"sum",c("<","=<"),c(20,30)), 
                 "Cannot have two 'less than' comparisonFuns or two")
    
    expect_error(comboGeneral(10,7,FALSE,NULL,NULL,NULL,"sum",c("<=","<"),c(20,20)), 
                 "The limitConstraints must be different")
    
    expect_error(comboGeneral(-100:100, 25, constraintFun = "prod",
                              comparisonFun = "<=", limitConstraints = 10, upper = 100), 
                 "The number of rows cannot exceed")
    
    expect_error(comboGeneral(5, 3, TRUE, constraintFun = "product", keepResults = TRUE), 
                 "contraintFun must be one of the following:")
    expect_error(comboGeneral((1:5)+.01, 3, TRUE, constraintFun = "product", keepResults = TRUE), 
                 "contraintFun must be one of the following:")
    expect_error(comboGeneral(5,3,freqs = c(1,2,3,"2",1)), "freqs must be of type numeric")
    expect_error(comboGeneral(5,3,freqs = c(1,2,3,-2,1)), "in freqs must be a positive")
    expect_error(comboGeneral(5,1000,freqs = rep(5000, 5)), "number of rows cannot exceed")
    expect_error(comboGeneral(5,freqs = rep(1,6)), "the length of freqs must equal the")
})
