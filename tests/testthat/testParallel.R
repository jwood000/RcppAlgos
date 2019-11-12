context("testing with Parallel enabled")

test_that("comboGeneral produces correct results with Parallel enabled and no constrainFun", {
    
    set.seed(12)
    ## N.B. Parallel has no effect when number of results is less than 20000
    expect_equal(comboGeneral(18, 9, nThreads = 2, upper = 5000), 
                 comboGeneral(18, 9, upper = 5000))
    
    ## For both of the usages below, only 2 threads will be spawned
    expect_equal(comboGeneral(18, 9, nThreads = 3, upper = 25000), 
                 comboGeneral(18, 9, Parallel = TRUE, upper = 25000))
    
    ######********************** All Results *******************#########
    #### NO Repetition
    ## comboCount(18, 9)
    ## [1] 48620 
    expect_equal(comboGeneral(18, 9, nThreads = 2), comboGeneral(18, 9))
    
    #### With Repetition
    ## comboCount(12, 8, T)
    ## [1] 75582
    numVec <- rnorm(12)
    expect_equal(comboGeneral(numVec, 8, TRUE, nThreads = 2), comboGeneral(numVec, 8, TRUE))
    
    #### Multisets
    ## comboCount(15, 7, freqs = rep(1:5, 3))
    ## [1] 73065
    expect_equal(comboGeneral(factor(1:15), 7, freqs = rep(1:5, 3), nThreads = 2), 
                 comboGeneral(factor(1:15), 7, freqs = rep(1:5, 3)))
    
    ######********************** Upper Only *******************#########
    #### NO Repetition
    expect_equal(comboGeneral(50, 10, nThreads = 2, upper = 30000), 
                 comboGeneral(50, 10, upper = 30000))
    
    ######********************** Lower Only *******************#########
    #### With Repetition
    numVec <- rnorm(25)
    total = comboCount(25, 8, TRUE)
    expect_equal(comboGeneral(numVec, 8, TRUE, nThreads = 2, lower = total - 30000),
                 comboGeneral(numVec, 8, TRUE, lower = total - 30000))
    
    ######********************** Upper & Lower *******************#########
    #### Multisets
    expect_equal(comboGeneral(factor(1:50), 7, freqs = rep(1:5, 10), nThreads = 2,
                              lower = 100000, upper = 130000), 
                 comboGeneral(factor(1:50), 7, freqs = rep(1:5, 10), lower = 100000, upper = 130000))
})

test_that("comboGeneral produces correct results with Parallel enabled and constrainFun", {
    
    set.seed(13)
    ## N.B. Parallel has no effect when number of results is less than 20000
    ## N.B. No longer need keepResult = TRUE if user only declares constraintFun
    
    ######********************** All Results *******************#########
    #### NO Repetition
    ## comboCount(18, 9)
    ## [1] 48620 
    expect_equal(comboGeneral(18, 9, nThreads = 2, constraintFun = "sum"), 
                 comboGeneral(18, 9, constraintFun = "sum", keepResults = TRUE))
    
    #### With Repetition
    ## comboCount(12, 8, T)
    ## [1] 75582
    numVec <- rnorm(12)
    expect_equal(comboGeneral(numVec, 8, TRUE, nThreads = 2, constraintFun = "prod"),
                 comboGeneral(numVec, 8, TRUE, constraintFun = "prod", keepResults = TRUE))
    
    #### Multisets
    ## comboCount(15, 7, freqs = rep(1:5, 3))
    ## [1] 73065
    numVec <- runif(15)
    expect_equal(comboGeneral(numVec, 7, freqs = rep(1:5, 3), nThreads = 2, constraintFun = "mean"), 
                 comboGeneral(numVec, 7, freqs = rep(1:5, 3), constraintFun = "mean", keepResults = TRUE))
    
    
    ######********************** Upper Only *******************#########
    #### NO Repetition
    expect_equal(comboGeneral(50, 10, nThreads = 2, upper = 30000, constraintFun = "max"), 
                 comboGeneral(50, 10, upper = 30000, constraintFun = "max", keepResults = TRUE))
    
    #### With Repetition
    numVec <- rnorm(25)
    expect_equal(comboGeneral(numVec, 8, TRUE, nThreads = 2, upper = 30000, constraintFun = "min"),
                 comboGeneral(numVec, 8, TRUE, upper = 30000, constraintFun = "min", keepResults = TRUE))
    
    #### Multisets
    numVec <- rnorm(50)
    expect_equal(comboGeneral(numVec, 7, freqs = rep(1:5, 10), nThreads = 2, 
                              upper = 30000, constraintFun = "sum"), 
                 comboGeneral(numVec, 7, freqs = rep(1:5, 10), upper = 30000, 
                              constraintFun = "sum", keepResults = TRUE))
    
    ######********************** Lower Only *******************#########
    #### NO Repetition
    total = comboCount(50, 10)
    expect_equal(comboGeneral(50, 10, nThreads = 2, lower = total - 30000, constraintFun = "prod"), 
                 comboGeneral(50, 10, lower = total - 30000, constraintFun = "prod", keepResults = TRUE))
    
    ######********************** Upper & Lower *******************#########
    #### Multisets
    expect_equal(comboGeneral(50, 7, freqs = rep(1:5, 10), nThreads = 2, 
                              lower = 100000, upper = 130000, constraintFun = "mean"), 
                 comboGeneral(50, 7, freqs = rep(1:5, 10), lower = 100000, 
                              upper = 130000, constraintFun = "mean", keepResults = TRUE))
})

test_that("permuteGeneral produces correct results with Parallel enabled and no constrainFun", {
    
    set.seed(14)
    ## N.B. Parallel has no effect when number of results is less than 20000
    
    ######********************** All Results *******************#########
    #### NO Repetition
    ## permuteCount(11, 5)
    ## [1] 55440
    expect_equal(permuteGeneral(11, 5, nThreads = 2), permuteGeneral(11, 5))
    
    #### With Repetition
    ## permuteCount(9, 5, T)
    ## [1] 59049
    numVec <- rnorm(9)
    expect_equal(permuteGeneral(numVec, 5, TRUE, nThreads = 2), permuteGeneral(numVec, 5, TRUE))
    
    #### Multisets
    ## permuteCount(9, 5, freqs = rep(1:3, 3))
    ## [1] 40260
    expect_equal(permuteGeneral(factor(1:9), 5, freqs = rep(1:3, 3), nThreads = 2), 
                 permuteGeneral(factor(1:9), 5, freqs = rep(1:3, 3)))
    
    
    ######********************** Upper Only *******************#########
    #### NO Repetition
    expect_equal(permuteGeneral(30, 10, nThreads = 2, upper = 30000), 
                 permuteGeneral(30, 10, upper = 30000))
    
    expect_equal(permuteGeneral(10, 6, nThreads = 2, upper = 30000), 
                 permuteGeneral(10, 6, upper = 30000))
    
    #### With Repetition
    expect_equal(permuteGeneral(10, 5, T, upper = 39999, nThreads = 2),
                 permuteGeneral(10, 5, T, upper = 39999))
    
    ######********************** Lower Only *******************#########
    #### With Repetition
    numVec <- rnorm(12)
    total = permuteCount(12, 8, TRUE)
    expect_equal(permuteGeneral(numVec, 8, TRUE, nThreads = 2, lower = total - 30000),
                 permuteGeneral(numVec, 8, TRUE, lower = total - 30000))
    
    ######********************** Upper & Lower *******************#########
    #### Multisets
    expect_equal(permuteGeneral(factor(1:30), 7, freqs = rep(1:5, 6), 
                                nThreads = 2, lower = 100000, upper = 130000), 
                 permuteGeneral(factor(1:30), 7, freqs = rep(1:5, 6), lower = 100000, upper = 130000))
})

test_that("permuteGeneral produces correct results with Parallel enabled with logical vector", {
    
    ## N.B. Parallel has no effect when number of results is less than 20000
    
    ######********************** All Results *******************#########
    #### With Repetition
    ## permuteCount(2, 15, T)
    ## [1] 32768
    expect_equal(permuteGeneral(c(TRUE, FALSE), 15, TRUE, nThreads = 2),
                 permuteGeneral(c(TRUE, FALSE), 15, TRUE))
    
    #### Multisets
    ## permuteCount(2, 15, freqs = c(9, 9))
    ## [1] 22880
    expect_equal(permuteGeneral(c(TRUE, FALSE), 15, freqs = c(9, 9), nThreads = 2),
                 permuteGeneral(c(TRUE, FALSE), 15, freqs = c(9, 9)))
    
    
    ######********************** Upper Only *******************#########
    #### With Repetition
    expect_equal(permuteGeneral(c(TRUE, FALSE), 30, TRUE, nThreads = 2, upper = 30000),
                 permuteGeneral(c(TRUE, FALSE), 30, TRUE, upper = 30000))
    
    ######********************** Lower Only *******************#########
    #### Multisets
    total = permuteCount(2, 30, freqs = c(20, 20))
    expect_equal(permuteGeneral(c(TRUE, FALSE), 30, freqs = c(20, 20), nThreads = 2, lower = total - 30000),
                 permuteGeneral(c(TRUE, FALSE), 30, freqs = c(20, 20), lower = total - 30000))
    
    ######********************** Upper & Lower *******************#########
    #### With Repetition
    expect_equal(permuteGeneral(c(TRUE, FALSE), 30, TRUE, nThreads = 2, lower = 100000, upper = 130000),
                 permuteGeneral(c(TRUE, FALSE), 30, TRUE, lower = 100000, upper = 130000))
})

test_that("permuteGeneral produces correct results with Parallel enabled and constrainFun", {
    
    set.seed(15)
    ## N.B. Parallel has no effect when number of results is less than 20000
    ## N.B. No longer need keepResult = TRUE if user only declares constraintFun
    
    ######********************** All Results *******************#########
    #### NO Repetition
    expect_equal(permuteGeneral(11, 5, nThreads = 2, constraintFun = "sum"), 
                 permuteGeneral(11, 5, constraintFun = "sum", keepResults = TRUE))
    
    #### With Repetition
    numVec <- rnorm(9)
    expect_equal(permuteGeneral(numVec, 5, TRUE, nThreads = 2, constraintFun = "prod"),
                 permuteGeneral(numVec, 5, TRUE, constraintFun = "prod", keepResults = TRUE))
    
    #### Multisets
    numVec <- runif(9)
    expect_equal(permuteGeneral(numVec, 5, freqs = rep(1:3, 3), nThreads = 2, constraintFun = "mean"), 
                 permuteGeneral(numVec, 5, freqs = rep(1:3, 3), constraintFun = "mean", keepResults = TRUE))
    
    
    ######********************** Upper Only *******************#########
    #### NO Repetition
    expect_equal(permuteGeneral(30, 10, nThreads = 2, upper = 30000, constraintFun = "max"), 
                 permuteGeneral(30, 10, upper = 30000, constraintFun = "max", keepResults = TRUE))
    
    
    ######********************** Lower Only *******************#########
    #### NO Repetition
    total = permuteCount(30, 10)
    expect_equal(permuteGeneral(30, 10, nThreads = 2, lower = total - 30000, constraintFun = "prod"), 
                 permuteGeneral(30, 10, lower = total - 30000, constraintFun = "prod", keepResults = TRUE))
    
    #### With Repetition
    numVec <- rnorm(12)
    total = permuteCount(12, 8, TRUE)
    expect_equal(permuteGeneral(numVec, 8, TRUE, nThreads = 2, lower = total - 30000, constraintFun = "mean"),
                 permuteGeneral(numVec, 8, TRUE, lower = total - 30000, constraintFun = "mean", keepResults = TRUE))
    
    ######********************** Upper & Lower *******************#########
    #### Multisets
    expect_equal(permuteGeneral(10, 7, freqs = rep(1:5, 2), nThreads = 2, 
                              lower = 100000, upper = 130000, constraintFun = "sum"), 
                 permuteGeneral(10, 7, freqs = rep(1:5, 2), lower = 100000, 
                              upper = 130000, constraintFun = "sum", keepResults = TRUE))
})

test_that("comboGeneral produces correct results with Parallel and GMP enabled", {
    
    set.seed(16)
    ## N.B. Parallel has no effect when number of results is less than 20000
    
    ######********************** lower Only *******************#########
    #### NO Repetition
    ## comboCount(100, 50) - 30000
    ## Big Integer ('bigz') :
    ## [1] 100891344545564193334812467256
    expect_equal(comboGeneral(100, 50, nThreads = 2, 
                              lower = "100891344545564193334812467256"),
                 comboGeneral(100, 50, lower = "100891344545564193334812467256"))
    
    ######********************** Upper & Lower *******************#########
    #### With Repetition
    numVec <- rnorm(100)
    expect_equal(comboGeneral(numVec, 20, TRUE, nThreads = 2,
                              lower = "24551856075980529000000",
                              upper = "24551856075980529030000"),
                 comboGeneral(numVec, 20, TRUE, 
                              lower = "24551856075980529000000",
                              upper = "24551856075980529030000"))
    
    #### Multisets
    expect_equal(comboGeneral(factor(1:100), 20, freqs = rep(1:5, 20), nThreads = 2,
                              lower = "12298205197060077000000", 
                              upper = "12298205197060077030000"), 
                 comboGeneral(factor(1:100), 20, freqs = rep(1:5, 20), 
                              lower = "12298205197060077000000", 
                              upper = "12298205197060077030000"))
})

test_that("comboGeneral produces correct results with Parallel, GMP, and constrainFun enabled", {
    
    set.seed(17)
    ## N.B. Parallel has no effect when number of results is less than 20000
    
    ######********************** lower Only *******************#########
    #### NO Repetition
    ## comboCount(100, 50) - 30000
    ## Big Integer ('bigz') :
    ## [1] 100891344545564193334812467256
    expect_equal(comboGeneral(100, 50, nThreads = 2, 
                              lower = "100891344545564193334812467256", 
                              constraintFun = "sum"),
                 comboGeneral(100, 50, lower = "100891344545564193334812467256",
                              constraintFun = "sum", keepResults = TRUE))
        
    #### Multisets
    ## comboCount(100, 20, freqs = rep(1:5, 20)) - 100000
    ## Big Integer ('bigz') :
    ## [1] 12298205197060077177152
    numVec <- rnorm(100)
    expect_equal(comboGeneral(numVec, 20, freqs = rep(1:5, 20), nThreads = 2,
                              lower = "12298205197060077177152", 
                              constraintFun = "mean"), 
                 comboGeneral(numVec, 20, freqs = rep(1:5, 20),
                              lower = "12298205197060077177152",
                              constraintFun = "mean", keepResults = TRUE))
    
    ######********************** Upper & Lower *******************#########
    #### With Repetition
    numVec <- rnorm(100)
    expect_equal(comboGeneral(numVec, 20, TRUE, nThreads = 2,
                              lower = "24551856075980529000000",
                              upper = "24551856075980529030000", 
                              constraintFun = "min"),
                 comboGeneral(numVec, 20, TRUE, 
                              lower = "24551856075980529000000",
                              upper = "24551856075980529030000",
                              constraintFun = "min", keepResults = TRUE))
})

test_that("permuteGeneral produces correct results with Parallel and GMP enabled", {
    
    set.seed(18)
    ## N.B. Parallel has no effect when number of results is less than 20000
    
    ######********************** Lower Only *******************#########
    #### NO Repetition
    # permuteCount(30, 15) - 30000
    # Big Integer ('bigz') :
    # [1] 202843204931727330000
    expect_equal(permuteGeneral(30, 15, nThreads = 2, lower = "202843204931727330000"), 
                 permuteGeneral(30, 15, lower = "202843204931727330000"))
    
    #### With Repetition
    # permuteCount(15, 15, TRUE) - 30000
    # Big Integer ('bigz') :
    # [1] 437893890380829375
    numVec <- rnorm(15)
    expect_equal(permuteGeneral(numVec, 15, TRUE, nThreads = 2, 
                                lower = "437893890380829375"),
                 permuteGeneral(numVec, 15, TRUE, 
                                lower = "437893890380829375"))
    
    ######********************** Upper & Lower *******************#########
    #### Multisets
    expect_equal(permuteGeneral(factor(1:30), 15, freqs = rep(1:5, 6), nThreads = 2,
                                lower = "7260161756223012000000",
                                upper = "7260161756223012030000"), 
                 permuteGeneral(factor(1:30), 15, freqs = rep(1:5, 6),
                                lower = "7260161756223012000000",
                                upper = "7260161756223012030000"))
})

test_that("permuteGeneral produces correct results with Parallel, GMP, and constrainFun", {
    
    set.seed(19)
    ## N.B. Parallel has no effect when number of results is less than 20000
    
    ######********************** Lower Only *******************#########
    #### NO Repetition
    # permuteCount(30, 15) - 30000
    # Big Integer ('bigz') :
    # [1] 202843204931727330000
    expect_equal(permuteGeneral(30, 15, nThreads = 2, 
                                lower = "202843204931727330000", constraintFun = "sum"), 
                 permuteGeneral(30, 15, lower = "202843204931727330000",
                                constraintFun = "sum", keepResults = TRUE))
    
    #### Multisets
    # permuteCount(30, 15, freqs = rep(1:5, 6)) - 30000
    # Big Integer ('bigz') :
    # [1] 7260161756223012226620
    numVec <- runif(30)
    expect_equal(permuteGeneral(numVec, 15, freqs = rep(1:5, 6), nThreads = 2,
                                lower = "7260161756223012226620", constraintFun = "mean"), 
                 permuteGeneral(numVec, 15, freqs = rep(1:5, 6), 
                                lower = "7260161756223012226620", 
                                constraintFun = "mean", keepResults = TRUE))
    
    ######********************** Upper & Lower *******************#########
    #### With Repetition
    numVec <- rnorm(15)
    expect_equal(permuteGeneral(numVec, 15, TRUE, nThreads = 2,
                                lower = "437893890380000000",
                                upper = "437893890380030000", constraintFun = "min"),
                 permuteGeneral(numVec, 15, TRUE,
                                lower = "437893890380000000",
                                upper = "437893890380030000", 
                                constraintFun = "min", keepResults = TRUE))
})
