context("testing with Parallel enabled")

test_that("comboGeneral produces correct results with Parallel enabled and no constrainFun", {
    
    set.seed(12)
    ## N.B. Parallel has no effect when number of results is less than 20000
    
    ######********************** All Results *******************#########
    #### NO Repetition
    ## comboCount(18, 9)
    ## [1] 48620 
    expect_equal(comboGeneral(18, 9, Parallel = TRUE), comboGeneral(18, 9))
    
    #### With Repetition
    ## comboCount(12, 8, T)
    ## [1] 75582
    numVec <- rnorm(12)
    expect_equal(comboGeneral(numVec, 8, TRUE, Parallel = TRUE), comboGeneral(numVec, 8, TRUE))
    
    #### Multisets
    ## comboCount(15, 7, freqs = rep(1:5, 3))
    ## [1] 73065
    expect_equal(comboGeneral(factor(1:15), 7, freqs = rep(1:5, 3), Parallel = TRUE), 
                 comboGeneral(factor(1:15), 7, freqs = rep(1:5, 3)))
    
    
    ######********************** Upper Only *******************#########
    #### NO Repetition
    expect_equal(comboGeneral(50, 10, Parallel = TRUE, upper = 30000), 
                 comboGeneral(50, 10, upper = 30000))
    
    #### With Repetition
    numVec <- rnorm(25)
    expect_equal(comboGeneral(numVec, 8, TRUE, Parallel = TRUE, upper = 30000),
                 comboGeneral(numVec, 8, TRUE, upper = 30000))
    
    #### Multisets
    expect_equal(comboGeneral(factor(1:50), 7, freqs = rep(1:5, 10), Parallel = TRUE, upper = 30000), 
                 comboGeneral(factor(1:50), 7, freqs = rep(1:5, 10), upper = 30000))
    
    
    ######********************** Lower Only *******************#########
    #### NO Repetition
    total = comboCount(50, 10)
    expect_equal(comboGeneral(50, 10, Parallel = TRUE, lower = total - 30000), 
                 comboGeneral(50, 10, lower = total - 30000))
    
    #### With Repetition
    numVec <- rnorm(25)
    total = comboCount(25, 8, TRUE)
    expect_equal(comboGeneral(numVec, 8, TRUE, Parallel = TRUE, lower = total - 30000),
                 comboGeneral(numVec, 8, TRUE, lower = total - 30000))
    
    #### Multisets
    total = comboCount(50, 7, freqs = rep(1:5, 10))
    expect_equal(comboGeneral(factor(1:50), 7, freqs = rep(1:5, 10), Parallel = TRUE, lower = total - 30000), 
                 comboGeneral(factor(1:50), 7, freqs = rep(1:5, 10), lower = total - 30000))
    
    ######********************** Upper & Lower *******************#########
    #### NO Repetition
    expect_equal(comboGeneral(50, 10, Parallel = TRUE, lower = 100000, upper = 130000), 
                 comboGeneral(50, 10, lower = 100000, upper = 130000))
    
    #### With Repetition
    numVec <- rnorm(25)
    expect_equal(comboGeneral(numVec, 8, TRUE, Parallel = TRUE, lower = 100000, upper = 130000),
                 comboGeneral(numVec, 8, TRUE, lower = 100000, upper = 130000))
    
    #### Multisets
    expect_equal(comboGeneral(factor(1:50), 7, freqs = rep(1:5, 10), Parallel = TRUE, lower = 100000, upper = 130000), 
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
    expect_equal(comboGeneral(18, 9, Parallel = TRUE, constraintFun = "sum"), 
                 comboGeneral(18, 9, constraintFun = "sum", keepResults = TRUE))
    
    #### With Repetition
    ## comboCount(12, 8, T)
    ## [1] 75582
    numVec <- rnorm(12)
    expect_equal(comboGeneral(numVec, 8, TRUE, Parallel = TRUE, constraintFun = "prod"),
                 comboGeneral(numVec, 8, TRUE, constraintFun = "prod", keepResults = TRUE))
    
    #### Multisets
    ## comboCount(15, 7, freqs = rep(1:5, 3))
    ## [1] 73065
    numVec <- runif(15)
    expect_equal(comboGeneral(numVec, 7, freqs = rep(1:5, 3), Parallel = TRUE, constraintFun = "mean"), 
                 comboGeneral(numVec, 7, freqs = rep(1:5, 3), constraintFun = "mean", keepResults = TRUE))
    
    
    ######********************** Upper Only *******************#########
    #### NO Repetition
    expect_equal(comboGeneral(50, 10, Parallel = TRUE, upper = 30000, constraintFun = "max"), 
                 comboGeneral(50, 10, upper = 30000, constraintFun = "max", keepResults = TRUE))
    
    #### With Repetition
    numVec <- rnorm(25)
    expect_equal(comboGeneral(numVec, 8, TRUE, Parallel = TRUE, upper = 30000, constraintFun = "min"),
                 comboGeneral(numVec, 8, TRUE, upper = 30000, constraintFun = "min", keepResults = TRUE))
    
    #### Multisets
    numVec <- rnorm(50)
    expect_equal(comboGeneral(numVec, 7, freqs = rep(1:5, 10), Parallel = TRUE, upper = 30000, constraintFun = "sum"), 
                 comboGeneral(numVec, 7, freqs = rep(1:5, 10), upper = 30000, constraintFun = "sum", keepResults = TRUE))
    
    
    ######********************** Lower Only *******************#########
    #### NO Repetition
    total = comboCount(50, 10)
    expect_equal(comboGeneral(50, 10, Parallel = TRUE, lower = total - 30000, constraintFun = "prod"), 
                 comboGeneral(50, 10, lower = total - 30000, constraintFun = "prod", keepResults = TRUE))
    
    #### With Repetition
    numVec <- rnorm(25)
    total = comboCount(25, 8, TRUE)
    expect_equal(comboGeneral(numVec, 8, TRUE, Parallel = TRUE, lower = total - 30000, constraintFun = "mean"),
                 comboGeneral(numVec, 8, TRUE, lower = total - 30000, constraintFun = "mean", keepResults = TRUE))
    
    #### Multisets
    total = comboCount(50, 7, freqs = rep(1:5, 10))
    expect_equal(comboGeneral(50, 7, freqs = rep(1:5, 10), Parallel = TRUE, 
                              lower = total - 30000, constraintFun = "max"), 
                 comboGeneral(50, 7, freqs = rep(1:5, 10), lower = total - 30000, 
                              constraintFun = "max", keepResults = TRUE))
    
    ######********************** Upper & Lower *******************#########
    #### NO Repetition
    expect_equal(comboGeneral(50, 10, Parallel = TRUE, lower = 100000, upper = 130000, constraintFun = "mean"), 
                 comboGeneral(50, 10, lower = 100000, upper = 130000, constraintFun = "mean", keepResults = TRUE))
    
    #### With Repetition
    numVec <- rnorm(25)
    expect_equal(comboGeneral(numVec, 8, TRUE, Parallel = TRUE, 
                              lower = 100000, upper = 130000, constraintFun = "prod"),
                 comboGeneral(numVec, 8, TRUE, 
                              lower = 100000, upper = 130000, constraintFun = "prod", keepResults = TRUE))
    
    #### Multisets
    expect_equal(comboGeneral(50, 7, freqs = rep(1:5, 10), Parallel = TRUE, 
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
    expect_equal(permuteGeneral(11, 5, Parallel = TRUE), permuteGeneral(11, 5))
    
    #### With Repetition
    ## permuteCount(9, 5, T)
    ## [1] 59049
    numVec <- rnorm(9)
    expect_equal(permuteGeneral(numVec, 5, TRUE, Parallel = TRUE), permuteGeneral(numVec, 5, TRUE))
    
    #### Multisets
    ## permuteCount(9, 5, freqs = rep(1:3, 3))
    ## [1] 40260
    expect_equal(permuteGeneral(factor(1:9), 5, freqs = rep(1:3, 3), Parallel = TRUE), 
                 permuteGeneral(factor(1:9), 5, freqs = rep(1:3, 3)))
    
    
    ######********************** Upper Only *******************#########
    #### NO Repetition
    expect_equal(permuteGeneral(30, 10, Parallel = TRUE, upper = 30000), 
                 permuteGeneral(30, 10, upper = 30000))
    
    #### With Repetition
    numVec <- rnorm(10)
    expect_equal(permuteGeneral(numVec, 8, TRUE, Parallel = TRUE, upper = 30000),
                 permuteGeneral(numVec, 8, TRUE, upper = 30000))
    
    #### Multisets
    expect_equal(permuteGeneral(factor(1:10), 7, freqs = rep(1:5, 2), Parallel = TRUE, upper = 30000), 
                 permuteGeneral(factor(1:10), 7, freqs = rep(1:5, 2), upper = 30000))
    
    
    ######********************** Lower Only *******************#########
    #### NO Repetition
    total = permuteCount(30, 10)
    expect_equal(permuteGeneral(30, 10, Parallel = TRUE, lower = total - 30000), 
                 permuteGeneral(30, 10, lower = total - 30000))
    
    #### With Repetition
    numVec <- rnorm(12)
    total = permuteCount(12, 8, TRUE)
    expect_equal(permuteGeneral(numVec, 8, TRUE, Parallel = TRUE, lower = total - 30000),
                 permuteGeneral(numVec, 8, TRUE, lower = total - 30000))
    
    #### Multisets
    total = permuteCount(30, 7, freqs = rep(1:5, 6))
    expect_equal(permuteGeneral(factor(1:30), 7, freqs = rep(1:5, 6), Parallel = TRUE, lower = total - 30000), 
                 permuteGeneral(factor(1:30), 7, freqs = rep(1:5, 6), lower = total - 30000))
    
    ######********************** Upper & Lower *******************#########
    #### NO Repetition
    expect_equal(permuteGeneral(30, 10, Parallel = TRUE, lower = 100000, upper = 130000), 
                 permuteGeneral(30, 10, lower = 100000, upper = 130000))
    
    #### With Repetition
    numVec <- rnorm(12)
    expect_equal(permuteGeneral(numVec, 8, TRUE, Parallel = TRUE, lower = 100000, upper = 130000),
                 permuteGeneral(numVec, 8, TRUE, lower = 100000, upper = 130000))
    
    #### Multisets
    expect_equal(permuteGeneral(factor(1:30), 7, freqs = rep(1:5, 6), Parallel = TRUE, lower = 100000, upper = 130000), 
                 permuteGeneral(factor(1:30), 7, freqs = rep(1:5, 6), lower = 100000, upper = 130000))
})

test_that("permuteGeneral produces correct results with Parallel enabled with logical vector", {
    
    ## N.B. Parallel has no effect when number of results is less than 20000
    
    ######********************** All Results *******************#########
    #### With Repetition
    ## permuteCount(2, 15, T)
    ## [1] 32768
    expect_equal(permuteGeneral(c(TRUE, FALSE), 15, TRUE, Parallel = TRUE), 
                 permuteGeneral(c(TRUE, FALSE), 15, TRUE))
    
    #### Multisets
    ## permuteCount(2, 15, freqs = c(9, 9))
    ## [1] 22880
    expect_equal(permuteGeneral(c(TRUE, FALSE), 15, freqs = c(9, 9), Parallel = TRUE), 
                 permuteGeneral(c(TRUE, FALSE), 15, freqs = c(9, 9)))
    
    
    ######********************** Upper Only *******************#########
    #### With Repetition
    expect_equal(permuteGeneral(c(TRUE, FALSE), 30, TRUE, Parallel = TRUE, upper = 30000), 
                 permuteGeneral(c(TRUE, FALSE), 30, TRUE, upper = 30000))
    
    #### Multisets
    expect_equal(permuteGeneral(c(TRUE, FALSE), 30, freqs = c(20, 20), Parallel = TRUE, upper = 30000), 
                 permuteGeneral(c(TRUE, FALSE), 30, freqs = c(20, 20), upper = 30000))
    
    ######********************** Lower Only *******************#########
    #### With Repetition
    total = permuteCount(2, 30, TRUE)
    expect_equal(permuteGeneral(c(TRUE, FALSE), 30, TRUE, Parallel = TRUE, lower = total - 30000), 
                 permuteGeneral(c(TRUE, FALSE), 30, TRUE, lower = total - 30000))
    
    #### Multisets
    total = permuteCount(2, 30, freqs = c(20, 20))
    expect_equal(permuteGeneral(c(TRUE, FALSE), 30, freqs = c(20, 20), Parallel = TRUE, lower = total - 30000), 
                 permuteGeneral(c(TRUE, FALSE), 30, freqs = c(20, 20), lower = total - 30000))
    
    ######********************** Upper & Lower *******************#########
    #### With Repetition
    expect_equal(permuteGeneral(c(TRUE, FALSE), 30, TRUE, Parallel = TRUE, lower = 100000, upper = 130000),
                 permuteGeneral(c(TRUE, FALSE), 30, TRUE, lower = 100000, upper = 130000))
    
    #### Multisets
    expect_equal(permuteGeneral(c(TRUE, FALSE), 30, freqs = c(20, 20),
                                Parallel = TRUE, lower = 100000, upper = 130000), 
                 permuteGeneral(c(TRUE, FALSE), 30, freqs = c(20, 20),
                                lower = 100000, upper = 130000))
})

test_that("permuteGeneral produces correct results with Parallel enabled and constrainFun", {
    
    set.seed(15)
    ## N.B. Parallel has no effect when number of results is less than 20000
    ## N.B. No longer need keepResult = TRUE if user only declares constraintFun
    
    ######********************** All Results *******************#########
    #### NO Repetition
    expect_equal(permuteGeneral(11, 5, Parallel = TRUE, constraintFun = "sum"), 
                 permuteGeneral(11, 5, constraintFun = "sum", keepResults = TRUE))
    
    #### With Repetition
    numVec <- rnorm(9)
    expect_equal(permuteGeneral(numVec, 5, TRUE, Parallel = TRUE, constraintFun = "prod"),
                 permuteGeneral(numVec, 5, TRUE, constraintFun = "prod", keepResults = TRUE))
    
    #### Multisets
    numVec <- runif(9)
    expect_equal(permuteGeneral(numVec, 5, freqs = rep(1:3, 3), Parallel = TRUE, constraintFun = "mean"), 
                 permuteGeneral(numVec, 5, freqs = rep(1:3, 3), constraintFun = "mean", keepResults = TRUE))
    
    
    ######********************** Upper Only *******************#########
    #### NO Repetition
    expect_equal(permuteGeneral(30, 10, Parallel = TRUE, upper = 30000, constraintFun = "max"), 
                 permuteGeneral(30, 10, upper = 30000, constraintFun = "max", keepResults = TRUE))
    
    #### With Repetition
    numVec <- rnorm(10)
    expect_equal(permuteGeneral(numVec, 8, TRUE, Parallel = TRUE, upper = 30000, constraintFun = "min"),
                 permuteGeneral(numVec, 8, TRUE, upper = 30000, constraintFun = "min", keepResults = TRUE))
    
    #### Multisets
    numVec <- rnorm(10)
    expect_equal(permuteGeneral(numVec, 7, freqs = rep(1:5, 2), Parallel = TRUE, upper = 30000, constraintFun = "sum"), 
                 permuteGeneral(numVec, 7, freqs = rep(1:5, 2), upper = 30000, constraintFun = "sum", keepResults = TRUE))
    
    
    ######********************** Lower Only *******************#########
    #### NO Repetition
    total = permuteCount(30, 10)
    expect_equal(permuteGeneral(30, 10, Parallel = TRUE, lower = total - 30000, constraintFun = "prod"), 
                 permuteGeneral(30, 10, lower = total - 30000, constraintFun = "prod", keepResults = TRUE))
    
    #### With Repetition
    numVec <- rnorm(12)
    total = permuteCount(12, 8, TRUE)
    expect_equal(permuteGeneral(numVec, 8, TRUE, Parallel = TRUE, lower = total - 30000, constraintFun = "mean"),
                 permuteGeneral(numVec, 8, TRUE, lower = total - 30000, constraintFun = "mean", keepResults = TRUE))
    
    #### Multisets
    total = permuteCount(10, 7, freqs = rep(1:5, 2))
    expect_equal(permuteGeneral(10, 7, freqs = rep(1:5, 2), Parallel = TRUE, 
                              lower = total - 30000, constraintFun = "max"), 
                 permuteGeneral(10, 7, freqs = rep(1:5, 2), lower = total - 30000, 
                              constraintFun = "max", keepResults = TRUE))
    
    ######********************** Upper & Lower *******************#########
    #### NO Repetition
    expect_equal(permuteGeneral(30, 10, Parallel = TRUE, lower = 100000, upper = 130000, constraintFun = "mean"), 
                 permuteGeneral(30, 10, lower = 100000, upper = 130000, constraintFun = "mean", keepResults = TRUE))
    
    #### With Repetition
    numVec <- rnorm(12)
    expect_equal(permuteGeneral(numVec, 8, TRUE, Parallel = TRUE, lower = 100000, 
                                upper = 130000, constraintFun = "prod"),
                 permuteGeneral(numVec, 8, TRUE, lower = 100000, upper = 130000, 
                                constraintFun = "prod", keepResults = TRUE))
    
    #### Multisets
    expect_equal(permuteGeneral(10, 7, freqs = rep(1:5, 2), Parallel = TRUE, 
                              lower = 100000, upper = 130000, constraintFun = "mean"), 
                 permuteGeneral(10, 7, freqs = rep(1:5, 2), lower = 100000, 
                              upper = 130000, constraintFun = "mean", keepResults = TRUE))
})

test_that("comboGeneral produces correct results with Parallel and GMP enabled", {
    
    set.seed(16)
    ## N.B. Parallel has no effect when number of results is less than 20000
    
    ######********************** lower Only *******************#########
    #### NO Repetition
    ## comboCount(100, 50) - 30000
    ## Big Integer ('bigz') :
    ## [1] 100891344545564193334812467256
    expect_equal(comboGeneral(100, 50, Parallel = TRUE, 
                              lower = "100891344545564193334812467256"),
                 comboGeneral(100, 50, lower = "100891344545564193334812467256"))
    
    #### With Repetition
    ## comboCount(100, 20, T) - 30000
    ## Big Integer ('bigz') :
    ## [1] 24551856075980529735105
    numVec <- rnorm(100)
    expect_equal(comboGeneral(numVec, 20, TRUE, Parallel = TRUE,
                              lower = "24551856075980529735105"),
                 comboGeneral(numVec, 20, TRUE, 
                              lower = "24551856075980529735105"))
    
    #### Multisets
    ## comboCount(100, 20, freqs = rep(1:5, 20)) - 30000
    ## Big Integer ('bigz') :
    ## [1] 12298205197060077177152
    expect_equal(comboGeneral(factor(1:100), 20, freqs = rep(1:5, 20), Parallel = TRUE,
                              lower = "12298205197060077177152"), 
                 comboGeneral(factor(1:100), 20, freqs = rep(1:5, 20),
                              lower = "12298205197060077177152"))
    
    ######********************** Upper & Lower *******************#########
    #### NO Repetition
    expect_equal(comboGeneral(100, 50, Parallel = TRUE, 
                              lower = "100891344545564193334812000000", 
                              upper = "100891344545564193334812030000"), 
                 comboGeneral(100, 50, 
                              lower = "100891344545564193334812000000", 
                              upper = "100891344545564193334812030000"))
    
    #### With Repetition
    numVec <- rnorm(100)
    expect_equal(comboGeneral(numVec, 20, TRUE, Parallel = TRUE,
                              lower = "24551856075980529000000",
                              upper = "24551856075980529030000"),
                 comboGeneral(numVec, 20, TRUE, 
                              lower = "24551856075980529000000",
                              upper = "24551856075980529030000"))
    
    #### Multisets
    expect_equal(comboGeneral(factor(1:100), 20, freqs = rep(1:5, 20), Parallel = TRUE,
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
    expect_equal(comboGeneral(100, 50, Parallel = TRUE, 
                              lower = "100891344545564193334812467256", 
                              constraintFun = "sum"),
                 comboGeneral(100, 50, lower = "100891344545564193334812467256",
                              constraintFun = "sum", keepResults = TRUE))
    
    expect_equal(rowSums(comboGeneral(100, 50, Parallel = TRUE, 
                              lower = "100891344545564193334812467256")), 
                 comboGeneral(100, 50, lower = "100891344545564193334812467256",
                              constraintFun = "sum", keepResults = TRUE)[,51])
    
    #### With Repetition
    ## comboCount(100, 20, T) - 30000
    ## Big Integer ('bigz') :
    ## [1] 24551856075980529735105
    numVec <- rnorm(100)
    expect_equal(comboGeneral(numVec, 20, TRUE, Parallel = TRUE,
                              lower = "24551856075980529735105", 
                              constraintFun = "prod"),
                 comboGeneral(numVec, 20, TRUE, 
                              lower = "24551856075980529735105",
                              constraintFun = "prod", keepResults = TRUE))
    
    ## double check the results with base function
    expect_equal(apply(comboGeneral(numVec, 20, TRUE, Parallel = TRUE,
                              lower = "24551856075980529735105"), 1, prod),
                 comboGeneral(numVec, 20, TRUE, 
                              lower = "24551856075980529735105",
                              constraintFun = "prod", keepResults = TRUE)[,21])
    
    #### Multisets
    ## comboCount(100, 20, freqs = rep(1:5, 20)) - 100000
    ## Big Integer ('bigz') :
    ## [1] 12298205197060077177152
    expect_equal(comboGeneral(numVec, 20, freqs = rep(1:5, 20), Parallel = TRUE,
                              lower = "12298205197060077177152", 
                              constraintFun = "mean"), 
                 comboGeneral(numVec, 20, freqs = rep(1:5, 20),
                              lower = "12298205197060077177152",
                              constraintFun = "mean", keepResults = TRUE))
    
    ## double check the results with base function
    expect_equal(rowMeans(comboGeneral(numVec, 20, freqs = rep(1:5, 20), Parallel = TRUE,
                              lower = "12298205197060077177152")), 
                 comboGeneral(numVec, 20, freqs = rep(1:5, 20),
                              lower = "12298205197060077177152",
                              constraintFun = "mean", keepResults = TRUE)[,21])
    
    ######********************** Upper & Lower *******************#########
    #### NO Repetition
    expect_equal(comboGeneral(100, 50, Parallel = TRUE, 
                              lower = "100891344545564193334812000000", 
                              upper = "100891344545564193334812030000", 
                              constraintFun = "max"), 
                 comboGeneral(100, 50, 
                              lower = "100891344545564193334812000000", 
                              upper = "100891344545564193334812030000",
                              constraintFun = "max", keepResults = TRUE))
    
    ## double check the results with base function
    expect_equal(apply(comboGeneral(100, 50, Parallel = TRUE, 
                              lower = "100891344545564193334812000000", 
                              upper = "100891344545564193334812030000"), 1, max), 
                 comboGeneral(100, 50, 
                              lower = "100891344545564193334812000000", 
                              upper = "100891344545564193334812030000",
                              constraintFun = "max", keepResults = TRUE)[,51])
    
    #### With Repetition
    numVec <- rnorm(100)
    expect_equal(comboGeneral(numVec, 20, TRUE, Parallel = TRUE,
                              lower = "24551856075980529000000",
                              upper = "24551856075980529030000", 
                              constraintFun = "min"),
                 comboGeneral(numVec, 20, TRUE, 
                              lower = "24551856075980529000000",
                              upper = "24551856075980529030000",
                              constraintFun = "min", keepResults = TRUE))
    
    ## double check the results with base function
    expect_equal(apply(comboGeneral(numVec, 20, TRUE, Parallel = TRUE,
                              lower = "24551856075980529000000",
                              upper = "24551856075980529030000"), 1, min),
                 comboGeneral(numVec, 20, TRUE, 
                              lower = "24551856075980529000000",
                              upper = "24551856075980529030000",
                              constraintFun = "min", keepResults = TRUE)[,21])
    
    #### Multisets
    expect_equal(comboGeneral(numVec, 20, freqs = rep(1:5, 20), Parallel = TRUE,
                              lower = "12298205197060077000000", 
                              upper = "12298205197060077030000", 
                              constraintFun = "sum"), 
                 comboGeneral(numVec, 20, freqs = rep(1:5, 20), 
                              lower = "12298205197060077000000", 
                              upper = "12298205197060077030000",
                              constraintFun = "sum", keepResults = TRUE))
    
    ## double check the results with base function
    expect_equal(rowSums(comboGeneral(numVec, 20, freqs = rep(1:5, 20), Parallel = TRUE,
                              lower = "12298205197060077000000", 
                              upper = "12298205197060077030000")),
                 comboGeneral(numVec, 20, freqs = rep(1:5, 20), 
                              lower = "12298205197060077000000", 
                              upper = "12298205197060077030000",
                              constraintFun = "sum", keepResults = TRUE)[,21])
})

test_that("permuteGeneral produces correct results with Parallel and GMP enabled", {
    
    set.seed(18)
    ## N.B. Parallel has no effect when number of results is less than 20000
    
    ######********************** Lower Only *******************#########
    #### NO Repetition
    # permuteCount(30, 15) - 30000
    # Big Integer ('bigz') :
    # [1] 202843204931727330000
    expect_equal(permuteGeneral(30, 15, Parallel = TRUE, lower = "202843204931727330000"), 
                 permuteGeneral(30, 15, lower = "202843204931727330000"))
    
    #### With Repetition
    # permuteCount(15, 15, TRUE) - 30000
    # Big Integer ('bigz') :
    # [1] 437893890380829375
    numVec <- rnorm(15)
    expect_equal(permuteGeneral(numVec, 15, TRUE, Parallel = TRUE, 
                                lower = "437893890380829375"),
                 permuteGeneral(numVec, 15, TRUE, 
                                lower = "437893890380829375"))
    
    #### Multisets
    # permuteCount(30, 15, freqs = rep(1:5, 6)) - 30000
    # Big Integer ('bigz') :
    # [1] 7260161756223012226620
    expect_equal(permuteGeneral(factor(1:30), 15, freqs = rep(1:5, 6), Parallel = TRUE,
                                lower = "7260161756223012226620"), 
                 permuteGeneral(factor(1:30), 15, freqs = rep(1:5, 6), 
                                lower = "7260161756223012226620"))
    
    ######********************** Upper & Lower *******************#########
    #### NO Repetition
    expect_equal(permuteGeneral(30, 15, Parallel = TRUE, 
                                lower = "202843204931727100000",
                                upper = "202843204931727130000"), 
                 permuteGeneral(30, 15, 
                                lower = "202843204931727100000",
                                upper = "202843204931727130000"))
    
    #### With Repetition
    numVec <- rnorm(15)
    expect_equal(permuteGeneral(numVec, 15, TRUE, Parallel = TRUE,
                                lower = "437893890380000000",
                                upper = "437893890380030000"),
                 permuteGeneral(numVec, 15, TRUE,
                                lower = "437893890380000000",
                                upper = "437893890380030000"))
    
    #### Multisets
    expect_equal(permuteGeneral(factor(1:30), 15, freqs = rep(1:5, 6), Parallel = TRUE,
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
    expect_equal(permuteGeneral(30, 15, Parallel = TRUE, 
                                lower = "202843204931727330000", constraintFun = "sum"), 
                 permuteGeneral(30, 15, lower = "202843204931727330000",
                                constraintFun = "sum", keepResults = TRUE))
    
    ## double check the results with base function
    expect_equal(rowSums(permuteGeneral(30, 15, Parallel = TRUE, 
                                lower = "202843204931727330000")), 
                 permuteGeneral(30, 15, lower = "202843204931727330000",
                                constraintFun = "sum", keepResults = TRUE)[,16])
    
    #### With Repetition
    # permuteCount(15, 15, TRUE) - 30000
    # Big Integer ('bigz') :
    # [1] 437893890380829375
    numVec <- rnorm(15)
    expect_equal(permuteGeneral(numVec, 15, TRUE, Parallel = TRUE, 
                                lower = "437893890380829375", constraintFun = "prod"),
                 permuteGeneral(numVec, 15, TRUE, 
                                lower = "437893890380829375", 
                                constraintFun = "prod", keepResults = TRUE))
    
    ## double check the results with base function
    expect_equal(apply(permuteGeneral(numVec, 15, TRUE, Parallel = TRUE, 
                                lower = "437893890380829375"), 1, prod),
                 permuteGeneral(numVec, 15, TRUE, 
                                lower = "437893890380829375", 
                                constraintFun = "prod", keepResults = TRUE)[,16])
    
    #### Multisets
    # permuteCount(30, 15, freqs = rep(1:5, 6)) - 30000
    # Big Integer ('bigz') :
    # [1] 7260161756223012226620
    numVec <- runif(30)
    expect_equal(permuteGeneral(numVec, 15, freqs = rep(1:5, 6), Parallel = TRUE,
                                lower = "7260161756223012226620", constraintFun = "mean"), 
                 permuteGeneral(numVec, 15, freqs = rep(1:5, 6), 
                                lower = "7260161756223012226620", 
                                constraintFun = "mean", keepResults = TRUE))
    
    ## double check the results with base function
    expect_equal(rowMeans(permuteGeneral(numVec, 15, freqs = rep(1:5, 6), Parallel = TRUE,
                                lower = "7260161756223012226620")), 
                 permuteGeneral(numVec, 15, freqs = rep(1:5, 6), 
                                lower = "7260161756223012226620", 
                                constraintFun = "mean", keepResults = TRUE)[,16])
    
    ######********************** Upper & Lower *******************#########
    #### NO Repetition
    expect_equal(permuteGeneral(30, 15, Parallel = TRUE, 
                                lower = "202843204931727100000",
                                upper = "202843204931727130000", constraintFun = "max"), 
                 permuteGeneral(30, 15, 
                                lower = "202843204931727100000",
                                upper = "202843204931727130000", 
                                constraintFun = "max", keepResults = TRUE))
    
    ## double check the results with base function
    expect_equal(apply(permuteGeneral(30, 15, Parallel = TRUE, 
                                lower = "202843204931727100000",
                                upper = "202843204931727130000"), 1, max), 
                 permuteGeneral(30, 15, 
                                lower = "202843204931727100000",
                                upper = "202843204931727130000", 
                                constraintFun = "max", keepResults = TRUE)[,16])
    
    #### With Repetition
    numVec <- rnorm(15)
    expect_equal(permuteGeneral(numVec, 15, TRUE, Parallel = TRUE,
                                lower = "437893890380000000",
                                upper = "437893890380030000", constraintFun = "min"),
                 permuteGeneral(numVec, 15, TRUE,
                                lower = "437893890380000000",
                                upper = "437893890380030000", 
                                constraintFun = "min", keepResults = TRUE))
    
    ## double check the results with base function
    expect_equal(apply(permuteGeneral(numVec, 15, TRUE, Parallel = TRUE,
                                lower = "437893890380000000",
                                upper = "437893890380030000"), 1, min),
                 permuteGeneral(numVec, 15, TRUE,
                                lower = "437893890380000000",
                                upper = "437893890380030000", 
                                constraintFun = "min", keepResults = TRUE)[,16])
    
    #### Multisets
    numVec <- runif(30)
    expect_equal(permuteGeneral(numVec, 15, freqs = rep(1:5, 6), Parallel = TRUE,
                                lower = "7260161756223012000000",
                                upper = "7260161756223012030000", constraintFun = "sum"), 
                 permuteGeneral(numVec, 15, freqs = rep(1:5, 6),
                                lower = "7260161756223012000000",
                                upper = "7260161756223012030000", 
                                constraintFun = "sum", keepResults = TRUE))
    
    ## double check the results with base function
    expect_equal(rowSums(permuteGeneral(numVec, 15, freqs = rep(1:5, 6), Parallel = TRUE,
                                lower = "7260161756223012000000",
                                upper = "7260161756223012030000")), 
                 permuteGeneral(numVec, 15, freqs = rep(1:5, 6),
                                lower = "7260161756223012000000",
                                upper = "7260161756223012030000", 
                                constraintFun = "sum", keepResults = TRUE)[,16])
})

test_that("combo/permuteGeneral produces correct error messages with Parallel", {
    expect_error(permuteGeneral(10, 5, Parallel = "TRUE"), "Not compatible with requested type")
})