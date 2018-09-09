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
    expect_equal(comboGeneral(50, 10, Parallel = TRUE, upper = 100000), 
                 comboGeneral(50, 10, upper = 100000))
    
    #### With Repetition
    numVec <- rnorm(25)
    expect_equal(comboGeneral(numVec, 8, TRUE, Parallel = TRUE, upper = 100000),
                 comboGeneral(numVec, 8, TRUE, upper = 100000))
    
    #### Multisets
    expect_equal(comboGeneral(factor(1:50), 7, freqs = rep(1:5, 10), Parallel = TRUE, upper = 100000), 
                 comboGeneral(factor(1:50), 7, freqs = rep(1:5, 10), upper = 100000))
    
    
    ######********************** Lower Only *******************#########
    #### NO Repetition
    total = comboCount(50, 10)
    expect_equal(comboGeneral(50, 10, Parallel = TRUE, lower = total - 100000), 
                 comboGeneral(50, 10, lower = total - 100000))
    
    #### With Repetition
    numVec <- rnorm(25)
    total = comboCount(25, 8, TRUE)
    expect_equal(comboGeneral(numVec, 8, TRUE, Parallel = TRUE, lower = total - 100000),
                 comboGeneral(numVec, 8, TRUE, lower = total - 100000))
    
    #### Multisets
    total = comboCount(50, 7, freqs = rep(1:5, 10))
    expect_equal(comboGeneral(factor(1:50), 7, freqs = rep(1:5, 10), Parallel = TRUE, lower = total - 100000), 
                 comboGeneral(factor(1:50), 7, freqs = rep(1:5, 10), lower = total - 100000))
    
    ######********************** Upper & Lower *******************#########
    #### NO Repetition
    expect_equal(comboGeneral(50, 10, Parallel = TRUE, lower = 100000, upper = 200000), 
                 comboGeneral(50, 10, lower = 100000, upper = 200000))
    
    #### With Repetition
    numVec <- rnorm(25)
    expect_equal(comboGeneral(numVec, 8, TRUE, Parallel = TRUE, lower = 100000, upper = 200000),
                 comboGeneral(numVec, 8, TRUE, lower = 100000, upper = 200000))
    
    #### Multisets
    expect_equal(comboGeneral(factor(1:50), 7, freqs = rep(1:5, 10), Parallel = TRUE, lower = 100000, upper = 200000), 
                 comboGeneral(factor(1:50), 7, freqs = rep(1:5, 10), lower = 100000, upper = 200000))
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
    expect_equal(comboGeneral(50, 10, Parallel = TRUE, upper = 100000, constraintFun = "max"), 
                 comboGeneral(50, 10, upper = 100000, constraintFun = "max", keepResults = TRUE))
    
    #### With Repetition
    numVec <- rnorm(25)
    expect_equal(comboGeneral(numVec, 8, TRUE, Parallel = TRUE, upper = 100000, constraintFun = "min"),
                 comboGeneral(numVec, 8, TRUE, upper = 100000, constraintFun = "min", keepResults = TRUE))
    
    #### Multisets
    numVec <- rnorm(50)
    expect_equal(comboGeneral(numVec, 7, freqs = rep(1:5, 10), Parallel = TRUE, upper = 100000, constraintFun = "sum"), 
                 comboGeneral(numVec, 7, freqs = rep(1:5, 10), upper = 100000, constraintFun = "sum", keepResults = TRUE))
    
    
    ######********************** Lower Only *******************#########
    #### NO Repetition
    total = comboCount(50, 10)
    expect_equal(comboGeneral(50, 10, Parallel = TRUE, lower = total - 100000, constraintFun = "prod"), 
                 comboGeneral(50, 10, lower = total - 100000, constraintFun = "prod", keepResults = TRUE))
    
    #### With Repetition
    numVec <- rnorm(25)
    total = comboCount(25, 8, TRUE)
    expect_equal(comboGeneral(numVec, 8, TRUE, Parallel = TRUE, lower = total - 100000, constraintFun = "mean"),
                 comboGeneral(numVec, 8, TRUE, lower = total - 100000, constraintFun = "mean", keepResults = TRUE))
    
    #### Multisets
    total = comboCount(50, 7, freqs = rep(1:5, 10))
    expect_equal(comboGeneral(50, 7, freqs = rep(1:5, 10), Parallel = TRUE, 
                              lower = total - 100000, constraintFun = "max"), 
                 comboGeneral(50, 7, freqs = rep(1:5, 10), lower = total - 100000, 
                              constraintFun = "max", keepResults = TRUE))
    
    ######********************** Upper & Lower *******************#########
    #### NO Repetition
    expect_equal(comboGeneral(50, 10, Parallel = TRUE, lower = 100000, upper = 200000, constraintFun = "mean"), 
                 comboGeneral(50, 10, lower = 100000, upper = 200000, constraintFun = "mean", keepResults = TRUE))
    
    #### With Repetition
    numVec <- rnorm(25)
    expect_equal(comboGeneral(numVec, 8, TRUE, Parallel = TRUE, lower = 100000, upper = 200000, constraintFun = "prod"),
                 comboGeneral(numVec, 8, TRUE, lower = 100000, upper = 200000, constraintFun = "prod", keepResults = TRUE))
    
    #### Multisets
    expect_equal(comboGeneral(50, 7, freqs = rep(1:5, 10), Parallel = TRUE, 
                              lower = 100000, upper = 200000, constraintFun = "mean"), 
                 comboGeneral(50, 7, freqs = rep(1:5, 10), lower = 100000, 
                              upper = 200000, constraintFun = "mean", keepResults = TRUE))
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
    ## permuteCount(10, 5, T)
    ## [1] 100000
    numVec <- rnorm(10)
    expect_equal(permuteGeneral(numVec, 5, TRUE, Parallel = TRUE), permuteGeneral(numVec, 5, TRUE))
    
    #### Multisets
    ## permuteCount(9, 5, freqs = rep(1:3, 3))
    ## [1] 40260
    expect_equal(permuteGeneral(factor(1:9), 5, freqs = rep(1:3, 3), Parallel = TRUE), 
                 permuteGeneral(factor(1:9), 5, freqs = rep(1:3, 3)))
    
    
    ######********************** Upper Only *******************#########
    #### NO Repetition
    expect_equal(permuteGeneral(30, 10, Parallel = TRUE, upper = 100000), 
                 permuteGeneral(30, 10, upper = 100000))
    
    #### With Repetition
    numVec <- rnorm(10)
    expect_equal(permuteGeneral(numVec, 8, TRUE, Parallel = TRUE, upper = 100000),
                 permuteGeneral(numVec, 8, TRUE, upper = 100000))
    
    #### Multisets
    expect_equal(permuteGeneral(factor(1:10), 7, freqs = rep(1:5, 2), Parallel = TRUE, upper = 100000), 
                 permuteGeneral(factor(1:10), 7, freqs = rep(1:5, 2), upper = 100000))
    
    
    ######********************** Lower Only *******************#########
    #### NO Repetition
    total = permuteCount(30, 10)
    expect_equal(permuteGeneral(30, 10, Parallel = TRUE, lower = total - 100000), 
                 permuteGeneral(30, 10, lower = total - 100000))
    
    #### With Repetition
    numVec <- rnorm(12)
    total = permuteCount(12, 8, TRUE)
    expect_equal(permuteGeneral(numVec, 8, TRUE, Parallel = TRUE, lower = total - 100000),
                 permuteGeneral(numVec, 8, TRUE, lower = total - 100000))
    
    #### Multisets
    total = permuteCount(30, 7, freqs = rep(1:5, 6))
    expect_equal(permuteGeneral(factor(1:30), 7, freqs = rep(1:5, 6), Parallel = TRUE, lower = total - 100000), 
                 permuteGeneral(factor(1:30), 7, freqs = rep(1:5, 6), lower = total - 100000))
    
    ######********************** Upper & Lower *******************#########
    #### NO Repetition
    expect_equal(permuteGeneral(30, 10, Parallel = TRUE, lower = 100000, upper = 200000), 
                 permuteGeneral(30, 10, lower = 100000, upper = 200000))
    
    #### With Repetition
    numVec <- rnorm(12)
    expect_equal(permuteGeneral(numVec, 8, TRUE, Parallel = TRUE, lower = 100000, upper = 200000),
                 permuteGeneral(numVec, 8, TRUE, lower = 100000, upper = 200000))
    
    #### Multisets
    expect_equal(permuteGeneral(factor(1:30), 7, freqs = rep(1:5, 6), Parallel = TRUE, lower = 100000, upper = 200000), 
                 permuteGeneral(factor(1:30), 7, freqs = rep(1:5, 6), lower = 100000, upper = 200000))
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
    numVec <- rnorm(10)
    expect_equal(permuteGeneral(numVec, 5, TRUE, Parallel = TRUE, constraintFun = "prod"),
                 permuteGeneral(numVec, 5, TRUE, constraintFun = "prod", keepResults = TRUE))
    
    #### Multisets
    numVec <- runif(9)
    expect_equal(permuteGeneral(numVec, 5, freqs = rep(1:3, 3), Parallel = TRUE, constraintFun = "mean"), 
                 permuteGeneral(numVec, 5, freqs = rep(1:3, 3), constraintFun = "mean", keepResults = TRUE))
    
    
    ######********************** Upper Only *******************#########
    #### NO Repetition
    expect_equal(permuteGeneral(30, 10, Parallel = TRUE, upper = 100000, constraintFun = "max"), 
                 permuteGeneral(30, 10, upper = 100000, constraintFun = "max", keepResults = TRUE))
    
    #### With Repetition
    numVec <- rnorm(10)
    expect_equal(permuteGeneral(numVec, 8, TRUE, Parallel = TRUE, upper = 100000, constraintFun = "min"),
                 permuteGeneral(numVec, 8, TRUE, upper = 100000, constraintFun = "min", keepResults = TRUE))
    
    #### Multisets
    numVec <- rnorm(10)
    expect_equal(permuteGeneral(numVec, 7, freqs = rep(1:5, 2), Parallel = TRUE, upper = 100000, constraintFun = "sum"), 
                 permuteGeneral(numVec, 7, freqs = rep(1:5, 2), upper = 100000, constraintFun = "sum", keepResults = TRUE))
    
    
    ######********************** Lower Only *******************#########
    #### NO Repetition
    total = permuteCount(30, 10)
    expect_equal(permuteGeneral(30, 10, Parallel = TRUE, lower = total - 100000, constraintFun = "prod"), 
                 permuteGeneral(30, 10, lower = total - 100000, constraintFun = "prod", keepResults = TRUE))
    
    #### With Repetition
    numVec <- rnorm(12)
    total = permuteCount(12, 8, TRUE)
    expect_equal(permuteGeneral(numVec, 8, TRUE, Parallel = TRUE, lower = total - 100000, constraintFun = "mean"),
                 permuteGeneral(numVec, 8, TRUE, lower = total - 100000, constraintFun = "mean", keepResults = TRUE))
    
    #### Multisets
    total = permuteCount(10, 7, freqs = rep(1:5, 2))
    expect_equal(permuteGeneral(10, 7, freqs = rep(1:5, 2), Parallel = TRUE, 
                              lower = total - 100000, constraintFun = "max"), 
                 permuteGeneral(10, 7, freqs = rep(1:5, 2), lower = total - 100000, 
                              constraintFun = "max", keepResults = TRUE))
    
    ######********************** Upper & Lower *******************#########
    #### NO Repetition
    expect_equal(permuteGeneral(30, 10, Parallel = TRUE, lower = 100000, upper = 200000, constraintFun = "mean"), 
                 permuteGeneral(30, 10, lower = 100000, upper = 200000, constraintFun = "mean", keepResults = TRUE))
    
    #### With Repetition
    numVec <- rnorm(12)
    expect_equal(permuteGeneral(numVec, 8, TRUE, Parallel = TRUE, lower = 100000, 
                                upper = 200000, constraintFun = "prod"),
                 permuteGeneral(numVec, 8, TRUE, lower = 100000, upper = 200000, 
                                constraintFun = "prod", keepResults = TRUE))
    
    #### Multisets
    expect_equal(permuteGeneral(10, 7, freqs = rep(1:5, 2), Parallel = TRUE, 
                              lower = 100000, upper = 200000, constraintFun = "mean"), 
                 permuteGeneral(10, 7, freqs = rep(1:5, 2), lower = 100000, 
                              upper = 200000, constraintFun = "mean", keepResults = TRUE))
})

test_that("comboGeneral produces correct results with Parallel and GMP enabled", {
    
    set.seed(16)
    ## N.B. Parallel has no effect when number of results is less than 20000
    
    ######********************** lower Only *******************#########
    #### NO Repetition
    ## comboCount(100, 50) - 100000
    ## Big Integer ('bigz') :
    ## [1] 100891344545564193334812397256
    expect_equal(comboGeneral(100, 50, Parallel = TRUE, 
                              lower = "100891344545564193334812397256"),
                 comboGeneral(100, 50, lower = "100891344545564193334812397256"))
    
    #### With Repetition
    ## comboCount(100, 20, T) - 100000
    ## Big Integer ('bigz') :
    ## [1] 24551856075980529665105
    numVec <- rnorm(100)
    expect_equal(comboGeneral(numVec, 20, TRUE, Parallel = TRUE,
                              lower = "24551856075980529665105"),
                 comboGeneral(numVec, 20, TRUE, 
                              lower = "24551856075980529665105"))
    
    #### Multisets
    ## comboCount(100, 20, freqs = rep(1:5, 20)) - 100000
    ## Big Integer ('bigz') :
    ## [1] 12298205197060077107152
    expect_equal(comboGeneral(factor(1:100), 20, freqs = rep(1:5, 20), Parallel = TRUE,
                              lower = "12298205197060077107152"), 
                 comboGeneral(factor(1:100), 20, freqs = rep(1:5, 20),
                              lower = "12298205197060077107152"))
    
    ######********************** Upper & Lower *******************#########
    #### NO Repetition
    expect_equal(comboGeneral(100, 50, Parallel = TRUE, 
                              lower = "100891344545564193334812000000", 
                              upper = "100891344545564193334812100000"), 
                 comboGeneral(100, 50, 
                              lower = "100891344545564193334812000000", 
                              upper = "100891344545564193334812100000"))
    
    #### With Repetition
    numVec <- rnorm(100)
    expect_equal(comboGeneral(numVec, 20, TRUE, Parallel = TRUE,
                              lower = "24551856075980529000000",
                              upper = "24551856075980529100000"),
                 comboGeneral(numVec, 20, TRUE, 
                              lower = "24551856075980529000000",
                              upper = "24551856075980529100000"))
    
    #### Multisets
    expect_equal(comboGeneral(factor(1:100), 20, freqs = rep(1:5, 20), Parallel = TRUE,
                              lower = "12298205197060077000000", 
                              upper = "12298205197060077100000"), 
                 comboGeneral(factor(1:100), 20, freqs = rep(1:5, 20), 
                              lower = "12298205197060077000000", 
                              upper = "12298205197060077100000"))
})

test_that("comboGeneral produces correct results with Parallel, GMP, and constrainFun enabled", {
    
    set.seed(17)
    ## N.B. Parallel has no effect when number of results is less than 20000
    
    ######********************** lower Only *******************#########
    #### NO Repetition
    ## comboCount(100, 50) - 100000
    ## Big Integer ('bigz') :
    ## [1] 100891344545564193334812397256
    expect_equal(comboGeneral(100, 50, Parallel = TRUE, 
                              lower = "100891344545564193334812397256", 
                              constraintFun = "sum"),
                 comboGeneral(100, 50, lower = "100891344545564193334812397256",
                              constraintFun = "sum", keepResults = TRUE))
    
    expect_equal(rowSums(comboGeneral(100, 50, Parallel = TRUE, 
                              lower = "100891344545564193334812397256")), 
                 comboGeneral(100, 50, lower = "100891344545564193334812397256",
                              constraintFun = "sum", keepResults = TRUE)[,51])
    
    #### With Repetition
    ## comboCount(100, 20, T) - 100000
    ## Big Integer ('bigz') :
    ## [1] 24551856075980529665105
    numVec <- rnorm(100)
    expect_equal(comboGeneral(numVec, 20, TRUE, Parallel = TRUE,
                              lower = "24551856075980529665105", 
                              constraintFun = "prod"),
                 comboGeneral(numVec, 20, TRUE, 
                              lower = "24551856075980529665105",
                              constraintFun = "prod", keepResults = TRUE))
    
    ## double check the results with base function
    expect_equal(apply(comboGeneral(numVec, 20, TRUE, Parallel = TRUE,
                              lower = "24551856075980529665105"), 1, prod),
                 comboGeneral(numVec, 20, TRUE, 
                              lower = "24551856075980529665105",
                              constraintFun = "prod", keepResults = TRUE)[,21])
    
    #### Multisets
    ## comboCount(100, 20, freqs = rep(1:5, 20)) - 100000
    ## Big Integer ('bigz') :
    ## [1] 12298205197060077107152
    expect_equal(comboGeneral(numVec, 20, freqs = rep(1:5, 20), Parallel = TRUE,
                              lower = "12298205197060077107152", 
                              constraintFun = "mean"), 
                 comboGeneral(numVec, 20, freqs = rep(1:5, 20),
                              lower = "12298205197060077107152",
                              constraintFun = "mean", keepResults = TRUE))
    
    ## double check the results with base function
    expect_equal(rowMeans(comboGeneral(numVec, 20, freqs = rep(1:5, 20), Parallel = TRUE,
                              lower = "12298205197060077107152")), 
                 comboGeneral(numVec, 20, freqs = rep(1:5, 20),
                              lower = "12298205197060077107152",
                              constraintFun = "mean", keepResults = TRUE)[,21])
    
    ######********************** Upper & Lower *******************#########
    #### NO Repetition
    expect_equal(comboGeneral(100, 50, Parallel = TRUE, 
                              lower = "100891344545564193334812000000", 
                              upper = "100891344545564193334812100000", 
                              constraintFun = "max"), 
                 comboGeneral(100, 50, 
                              lower = "100891344545564193334812000000", 
                              upper = "100891344545564193334812100000",
                              constraintFun = "max", keepResults = TRUE))
    
    ## double check the results with base function
    expect_equal(apply(comboGeneral(100, 50, Parallel = TRUE, 
                              lower = "100891344545564193334812000000", 
                              upper = "100891344545564193334812100000"), 1, max), 
                 comboGeneral(100, 50, 
                              lower = "100891344545564193334812000000", 
                              upper = "100891344545564193334812100000",
                              constraintFun = "max", keepResults = TRUE)[,51])
    
    #### With Repetition
    numVec <- rnorm(100)
    expect_equal(comboGeneral(numVec, 20, TRUE, Parallel = TRUE,
                              lower = "24551856075980529000000",
                              upper = "24551856075980529100000", 
                              constraintFun = "min"),
                 comboGeneral(numVec, 20, TRUE, 
                              lower = "24551856075980529000000",
                              upper = "24551856075980529100000",
                              constraintFun = "min", keepResults = TRUE))
    
    ## double check the results with base function
    expect_equal(apply(comboGeneral(numVec, 20, TRUE, Parallel = TRUE,
                              lower = "24551856075980529000000",
                              upper = "24551856075980529100000"), 1, min),
                 comboGeneral(numVec, 20, TRUE, 
                              lower = "24551856075980529000000",
                              upper = "24551856075980529100000",
                              constraintFun = "min", keepResults = TRUE)[,21])
    
    #### Multisets
    expect_equal(comboGeneral(numVec, 20, freqs = rep(1:5, 20), Parallel = TRUE,
                              lower = "12298205197060077000000", 
                              upper = "12298205197060077100000", 
                              constraintFun = "sum"), 
                 comboGeneral(numVec, 20, freqs = rep(1:5, 20), 
                              lower = "12298205197060077000000", 
                              upper = "12298205197060077100000",
                              constraintFun = "sum", keepResults = TRUE))
    
    ## double check the results with base function
    expect_equal(rowSums(comboGeneral(numVec, 20, freqs = rep(1:5, 20), Parallel = TRUE,
                              lower = "12298205197060077000000", 
                              upper = "12298205197060077100000")),
                 comboGeneral(numVec, 20, freqs = rep(1:5, 20), 
                              lower = "12298205197060077000000", 
                              upper = "12298205197060077100000",
                              constraintFun = "sum", keepResults = TRUE)[,21])
})

test_that("permuteGeneral produces correct results with Parallel and GMP enabled", {
    
    set.seed(18)
    ## N.B. Parallel has no effect when number of results is less than 20000
    
    ######********************** Lower Only *******************#########
    #### NO Repetition
    # permuteCount(30, 15) - 100000
    # Big Integer ('bigz') :
    # [1] 202843204931727260000
    expect_equal(permuteGeneral(30, 15, Parallel = TRUE, lower = "202843204931727260000"), 
                 permuteGeneral(30, 15, lower = "202843204931727260000"))
    
    #### With Repetition
    # permuteCount(15, 15, TRUE) - 100000
    # Big Integer ('bigz') :
    # [1] 437893890380759375
    numVec <- rnorm(15)
    expect_equal(permuteGeneral(numVec, 15, TRUE, Parallel = TRUE, 
                                lower = "437893890380759375"),
                 permuteGeneral(numVec, 15, TRUE, 
                                lower = "437893890380759375"))
    
    #### Multisets
    # permuteCount(30, 15, freqs = rep(1:5, 6)) - 100000
    # Big Integer ('bigz') :
    # [1] 7260161756223012156620
    expect_equal(permuteGeneral(factor(1:30), 15, freqs = rep(1:5, 6), Parallel = TRUE,
                                lower = "7260161756223012156620"), 
                 permuteGeneral(factor(1:30), 15, freqs = rep(1:5, 6), 
                                lower = "7260161756223012156620"))
    
    ######********************** Upper & Lower *******************#########
    #### NO Repetition
    expect_equal(permuteGeneral(30, 15, Parallel = TRUE, 
                                lower = "202843204931727100000",
                                upper = "202843204931727200000"), 
                 permuteGeneral(30, 15, 
                                lower = "202843204931727100000",
                                upper = "202843204931727200000"))
    
    #### With Repetition
    numVec <- rnorm(15)
    expect_equal(permuteGeneral(numVec, 15, TRUE, Parallel = TRUE,
                                lower = "437893890380000000",
                                upper = "437893890380100000"),
                 permuteGeneral(numVec, 15, TRUE,
                                lower = "437893890380000000",
                                upper = "437893890380100000"))
    
    #### Multisets
    expect_equal(permuteGeneral(factor(1:30), 15, freqs = rep(1:5, 6), Parallel = TRUE,
                                lower = "7260161756223012000000",
                                upper = "7260161756223012100000"), 
                 permuteGeneral(factor(1:30), 15, freqs = rep(1:5, 6),
                                lower = "7260161756223012000000",
                                upper = "7260161756223012100000"))
})

test_that("permuteGeneral produces correct results with Parallel, GMP, and constrainFun", {
    
    set.seed(19)
    ## N.B. Parallel has no effect when number of results is less than 20000
    
    ######********************** Lower Only *******************#########
    #### NO Repetition
    # permuteCount(30, 15) - 100000
    # Big Integer ('bigz') :
    # [1] 202843204931727260000
    expect_equal(permuteGeneral(30, 15, Parallel = TRUE, 
                                lower = "202843204931727260000", constraintFun = "sum"), 
                 permuteGeneral(30, 15, lower = "202843204931727260000",
                                constraintFun = "sum", keepResults = TRUE))
    
    ## double check the results with base function
    expect_equal(rowSums(permuteGeneral(30, 15, Parallel = TRUE, 
                                lower = "202843204931727260000")), 
                 permuteGeneral(30, 15, lower = "202843204931727260000",
                                constraintFun = "sum", keepResults = TRUE)[,16])
    
    #### With Repetition
    # permuteCount(15, 15, TRUE) - 100000
    # Big Integer ('bigz') :
    # [1] 437893890380759375
    numVec <- rnorm(15)
    expect_equal(permuteGeneral(numVec, 15, TRUE, Parallel = TRUE, 
                                lower = "437893890380759375", constraintFun = "prod"),
                 permuteGeneral(numVec, 15, TRUE, 
                                lower = "437893890380759375", 
                                constraintFun = "prod", keepResults = TRUE))
    
    ## double check the results with base function
    expect_equal(apply(permuteGeneral(numVec, 15, TRUE, Parallel = TRUE, 
                                lower = "437893890380759375"), 1, prod),
                 permuteGeneral(numVec, 15, TRUE, 
                                lower = "437893890380759375", 
                                constraintFun = "prod", keepResults = TRUE)[,16])
    
    #### Multisets
    # permuteCount(30, 15, freqs = rep(1:5, 6)) - 100000
    # Big Integer ('bigz') :
    # [1] 7260161756223012156620
    numVec <- runif(30)
    expect_equal(permuteGeneral(numVec, 15, freqs = rep(1:5, 6), Parallel = TRUE,
                                lower = "7260161756223012156620", constraintFun = "mean"), 
                 permuteGeneral(numVec, 15, freqs = rep(1:5, 6), 
                                lower = "7260161756223012156620", 
                                constraintFun = "mean", keepResults = TRUE))
    
    ## double check the results with base function
    expect_equal(rowMeans(permuteGeneral(numVec, 15, freqs = rep(1:5, 6), Parallel = TRUE,
                                lower = "7260161756223012156620")), 
                 permuteGeneral(numVec, 15, freqs = rep(1:5, 6), 
                                lower = "7260161756223012156620", 
                                constraintFun = "mean", keepResults = TRUE)[,16])
    
    ######********************** Upper & Lower *******************#########
    #### NO Repetition
    expect_equal(permuteGeneral(30, 15, Parallel = TRUE, 
                                lower = "202843204931727100000",
                                upper = "202843204931727200000", constraintFun = "max"), 
                 permuteGeneral(30, 15, 
                                lower = "202843204931727100000",
                                upper = "202843204931727200000", 
                                constraintFun = "max", keepResults = TRUE))
    
    ## double check the results with base function
    expect_equal(apply(permuteGeneral(30, 15, Parallel = TRUE, 
                                lower = "202843204931727100000",
                                upper = "202843204931727200000"), 1, max), 
                 permuteGeneral(30, 15, 
                                lower = "202843204931727100000",
                                upper = "202843204931727200000", 
                                constraintFun = "max", keepResults = TRUE)[,16])
    
    #### With Repetition
    numVec <- rnorm(15)
    expect_equal(permuteGeneral(numVec, 15, TRUE, Parallel = TRUE,
                                lower = "437893890380000000",
                                upper = "437893890380100000", constraintFun = "min"),
                 permuteGeneral(numVec, 15, TRUE,
                                lower = "437893890380000000",
                                upper = "437893890380100000", 
                                constraintFun = "min", keepResults = TRUE))
    
    ## double check the results with base function
    expect_equal(apply(permuteGeneral(numVec, 15, TRUE, Parallel = TRUE,
                                lower = "437893890380000000",
                                upper = "437893890380100000"), 1, min),
                 permuteGeneral(numVec, 15, TRUE,
                                lower = "437893890380000000",
                                upper = "437893890380100000", 
                                constraintFun = "min", keepResults = TRUE)[,16])
    
    #### Multisets
    numVec <- runif(30)
    expect_equal(permuteGeneral(numVec, 15, freqs = rep(1:5, 6), Parallel = TRUE,
                                lower = "7260161756223012000000",
                                upper = "7260161756223012100000", constraintFun = "sum"), 
                 permuteGeneral(numVec, 15, freqs = rep(1:5, 6),
                                lower = "7260161756223012000000",
                                upper = "7260161756223012100000", 
                                constraintFun = "sum", keepResults = TRUE))
    
    ## double check the results with base function
    expect_equal(rowSums(permuteGeneral(numVec, 15, freqs = rep(1:5, 6), Parallel = TRUE,
                                lower = "7260161756223012000000",
                                upper = "7260161756223012100000")), 
                 permuteGeneral(numVec, 15, freqs = rep(1:5, 6),
                                lower = "7260161756223012000000",
                                upper = "7260161756223012100000", 
                                constraintFun = "sum", keepResults = TRUE)[,16])
})

test_that("comboSample produces appropriate error messages", {
    expect_error(comboSample(5, 3), "n and sampleVec cannot both be NULL")
    expect_error(comboSample(5,3,freqs = c(1,2,3,-2,1)), "in freqs must be a positive")
    expect_error(comboSample(5,3, n = 100), "n exceeds the maximum number of possible results")
    expect_error(comboSample(5,3, comboSample(5,3, sampleVec = 1:200)), "exceeds the maximum number of possible results")
    expect_error(comboSample(5,freqs = rep(1,6)), "the length of freqs must equal the")
    
    expect_error(comboSample(5,3, n = "5"), 
                 "n must be a number")
    expect_error(comboSample(5,3, n = 1:5), 
                 "length of n must be 1")
})

test_that("permuteSample produces appropriate error messages", {
    expect_error(permuteSample(5, 3), 
                 "n and sampleVec cannot both be NULL")
    expect_error(permuteSample(5,3,freqs = c(1,2,3,-2,1)), 
                 "in freqs must be a positive")
    expect_error(permuteSample(5,3, n = 100), 
                 "n exceeds the maximum number of possible results")
    expect_error(permuteSample(5,3, n = "5"), 
                 "n must be a number")
    expect_error(permuteSample(5,3, n = 1:5), 
                 "length of n must be 1")
    expect_error(permuteSample(5,3, permuteSample(5,3, sampleVec = 1:200)), 
                 "exceeds the maximum number of possible results")
    expect_error(permuteSample(5,freqs = rep(1,6)), 
                 "the length of freqs must equal the")
    expect_error(permuteSample(5, 4, sampleVec = "adfs"), 
                 "sampleVec must be convertible to a real number")
    expect_error(permuteSample(100, 10, n = 2^50),
                 "The number of rows cannot exceed")
    expect_error(permuteSample(100, 10, n = -200),
                 "n must be a positive number")
    expect_error(permuteSample(100, 10, sampleVec = c("62815650955529472001")),
                 "One or more of the requested values in sampleVec")
})