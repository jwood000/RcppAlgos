context("testing comboGeneral and permuteGeneral results")

test_that("comboGeneral produces correct results with constraintFun only", {

    set.seed(13)
    ## N.B. No longer need keepResult = TRUE if user only declares constraintFun

    ######********************** All Results *******************#########
    #### NO Repetition
    ## comboCount(15, 8)
    ## [1] 6435
    expect_equal(rowSums(comboGeneral(15, 8)),
                 comboGeneral(15, 8, constraintFun = "sum")[, 9])

    #### With Repetition
    ## comboCount(10, 5, T)
    ## [1] 2002
    numVec <- rnorm(10)
    expect_equal(apply(comboGeneral(numVec, 5, TRUE), 1, prod),
                 comboGeneral(numVec, 5, TRUE, constraintFun = "prod")[, 6])

    #### Multisets
    ## comboCount(10, 7, freqs = rep(1:5, 2))
    ## [1] 5932
    numVec <- runif(10)
    expect_equal(rowMeans(comboGeneral(numVec, 7, freqs = rep(1:5, 2))),
                 comboGeneral(numVec, 7, freqs = rep(1:5, 2), constraintFun = "mean")[,8])


    ######********************** Upper Only *******************#########
    #### NO Repetition
    expect_equal(apply(comboGeneral(50, 10, upper = 3000), 1, max),
                 comboGeneral(50, 10, upper = 3000, constraintFun = "max")[,11])

    #### With Repetition
    numVec <- rnorm(25)
    expect_equal(apply(comboGeneral(numVec, 8, TRUE, upper = 3000), 1, min),
                 comboGeneral(numVec, 8, TRUE, upper = 3000, constraintFun = "min")[,9])

    #### Multisets
    numVec <- rnorm(50)
    expect_equal(rowSums(comboGeneral(numVec, 7, freqs = rep(1:5, 10), upper = 3000)),
                 comboGeneral(numVec, 7, freqs = rep(1:5, 10), upper = 3000, constraintFun = "sum")[,8])


    ######********************** Lower Only *******************#########
    #### NO Repetition
    total = comboCount(50, 10)
    expect_equal(apply(comboGeneral(50, 10, lower = total - 3000), 1, prod),
                 comboGeneral(50, 10, lower = total - 3000, constraintFun = "prod")[,11])

    #### With Repetition
    numVec <- rnorm(25)
    total = comboCount(25, 8, TRUE)
    expect_equal(rowMeans(comboGeneral(numVec, 8, TRUE, lower = total - 3000)),
                 comboGeneral(numVec, 8, TRUE, lower = total - 3000, constraintFun = "mean")[,9])

    #### Multisets
    total = comboCount(50, 7, freqs = rep(1:5, 10))
    expect_equal(apply(comboGeneral(50, 7, freqs = rep(1:5, 10),
                              lower = total - 3000), 1, max),
                 comboGeneral(50, 7, freqs = rep(1:5, 10), lower = total - 3000,
                              constraintFun = "max")[,8])

    ######********************** Upper & Lower *******************#########
    #### NO Repetition
    expect_equal(rowMeans(comboGeneral(50, 10, lower = 100000, upper = 103000)),
                 comboGeneral(50, 10, lower = 100000, upper = 103000, constraintFun = "mean")[,11])

    #### With Repetition
    numVec <- rnorm(25)
    expect_equal(apply(comboGeneral(numVec, 8, TRUE,
                              lower = 100000, upper = 103000), 1, prod),
                 comboGeneral(numVec, 8, TRUE,
                              lower = 100000, upper = 103000, constraintFun = "prod")[,9])

    #### Multisets
    expect_equal(rowMeans(comboGeneral(50, 7, freqs = rep(1:5, 10),
                              lower = 100000, upper = 103000)),
                 comboGeneral(50, 7, freqs = rep(1:5, 10), lower = 100000,
                              upper = 103000, constraintFun = "mean")[,8])
    
    set.seed(12345)
    myNums = rnorm(15)
    expect_equal(comboGeneral(myNums, 8, constraintFun = "sum", keepResults = TRUE)[, 9], 
                 rowSums(comboGeneral(myNums, 8)))
    
    expect_equal(comboGeneral(myNums, 5, TRUE, constraintFun = "sum", keepResults = TRUE)[, 6], 
                 rowSums(comboGeneral(myNums, 5, TRUE)))
    
    expect_equal(comboGeneral(myNums, 5, freqs = rep(2, 15), 
                              constraintFun = "sum", keepResults = TRUE)[, 6], 
                 rowSums(comboGeneral(myNums, 5, freqs = rep(2, 15))))
    
    ## test product
    expect_equal(comboGeneral(myNums, 8, constraintFun = "prod", keepResults = TRUE)[, 9], 
                 apply(comboGeneral(myNums, 8), 1, prod))
    
    ## test mean
    expect_equal(comboGeneral(myNums, 5, freqs = rep(2, 15), 
                              constraintFun = "mean", keepResults = TRUE)[, 6], 
                 rowMeans(comboGeneral(myNums, 5, freqs = rep(2, 15))))
    
    ## test max
    expect_equal(comboGeneral(myNums, 8, constraintFun = "max", keepResults = TRUE)[, 9], 
                 apply(comboGeneral(myNums, 8), 1, max))
    
    ## test min
    expect_equal(comboGeneral(myNums, 5, TRUE, constraintFun = "min", keepResults = TRUE)[, 6], 
                 apply(comboGeneral(myNums, 5, TRUE), 1, min))
})

test_that("comboGeneral produces correct results with GMP and constrainFun enabled", {

    set.seed(17)
    ######********************** lower Only *******************#########
    #### NO Repetition
    ## comboCount(100, 50) - 30000
    ## Big Integer ('bigz') :
    ## [1] 100891344545564193334812467256
    expect_equal(rowSums(comboGeneral(100, 50,
                              lower = "100891344545564193334812467256")),
                 comboGeneral(100, 50, lower = "100891344545564193334812467256",
                              constraintFun = "sum")[,51])

    #### With Repetition
    ## comboCount(100, 20, T) - 30000
    ## Big Integer ('bigz') :
    ## [1] 24551856075980529735105
    numVec <- rnorm(100)
    expect_equal(apply(comboGeneral(numVec, 20, TRUE,
                              lower = "24551856075980529735105"), 1, prod),
                 comboGeneral(numVec, 20, TRUE,
                              lower = "24551856075980529735105",
                              constraintFun = "prod")[,21])

    #### Multisets
    ## comboCount(100, 20, freqs = rep(1:5, 20)) - 100000
    ## Big Integer ('bigz') :
    ## [1] 12298205197060077177152
    expect_equal(rowMeans(comboGeneral(numVec, 20, freqs = rep(1:5, 20),
                              lower = "12298205197060077177152")),
                 comboGeneral(numVec, 20, freqs = rep(1:5, 20),
                              lower = "12298205197060077177152",
                              constraintFun = "mean")[,21])

    ######********************** Upper & Lower *******************#########
    #### NO Repetition
    expect_equal(apply(comboGeneral(100, 50,
                              lower = "100891344545564193334812000000",
                              upper = "100891344545564193334812003000"), 1, max),
                 comboGeneral(100, 50,
                              lower = "100891344545564193334812000000",
                              upper = "100891344545564193334812003000",
                              constraintFun = "max")[,51])

    #### With Repetition
    numVec <- rnorm(100)
    expect_equal(apply(comboGeneral(numVec, 20, TRUE,
                              lower = "24551856075980529000000",
                              upper = "24551856075980529003000"), 1, min),
                 comboGeneral(numVec, 20, TRUE,
                              lower = "24551856075980529000000",
                              upper = "24551856075980529003000",
                              constraintFun = "min")[,21])

    #### Multisets
    expect_equal(rowSums(comboGeneral(numVec, 20, freqs = rep(1:5, 20),
                              lower = "12298205197060077000000",
                              upper = "12298205197060077003000")),
                 comboGeneral(numVec, 20, freqs = rep(1:5, 20),
                              lower = "12298205197060077000000",
                              upper = "12298205197060077003000",
                              constraintFun = "sum")[,21])
})

test_that("permuteGeneral produces correct results with constrainFun only", {

    set.seed(15)
    ## N.B. No longer need keepResult = TRUE if user only declares constraintFun

    ######********************** All Results *******************#########
    #### NO Repetition
    expect_equal(rowSums(permuteGeneral(11, 5)),
                 permuteGeneral(11, 5, constraintFun = "sum")[,6])

    #### With Repetition
    numVec <- rnorm(9)
    expect_equal(apply(permuteGeneral(numVec, 5, TRUE), 1, prod),
                 permuteGeneral(numVec, 5, TRUE, constraintFun = "prod")[,6])

    #### Multisets
    numVec <- runif(9)
    expect_equal(rowMeans(permuteGeneral(numVec, 5, freqs = rep(1:3, 3))),
                 permuteGeneral(numVec, 5, freqs = rep(1:3, 3), constraintFun = "mean")[,6])


    ######********************** Upper Only *******************#########
    #### NO Repetition
    expect_equal(apply(permuteGeneral(30, 10, upper = 3000), 1, max),
                 permuteGeneral(30, 10, upper = 3000, constraintFun = "max")[,11])

    #### With Repetition
    numVec <- rnorm(10)
    expect_equal(apply(permuteGeneral(numVec, 8, TRUE, upper = 3000), 1, min),
                 permuteGeneral(numVec, 8, TRUE, upper = 3000, constraintFun = "min")[,9])

    #### Multisets
    numVec <- rnorm(10)
    expect_equal(rowSums(permuteGeneral(numVec, 7, freqs = rep(1:5, 2), upper = 3000)),
                 permuteGeneral(numVec, 7, freqs = rep(1:5, 2), upper = 3000, constraintFun = "sum")[,8])


    ######********************** Lower Only *******************#########
    #### NO Repetition
    total = permuteCount(30, 10)
    expect_equal(apply(permuteGeneral(30, 10, lower = total - 3000), 1, prod),
                 permuteGeneral(30, 10, lower = total - 3000, constraintFun = "prod")[,11])

    #### With Repetition
    numVec <- rnorm(12)
    total = permuteCount(12, 8, TRUE)
    expect_equal(rowMeans(permuteGeneral(numVec, 8, TRUE, lower = total - 3000)),
                 permuteGeneral(numVec, 8, TRUE, lower = total - 3000, constraintFun = "mean")[,9])

    #### Multisets
    total = permuteCount(10, 7, freqs = rep(1:5, 2))
    expect_equal(apply(permuteGeneral(10, 7, freqs = rep(1:5, 2),
                              lower = total - 3000), 1, max),
                 permuteGeneral(10, 7, freqs = rep(1:5, 2), lower = total - 3000,
                              constraintFun = "max")[,8])

    ######********************** Upper & Lower *******************#########
    #### NO Repetition
    expect_equal(rowMeans(permuteGeneral(30, 10, lower = 100000, upper = 103000)),
                 permuteGeneral(30, 10, lower = 100000, upper = 103000, constraintFun = "mean")[,11])

    #### With Repetition
    numVec <- rnorm(12)
    expect_equal(apply(permuteGeneral(numVec, 8, TRUE, lower = 100000,
                                upper = 103000), 1, prod),
                 permuteGeneral(numVec, 8, TRUE, lower = 100000, upper = 103000,
                                constraintFun = "prod")[,9])

    #### Multisets
    expect_equal(rowMeans(permuteGeneral(10, 7, freqs = rep(1:5, 2),
                              lower = 100000, upper = 103000)),
                 permuteGeneral(10, 7, freqs = rep(1:5, 2), lower = 100000,
                              upper = 103000, constraintFun = "mean")[,8])
})

test_that("permuteGeneral produces correct results with GMP and constrainFun", {

    set.seed(19)
    ######********************** Lower Only *******************#########
    #### NO Repetition
    # permuteCount(30, 15) - 30000
    # Big Integer ('bigz') :
    # [1] 202843204931727330000
    expect_equal(rowSums(permuteGeneral(30, 15,
                                lower = "202843204931727330000")),
                 permuteGeneral(30, 15, lower = "202843204931727330000",
                                constraintFun = "sum")[,16])

    #### With Repetition
    # permuteCount(15, 15, TRUE) - 30000
    # Big Integer ('bigz') :
    # [1] 437893890380829375
    numVec <- rnorm(15)
    expect_equal(apply(permuteGeneral(numVec, 15, TRUE,
                                lower = "437893890380829375"), 1, prod),
                 permuteGeneral(numVec, 15, TRUE,
                                lower = "437893890380829375",
                                constraintFun = "prod")[,16])

    #### Multisets
    # permuteCount(30, 15, freqs = rep(1:5, 6)) - 30000
    # Big Integer ('bigz') :
    # [1] 7260161756223012226620
    numVec <- runif(30)
    expect_equal(rowMeans(permuteGeneral(numVec, 15, freqs = rep(1:5, 6),
                                lower = "7260161756223012226620")),
                 permuteGeneral(numVec, 15, freqs = rep(1:5, 6),
                                lower = "7260161756223012226620",
                                constraintFun = "mean")[,16])

    ######********************** Upper & Lower *******************#########
    #### NO Repetition
    expect_equal(apply(permuteGeneral(30, 15,
                                lower = "202843204931727100000",
                                upper = "202843204931727103000"), 1, max),
                 permuteGeneral(30, 15,
                                lower = "202843204931727100000",
                                upper = "202843204931727103000",
                                constraintFun = "max")[,16])

    #### With Repetition
    numVec <- rnorm(15)
    expect_equal(apply(permuteGeneral(numVec, 15, TRUE,
                                lower = "437893890380000000",
                                upper = "437893890380003000"), 1, min),
                 permuteGeneral(numVec, 15, TRUE,
                                lower = "437893890380000000",
                                upper = "437893890380003000",
                                constraintFun = "min")[,16])

    #### Multisets
    numVec <- runif(30)
    expect_equal(rowSums(permuteGeneral(numVec, 15, freqs = rep(1:5, 6),
                                lower = "7260161756223012000000",
                                upper = "7260161756223012003000")),
                 permuteGeneral(numVec, 15, freqs = rep(1:5, 6),
                                lower = "7260161756223012000000",
                                upper = "7260161756223012003000",
                                constraintFun = "sum")[,16])
})
