reprex::reprex({
    #' Here is an answer that is much faster than the OP's solution on large test cases. It doesn't rely on `paste`, but rather we take advantage of properties of numbers and vectorized operations.  We also use `comboGeneral` from the `RcppAlgos` package (I am the author) which is much faster than `combn` and `combn_prim` from the linked answer for generating combinations of a vector. First we show the efficiency gains of `comboGeneral` over the other functions:

    ## library(gRbase)
    library(RcppAlgos)
    library(microbenchmark)
    options(digits = 4)
    options(width = 90)
    options(scipen = 99)

    microbenchmark(
        gRBase = gRbase::combn_prim(300, 2),
        utils  = combn(300, 2),
        RcppAlgos = comboGeneral(300, 2),
        unit = "relative"
    )

    #' Now, we create a function to create some random reproducible data that will be passed to our test functions:

    makeTestSet <- function(vectorSize, elementSize, mySeed = 42, withRep = FALSE) {
        set.seed(mySeed)
        sapply(1:vectorSize, function(x) {
            paste(sample(10^6, s1 <- sample(2:elementSize, 1),
                         replace = withRep), collapse = " ")
        })
    }

    makeTestSet(5, 3)

    #' That looks good. Now, lets see if setting `fixed = TRUE` gets us any gains (as suggested above by @MichaelChirico):

    bigVec <- makeTestSet(10, 100000)

    microbenchmark(standard = strsplit(bigVec, " "),
                   withFixed = strsplit(bigVec, " ", fixed = TRUE),
                   times = 15, unit = "relative")

    #' @MichaelChirico was spot on. Putting it all together we get:

    combPairFast <- function(testVec) {
        lapply(strsplit(testVec, " ", fixed = TRUE), function(x) {
            combs <- RcppAlgos::comboGeneral(as.numeric(x), 2)
            unique(combs[,1] * (10)^(as.integer(log10(combs[,2])) + 1L) + combs[,2])
        })
    }

    ## test.vector defined above by OP
    test.vector <- c(
        "335261 344015 537633", "22404 132858",
        "254654 355860 488288", "219943 373817",
        "331839 404477"
    )

    combPairFast(test.vector)

    ## OP original code
    combPairOP <- function(testVec) {
        lapply(strsplit(testVec, " "), function(x) unique(
            apply(combn(x, 2), 2, function(y) paste0(y, collapse = "")))
        )
    }

    #' As stated in the comments by the OP, the maximum number is less than a million (600000 to be exact), which means that after we multiply one of the numbers by at most 10^6 and add it to another 6 digit number (equivalent to simply concatenating two strings of numbers), we are guaranteed to be within the numerical precision of base R (i.e. `2^53 - 1`). This is good because arithmetic operations on numerical numbers is much more efficient than strings operations.
    #'
    #' All that is left is to benchmark:

    test.vector <- makeTestSet(100, 50)

    microbenchmark(
        original = combPairOP(test.vector),
        new_impl = combPairFast(test.vector),
        times = 20,
        unit = "relative"
    )

    #' And on larger vectors:

    bigTest.vector <- makeTestSet(1000, 100, mySeed = 22, withRep = TRUE)

    ## Duplicate values exist
    any(sapply(strsplit(bigTest.vector, " ", fixed = TRUE), function(x) {
        any(duplicated(x))
    }))

    system.time(t1 <- combPairFast(bigTest.vector))

    system.time(t2 <- combPairOP(bigTest.vector))

    ## results are the same
    all.equal(t1, lapply(t2, as.numeric))
}, advertise = FALSE, venue = "r", html_preview = FALSE, wd = ".")
