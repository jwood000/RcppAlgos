reprex::reprex({
    #' Update
    #' ---
    #' Here is a solution that is a good starting place, but could probably be improved:

    library(RcppAlgos)
    getCombs <- function(myMat, upper = NULL, minYears = NULL) {

        numRows <- nrow(myMat)
        myColNames <- colnames(myMat)

        if (is.null(minYears))  ## set default
            repZero <- numRows - 1
        else if (minYears >= numRows || minYears < 1)  ## check for extreme cases
            repZero <- numRows - 1
        else
            repZero <- numRows - minYears

        combs <- comboGeneral(v = 0:numRows, m = numRows,
                              freqs = c(repZero, rep(1, numRows)),
                              upper = upper)

        ## I think this part could be improved
        out <- lapply(1:nrow(combs), \(x) {
            myRows <- myMat[combs[x,],]

            if (is.null(nrow(myRows)))
                result <- !is.na(myRows)
            else
                result <- complete.cases(t(myRows))

            myColNames[result]
        })

        myRowNames <- rownames(myMat)
        names(out) <- lapply(1:nrow(combs), \(x) {
            myRowNames[combs[x, combs[x, ] > 0]]
        })
        out
    }

    #' Here is the output for the OP's example. (The OP is missing 5 of the outcomes below):

    testMat <- matrix(
        c(NA, 50, NA, 85,
          110, 75, 76, 86,
          120, NA, 78, 87,
          130, 100, 80, 88),
        nrow = 4, byrow = TRUE
    )

    row.names(testMat) <- 2000:2003
    colnames(testMat) <- LETTERS[1:4]
    getCombs(testMat)

    #' However, this answer, or any future answer for that matter, won't get you every combination as you have 144 countries and 47 years of data.  That produces a very _very_ large number. Every combination of any length up to _n_ is equivalent to the [power set](https://en.wikipedia.org/wiki/Power_set). The number of elements in the power set is simply `2^n`. Since we are not counting the equivalent of the empty set, we need to subtract one, thus:

    suppressWarnings(suppressMessages(library(gmp)))
    sub.bigz(pow.bigz(2, 47), 1)

    #' Yes, that is over one hundred trillion!!! You will probably need to rethink your approach as there are simply too many outcomes.
    #'
    #' All is not lost! You can make use of the `upper` argument, to limit the number of outcomes, so you can still investigate possible combinations. Observe:

    set.seed(11111)
    biggerTest <- matrix(sample(100, 20*20, replace = TRUE), nrow = 20)

    colnames(biggerTest) <- LETTERS[1:20]
    rownames(biggerTest) <- 1988:2007

    ## set 10% of values to NA
    myNAs <- sample(400, 400 / 10)
    biggerTest[myNAs] <- NA

    biggerTest[1:6, 1:10]

    ## Getting all 1,048,575 results takes a good bit of time
    system.time(allResults <- getCombs(biggerTest))

    ## Using upper greatly reduces the amount of time
    system.time(smallSampTest <- getCombs(biggerTest, upper = 10000))

    #' Alternatively, you can use the `minYears` argument to only return results with a minimum number of combinations of years. For example, per the OP's comments to @CPak's answer, if you only want to see results with 15 or more years of combinations, we have:

    system.time(minYearTest <- getCombs(biggerTest, minYears = 15))

    set.seed(123)
    minYearTest[sample(length(minYearTest), 5)]

    #' Or use both arguments together:

    system.time(
        bothConstraintsTest <- getCombs(biggerTest, 10000, minYears = 10)
    )

    bothConstraintsTest[1:5]

    #' Explanation
    #' ===
    #' The first thing we need to do is to determine every combination of _n_ years. This boils down to finding all _n_-tuples of the [multiset](https://en.wikipedia.org/wiki/Multiset) `c(rep(0, n-1), 1:n)` or equivalently, the power set of an _n_ element set minus the empty set. For example, for the years `2000:2003` (4 year span), the possible combinations are given by:

    comboGeneral(v = 0:4, m = 4, freqs = c(3, rep(1, 4)))

    #' Now, we iterate over each row of our combinations where each row tells us which combination of rows from the original matrix to test for `NAs`. If the particular combination only contains one result, we determine which indices are not `NA`. That is easily carried out by `!is.na(`. If we have more than one row, we employ `complete.cases(t` to obtain the columns that have only numbers (i.e. no occurrences of `NA`).
    #'
    #' After this we are just using indexing to obtain names for our outcomes and Voila, we have our desired results.
}, advertise = FALSE, venue = "r", html_preview = FALSE, wd = ".")
