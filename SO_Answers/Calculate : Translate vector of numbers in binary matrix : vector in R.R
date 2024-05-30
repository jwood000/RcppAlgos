reprex::reprex({
    #' Here is a one-liner using `RcppAlgos` (I am the author):

    vector <- c(1,2,3)

    library(RcppAlgos)
    toBinary <- function(v) permuteGeneral(0:1, length(v), TRUE)[-1,]

    toBinary(vector)

    #' The `[-1, ]` is to remove the row of all zeros.  This row would represent the empty set in a [power set](https://en.wikipedia.org/wiki/Power_set). In fact, what you are asking for is technically a mapping from the power set of a vector (minus the empty set of course) to a binary matrix.
    #'
    #' If you really want the `row.names` to be the actual permutations, you can use the `powerSet` function from the `rje` package. Observe:

    library(rje)
    nameTest <- toBinary(vector)
    row.names(nameTest) <- lapply(powerSet(rev(vector))[-1], sort)

    nameTest
}, advertise = FALSE, venue = "r", html_preview = FALSE, wd = ".")