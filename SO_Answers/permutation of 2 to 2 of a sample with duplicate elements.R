reprex::reprex({
    #' This is easily achieved with one of the many packages for generating permutations with repetition.

    library(gtools)
    gtools::permutations(3, 2, c(1, 2, 2), set = FALSE, repeats.allowed = TRUE)

    library(arrangements)
    ## output same as above
    arrangements::permutations(x = c(1,2,2), k = 2, replace = TRUE)


    library(RcppAlgos) ### I am the author
    ## output same as above
    RcppAlgos::permuteGeneral(c(1,2,2), 2, TRUE)
})