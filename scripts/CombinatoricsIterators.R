reprex::reprex({
    #'
    #' This document covers working with combinatorial iterators in `RcppAlgos`. Combinatorial iterators in `RcppAlgos` are memory efficient like traditional iterator objects. They allow traversal of combinations/permutations/partitions/compositions one by one without the necessity for storing all results in memory.
    #'
    #' Unlike traditional combinatorial iterators, the iterators in `RcppAlgos` offers random access via the `[[` operator. This means, we can access the _n<sup>th</sup>_ [lexicographical order](<https://en.wikipedia.org/wiki/Lexicographical_order>) result on demand without having to first iterate over the previous _n - 1_ results.
    #'
    #' ***
    #'
    #' ## Iterating over Combinations and Permutations
    #'
    #' In order to iterate, we must initialize an iterator via `comboIter` or `permuteIter`. The interface is very similar to `comboGeneral` and `permuteGeneral`.
    #'

    library(RcppAlgos)
    options(width = 90)

    ## Initialize the iterator
    a = comboIter(5, 3)

    ## Get the first combination
    a$nextIter()

    ## And the next
    a$nextIter()

    ## Set the current iterator to a variable
    iter = a$currIter()
    i = 1

    ## Iterate until there are no more
    while (!is.null(iter)) {
        cat(i, " ", iter, "\n")
        iter = a$nextIter()
        i = i + 1
    }

    ## See the output of comboGeneral for comparison
    comboGeneral(5, 3, lower = 2)

    ## Call the summary method to see information about our iterator
    a$summary()

    #'
    #' ### Bidirectional Iterators
    #'
    #' Some of the combinatorial iterators in `RcppAlgos` are bidirectional iterators. This means that not only can we iterate in a forward manner (i.e. lexicographically), but we can also iterate backwards (i.e. [Reverse Lexicographical Order](<https://oeis.org/wiki/Orderings#Reverse_lexicographic_order>)) via the `prevIter` method(s).
    #'

    ## Using the same iterable from the previous section
    a$currIter()

    ## As the comment says, we call the prevIter method to see the last result
    a$prevIter()

    ## Get the previous result
    a$prevIter()

    ## As in the previous example, we set the current iterator to a variable
    iter = a$currIter()

    ## Defined above
    print(i)

    ## Iterate until we are at the very beginning. Note that the
    ## output is exactly the same as above, but in reverse order
    while (!is.null(iter)) {
        i = i - 1
        cat(i, " ", iter, "\n")
        iter = a$prevIter()
    }

    ## Call the summary method to see information about our iterator
    a$summary()

    #'
    #' ### Retrieving More than One Result at a Time
    #'
    #' There are four methods which allow for obtaining more than one result at a time: `nextNIter`, `prevNIter`, `nextRemaining`, and `prevRemaining`.
    #'

    ## Reset the iterator
    a$startOver()

    ## Get the next 4 combinations
    a$nextNIter(4)

    ## Get the summary. Note that the index has been updated
    a$summary()

    ## View the current combination
    a$currIter()

    ## Get the remaining combinations with nextRemaining
    a$nextRemaining()

    a$summary()

    #'
    #' Now, we look at the opposite direction.
    #'

    ## Get the previous 4 combinations
    a$prevNIter(4)

    ## Get the summary. Note that the index has been updated
    a$summary()

    ## View the current combination
    a$currIter()

    ## Get the remaining previous combinations with prevRemaining
    a$prevRemaining()

    a$summary()

    #'
    #' ### Random Access Iterator
    #'
    #' As with the bidirectional iterators, with some of the combinatorial iterators in `RcppAlgos`, we can jump to the _n<sup>th</sup>_ result without the need for iterating over the first _n - 1_ results.
    #'

    ## Reset the iterator
    a$startOver()

    ## How many total combinations do we have?
    a$summary()$totalResults

    ## Let's get the 3rd combination
    a[[3]]

    ## See the summary. Note that the index has been updated
    a$summary()

    ## Let's see the 9th combination
    a[[9]]

    ## What about the first and last combination?
    a$front()

    a$back()

    ## Again the index has been updated
    a$summary()

    a$currIter()

    #'
    #' We can also easily return a random sample of combinations with the `[[` operator by passing a vector of indices. In these cases, it should be noted that the current index will not be updated.
    #'

    ## Set the current index to the second combination
    a[[2]]

    set.seed(121)
    samp = sample(a$summary()$totalResults, 4)

    samp

    a[[samp]]

    ## Note that the current index remains unchanged
    a$summary()

    #'
    #' ### User Defined Functions
    #'
    #' Just as with `comboGeneral` and `permuteGeneral`, we can pass a user defined function to `comboIter` and `permuteIter`.
    #'

    ## Initialize the iterator
    b = permuteIter(LETTERS[1:4], 3, FUN = function(p) paste(p, collapse = ""),
                    FUN.VALUE = "a")
    b$nextIter()

    b$nextNIter(5)

    b$back()

    b$prevIter()

    b$prevNIter(5)

    b$nextRemaining()

    ## Random access
    b[[5]]

    b$prevRemaining()

    ## View the source vector
    b$sourceVector()

    #'
    #' ## New in Verison `2.5.0`
    #'
    #' As of version `2.5.0`, we no longer rely on `Rcpp` as a dependency, which means that we do not utilize `Rcpp` modules for exposing C++ classes. This is now carried out using external pointers (See [External pointers and weak references](<https://cran.r-project.org/doc/manuals/r-release/R-exts.html#External-pointers-and-weak-references>)) along with [S4 Classes](<http://adv-r.had.co.nz/S4.html>). We use the slots of `S4` classes for exposing each method so access is carried out with the "at sign", `@`. We have also added the ability to access each method with the "dollar sign", `$`, for backwards compatibility.
    #'
    #' ### Access Efficiency in `2.5.0+`
    #'
    #' Our tests show that accessing methods is much more efficient in `2.5.0+` compared to prior versions. In the below tests, we measure excecution time of calling `nextIter` multiple times in different versions. We will use the function `test_nextIter` for our testing. If one needs to reproduce, simply download the `2.4.3` tar here: https://cran.r-project.org/src/contrib/Archive/RcppAlgos/, change `RcppAlgos` to `RcppAlgos243` in a few place (e.g. `DESCRIPTION`, `NAMESPACE`, etc.), and rebuild.
    #'

    test_nextIter <- function(n, m, get_val = FALSE, v = 243) {
        a <- if (v == 243) {
            RcppAlgos243::comboIter(n, m)
        } else {
            comboIter(n, m)
        }

        total <- comboCount(n, m)

        if (get_val) {
            mat <- matrix(0L, nrow = total, ncol = m)
            for (i in 1:total) mat[i, ] <- a$nextIter()
            return(mat)
        } else {
            if (v == 243) {
                for (i in 1:total) a$nextIter()
            } else {
                for (i in 1:total) a@nextIter()
            }

            invisible(NULL)
        }
    }

    #'
    #' #### Version `2.4.3` Using `Rcpp`
    #'

    library(microbenchmark)
    ## Using R version 4.1.3
    comboCount(15, 8)

    microbenchmark(v243 = test_nextIter(15, 8))

    identical(test_nextIter(15, 8, get_val = TRUE),
              comboGeneral(15, 8))

    comboCount(25, 10)

    system.time(test_nextIter(25, 10))

    Rprof("Version243.out", memory.profiling = TRUE)
    test_nextIter(25, 10)
    Rprof(NULL)
    lapply(summaryRprof("Version243.out", memory = "both"), head)

    #'
    #' #### Version ``r packageVersion("RcppAlgos")`` (No `Rcpp`)
    #'

    curr_version <- as.integer(gsub("\\.", "", packageVersion("RcppAlgos")))
    curr_version
    microbenchmark(curr_v = test_nextIter(15, 8, v = curr_version))

    system.time(test_nextIter(25, 10, v = curr_version))

    identical(test_nextIter(15, 8, get_val = TRUE, v = curr_version),
              comboGeneral(15, 8))

    Rprof("Version250.out", memory.profiling = TRUE)
    test_nextIter(25, 10, v = curr_version)
    Rprof(NULL)
    lapply(summaryRprof("Version250.out", memory = "both"), head)

    #'
    #' #### Conclusions
    #'
    #' It appears that memory is the issue in previous versions. Indeed, if we look at [Memory statistics from Rprof](<https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Memory-statistics-from-Rprof>), and view both files with `memory = "stats"` we see that the C funciton, `duplicate`, appears to be the main culprit.
    #'

    ## We set index = 1 to ensure we get the very bottom of the stack

    ### Verison 2.4.3
    v243 <- summaryRprof("Version243.out", memory = "stats", index = 1)
    v243

    ## Version 2.5.0
    v250 <- summaryRprof("Version250.out", memory = "stats", index = 1)
    v250

    #'
    #' With verison `2.5.0+` there are only `r v250[[1]][["tot.duplications"]]` `tot.duplications` whereas with version `2.4.3` there are millions of `tot.duplications`. In fact, there are a total of `r format(v243[[1]][["tot.duplications"]], scientific=FALSE)` duplications with version `2.4.3`. This together with `comboCount(25, 10) = 3,268,760` implies that the C funciton, `duplicate`, is called about 3 times per iteration with older versions (i.e. `r format(v243[[1]][["tot.duplications"]], scientific=FALSE)` ` / 3268760 ~= ` `r round(v243[[1]][["tot.duplications"]] / 3268760, 4)`).
    #'
    #' ### Iterating over Partitions and Compositions of a Number
    #'
    #' For most partition cases, we have all of the capabilities of the standard `comboIter` and `permuteIter` except for bidirectionality (i.e. the `prevIter` methods). For cases involving standard multisets we also don't have random access methods.
    #'

    ## Similar illustration of comboIter(5, 3) at the top
    p = partitionsIter(16, 4)
    p@nextIter()

    p@nextIter()

    iter = p@currIter()
    i = 1

    while (!is.null(iter)) {
        cat(i, " ", iter, "\n")
        iter = p@nextIter()
        i = i + 1
    }

    partitionsGeneral(16, 4, lower = 2)

    p@summary()

    ## Using random access
    p[[7]]

    ## No previous iterators
    p@prevIter()

    #'
    #' For compositions, the options are limited to a subset of compositions with repetition.
    #'

    ## Similar illustration of comboIter(5, 3) at the top
    p = compositionsIter(6, 3, TRUE)
    p@nextIter()

    p@nextIter()

    iter = p@currIter()
    i = 1

    while (!is.null(iter)) {
        cat(i, " ", iter, "\n")
        iter = p@nextIter()
        i = i + 1
    }

    compositionsGeneral(6, 3, TRUE, lower = 2)

    p@summary()

    ## Using random access
    p[[7]]

    ## No previous iterators
    p@prevIter()

    #'
    #' ### Iterating over Constrained Combinations/Permutations
    #'
    #' Now, the combinatorial iterators have all of the features of their "general" analogs (I.e. `{combo|permute|partitions|compositions}General`), which includes constrained results.
    #'
    #' For general constrained cases, these iterators offer huge advantages over their "general" counterparts. Previously, one had to guess how many results there would be using the `upper` parameter as executing the function with no constraints meant the user could be waiting for a while or consume a large amount of resources.
    #'
    #' Another drawback is that it difficult to start generating from a particular point. With the "general" functions, if the `lower` parameter is used, we have to make a decision in order to disambiguate the use. Without constraints, using `lower` is easy to understand. It simply means to start generating results starting at a particular lexicographical result, which we can do efficiently (i.e. no need to generate the first `lower - 1` results). With constraints, it could mean one of two things:
    #'
    #'   1. Start checking from a particular lexicographical result without considering the constraint (as we do normally).
    #'   2. Start generating results from a particular result with regards to the final constrained output.
    #'
    #' In `RcppAlgos` we have always used the first interpretation. A big downside for the second point is that we don't have any fast algorithms for enumerating the total number of results, which reduces determining the _n<sup>th</sup>_ result to a brute force approach.
    #'
    #' With iterators, we can generate _n_ results with `nextNIter(n)` or calling `nextIter()` _n_ times (or some combination of the two). Then, if we want to continue iterating, we pick up where we left off fetching the _(n + 1)<sup>th</sup>_ result and beyond (if there are any results left). This allows us to keep memory low without sacrificing our current state.
    #'

    set.seed(55)
    s = runif(10, -5, 5)

    print(s)

    ## Using comboGeneral to retrieve all results
    comboGeneral(s, 5, constraintFun = "mean",
                 comparisonFun = "<", limitConstraints = -3)


    ## Using comboIter
    a = comboIter(s, 5, constraintFun = "mean",
                  comparisonFun = "<", limitConstraints = -3)

    ## See the first result
    a@nextIter()

    ## Get the next three
    a@nextNIter(3)

    ## See the summary... Note the totalResults and totalRemaining
    ## fields are NA as we are not able to calculate this upfront.
    a@summary()


    a@nextNIter(3)

    ## Get the rest
    a@nextRemaining()

    #'
    #' They are very efficient as well. Consider the example below where we use `comboGeneral` to generate all results without capping the output. Again, we are in a situation where we don't know _a priori_ how many results we will obtain.
    #'

    set.seed(77)
    s = runif(50, 20, 100)

    ## Over one trillion results to sift through
    comboCount(s, 15)

    time_all <- system.time({
        print(
            nrow(
                comboGeneral(s, 15,
                             constraintFun = "mean",
                             comparisonFun = ">",
                             limitConstraints = 83)
            )
        )
    })
    time_all

    ## Over 4 GBs of results
    (38935252 * 15 * 8) / 2^30

    #'
    #' Just over `r trunc(time_all[["elapsed"]])` seconds isn't bad, however 4 GBs could put a strain on your computer.
    #'
    #' Let's use iterators instead and only generate ten thousand at a time to keep memory low. We should mention here that the iterators are "smart" in that there is no fear in requesting more results than what is actually left. For example, if in the problem above, we had iterated to the 38<sup>th</sup> million result and requested 10 million more, we would only obtain 935,252 results.
    #'

    invisible(gc())
    time_iter <- system.time({
        a = comboIter(s, 15,
                      constraintFun = "mean",
                      comparisonFun = ">",
                      limitConstraints = 83)
        while (!is.null(a@nextNIter(1e4))) {}
        print(a@summary())
    })
    time_iter

    ## Only 1 MBs per iteration
    (1e4 * 15 * 8) / 2^20

    #' Wow! Using the iterator approach is not only easier on your RAM, but faster as well (`r time_all[["elapsed"]]` ` / ` `r time_iter[["elapsed"]]` ` ~= ` `r round(time_all[["elapsed"]] / time_iter[["elapsed"]], 4)`)! Our gains came strictly from memory efficiency (From over 4 GBs to just over 1 MB) as the underlying algorithm is exactly the same.
    #'
    #' Lastly, using iterators make some problems possible that would otherwise be intractable because of hardware. For instance, using the example above, if we changed the `limitConstraints` from 83 to 81 and tried the first approach your computer will most certainly become unusable (at least mine did). My memory usage shot up to over 30 GB and R became unresponsive. After a restart, I tried the second approach and obtained my result in just over 10 seconds barely noticing any jumps in memory:
    #'

    ## Don't run... consumes a huge chunk of memory
    # time_all <- system.time({
    #     print(
    #         nrow(
    #             comboGeneral(s, 15,
    #                          constraintFun = "mean",
    #                          comparisonFun = ">",
    #                          limitConstraints = 81)  ## 83 -->> 81
    #         )
    #     )
    # })

    ## No problem with iterators
    invisible(gc())
    system.time({
        a = comboIter(s, 15,
                      constraintFun = "mean",
                      comparisonFun = ">",
                      limitConstraints = 81)   ## 83 -->> 81
        while (!is.null(a@nextNIter(1e4))) {}
        print(a@summary())
    })

}, advertise = FALSE, venue = "r", html_preview = FALSE, wd = ".")