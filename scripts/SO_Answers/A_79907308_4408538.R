reprex::reprex({

    #' The accepted answer provides a good workaround for older versions of `RcppAlgos`, using `tryCatch()` to fall back to `permuteGeneral()`.
    #'
    #' However, as of version `2.10.0`, this limitation has been removed. The underlying composition algorithms were extended, so this case is now supported directly:

    library(RcppAlgos)
    packageVersion("RcppAlgos")

    head(compositionsGeneral(3, 6, repetition = TRUE, target = 10))

    tail(compositionsGeneral(3, 6, repetition = TRUE, target = 10))

    #'
    #' ### Memory Efficient Iterators
    #'
    #' Since the question is tagged **memory-efficient**, it is also worth noting that `RcppAlgos` provides memory efficient iterators.
    #'
    #' For larger problems, the number of results can be enormous, so generating them all at once may be impractical. This is where iterators really start to shine.
    #'

    ## Huge example
    compositionsCount(100, 10, repetition = TRUE, target = 500)

    it <- compositionsIter(100, 10, repetition = TRUE, target = 500)

    ## See the first few results
    replicate(3, it@nextIter(), simplify = FALSE)

    ## Go to the nth result instantly
    it[[1e12]]

    ## Continuing iterating, producing the next 3 results
    replicate(3, it@nextIter(), simplify = FALSE)

    ## Note the state
    it@summary()

    #'
    #' And to show that memory doesn't explode, we will use the library `profmem` to demonstrate. In the example below, we will use a simple for loop instead of `replicate` so we don't incur a penalty for the returned list and compare this to generating chunks via `nextNIter()`.
    #'

    library(profmem)

    ## Very small threshold... only 2000 bytes
    options(profmem.threshold = 2000)

    ## Warm-up (see explanation below)
    it <- compositionsIter(100, 10, repetition = TRUE, target = 500)
    for (i in seq_len(1000)) it@nextIter()

    ## Fresh iterator, profile after warm-up
    profmem({
        it <- compositionsIter(100, 10, repetition = TRUE, target = 500)
        for (i in seq_len(1e6)) it@nextIter()
    })

    ## This will trigger a large allocation... a matrix of 1e6 x 10
    profmem({
        it <- compositionsIter(100, 10, repetition = TRUE, target = 500)
        it@nextNIter(1e+06)
    })

    #'
    #' Note that there were no allocations above the 2000 byte threshold when advancing results one at a time with `nextIter()`!
    #'
    #' ### Warm-up Explanation
    #'
    #' When profiling memory in R it is helpful to perform a short warm-up run first. From the vignette [Simple Memory Profiling in R](<https://cran.r-project.org/web/packages/profmem/vignettes/profmem.html>):
    #'
    #' > "The `profmem()` function builds upon existing memory profiling features available in R. It logs *every* memory allocation done by plain R code as well as those done by native code such as C and Fortran. For each entry, it records the size (in bytes) and the name of the functions on the call stack."
    #'
    #' If you don't perform the initial warm-up run, you will see many results like this:
    #'
    #' > `tryCatch() -> tryCatchList() -> tryCatchOne() -> doTryCatch() -> compile() -> genCode() -> cmp() -> cmpCall() -> tryInline() -> h()`
    #'
    #' Which really has nothing to with the iterator itself.
    #'
    #' *Disclosure: Joseph Wood is the author of the `RcppAlgos` package.*
    #'

}, advertise = FALSE, venue = "r", html_preview = FALSE, wd = ".")
