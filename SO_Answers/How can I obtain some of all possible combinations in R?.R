reprex::reprex({
    #' It sounds like the OP is describing an [iterator](https://en.wikipedia.org/wiki/Iterator#:~:text=In%20computer%20programming%2C%20an%20iterator,via%20a%20container's%20interface.). There is a package called `iterators` that seems promising, however a major drawback is that one needs to generate the object upfront in order to iterate. That does not help us in the OP's case where we are trying to avoid generating the entire object upfront.
    #'
    #' There is a package `RcppAlgos` (I am the author) that provides flexible combinatorial iterators that allow one to traverse forwards, backwards, and even allows random access via `[[`. The underlying algorithms are written in `C++` for maximum efficiency. Results are produces on the fly, which keeps memory usage in check, all while preserving the state.

    library(RcppAlgos)
    it <- comboIter(25, 5)

    ## Get the first iteration
    it@nextIter()

    ## See the current state
    it@summary()

    ## Get the next 7 iterations
    it@nextNIter(n = 7)

    ## See the current state
    it@summary()

    ## Go back to the previous iteration
    it@prevIter()

    ## Skip ahead to the 20000th iteration instantly
    it[[20000]]

    ## See the current state. Notice we are at 20000
    it@summary()

    ## Start iterating from there
    it@nextIter()

    ## Reset the iterator
    it@startOver()

    ## Results are identical
    identical(
        replicate(choose(25, 5), it@nextIter(), simplify = "array"),
        combn(25, 5)
    )

    #' These iterators are very efficient. Even with all of the communication back and forth between `R` and `C++`, iterating "one at a time" is almost as fast as `combn` generating all of them upfront.

    library(microbenchmark)
    options(digits = 4)
    options(width = 90)

    ## helper functions for resetting the iterator
    one_at_a_time <- function() {
        it@startOver()
        replicate(choose(25, 5), it@nextIter())
    }

    microbenchmark(
        RcppAlgos_single = one_at_a_time(),
        combn_all = combn(25, 5),
        unit = "relative"
    )

    #' Even better, when we shift from "one at a time" to just a few at a time, the speed up is substantial. Observe:

    multiple_at_a_time <- function(n) {
        it@startOver()
        replicate(choose(25, 5) / n, it@nextNIter(n = n))
    }

    microbenchmark(
        RcppAlgos_chunks = multiple_at_a_time(30),
        combn_all = combn(25, 5),
        unit = "relative"
    )

    #' We can even evaluate each combination on the fly using the `FUN` parameter:

    ## the FUN.VALUE parameter is optional. If NULL (the default), when
    ## multiple results are requested via nextNIter, prevNIter, nextRemaining,
    ## and prevRemaining, a list will be returned. This parameter is modeled
    ## after the usage in vapply
    it_fun <- comboIter(25, 5, FUN = cumprod, FUN.VALUE = cumprod(1:5))

    it_fun@nextIter()

    it_fun@nextNIter(n = 2)

    ## See the previous results in reverse order
    it_fun@prevRemaining()

    #' Finally, the `[[` method allows for random access of a single result or multiple results:

    ## As seen above
    it[[20000]]

    ## Pass a random sample. N.B. In this case the state is unaffected. That
    ## is, it will remain whatever it was prior to passing the vector. In our
    ## case, it will still be on the 20000 index.
    set.seed(42)
    it[[sample(choose(25, 5), 10)]]

    ## State unaffected
    it@summary()

    ## Sample with replacement
    it[[sample(choose(25, 5), 5, replace = TRUE)]]

    #' There are also combinatorial sampling functions in `RcppAlgos`. For our current case, we would call upon `comboSample`.

    ## Same as above:
    ##     set.seed(42)
    ##     it[[sample(choose(25, 5), 10)]]
    comboSample(25, 5, seed = 42, n = 10)

    #' See my answer to [How to resample in R without repeating permutations?](https://stats.stackexchange.com/a/647976/139295) for more info.

}, advertise = FALSE, venue = "r", html_preview = FALSE, wd = ".")