reprex::reprex({
    #' I was re-directed here from [List All Combinations With combn](https://stackoverflow.com/q/49570793/4408538) as this was one of the dupe targets. This is an old question and the answer provided by @RichScriven is very nice, but I wanted to give the community a few more options that are arguably more natural and more efficient (the last two).
    #'
    #' We first note that the output is very similar to the [Power Set](https://en.wikipedia.org/wiki/Power_set). Calling `powerSet` from the `rje` package, we see that indeed our output matches every element from the power set except the first element which is equivalent to the [Empty Set](https://en.wikipedia.org/wiki/Empty_set):

    x <- c("red", "blue", "black")
    rje::powerSet(x)

    #' If you don't want the first element, you can easily add a `[-1]` to the end of your function call like so : `rje::powerSet(x)[-1]`.
    #'
    #' The next two solutions are from the newer packages `arrangements` and `RcppAlgos` (I am the author), that will offer the user great gains in efficiency. Both of these packages are capable of generating combinations of [Multisets](https://en.wikipedia.org/wiki/Multiset).
    #'
    #' > **Why is this important?**
    #'
    #' It can be shown that there is a [one-to-one mapping](https://en.wikipedia.org/wiki/Injective_function) from the power set of `A` to all combinations of the multiset `c(rep(emptyElement, length(A)), A)` choose `length(A)`, where `emptyElement` is a representation of the empty set (like zero or a blank). With this in mind, observe:

    ## There is also a function called combinations in the
    ## rje package, so we fully declare the function with scope operator
    ##
    ## library(arrangements)
    arrangements::combinations(x = c("",x), k = 3, freq = c(2, rep(1, 3)))

    library(RcppAlgos)
    comboGeneral(c("",x), 3, freqs = c(2, rep(1, 3)))

    #' If you don't like dealing with blank elements and/or matrices, you can also return a list making use of `lapply`.

    lapply(seq_along(x), comboGeneral, v = x)

    lapply(seq_along(x), function(y) arrangements::combinations(x, y))

    #' Now we show that the last two methods are much more efficient (N.B. I removed `do.call(c, ` and `simplify = FALSE` from the answer provided by @RichSciven in order to compare generation of similar outputs. I also included `rje::powerSet` for good measure):

    library(microbenchmark)
    options(digits = 4)
    options(width = 90)

    set.seed(8128)
    bigX <- sort(sample(10^6, 20)) ## With this as an input, we will get 2^20 - 1 results.. i.e. 1,048,575

    microbenchmark(
        powSetRje = rje::powerSet(bigX),
        powSetRich = lapply(seq_along(bigX), combn, x = bigX),
        powSetArrange = lapply(seq_along(bigX), \(y) {
            arrangements::combinations(x = bigX, k = y)
        }),
        powSetAlgos = lapply(seq_along(bigX), comboGeneral, v = bigX),
        unit = "relative"
    )

    #' Even further, `arrangements` comes equipped with an argument called `layout` which lets the user choose a particular format for their output. One of those is `layout = "l"` for list. It is similar to setting `simplify = FALSE` in `combn` and allows us to obtain output like that of `powerSet`. Observe:

    do.call(c, lapply(seq_along(x), function(y) {
        arrangements::combinations(x, y, layout = "l")
    }))

    #' And the benchmarks:

    microbenchmark(
        powSetRje = rje::powerSet(bigX)[-1],
        powSetRich = do.call(c, lapply(seq_along(bigX), combn,
                                       x = bigX, simplify = FALSE)),
        powSetArrange = do.call(c, lapply(seq_along(bigX), \(y) {
            arrangements::combinations(bigX, y, layout = "l")}
        )),
        unit = "relative",
        times = 15
    )
}, advertise = FALSE, venue = "r", html_preview = FALSE, wd = ".")
