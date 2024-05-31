reprex::reprex({
    #' There are a few packages specifically built for this. The basic premise is that we need combinations with repetition of length `m` where `m` could be larger than the input vector. We start with the classic `gtools`:

    library(gtools)
    combinations(2, 3, letters[1:2], repeats.allowed = TRUE)

    #' And then there is `arrangements` which is a replacement for `iterpc` (the package linked by @Gregor in the comments above):

    ## library(arrangements)
    arrangements::combinations(letters[1:2], 3, replace = TRUE)

    #' And finally there is `RcppAlgos`, which I authored:

    library(RcppAlgos)
    comboGeneral(letters[1:2], 3, TRUE)

    #' `combn` is an awesome function that ships as one of the base packages with `R`, however one of its shortcomings is that it doesn't allow repetition (which is what is required here). I wrote a comprehensive overview for problems exactly like this one that can be found here: [A Walk Through a Slice of Combinatorics in R](https://stackoverflow.com/a/47983855/4408538).
}, advertise = FALSE, venue = "r", html_preview = FALSE, wd = ".")