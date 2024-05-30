reprex::reprex({
    #' There are several packages that will produce exactly what you need. Let's start with the classic `gtools`. Also, from the looks of the example provided by the OP, we are looking for permutations without repetition, not combinations with repetition.

    wordsList <- c("alice", "moon", "walks", "mars", "sings", "guitar", "bravo")

    library(gtools)
    attempt1 <- permutations(length(wordsList), 3, wordsList)
    head(attempt1)

    #' Then there is `arrangements`:

    library(arrangements)
    attempt2 <- permutations(wordsList, 3)
    head(attempt2)

    #' And finally, `RcppAlgos` (I am the author):

    library(RcppAlgos)
    attempt3 <- permuteGeneral(wordsList, 3)

    #' They are all fairly efficient and produce similar outcomes (different orderings)

    identical(attempt1[do.call(order,as.data.frame(attempt1)),],
              attempt2[do.call(order,as.data.frame(attempt3)),])

    identical(attempt2, attempt3)

    #' If you really want permutations with repetition, each function provides an argument for carrying out that function. For `gtools` set `repeats.allowed = TRUE`, for `arrangments` set `replace = TRUE`, and finally for `RcppAlgos` set `repetition = TRUE`.
    #'
    #' Since the OP is working with a `wordsList` with many words and is looking for all permutations chosen 15 at a time, the aforementioned methods will fail. However, there are some alternatives from `arrangements` as well as `RcppAlgos` that _may_ help.
    #'
    #' With `arrangements` you can use the function `getnext` and produce successive permutations. Generating them all will still take quite a long time.
    #'
    #' As with `arrangements`, `RcppAlgos` provides access to memory efficient iterators. They are quite flexible and very fast.

    it <- permuteIter(wordsList, 3)

    ## start iterating
    it@nextIter()

    ## draw n at a time
    it@nextNIter(n = 3)

    ## get the last one
    it@back()

    ## iterate in reverse
    it@prevIter()
}, advertise = FALSE, venue = "r", html_preview = FALSE, wd = ".")