reprex::reprex({

    #' Several answers have already demonstrated how to obtain the lexicographic rank of UDAYPUR.
    #'
    #' If you have several problems like this with varying number of repeated letters, it can be very helpful if you can quickly verify your result.
    #'
    #' Here is a simple R script using the library `RcppAlgos` (I am the author), which quickly computes the rank without generating all permutations.
    #'
    #' The math/logic underneath the hood is exactly that mentioned in other posts.

    library(RcppAlgos)

    getWordRank <- function(word) {
        x   <- strsplit(word, "", fixed = TRUE)[[1]]
        tbl <- table(x)
        v   <- sort(names(tbl))
        fr  <- as.integer(tbl[v])
        permuteRank(x, v = v, freqs = fr)
    }

    ## O.P. example
    getWordRank("UDAYPUR")

    ## What about other examples??
    set.seed(42)

    random_words <- replicate(10, {
        paste(sample(LETTERS, sample(3:15, 1), replace = TRUE), collapse = "")
    })

    sapply(random_words, getWordRank)

}, advertise = FALSE, venue = "r", html_preview = FALSE, wd = ".")