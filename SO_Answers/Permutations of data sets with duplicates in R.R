reprex::reprex({
    #' What you are looking for is permutations of [multisets](https://en.wikipedia.org/wiki/Multiset).

    library(RcppAlgos)
    multiPerm <- permuteGeneral(1:4, freqs = rep(2,4))

    head(multiPerm)

    v <- rep(1:4, each = 2)
    v

    ## Also, utilize S3 'table' method
    head(permuteGeneral(table(v)))

    #' Or using the package `arrangements`:

    multiPerm2 <- arrangements::permutations(1:4, freq = rep(2, 4))

    #' Sanity check:

    identical(multiPerm, multiPerm2)

    suppressWarnings(suppressMessages(library(combinat)))
    library(gtools)

    OPTestOne <- unlist(unique(permn(v, paste0, collapse = "")))
    all.equal(sort(apply(multiPerm, 1, paste, collapse = "")),
              sort(OPTestOne))

    OPTestTwo <- unique(permutations(8, 8, v, set = FALSE))
    all.equal(OPTestTwo, multiPerm)

    #' Here are some benchmarks:

    library(microbenchmark)
    options(digits = 4)

    microbenchmark(OP_One = unique(permn(v, paste0, collapse = "")),
                   Arnge = arrangements::permutations(1:4, freq = rep(2, 4)),
                   OP_Two = unique(permutations(8, 8, v, set = FALSE)),
                   times = 5, unit = "relative")

    #' Finding permutations of multisets choose _m_ is no problem either.

    microbenchmark(
        algos = permuteGeneral(1:4, m = 11, freqs = rep(4, 4)),
        arnge = arrangements::permutations(1:4, 11, freq = rep(4, 4))
    )

    #' The OP's guess at how many permutations (525,525) for this latter example is not correct. Finding this is a [little more involved](https://github.com/jwood000/RcppAlgos/blob/2c2a6bfd9ea56488a718addb4aa30df3d9a86ee1/src/CombinatoricsContainer.cpp#L622) than the one liner offered.

    permuteCount(table(rep(1:4, each = 4)), m = 11)

    arrangements::npermutations(1:4, 11, freq = rep(4, 4))

    #' And just to show that these are all unique:

    length(
        unique(
            apply(
                permuteGeneral(1:4, m = 11, freqs = rep(4, 4)),
                MARGIN = 1, paste, collapse = ""
            )
        )
    )

    #' For more information regarding problems like this in R, I wrote a [thorough overview](https://stackoverflow.com/a/47983855/4408538) to the question: [How to generate permutations or combinations of object in R?](https://stackoverflow.com/q/22569176/4408538) by @RandyLai (author of `arrangements`)
}, advertise = FALSE, venue = "r", html_preview = FALSE, wd = ".")