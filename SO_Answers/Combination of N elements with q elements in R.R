reprex::reprex({
    #' These are permutations of multisets. Currently there are 4 packages on CRAN that can handle such tasks. They are: `partitions`, `multicool`, `arrangements` (which as @StéphaneLaurent points out is a replacement for `iterpc`), and `RcppAlgos` (I am the author).
    #'
    #' Both `arrangements` and `RcppAlgos` produce results in lexicographical order, producing a `matrix` where each row is a different permutation. Both of these packages offer access to memory efficient iterators for attacking larger cases.

    RcppAlgos::permuteGeneral(3, freqs = c(1, 2, 1))

    ## Or using S3 method for class 'table'
    RcppAlgos::permuteGeneral(table(c(1, 2, 2, 3)))

    #' The package `multicool` produces permutations row by row in a matrix in coolex order which is similar to [colexicographical order](https://en.wikipedia.org/wiki/Colexicographical_order). This package also offers iterators via `nextPerm`.

    multicool::allPerm(multicool::initMC(c(1, 2, 2, 3)))

    #' Finally, `partitions::multiset` outputs a `partitions` objects whereby each column represents a new permutation. It is based off of Knuth's algorithm as outlined in the [The Art of Computer Programming](https://en.wikipedia.org/wiki/The_Art_of_Computer_Programming).

    partitions::multiset(c(1, 2, 2, 3))

    #' Here are some benchmarks:

    library(microbenchmark)
    options(digits = 4)

    microbenchmark(
        arrangements = arrangements::permutations(x = 0:3, freq = 5:2),
        multicool    = multicool::allPerm(multicool::initMC(rep(0:3, times = 5:2))),
        partitions   = partitions::multiset(rep(0:3, times = 5:2)),
        RcppAlgos    = RcppAlgos::permuteGeneral(0:3, freqs = 5:2),
        unit = "relative",
        times = 40
    )

    #' For more information regarding problems like this in `R`, I wrote an [extensive overview](https://stackoverflow.com/a/47983855/4408538) to the question : [How to generate permutations or combinations of object in R?](https://stackoverflow.com/q/22569176/4408538).
}, advertise = FALSE, venue = "r", html_preview = FALSE, wd = ".")
