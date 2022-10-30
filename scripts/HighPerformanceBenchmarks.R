reprex::reprex({
    #' This document serves as an overview for measuring the performance of `RcppAlgos` against other tools for generating combinations, permutations, and partitions. This stackoverflow post: [How to generate permutations or combinations of object in R?](<https://stackoverflow.com/a/47983855/4408538>) has some benchmarks. You will note that the examples in that post were relatively small. The benchmarks below will focus on larger examples where performance really matters and for this reason we only consider the packages [arrangements](<https://cran.r-project.org/package=arrangements>), [partitions](<https://cran.r-project.org/package=partitions>), and [RcppAlgos](<https://cran.r-project.org/package=RcppAlgos>).
    #'
    #' ## Setup Information
    #'
    #' For the benchmarks below, we used a `2022 Macbook Air Apple M2 24 GB` machine. We also tested on a Windows and Linux machine with similar specs and obtained similar results.
    #'

    library(RcppAlgos)
    library(partitions)
    library(arrangements)
    library(microbenchmark)

    options(digits = 4)
    options(width = 90)

    pertinent_output <- capture.output(sessionInfo())
    cat(paste(pertinent_output[1:3], collapse = "\n"))

    pkgs <- c("RcppAlgos", "arrangements", "partitions", "microbenchmark")
    sapply(pkgs, packageVersion, simplify = FALSE)

    numThreads <- min(as.integer(RcppAlgos::stdThreadMax() / 2), 6)
    numThreads

    #'
    #' ## Combinations
    #'
    #' ### Combinations - Distinct
    #'

    set.seed(13)
    v1 <- sort(sample(100, 30))
    m <- 21
    t1 <- comboGeneral(v1, m, Parallel = T)
    t2 <- combinations(v1, m)
    stopifnot(identical(t1, t2))
    dim(t1)
    rm(t1, t2)
    invisible(gc())
    microbenchmark(cbRcppAlgosPar = comboGeneral(v1, m, nThreads = numThreads),
                   cbRcppAlgosSer = comboGeneral(v1, m),
            	   cbArrangements = combinations(v1, m),
            	   times = 15, unit = "relative")
    #'
    #'
    #' ### Combinations - Repetition
    #'

    v2 <- v1[1:10]
    m <- 20
    t1 <- comboGeneral(v2, m, repetition = TRUE, nThreads = numThreads)
    t2 <- combinations(v2, m, replace = TRUE)
    stopifnot(identical(t1, t2))
    dim(t1)
    rm(t1, t2)
    invisible(gc())
    microbenchmark(cbRcppAlgosPar = comboGeneral(v2, m, TRUE, nThreads = numThreads),
                   cbRcppAlgosSer = comboGeneral(v2, m, TRUE),
                   cbArrangements = combinations(v2, m, replace = TRUE),
                   times = 15, unit = "relative")
    #'
    #'
    #' ### Combinations - Multisets
    #'

    myFreqs <- c(2, 4, 4, 5, 3, 2, 2, 2, 3, 4, 1, 4, 2, 5)
    v3 <- as.integer(c(1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610))
    t1 <- comboGeneral(v3, 20, freqs = myFreqs, nThreads = numThreads)
    t2 <- combinations(freq = myFreqs, k = 20, x = v3)
    stopifnot(identical(t1, t2))
    dim(t1)
    rm(t1, t2)
    invisible(gc())
    microbenchmark(cbRcppAlgosPar = comboGeneral(v3, 20, freqs = myFreqs, nThreads = numThreads),
                   cbRcppAlgosSer = comboGeneral(v3, 20, freqs = myFreqs),
                   cbArrangements = combinations(freq = myFreqs, k = 20, x = v3),
                   times = 10, unit = "relative")
    #'
    #'
    #' ## Permutations
    #'
    #' ### Permutations - Distinct
    #'

    v4 <- as.integer(c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59))
    t1 <- permuteGeneral(v4, 6, nThreads = numThreads)
    t2 <- permutations(v4, 6)
    stopifnot(identical(t1, t2))
    dim(t1)
    rm(t1, t2)
    invisible(gc())
    microbenchmark(cbRcppAlgosPar = permuteGeneral(v4, 6, nThreads = numThreads),
                   cbRcppAlgosSer = permuteGeneral(v4, 6),
                   cbArrangements = permutations(v4, 6),
                   times = 15, unit = "relative")


    ## Indexing permutation example with the partitions package
    t1 <- permuteGeneral(11, nThreads = 4)
    t2 <- permutations(11)
    t3 <- perms(11)

    dim(t1)

    stopifnot(identical(t1, t2), identical(t1, t(as.matrix(t3))))
    rm(t1, t2, t3)
    invisible(gc())

    microbenchmark(cbRcppAlgosPar = permuteGeneral(11, nThreads = 4),
                   cbRcppAlgosSer = permuteGeneral(11),
                   cbArrangements = permutations(11),
                   cbPartitions   = perms(11),
                   times = 5, unit = "relative")
    #'
    #'
    #' ### Permutations - Repetition
    #'

    v5 <- v3[1:5]
    t1 <- permuteGeneral(v5, 10, repetition = TRUE, nThreads = numThreads)
    t2 <- permutations(v5, 10, replace = TRUE)
    stopifnot(identical(t1, t2))
    dim(t1)
    rm(t1, t2)
    invisible(gc())
    microbenchmark(cbRcppAlgosPar = permuteGeneral(v5, 10, TRUE, nThreads = numThreads),
                   cbRcppAlgosSer = permuteGeneral(v5, 10, TRUE),
                   cbArrangements = permutations(x = v5, k = 10, replace = TRUE),
                   times = 10, unit = "relative")
    #'
    #'
    #' ### Permutations - Multisets
    #'

    v6 <- sort(runif(12))
    t1 <- permuteGeneral(v6, 7, freqs = rep(1:3, 4), nThreads = numThreads)
    t2 <- permutations(freq = rep(1:3, 4), k = 7, x = v6)
    stopifnot(identical(t1, t2))
    dim(t1)
    rm(t1, t2)
    invisible(gc())
    microbenchmark(cbRcppAlgosPar = permuteGeneral(v6, 7, freqs = rep(1:3, 4), nThreads = numThreads),
                   cbRcppAlgosSer = permuteGeneral(v6, 7, freqs = rep(1:3, 4)),
                   cbArrangements = permutations(freq = rep(1:3, 4), k = 7, x = v6),
                   times = 10, unit = "relative")
    #'
    #'
    #' ## Partitions
    #'
    #' ### Partitions - Distinct
    #'
    #' #### All Distinct Partitions
    #'

    t1 <- comboGeneral(0:140, freqs=c(140, rep(1, 140)),
    				   constraintFun = "sum", comparisonFun = "==",
    				   limitConstraints = 140)
    t2 <- partitions(140, distinct = TRUE)
    t3 <- diffparts(140)

    # Each package has different output formats... we only examine dimensions
    # and that each result is a partition of 140
    stopifnot(identical(dim(t1), dim(t2)), identical(dim(t1), dim(t(t3))),
                        all(rowSums(t1) == 140), all(rowSums(t2) == 140),
                        all(colSums(t3) == 140))
    dim(t1)
    rm(t1, t2, t3)
    invisible(gc())
    microbenchmark(cbRcppAlgosPar = partitionsGeneral(0:140, freqs=c(140, rep(1, 140)), nThreads = numThreads),
                   cbRcppAlgosSer = partitionsGeneral(0:140, freqs=c(140, rep(1, 140))),
                   cbArrangements = partitions(140, distinct = TRUE),
                   cbPartitions   = diffparts(140),
                   times = 10, unit = "relative")
    #'
    #'
    #' #### Restricted Distinct Partitions
    #'

    t1 <- comboGeneral(160, 10,
    				   constraintFun = "sum", comparisonFun = "==",
    				   limitConstraints = 160)
    t2 <- partitions(160, 10, distinct = TRUE)
    stopifnot(identical(t1, t2))
    dim(t1)
    rm(t1, t2)
    invisible(gc())
    microbenchmark(cbRcppAlgosPar = partitionsGeneral(160, 10, nThreads = numThreads),
                   cbRcppAlgosSer = partitionsGeneral(160, 10),
                   cbArrangements = partitions(160, 10, distinct = TRUE),
                   times = 10, unit = "relative")
    #'
    #'
    #' ### Partitions - Repetition
    #'
    #' #### All Partitions
    #'

    t1 <- comboGeneral(0:65, repetition = TRUE, constraintFun = "sum",
     		           comparisonFun = "==", limitConstraints = 65)
    t2 <- partitions(65)
    t3 <- parts(65)

    # Each package has different output formats... we only examine dimensions
    # and that each result is a partition of 65
    stopifnot(identical(dim(t1), dim(t2)), identical(dim(t1), dim(t(t3))),
              all(rowSums(t1) == 65), all(rowSums(t2) == 65),
              all(colSums(t3) == 65))
    dim(t1)
    rm(t1, t2, t3)
    invisible(gc())
    microbenchmark(cbRcppAlgosPar = partitionsGeneral(0:65, repetition = TRUE,
                                                      nThreads = numThreads),
                   cbRcppAlgosSer = partitionsGeneral(0:65, repetition = TRUE),
     			   cbArrangements = partitions(65),
     			   cbPartitions   = parts(65),
     			   times = 20, unit = "relative")
    #'
    #'
    #' #### Restricted Partitions
    #'

    t1 <- comboGeneral(100, 15, TRUE, constraintFun = "sum",
                       comparisonFun = "==", limitConstraints = 100)
    t2 <- partitions(100, 15)
    stopifnot(identical(t1, t2))
    dim(t1)
    rm(t1, t2)

    # This takes a really long time... not because of restrictedparts,
    # but because apply is not that fast. This transformation is
    # needed for proper comparisons. As a result, we will compare
    # a smaller example
    # t3 <- t(apply(as.matrix(restrictedparts(100, 15, include.zero = F)), 2, sort))
    t3 <- t(apply(as.matrix(restrictedparts(50, 15, include.zero = F)), 2, sort))
    stopifnot(identical(partitions(50, 15), t3))
    rm(t3)
    invisible(gc())
    microbenchmark(cbRcppAlgosPar = partitionsGeneral(100, 15, TRUE,
                                                      nThreads = numThreads),
                   cbRcppAlgosSer = partitionsGeneral(100, 15, TRUE),
    		       cbArrangements = partitions(100, 15),
    		       cbPartitions   = restrictedparts(100, 15,
    		                                        include.zero = FALSE),
    		       times = 10, unit = "relative")
    #'
    #'
    #' ### Partitions - Multisets
    #'
    #' Currenlty, `RcppAlgos` is the only package capable of efficiently generating partitions of multisets. Therefore, we will only time `RcppAlgos` and use this as a reference for future improvements.
    #'

    t1 <- comboGeneral(120, 10, freqs=rep(1:8, 15),
    				   constraintFun = "sum", comparisonFun = "==",
    				   limitConstraints = 120)
    dim(t1)
    stopifnot(all(rowSums(t1) == 120))
    microbenchmark(cbRcppAlgos = partitionsGeneral(120, 10, freqs=rep(1:8, 15)),
                   times = 10)

    ### In RcppAlgos 2.3.6 - 2.4.3
    #'
    #'
    #' ## Compositions
    #'
    #' ### Compositions - Repetition
    #'
    #' #### All Compositions (Small case)
    #'

    t1 <- compositionsGeneral(0:15, repetition = TRUE)
    t2 <- arrangements::compositions(15)
    t3 <- partitions::compositions(15)

    # Each package has different output formats... we only examine dimensions
    # and that each result is a partition of 15
    stopifnot(identical(dim(t1), dim(t2)), identical(dim(t1), dim(t(t3))),
              all(rowSums(t1) == 15), all(rowSums(t2) == 15),
              all(colSums(t3) == 15))
    dim(t1)
    rm(t1, t2, t3)
    invisible(gc())
    microbenchmark(cbRcppAlgosSer = compositionsGeneral(0:15, repetition = TRUE),
     			   cbArrangements = arrangements::compositions(15),
     			   cbPartitions   = partitions::compositions(15),
     			   times = 20, unit = "relative")
    #'
    #'
    #' For the next two examples, we will exclude the `partitions` package for efficiency reasons.
    #'
    #' #### All Compositions (Larger case)
    #'

    t1 <- compositionsGeneral(0:23, repetition = TRUE)
    t2 <- arrangements::compositions(23)

    # Each package has different output formats... we only examine dimensions
    # and that each result is a partition of 23
    stopifnot(identical(dim(t1), dim(t2)), all(rowSums(t1) == 23),
              all(rowSums(t2) == 23))
    dim(t1)
    rm(t1, t2)
    invisible(gc())
    microbenchmark(cbRcppAlgosPar = compositionsGeneral(0:23, repetition = TRUE,
                                                        nThreads = numThreads),
                   cbRcppAlgosSer = compositionsGeneral(0:23, repetition = TRUE),
     			   cbArrangements = arrangements::compositions(23),
     			   times = 20, unit = "relative")
    #'
    #'
    #' #### Restricted Compositions
    #'

    t1 <- compositionsGeneral(30, 10, repetition = TRUE)
    t2 <- arrangements::compositions(30, 10)

    stopifnot(identical(t1, t2), all(rowSums(t1) == 30))
    dim(t1)
    rm(t1, t2)
    invisible(gc())
    microbenchmark(cbRcppAlgosPar = compositionsGeneral(30, 10, repetition = TRUE,
                                                        nThreads = numThreads),
                   cbRcppAlgosSer = compositionsGeneral(30, 10, repetition = TRUE),
     			   cbArrangements = arrangements::compositions(30, 10),
     			   times = 20, unit = "relative")
    #'
    #'
    #' ## Iterators
    #'
    #' We will show one example from each category to demonstrate the efficiency of the iterators in `RcppAlgos`. The results are similar for the rest of the cases not shown.
    #'
    #' ### Combinations
    #'

    pkg_arrangements <- function(n, total) {
        a <- icombinations(n, as.integer(n / 2))
        for (i in 1:total) a$getnext()
    }

    pkg_RcppAlgos <- function(n, total) {
        a <- comboIter(n, as.integer(n / 2))
        for (i in 1:total) a@nextIter()
    }

    total <- comboCount(18, 9)
    total

    microbenchmark(cbRcppAlgos    = pkg_RcppAlgos(18, total),
                   cbArrangements = pkg_arrangements(18, total),
                   times = 15, unit = "relative")
    #'
    #'
    #' ### Permutations
    #'

    pkg_arrangements <- function(n, total) {
        a <- ipermutations(n)
        for (i in 1:total) a$getnext()
    }

    pkg_RcppAlgos <- function(n, total) {
        a <- permuteIter(n)
        for (i in 1:total) a@nextIter()
    }

    total <- permuteCount(8)
    total

    microbenchmark(cbRcppAlgos    = pkg_RcppAlgos(8, total),
                   cbArrangements = pkg_arrangements(8, total),
                   times = 15, unit = "relative")
    #'
    #'
    #' ### Partitions
    #'

    pkg_partitions <- function(n, total) {
        a <- firstpart(n)
        for (i in 1:(total - 1)) a <- nextpart(a)
    }

    pkg_arrangements <- function(n, total) {
        a <- ipartitions(n)
        for (i in 1:total) a$getnext()
    }

    pkg_RcppAlgos <- function(n, total) {
        a <- partitionsIter(0:n, repetition = TRUE)
        for (i in 1:total) a@nextIter()
    }

    total <- partitionsCount(0:40, repetition = TRUE)
    total

    microbenchmark(cbRcppAlgos    = pkg_RcppAlgos(40, total),
                   cbArrangements = pkg_arrangements(40, total),
                   cbPartitions   = pkg_partitions(40, total),
                   times = 15, unit = "relative")
    #'
    #'
    #' ### Compositions
    #'

    pkg_partitions <- function(n, total) {
        a <- firstcomposition(n)
        for (i in 1:(total - 1)) a <- nextcomposition(a, FALSE)
    }

    pkg_arrangements <- function(n, total) {
        a <- icompositions(n)
        for (i in 1:total) a$getnext()
    }

    pkg_RcppAlgos <- function(n, total) {
        a <- compositionsIter(0:n, repetition = TRUE)
        for (i in 1:total) a@nextIter()
    }

    total <- compositionsCount(0:15, repetition = TRUE)
    total

    microbenchmark(cbRcppAlgos    = pkg_RcppAlgos(15, total),
                   cbArrangements = pkg_arrangements(15, total),
                   cbPartitions   = pkg_partitions(15, total),
                   times = 15, unit = "relative")
}, advertise = FALSE, venue = "r", html_preview = FALSE, wd = ".")