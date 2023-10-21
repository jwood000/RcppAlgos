permuteIter <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, constraintFun = NULL,
    comparisonFun = NULL, limitConstraints = NULL, keepResults = NULL,
    FUN = NULL, Parallel = FALSE, nThreads = NULL, tolerance = NULL,
    FUN.VALUE = NULL
) {
    UseMethod("permuteIter")
}

permuteIter.default <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, constraintFun = NULL,
    comparisonFun = NULL, limitConstraints = NULL, keepResults = NULL,
    FUN = NULL, Parallel = FALSE, nThreads = NULL, tolerance = NULL,
    FUN.VALUE = NULL
) {
    ComboPermIter(
        v, m, repetition, freqs, constraintFun, comparisonFun,
        limitConstraints, keepResults, FUN, Parallel, nThreads,
        tolerance, FUN.VALUE, FALSE
    )
}

permuteIter.table <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, constraintFun = NULL,
    comparisonFun = NULL, limitConstraints = NULL, keepResults = NULL,
    FUN = NULL, Parallel = FALSE, nThreads = NULL, tolerance = NULL,
    FUN.VALUE = NULL
) {
    clean <- ResolveVFreqs(v, freqs)
    ComboPermIter(
        clean$v, m, repetition, clean$freqs, constraintFun,
        comparisonFun, limitConstraints, keepResults, FUN, Parallel,
        nThreads, tolerance, FUN.VALUE, FALSE
    )
}

permuteIter.list <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, constraintFun = NULL,
    comparisonFun = NULL, limitConstraints = NULL, keepResults = NULL,
    FUN = NULL, Parallel = FALSE, nThreads = NULL, tolerance = NULL,
    FUN.VALUE = NULL
) {
    ComboPermIter(
        seq_along(v), m, repetition, freqs, constraintFun, comparisonFun,
        limitConstraints, keepResults, FUN = function(x) v[x], Parallel,
        nThreads, tolerance, FUN.VALUE = NULL, FALSE
    )
}
