permuteIter <- function(v, m = NULL, ...) {
    UseMethod("permuteIter")
}

permuteIter.integer <-
permuteIter.numeric <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, constraintFun = NULL,
    comparisonFun = NULL, limitConstraints = NULL, keepResults = NULL,
    FUN = NULL, Parallel = FALSE, nThreads = NULL, tolerance = NULL,
    FUN.VALUE = NULL, ...
) {
    ComboPermuteIter(
        v, m, repetition, freqs, constraintFun, comparisonFun,
        limitConstraints, keepResults, FUN, Parallel, nThreads,
        tolerance, FUN.VALUE, FALSE
    )
}

permuteIter.factor <-
permuteIter.logical <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, FUN = NULL,
    Parallel = FALSE, nThreads = NULL, FUN.VALUE = NULL, ...
) {
    ComboPermuteIter(
        v, m, repetition, freqs, NULL, NULL, NULL, NULL,
        FUN, Parallel, nThreads, NULL, FUN.VALUE, FALSE
    )
}

permuteIter.default <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL,
    FUN = NULL, FUN.VALUE = NULL, ...
) {
    ComboPermuteIter(
        v, m, repetition, freqs, NULL, NULL, NULL, NULL,
        FUN, FALSE, NULL, NULL, FUN.VALUE, FALSE
    )
}

permuteIter.table <- function(
    v, m = NULL, constraintFun = NULL, comparisonFun = NULL,
    limitConstraints = NULL, keepResults = NULL, FUN = NULL,
    Parallel = FALSE, nThreads = NULL, tolerance = NULL, FUN.VALUE = NULL, ...
) {
    clean <- ResolveVFreqs(v)
    ComboPermuteIter(
        clean$v, m, FALSE, clean$freqs, constraintFun,
        comparisonFun, limitConstraints, keepResults, FUN,
        Parallel, nThreads, tolerance, FUN.VALUE, FALSE
    )
}

permuteIter.list <- function(v, m = NULL, repetition = FALSE,
                             freqs = NULL, ...) {
    ComboPermuteIter(
        seq_along(v), m, repetition, freqs, NULL, NULL, NULL, NULL,
        FUN = function(x) v[x], FALSE, NULL, NULL, NULL, FALSE
    )
}
