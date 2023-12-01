comboIter <- function(v, m = NULL, ...) {
    UseMethod("comboIter")
}

comboIter.integer <-
comboIter.numeric <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, constraintFun = NULL,
    comparisonFun = NULL, limitConstraints = NULL, keepResults = NULL,
    FUN = NULL, Parallel = FALSE, nThreads = NULL, tolerance = NULL,
    FUN.VALUE = NULL, ...
) {
    ComboPermuteIter(
        v, m, repetition, freqs, constraintFun, comparisonFun,
        limitConstraints, keepResults, FUN, Parallel, nThreads,
        tolerance, FUN.VALUE, TRUE
    )
}

comboIter.factor <-
comboIter.logical <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, FUN = NULL,
    Parallel = FALSE, nThreads = NULL, FUN.VALUE = NULL, ...
) {
    ComboPermuteIter(
        v, m, repetition, freqs, NULL, NULL, NULL, NULL,
        FUN, Parallel, nThreads, NULL, FUN.VALUE, TRUE
    )
}

comboIter.default <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL,
    FUN = NULL, FUN.VALUE = NULL, ...
) {
    ComboPermuteIter(
        v, m, repetition, freqs, NULL, NULL, NULL, NULL,
        FUN, FALSE, NULL, NULL, FUN.VALUE, TRUE
    )
}

comboIter.table <- function(
    v, m = NULL, constraintFun = NULL, comparisonFun = NULL,
    limitConstraints = NULL, keepResults = NULL, FUN = NULL,
    Parallel = FALSE, nThreads = NULL, tolerance = NULL, FUN.VALUE = NULL, ...
) {
    clean <- ResolveVFreqs(v)
    ComboPermuteIter(
        clean$v, m, FALSE, clean$freqs, constraintFun,
        comparisonFun, limitConstraints, keepResults, FUN,
        Parallel, nThreads, tolerance, FUN.VALUE, TRUE
    )
}

comboIter.list <- function(v, m = NULL, repetition = FALSE,
                           freqs = NULL, ...) {
    ComboPermuteIter(
        seq_along(v), m, repetition, freqs, NULL, NULL, NULL, NULL,
        FUN = function(x) v[x], FALSE, NULL, NULL, NULL, TRUE
    )
}
