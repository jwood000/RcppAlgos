comboGeneral <- function(v, m = NULL, ...) {
    UseMethod("comboGeneral")
}

comboGeneral.integer <-
comboGeneral.numeric <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, lower = NULL, upper = NULL,
    constraintFun = NULL, comparisonFun = NULL, limitConstraints = NULL,
    keepResults = NULL, FUN = NULL, Parallel = FALSE, nThreads = NULL,
    tolerance = NULL, FUN.VALUE = NULL, ...
) {
    ComboPermuteGen(
        v, m, repetition, freqs, lower, upper, constraintFun,
        comparisonFun, limitConstraints, keepResults, FUN, Parallel,
        nThreads, tolerance, FUN.VALUE, TRUE, ...
    )
}

comboGeneral.factor <-
comboGeneral.logical <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, lower = NULL, upper = NULL,
    FUN = NULL, Parallel = FALSE, nThreads = NULL, FUN.VALUE = NULL, ...
) {
    ComboPermuteGen(v, m, repetition, freqs, lower, upper, NULL, NULL, NULL,
                    NULL, FUN, Parallel, nThreads, NULL, FUN.VALUE, TRUE, ...)
}

comboGeneral.default <- function(v, m = NULL, repetition = FALSE,
                                 freqs = NULL, lower = NULL, upper = NULL,
                                 FUN = NULL, FUN.VALUE = NULL, ...) {
    ComboPermuteGen(v, m, repetition, freqs, lower, upper, NULL, NULL, NULL,
                    NULL, FUN, FALSE, NULL, NULL, FUN.VALUE, TRUE, ...)
}

comboGeneral.table <- function(
    v, m = NULL, lower = NULL, upper = NULL, constraintFun = NULL,
    comparisonFun = NULL, limitConstraints = NULL, keepResults = NULL,
    FUN = NULL, Parallel = FALSE, nThreads = NULL, tolerance = NULL,
    FUN.VALUE = NULL, ...
) {
    clean <- ResolveVFreqs(v)
    ComboPermuteGen(
        clean$v, m, FALSE, clean$freqs, lower, upper, constraintFun,
        comparisonFun, limitConstraints, keepResults, FUN, Parallel,
        nThreads, tolerance, FUN.VALUE, TRUE, ...
    )
}

comboGeneral.list <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                              lower = NULL, upper = NULL, ...) {
    ComboPermuteGen(
        seq_along(v), m, repetition, freqs, lower, upper, NULL, NULL, NULL,
        NULL, FUN = function(x) v[x], FALSE, NULL, NULL, NULL, TRUE, ...
    )
}
