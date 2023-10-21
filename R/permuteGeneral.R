permuteGeneral <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, lower = NULL,
    upper = NULL, constraintFun = NULL, comparisonFun = NULL,
    limitConstraints = NULL, keepResults = NULL, FUN = NULL,
    Parallel = FALSE, nThreads = NULL, tolerance = NULL, FUN.VALUE = NULL
) {
    UseMethod("permuteGeneral")
}

permuteGeneral.default <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, lower = NULL,
    upper = NULL, constraintFun = NULL, comparisonFun = NULL,
    limitConstraints = NULL, keepResults = NULL, FUN = NULL,
    Parallel = FALSE, nThreads = NULL, tolerance = NULL, FUN.VALUE = NULL
) {
    ComboPermGen(
        v, m, repetition, freqs, lower, upper, constraintFun,
        comparisonFun, limitConstraints, keepResults, FUN, Parallel,
        nThreads, tolerance, FUN.VALUE, FALSE
    )
}

permuteGeneral.table <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, lower = NULL,
    upper = NULL, constraintFun = NULL, comparisonFun = NULL,
    limitConstraints = NULL, keepResults = NULL, FUN = NULL,
    Parallel = FALSE, nThreads = NULL, tolerance = NULL, FUN.VALUE = NULL
) {
    clean <- ResolveVFreqs(v, freqs)
    ComboPermGen(
        clean$v, m, repetition, clean$freqs, lower, upper, constraintFun,
        comparisonFun, limitConstraints, keepResults, FUN, Parallel,
        nThreads, tolerance, FUN.VALUE, FALSE
    )
}

permuteGeneral.list <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, lower = NULL,
    upper = NULL, constraintFun = NULL, comparisonFun = NULL,
    limitConstraints = NULL, keepResults = NULL, FUN = NULL,
    Parallel = FALSE, nThreads = NULL, tolerance = NULL, FUN.VALUE = NULL
) {
    ComboPermGen(
        seq_along(v), m, repetition, freqs, lower, upper, constraintFun,
        comparisonFun, limitConstraints, keepResults, FUN = function(x) {
            v[x]
        }, Parallel, nThreads, tolerance, FUN.VALUE = NULL, FALSE
    )
}
