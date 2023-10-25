compositionsGeneral <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, target = NULL, weak = FALSE,
    lower = NULL, upper = NULL, nThreads = NULL, tolerance = NULL
) {
    stopifnot(is.numeric(v))
    UseMethod("compositionsGeneral")
}

compositionsGeneral.default <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, target = NULL, weak = FALSE,
    lower = NULL, upper = NULL, nThreads = NULL, tolerance = NULL
) {
    return(.Call(`_RcppAlgos_CombinatoricsCnstrt`, v, m, repetition,
                 freqs, lower, upper, "sum", "==", GetTarget(v, target),
                 FALSE, FALSE, FALSE, nThreads, pkgEnv$nThreads, tolerance,
                 TRUE, weak))
}

compositionsGeneral.table <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, target = NULL, weak = FALSE,
    lower = NULL, upper = NULL, nThreads = NULL, tolerance = NULL
) {
    clean <- ResolveVFreqs(v, freqs)
    return(.Call(
        `_RcppAlgos_CombinatoricsCnstrt`, clean$v, m, repetition, clean$freqs,
        lower, upper, "sum", "==", GetTarget(clean$v, target), FALSE, FALSE,
        FALSE, nThreads, pkgEnv$nThreads, tolerance, TRUE, weak
    ))
}
