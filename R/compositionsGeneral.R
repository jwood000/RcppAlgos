compositionsGeneral <- function(v, m = NULL, ...) {
    stopifnot(is.numeric(v))
    UseMethod("compositionsGeneral")
}

compositionsGeneral.default <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, target = NULL, weak = FALSE,
    lower = NULL, upper = NULL, nThreads = NULL, tolerance = NULL, ...
) {
    return(.Call(`_RcppAlgos_CombinatoricsCnstrt`, v, m, repetition,
                 freqs, lower, upper, "sum", "==", GetTarget(v, target),
                 FALSE, FALSE, FALSE, nThreads, pkgEnv$nThreads, tolerance,
                 TRUE, weak))
}

compositionsGeneral.table <- function(
    v, m = NULL, target = NULL, weak = FALSE, lower = NULL,
    upper = NULL, nThreads = NULL, tolerance = NULL, ...
) {
    clean <- ResolveVFreqs(v)
    return(.Call(
        `_RcppAlgos_CombinatoricsCnstrt`, clean$v, m, FALSE, clean$freqs,
        lower, upper, "sum", "==", GetTarget(clean$v, target), FALSE, FALSE,
        FALSE, nThreads, pkgEnv$nThreads, tolerance, TRUE, weak
    ))
}
