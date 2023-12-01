partitionsGeneral <- function(v, m = NULL, ...) {
    stopifnot(is.numeric(v))
    UseMethod("partitionsGeneral")
}

partitionsGeneral.default <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, target = NULL,
    lower = NULL, upper = NULL, nThreads = NULL, tolerance = NULL, ...
) {
    return(.Call(`_RcppAlgos_CombinatoricsCnstrt`, v, m, repetition,
                 freqs, lower, upper, "sum", "==", GetTarget(v, target),
                 TRUE, FALSE, FALSE, nThreads, pkgEnv$nThreads, tolerance,
                 FALSE, FALSE))
}

partitionsGeneral.table <- function(
    v, m = NULL, target = NULL, lower = NULL,
    upper = NULL, nThreads = NULL, tolerance = NULL, ...
) {
    clean <- ResolveVFreqs(v)
    return(.Call(
        `_RcppAlgos_CombinatoricsCnstrt`, clean$v, m, FALSE, clean$freqs,
        lower, upper, "sum", "==", GetTarget(clean$v, target), TRUE, FALSE,
        FALSE, nThreads, pkgEnv$nThreads, tolerance, FALSE, FALSE
    ))
}
