compositionsIter <- function(v, m = NULL, ...) {
    stopifnot(is.numeric(v))
    UseMethod("compositionsIter")
}

compositionsIter.default <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, target = NULL,
    weak = FALSE, nThreads = NULL, tolerance = NULL, ...
) {

    InitVals <- .Call(`_RcppAlgos_GetClassVals`, v, m, repetition,
                      freqs, FALSE, NULL, nThreads, pkgEnv$nThreads,
                      TRUE, TRUE, weak, NULL, NULL, NULL)

    new("Partitions", InitVals, FALSE, "sum", "==",
        GetTarget(v, target), FALSE, tolerance, is.null(m))
}

compositionsIter.table <- function(
    v, m = NULL, target = NULL, weak = FALSE,
    nThreads = NULL, tolerance = NULL, ...
) {

    clean <- ResolveVFreqs(v)
    InitVals <- .Call(`_RcppAlgos_GetClassVals`, clean$v, m, FALSE,
                      clean$freqs, FALSE, NULL, nThreads, pkgEnv$nThreads,
                      TRUE, TRUE, weak, NULL, NULL, NULL)

    new("Partitions", InitVals, FALSE, "sum", "==",
        GetTarget(clean$v, target), FALSE, tolerance, is.null(m))
}
