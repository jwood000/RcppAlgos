partitionsIter <- function(v, m = NULL, repetition = FALSE,
                           freqs = NULL, target = NULL,
                           nThreads = NULL, tolerance = NULL) {

    stopifnot(is.numeric(v))
    UseMethod("partitionsIter")
}

partitionsIter.default <- function(v, m = NULL, repetition = FALSE,
                                   freqs = NULL, target = NULL,
                                   nThreads = NULL, tolerance = NULL) {

    InitVals <- .Call(`_RcppAlgos_GetClassVals`, v, m, repetition,
                      freqs, TRUE, NULL, nThreads, pkgEnv$nThreads,
                      TRUE, FALSE, FALSE, NULL, NULL, NULL)

    new("Partitions", InitVals, FALSE, "sum", "==",
        GetTarget(v, target), FALSE, tolerance, is.null(m))
}

partitionsIter.table <- function(v, m = NULL, repetition = FALSE,
                                 freqs = NULL, target = NULL,
                                 nThreads = NULL, tolerance = NULL) {

    clean <- ResolveVFreqs(v, freqs)
    InitVals <- .Call(`_RcppAlgos_GetClassVals`, clean$v, m, repetition,
                      clean$freqs, TRUE, NULL, nThreads, pkgEnv$nThreads,
                      TRUE, FALSE, FALSE, NULL, NULL, NULL)

    new("Partitions", InitVals, FALSE, "sum", "==",
        GetTarget(clean$v, target), FALSE, tolerance, is.null(m))
}
