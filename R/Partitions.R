partitionsGeneral <- function(v, m = NULL, repetition = FALSE,
                              freqs = NULL, target = NULL, lower = NULL,
                              upper = NULL, nThreads = NULL,
                              tolerance = NULL) {

    stopifnot(is.numeric(v))
    return(.Call(`_RcppAlgos_CombinatoricsCnstrt`, v, m, repetition,
                 freqs, lower, upper, "sum", "==", GetTarget(v, target),
                 TRUE, FALSE, FALSE, nThreads, pkgEnv$nThreads, tolerance,
                 FALSE, FALSE))
}

partitionsCount <- function(v, m = NULL, repetition = FALSE,
                            freqs = NULL, target = NULL) {

    stopifnot(is.numeric(v))
    return(.Call(`_RcppAlgos_PartitionsCount`, GetTarget(v, target),
                 v, m, repetition, freqs, "==", NULL, NULL, FALSE,
                 FALSE, FALSE, FALSE))
}

partitionsDesign <- function(v, m = NULL, repetition = FALSE,
                             freqs = NULL, target = NULL,
                             showDesign = FALSE) {

    stopifnot(is.numeric(v))
    return(.Call(`_RcppAlgos_PartitionsCount`, GetTarget(v, target),
                 v, m, repetition, freqs, "==", NULL, NULL, TRUE,
                 showDesign, FALSE, FALSE))
}

partitionsRank <- function(..., v, repetition = FALSE,
                           freqs = NULL, target = NULL) {
    GetRankPart(..., v = v, repetition = repetition, freqs = freqs,
                target = target, IsComposition = FALSE, weak = FALSE)
}

partitionsSample <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                             target = NULL, n = NULL, sampleVec = NULL,
                             seed = NULL, nThreads = NULL,
                             namedSample = FALSE) {

    stopifnot(is.numeric(v))

    if (!is.null(seed)) {
        set.seed(seed)
    }

    return(.Call(`_RcppAlgos_SamplePartitions`, v, m, repetition, freqs,
                 sampleVec, seed, n, sample, FALSE, nThreads,
                 pkgEnv$nThreads, namedSample, "==",
                 GetTarget(v, target), NULL, new.env(), FALSE, FALSE))
}

partitionsIter <- function(v, m = NULL, repetition = FALSE,
                           freqs = NULL, target = NULL,
                           nThreads = NULL, tolerance = NULL) {

    stopifnot(is.numeric(v))
    InitVals <- .Call(`_RcppAlgos_GetClassVals`, v, m, repetition,
                      freqs, TRUE, NULL, nThreads, pkgEnv$nThreads,
                      TRUE, FALSE, FALSE, NULL, NULL, NULL)

    new("Partitions", InitVals, FALSE, "sum", "==",
        GetTarget(v, target), FALSE, tolerance, is.null(m))
}
