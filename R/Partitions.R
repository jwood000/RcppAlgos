GetTarget <- function(v, target) {
    if (is.null(target)) {target = max(v, na.rm = TRUE)}
    return(target)
}

partitionsGeneral <- function(v, m = NULL, repetition = FALSE,
                              freqs = NULL, target = NULL, lower = NULL,
                              upper = NULL, nThreads = NULL,
                              tolerance = NULL) {

    return(.Call(`_RcppAlgos_CombinatoricsCnstrt`, v, m, repetition, freqs,
                 lower, upper, "sum", "==", GetTarget(v, target), TRUE,
                 FALSE, FALSE, nThreads, pkgEnv$nThreads, tolerance))
}

partitionsCount <- function(v, m = NULL, repetition = FALSE,
                            freqs = NULL, target = NULL) {
    return(.Call(`_RcppAlgos_PartitionsCount`, GetTarget(v, target), v, m,
                 repetition, freqs, "==", NULL, NULL, FALSE, FALSE))
}

partitionsDesign <- function(v, m = NULL, repetition = FALSE,
                             freqs = NULL, target = NULL,
                             showDesign = FALSE) {

    return(.Call(`_RcppAlgos_PartitionsCount`, GetTarget(v, target), v, m,
                 repetition, freqs, "==", NULL, NULL, TRUE, showDesign))
}

partitionsSample <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                             target = NULL, n = NULL, sampleVec = NULL,
                             seed = NULL, nThreads = NULL,
                             namedSample = FALSE) {

    if (!is.null(seed)) {
        set.seed(seed)
    }

    return(.Call(`_RcppAlgos_SamplePartitions`, v, m, repetition, freqs,
                 sampleVec, seed, n, sample, FALSE, nThreads,
                 pkgEnv$nThreads, namedSample, "==",
                 GetTarget(v, target), NULL, new.env()))
}

partitionsIter <- function(v, m = NULL, repetition = FALSE,
                           freqs = NULL, target = NULL,
                           nThreads = NULL, tolerance = NULL) {

    InitVals <- .Call(`_RcppAlgos_GetClassVals`, v, m, repetition, freqs,
                      TRUE, NULL, nThreads, pkgEnv$nThreads, TRUE)

    new("Partitions", InitVals, FALSE, "sum", "==",
        GetTarget(v, target), FALSE, tolerance, is.null(m))
}
