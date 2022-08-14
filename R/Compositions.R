compositionsGeneral <- function(v, m = NULL, repetition = FALSE,
                                freqs = NULL, target = NULL, weak = FALSE,
                                lower = NULL, upper = NULL,
                                nThreads = NULL, tolerance = NULL) {

    stopifnot(is.numeric(v))
    return(.Call(`_RcppAlgos_CombinatoricsCnstrt`, v, m, repetition,
                 freqs, lower, upper, "sum", "==", GetTarget(v, target),
                 FALSE, FALSE, FALSE, nThreads, pkgEnv$nThreads, tolerance,
                 TRUE, weak))
}

compositionsCount <- function(v, m = NULL, repetition = FALSE,
                              freqs = NULL, target = NULL, weak = FALSE) {

    stopifnot(is.numeric(v))
    return(.Call(`_RcppAlgos_PartitionsCount`, GetTarget(v, target),
                 v, m, repetition, freqs, "==", NULL, NULL, FALSE,
                 FALSE, TRUE, weak))
}

compositionsDesign <- function(v, m = NULL, repetition = FALSE,
                               freqs = NULL, target = NULL,
                               weak = FALSE, showDesign = FALSE) {

    stopifnot(is.numeric(v))
    return(.Call(`_RcppAlgos_PartitionsCount`, GetTarget(v, target),
                 v, m, repetition, freqs, "==", NULL, NULL, TRUE,
                 showDesign, TRUE, weak))
}

compositionsRank <- function(..., v, repetition = FALSE,
                             freqs = NULL, target = NULL, weak = FALSE) {
    GetRankPart(..., v = v, repetition = repetition, freqs = freqs,
                target = target, IsComposition = TRUE, weak = weak)
}

compositionsSample <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                               target = NULL, weak = FALSE, n = NULL,
                               sampleVec = NULL, seed = NULL,
                               nThreads = NULL, namedSample = FALSE) {

    stopifnot(is.numeric(v))

    if (!is.null(seed)) {
        set.seed(seed)
    }

    return(.Call(`_RcppAlgos_SamplePartitions`, v, m, repetition, freqs,
                 sampleVec, seed, n, sample, FALSE, nThreads,
                 pkgEnv$nThreads, namedSample, "==",
                 GetTarget(v, target), NULL, new.env(), TRUE, weak))
}

compositionsIter <- function(v, m = NULL, repetition = FALSE,
                             freqs = NULL, target = NULL, weak = FALSE,
                             nThreads = NULL, tolerance = NULL) {

    stopifnot(is.numeric(v))
    InitVals <- .Call(`_RcppAlgos_GetClassVals`, v, m, repetition,
                      freqs, TRUE, NULL, nThreads, pkgEnv$nThreads,
                      TRUE, TRUE, weak)

    new("Partitions", InitVals, FALSE, "sum", "==",
        GetTarget(v, target), FALSE, tolerance, is.null(m))
}
