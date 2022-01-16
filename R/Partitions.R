GetTarget <- function(v, target) {
    if (is.null(target)) {target = max(v)}
    return(target)
}

partitionsGeneral <- function(v, m = NULL, repetition = FALSE,
                              freqs = NULL, target = NULL, lower = NULL,
                              upper = NULL, nThreads = NULL,
                              tolerance = NULL) {

    return(.Call(Algos_CombinatoricsCnstrt, v, m, repetition, freqs,
                 lower, upper, "sum", "==", GetTarget(v, target), TRUE,
                 FALSE, FALSE, nThreads, pkgEnv$nThreads, tolerance))
}

partitionsCount <- function(v, m = NULL, repetition = FALSE,
                            freqs = NULL, target = NULL) {
    return(.Call(Algos_PartitionsCount, GetTarget(v, target), v, m,
                 repetition, freqs, "==", NULL, NULL, FALSE, FALSE))
}

partitionsDesign <- function(v, m = NULL, repetition = FALSE,
                             freqs = NULL, target = NULL,
                             showDesign = FALSE) {

    return(.Call(Algos_PartitionsCount, GetTarget(v, target), v, m,
                 repetition, freqs, "==", NULL, NULL, TRUE, showDesign))
}

partitionsSample <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                             target = NULL, n = NULL, sampleVec = NULL,
                             seed = NULL, Parallel = FALSE,
                             nThreads = NULL, namedSample = FALSE,
                             tolerance = NULL) {

    if (!is.null(seed)) {
        set.seed(seed)
    }

    return(.Call(Algos_SamplePartitions, v, m, repetition, freqs,
                 sampleVec, seed, n, sample, Parallel, nThreads,
                 pkgEnv$nThreads, namedSample, "==", GetTarget(v, target),
                 tolerance, new.env()))
}

partitionsIter <- function(v, m = NULL, repetition = FALSE,
                           freqs = NULL, target = NULL, lower = NULL,
                           upper = NULL, nThreads = NULL,
                           tolerance = NULL) {

    InitVals <- .Call(Algos_GetClassVals, v, m, repetition, freqs,
                      TRUE, NULL, nThreads, pkgEnv$nThreads, TRUE)

    new("Partitions", InitVals, FALSE, "sum", "==",
        GetTarget(v, target), FALSE, tolerance, is.null(m))
}

partLen <- function(tar, m, v, rep = FALSE, fr = NULL, comb = TRUE) {
    if (comb) {
        RcppAlgos243::comboGeneral(v, m, rep, fr, constraintFun = "sum",
                                   comparisonFun = "==", limitConstraints = tar)
    } else {
        RcppAlgos243::permuteGeneral(v, m, rep, fr, constraintFun = "sum",
                                     comparisonFun = "==", limitConstraints = tar)
    }
}
