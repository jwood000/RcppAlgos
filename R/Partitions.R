GetTarget <- function(v, target) {
    if (is.null(target)) {target = max(v)}
    return(target)
}

partitionsGeneral <- function(v, m = NULL, repetition = FALSE,
                              freqs = NULL, target = NULL, lower = NULL,
                              upper = NULL, nThreads = NULL,
                              tolerance = NULL) {
    
    return(.Call(CombinatoricsCnstrt, v, m, repetition,
                 freqs, lower, upper, "sum", "==",
                 GetTarget(v, target), TRUE, FALSE, FALSE,
                 nThreads, pkgEnv$nThreads, tolerance,
                 PACKAGE = "RcppAlgos"))
}

partitionsCount <- function(v, m = NULL, repetition = FALSE,
                            freqs = NULL, target = NULL) {
    .Call(PartitionsCount, GetTarget(v, target), v, m, repetition,
          freqs, "==", NULL, NULL, FALSE, FALSE, PACKAGE = "RcppAlgos")
}

partitionsDesign <- function(v, m = NULL, repetition = FALSE,
                             freqs = NULL, target = NULL,
                             showDesign = FALSE) {
    
    .Call(PartitionsCount, GetTarget(v, target),
          v, m, repetition, freqs, "==", NULL, NULL,
          TRUE, showDesign, PACKAGE = "RcppAlgos")
}

partitionsSample <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                             target = NULL, n = NULL, sampleVec = NULL,
                             seed = NULL, Parallel = FALSE,
                             nThreads = NULL, namedSample = FALSE,
                             tolerance = NULL) {

    if (!is.null(seed)) {
        set.seed(seed)
    }
    
    return(.Call(SamplePartitions, v, m, repetition, freqs, sampleVec, seed,
                 n, sample, Parallel, nThreads, pkgEnv$nThreads, namedSample,
                 "==", GetTarget(v, target), tolerance, new.env(),
                 PACKAGE = "RcppAlgos"))
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
