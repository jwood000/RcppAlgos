GetTarget <- function(v, target) {
    if (is.null(target)) {target = max(v)}
    return(target)
}

partitionGeneral <- function(v, m = NULL, repetition = FALSE,
                             freqs = NULL, target = NULL, lower = NULL,
                             upper = NULL, nThreads = NULL,
                             tolerance = NULL) {
    
    return(.Call(CombinatoricsCnstrt, v, m, repetition,
                 freqs, lower, upper, "sum", "==",
                 GetTarget(v, target), TRUE, FALSE, FALSE,
                 nThreads, pkgEnv$nThreads, tolerance,
                 PACKAGE = "RcppAlgos"))
}

partitionCount <- function(v, m = NULL, repetition = FALSE,
                           freqs = NULL, target = NULL) {
    .Call(PartitionsCount, GetTarget(v, target), v, m, repetition,
          freqs, "==", NULL, NULL, FALSE, FALSE, PACKAGE = "RcppAlgos");
}

partitionDesign <- function(v, m = NULL, repetition = FALSE,
                            freqs = NULL, target = NULL,
                            showDesign = FALSE) {
    
    .Call(PartitionsCount, GetTarget(v, target),
          v, m, repetition, freqs, "==", NULL, NULL,
          TRUE, showDesign, PACKAGE = "RcppAlgos");
}

partLen <- function(tar, m, v, rep = FALSE, fr = NULL) {
    RcppAlgos241::comboGeneral(v, m, rep, fr, constraintFun = "sum", comparisonFun = "==", limitConstraints = tar)
}
