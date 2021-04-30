partitionGeneral <- function(v, m = NULL, repetition = FALSE,
                             freqs = NULL, target, lower = NULL,
                             upper = NULL, nThreads = NULL,
                             tolerance = NULL) {
    
    return(.Call(CombinatoricsCnstrt, v, m, repetition,
                 freqs, lower, upper, "sum", "==",
                 target, TRUE, FALSE, FALSE, nThreads,
                 pkgEnv$nThreads, tolerance,
                 PACKAGE = "RcppAlgos"))
}

partitionCount <-  function(v, m = NULL, repetition = FALSE,
                            freqs = NULL, target) {
    .Call(PartitionsCount, target, v, m, repetition, freqs,
          "==", NULL, NULL, FALSE, PACKAGE = "RcppAlgos");
}

partitionDesign <-  function(v, m = NULL, repetition = FALSE,
                             freqs = NULL, target) {
    .Call(PartitionsCount, target, v, m, repetition, freqs,
          "==", NULL, NULL, TRUE, PACKAGE = "RcppAlgos");
}

partLen <- function(tar, m, v, rep = FALSE, fr = NULL) {
    RcppAlgos241::comboGeneral(v, m, rep, fr, constraintFun = "sum", comparisonFun = "==", limitConstraints = tar)
}
