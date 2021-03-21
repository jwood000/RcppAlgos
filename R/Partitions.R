partitionGeneral <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                             lower = NULL, upper = NULL, target = NULL, 
                             nThreads = NULL, tolerance = NULL) {
    
    return(.Call(CombinatoricsCnstrt, v, m, repetition,
                 freqs, lower, upper, "sum", "==",
                 target, TRUE, FALSE, FALSE,
                 nThreads, pkgEnv$nThreads, tolerance,
                 PACKAGE = "RcppAlgos"))
}

# comboSample <- function(v, m = NULL, repetition = FALSE, freqs = NULL, n = NULL, 
#                         sampleVec = NULL, seed = NULL, FUN = NULL, Parallel = FALSE,
#                         nThreads = NULL, namedSample = FALSE) {
#     
#     if (!is.null(seed)) {set.seed(seed)}
#     SampleRcpp(v, m, repetition, freqs, sampleVec, TRUE, is.factor(v), seed, n,
#                sample, FUN, new.env(), Parallel, nThreads, pkgEnv$nThreads, namedSample)
# }

partitionCount <-  function(v, m = NULL, repetition = FALSE, freqs = NULL) {
    .Call(CombinatoricsCount, v, m,
          repetition, freqs, TRUE,
          PACKAGE = "RcppAlgos");
}
