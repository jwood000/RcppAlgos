permuteGeneral <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                           lower = NULL, upper = NULL, constraintFun = NULL,
                           comparisonFun = NULL, limitConstraints = NULL,
                           keepResults = NULL, FUN = NULL, Parallel = FALSE,
                           nThreads = NULL, tolerance = NULL, FUN.VALUE = NULL) {

    RetValue <- .Call(CheckReturn, v, constraintFun,
                      comparisonFun, limitConstraints,
                      keepResults, FUN, PACKAGE = "RcppAglos")
    
    if (RetValue == 1) {
        return(.Call(CombinatoricsStndrd, v, m, repetition,
                     freqs, lower, upper, Parallel, nThreads,
                     pkgEnv$nThreads, FALSE, PACKAGE = "RcppAlgos"))
    } else if (RetValue == 2) {
        return(.Call(CombinatoricsApply, v, m,
                     repetition, freqs, lower, upper,
                     FUN, new.env(), FUN.VALUE, FALSE,
                     PACKAGE = "RcppAlgos"))
    } else {
        return(.Call(CombinatoricsCnstrt, v, m, repetition,
                     freqs, lower, upper, constraintFun, comparisonFun,
                     limitConstraints, FALSE, keepResults, Parallel,
                     nThreads, pkgEnv$nThreads, tolerance,
                     PACKAGE = "RcppAlgos"))
    }
}

# 
# permuteSample <- function(v, m = NULL, repetition = FALSE, freqs = NULL, n = NULL,
#                           sampleVec = NULL, seed = NULL, FUN = NULL, Parallel = FALSE,
#                           nThreads = NULL, namedSample = FALSE) {
#     
#     if (!is.null(seed)) {set.seed(seed)}
#     SampleRcpp(v, m, repetition, freqs, sampleVec, FALSE, is.factor(v), seed, n,
#                sample, FUN, new.env(), Parallel, nThreads, pkgEnv$nThreads, namedSample)
# }

permuteCount <- function(v, m = NULL, repetition = FALSE, freqs = NULL) {
    .Call(CombinatoricsCount, v, m,
          repetition, freqs, FALSE, PACKAGE = "RcppAlgos");
}
