comboGeneral <- function(v, m = NULL, repetition = FALSE, freqs = NULL, lower = NULL,
                         upper = NULL, constraintFun = NULL, comparisonFun = NULL,
                         limitConstraints = NULL, keepResults = NULL, FUN = NULL,
                         Parallel = FALSE, nThreads = NULL, tolerance = NULL,
                         FUN.VALUE = NULL, simplify = FALSE) {
    
    RetValue <- .Call(CheckReturn, v, constraintFun,
                      comparisonFun, limitConstraints,
                      keepResults, FUN, PACKAGE = "RcppAglos")

    if (RetValue == 1) {
        return(.Call(CombinatoricsStndrd, v, m, repetition,
                     freqs, lower, upper, Parallel, nThreads,
                     pkgEnv$nThreads, TRUE, PACKAGE = "RcppAlgos"))
    } else if (RetValue == 2) {
        if (simplify && is.null(FUN.VALUE)) {
            res <- .Call(CombinatoricsApply, v, m,
                         repetition, freqs, lower, upper,
                         FUN, new.env(), FUN.VALUE, TRUE,
                         FALSE, PACKAGE = "RcppAlgos")
            lens <- lengths(res)

            if (max(lens) == min(lens)) {
                mat <- matrix(unlist(res, use.names = FALSE),
                              nrow = lens[1L], ncol = length(res))
                dim(mat) <- c(dim(res[[1]]), length(res))
                return(mat)
            } else {
                return(res)
            }
        } else if (!is.null(FUN.VALUE)) {
            return(.Call(CombinatoricsApply, v, m,
                         repetition, freqs, lower, upper,
                         FUN, new.env(), FUN.VALUE, TRUE,
                         TRUE, PACKAGE = "RcppAlgos"))
        } else {
            return(.Call(CombinatoricsApply, v, m,
                         repetition, freqs, lower, upper,
                         FUN, new.env(), FUN.VALUE, TRUE,
                         FALSE, PACKAGE = "RcppAlgos"))
        }
    } else {
        return(.Call(CombinatoricsCnstrt, v, m, repetition,
                     freqs, lower, upper, constraintFun, comparisonFun,
                     limitConstraints, TRUE, keepResults, Parallel,
                     nThreads, pkgEnv$nThreads, tolerance,
                     PACKAGE = "RcppAlgos"))
    }
}

# comboSample <- function(v, m = NULL, repetition = FALSE, freqs = NULL, n = NULL, 
#                         sampleVec = NULL, seed = NULL, FUN = NULL, Parallel = FALSE,
#                         nThreads = NULL, namedSample = FALSE) {
#     
#     if (!is.null(seed)) {set.seed(seed)}
#     SampleRcpp(v, m, repetition, freqs, sampleVec, TRUE, is.factor(v), seed, n,
#                sample, FUN, new.env(), Parallel, nThreads, pkgEnv$nThreads, namedSample)
# }
# 
# comboCount <-  function(v, m = NULL, repetition = FALSE, freqs = NULL) {
#     CombinatoricsCount(v, m, repetition, freqs, TRUE);
# }

stdThreadMax <- function() {
    nThreads <- .Call(cpp11GetNumThreads)
    if (nThreads < 1L || is.na(nThreads) || is.null(nThreads)) nThreads = 1L
    nThreads
}
