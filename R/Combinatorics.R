comboGeneral <- function(v, m = NULL, repetition = FALSE, freqs = NULL, 
                         lower = NULL, upper = NULL, constraintFun = NULL,
                         comparisonFun = NULL, limitConstraints = NULL, 
                         keepResults = NULL, FUN = NULL, Parallel = FALSE,
                         nThreads = NULL, tolerance = NULL) {
    
    isFactor <- is.factor(v)
    CombinatoricsRcpp(v, m, repetition, freqs, lower, upper, constraintFun, 
                      comparisonFun, limitConstraints, TRUE, keepResults, isFactor,
                      FALSE, FUN, new.env(), Parallel, nThreads, pkgEnv$nThreads, tolerance)
}

permuteGeneral <- function(v, m = NULL, repetition = FALSE, freqs = NULL, 
                           lower = NULL, upper = NULL, constraintFun = NULL,
                           comparisonFun = NULL, limitConstraints = NULL,
                           keepResults = NULL, FUN = NULL, Parallel = FALSE,
                           nThreads = NULL, tolerance = NULL) {
    
    isFactor <- is.factor(v)
    CombinatoricsRcpp(v, m, repetition, freqs, lower, upper, constraintFun, 
                      comparisonFun, limitConstraints, FALSE, keepResults, isFactor,
                      FALSE, FUN, new.env(), Parallel, nThreads, pkgEnv$nThreads, tolerance)
}

comboCount <-  function(v, m = NULL, repetition = FALSE, freqs = NULL) {
    
    CombinatoricsRcpp(v, m, repetition, freqs, NULL, NULL,
                      NULL, NULL, NULL, TRUE, FALSE, FALSE,
                      TRUE, NULL, NULL, FALSE, NULL, 0, NULL)
}

permuteCount <- function(v, m = NULL, repetition = FALSE, freqs = NULL) {
    
    CombinatoricsRcpp(v, m, repetition, freqs, NULL, NULL,
                      NULL, NULL, NULL, FALSE, FALSE, FALSE,
                      TRUE, NULL, NULL, FALSE, NULL, 0, NULL)
}

comboSample <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                        n = NULL, sampleVec = NULL, seed = NULL, FUN = NULL,
                        Parallel = FALSE, nThreads = NULL, namedSample = FALSE) {
    
    isFactor <- is.factor(v)
    if (!is.null(seed)) {set.seed(seed)}
    SampleRcpp(v, m, repetition, freqs, sampleVec, TRUE, isFactor,
               seed, n, sample, FUN, new.env(), Parallel, nThreads,
               pkgEnv$nThreads, namedSample)
}

permuteSample <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                          n = NULL, sampleVec = NULL, seed = NULL, FUN = NULL,
                          Parallel = FALSE, nThreads = NULL, namedSample = FALSE) {
    
    isFactor <- is.factor(v)
    if (!is.null(seed)) {set.seed(seed)}
    SampleRcpp(v, m, repetition, freqs, sampleVec, FALSE, isFactor,
               seed, n, sample, FUN, new.env(), Parallel, nThreads,
               pkgEnv$nThreads, namedSample)
}

comboGroups <- function(v, numGroups, retType = "matrix", lower = NULL,
                       upper = NULL, Parallel = FALSE, nThreads = NULL) {
    isFactor <- is.factor(v)
    ComboGroupsRcpp(v, numGroups, retType, lower, upper, 
                    isFactor, FALSE, Parallel, nThreads, pkgEnv$nThreads,
                    FALSE, NULL, NULL, NULL, sample, FALSE)
}

comboGroupsSample <- function(v, numGroups, retType = "matrix", n = NULL,
                              sampleVec = NULL, seed = NULL, Parallel = FALSE, 
                              nThreads = NULL, namedSample = FALSE) {
    isFactor <- is.factor(v)
    if (!is.null(seed)) {set.seed(seed)}
    ComboGroupsRcpp(v, numGroups, retType, NULL, NULL, 
                    isFactor, FALSE, Parallel, nThreads, pkgEnv$nThreads,
                    TRUE, sampleVec, seed, n, sample, namedSample)
}

comboGroupsCount <- function(v, numGroups) {
    ComboGroupsRcpp(v, numGroups, NULL, NULL, NULL, FALSE, TRUE, 
                    FALSE, FALSE, 0, FALSE, NULL, NULL, NULL, sample, FALSE)
}

