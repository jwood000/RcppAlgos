comboGeneral <- function(v, m = NULL, repetition = FALSE, freqs = NULL, lower = NULL,
                         upper = NULL, constraintFun = NULL, comparisonFun = NULL,
                         limitConstraints = NULL, keepResults = NULL, FUN = NULL,
                         Parallel = FALSE, nThreads = NULL, tolerance = NULL) {
    
    CombinatoricsRcpp(v, m, repetition, freqs, lower, upper, constraintFun, 
                      comparisonFun, limitConstraints, TRUE, keepResults, is.factor(v), 
                      FALSE, FUN, new.env(), Parallel, nThreads, pkgEnv$nThreads, tolerance)
}

permuteGeneral <- function(v, m = NULL, repetition = FALSE, freqs = NULL, lower = NULL,
                           upper = NULL, constraintFun = NULL, comparisonFun = NULL, 
                           limitConstraints = NULL, keepResults = NULL, FUN = NULL,
                           Parallel = FALSE, nThreads = NULL, tolerance = NULL) {

    CombinatoricsRcpp(v, m, repetition, freqs, lower, upper, constraintFun, 
                      comparisonFun, limitConstraints, FALSE, keepResults, is.factor(v), 
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

comboSample <- function(v, m = NULL, repetition = FALSE, freqs = NULL, n = NULL, 
                        sampleVec = NULL, seed = NULL, FUN = NULL, Parallel = FALSE,
                        nThreads = NULL, namedSample = FALSE) {
    
    if (!is.null(seed)) {set.seed(seed)}
    SampleRcpp(v, m, repetition, freqs, sampleVec, TRUE, is.factor(v), seed, n,
               sample, FUN, new.env(), Parallel, nThreads, pkgEnv$nThreads, namedSample)
}

permuteSample <- function(v, m = NULL, repetition = FALSE, freqs = NULL, n = NULL,
                          sampleVec = NULL, seed = NULL, FUN = NULL, Parallel = FALSE,
                          nThreads = NULL, namedSample = FALSE) {
    
    if (!is.null(seed)) {set.seed(seed)}
    SampleRcpp(v, m, repetition, freqs, sampleVec, FALSE, is.factor(v), seed, n,
               sample, FUN, new.env(), Parallel, nThreads, pkgEnv$nThreads, namedSample)
}

comboGroups <- function(v, numGroups, retType = "matrix", lower = NULL,
                       upper = NULL, Parallel = FALSE, nThreads = NULL) {
    
    ComboGroupsRcpp(v, numGroups, retType, lower, upper, is.factor(v), FALSE, 
                    Parallel, nThreads, pkgEnv$nThreads, FALSE, NULL, NULL,
                    NULL, sample, FALSE)
}

comboGroupsSample <- function(v, numGroups, retType = "matrix", n = NULL,
                              sampleVec = NULL, seed = NULL, Parallel = FALSE, 
                              nThreads = NULL, namedSample = FALSE) {
    
    if (!is.null(seed)) {set.seed(seed)}
    ComboGroupsRcpp(v, numGroups, retType, NULL, NULL, is.factor(v), FALSE,
                    Parallel, nThreads, pkgEnv$nThreads, TRUE, sampleVec, seed,
                    n, sample, namedSample)
}

comboGroupsCount <- function(v, numGroups) {
    ComboGroupsRcpp(v, numGroups, NULL, NULL, NULL, FALSE, TRUE, 
                    FALSE, FALSE, 0, FALSE, NULL, NULL, NULL, sample, FALSE)
}

