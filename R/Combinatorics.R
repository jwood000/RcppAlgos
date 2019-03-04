comboGeneral <- function(v, m = NULL, repetition = FALSE, freqs = NULL, 
                         lower = NULL, upper = NULL, constraintFun = NULL,
                         comparisonFun = NULL, limitConstraints = NULL, 
                         keepResults = NULL, FUN = NULL, 
                         Parallel = FALSE, nThreads = NULL) {
    
    isFactor <- is.factor(v)
    CombinatoricsRcpp(v, m, repetition, freqs, lower, upper, constraintFun, 
                      comparisonFun, limitConstraints, TRUE, keepResults, isFactor,
                      FALSE, FUN, new.env(), Parallel, nThreads, pkgEnv$nThreads)
}

permuteGeneral <- function(v, m = NULL, repetition = FALSE, freqs = NULL, 
                           lower = NULL, upper = NULL, constraintFun = NULL,
                           comparisonFun = NULL, limitConstraints = NULL,
                           keepResults = NULL, FUN = NULL, 
                           Parallel = FALSE, nThreads = NULL) {
    
    isFactor <- is.factor(v)
    CombinatoricsRcpp(v, m, repetition, freqs, lower, upper, constraintFun, 
                      comparisonFun, limitConstraints, FALSE, keepResults, isFactor,
                      FALSE, FUN, new.env(), Parallel, nThreads, pkgEnv$nThreads)
}

comboCount <-  function(v, m = NULL, repetition = FALSE, freqs = NULL) {
    
    CombinatoricsRcpp(v, m, repetition, freqs, NULL, NULL,
                      NULL, NULL, NULL, TRUE, FALSE, FALSE,
                      TRUE, NULL, NULL, FALSE, NULL, 0)
}

permuteCount <- function(v, m = NULL, repetition = FALSE, freqs = NULL) {
    
    CombinatoricsRcpp(v, m, repetition, freqs, NULL, NULL,
                      NULL, NULL, NULL, FALSE, FALSE, FALSE,
                      TRUE, NULL, NULL, FALSE, NULL, 0)
}

comboSample <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                        n = NULL, sampleVec = NULL, seed = NULL, FUN = NULL,
                        Parallel = FALSE, nThreads = NULL) {
    
    isFactor <- is.factor(v)
    if (!is.null(seed)) {set.seed(seed)}
    SampleRcpp(v, m, repetition, freqs, sampleVec, TRUE, isFactor,
               seed, n, sample, FUN, new.env(), Parallel, nThreads, pkgEnv$nThreads)
}

permuteSample <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                          n = NULL, sampleVec = NULL, seed = NULL, FUN = NULL,
                          Parallel = FALSE, nThreads = NULL) {
    
    isFactor <- is.factor(v)
    if (!is.null(seed)) {set.seed(seed)}
    SampleRcpp(v, m, repetition, freqs, sampleVec, FALSE, isFactor,
               seed, n, sample, FUN, new.env(), Parallel, nThreads, pkgEnv$nThreads)
}

