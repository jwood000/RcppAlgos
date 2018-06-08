comboGeneral <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                         lower = NULL, upper = NULL, constraintFun = NULL,
                         comparisonFun = NULL, limitConstraints = NULL, 
                         keepResults = FALSE) {
    
    isFactor <- is.factor(v)
    CombinatoricsRcpp(v, m, repetition, freqs, lower,
                      upper, constraintFun, 
                      comparisonFun, limitConstraints,
                      TRUE, keepResults, isFactor, FALSE)
}

permuteGeneral <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                           lower = NULL, upper = NULL, constraintFun = NULL,
                           comparisonFun = NULL, limitConstraints = NULL,
                           keepResults = FALSE) {
    
    isFactor <- is.factor(v)
    CombinatoricsRcpp(v, m, repetition, freqs, lower,
                      upper, constraintFun, 
                      comparisonFun, limitConstraints,
                      FALSE, keepResults, isFactor, FALSE)
}

comboCount <-  function(v, m = NULL, repetition = FALSE, freqs = NULL) {
    
    CombinatoricsRcpp(v, m, repetition, freqs, NULL,
                      NULL, NULL, NULL, NULL,
                      TRUE, FALSE, FALSE, TRUE)
}

permuteCount <- function(v, m = NULL, repetition = FALSE, freqs = NULL) {
    
    CombinatoricsRcpp(v, m, repetition, freqs, NULL,
                      NULL, NULL, NULL, NULL,
                      FALSE, FALSE, FALSE, TRUE)
}

comboSample <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                         n = NULL, sampleVec = NULL) {
    
    myCount <- comboCount(v, m, repetition, freqs)
    
    if (is.null(sampleVec)) {
        if (is.null(n)){
            stop("n and sampleVec cannot both be NULL")
        } else {
            if (!is.numeric(n))
                stop("n must be a number")
            else if (length(n) > 1)
                stop("length of n must be 1. For specific combinations, use sampleVec.")
            else if (n > myCount)
                stop("n exceeds the maximum number of possible results")
        } 
        sampleVec = sample(myCount, n)
    }
    
    isFactor <- is.factor(v)
    n <- length(sampleVec)
    
    SampleRcpp(v, m, repetition, freqs, sampleVec, TRUE, isFactor, myCount)
}

permuteSample <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                          n = NULL, sampleVec = NULL) {
    
    myCount <- permuteCount(v, m, repetition, freqs)
    
    if (is.null(sampleVec)) {
        if (is.null(n)){
            stop("n and sampleVec cannot both be NULL")
        } else {
            if (!is.numeric(n))
                stop("n must be a number")
            else if (length(n) > 1)
                stop("length of n must be 1. For specific permutation, use sampleVec.")
            else if (n > myCount)
                stop("n exceeds the maximum number of possible results")
        } 
        sampleVec = sample(myCount, n)
    }
    
    isFactor <- is.factor(v)
    n <- length(sampleVec)
    
    SampleRcpp(v, m, repetition, freqs, sampleVec, FALSE, isFactor, myCount)
}

