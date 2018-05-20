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