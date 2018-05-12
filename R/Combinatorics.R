comboGeneral <- function(v, m=NULL, repetition=FALSE, constraintFun=NULL,
                         comparisonFun=NULL, limitConstraints=NULL,
                         rowCap=NULL, keepResults=FALSE, freqs=NULL) {
    isFactor <- is.factor(v)
    CombinatoricsRcpp(v, m, repetition, constraintFun, 
                      comparisonFun, limitConstraints,
                      rowCap, TRUE, isFactor, keepResults, freqs)
}

permuteGeneral <- function(v, m=NULL, repetition=FALSE, constraintFun=NULL,
                           comparisonFun=NULL, limitConstraints=NULL,
                           rowCap=NULL, keepResults=FALSE, freqs=NULL) {
    isFactor <- is.factor(v)
    CombinatoricsRcpp(v, m, repetition, constraintFun, 
                      comparisonFun, limitConstraints,
                      rowCap, FALSE, isFactor, keepResults, freqs)
}

nthCombination <- function(v, m=NULL, index, repetition=FALSE, freqs=NULL) {
    isFactor <- is.factor(v)
    NthResultRcpp(v, m, index, repetition, TRUE, isFactor, freqs)
}

nthPermutation <- function(v, m=NULL, index, repetition=FALSE, freqs=NULL) {
    isFactor <- is.factor(v)
    NthResultRcpp(v, m, index, repetition, FALSE, isFactor, freqs)
}