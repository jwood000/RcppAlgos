comboGeneral <- function(v, m, repetition=FALSE, constraintFun=NULL,
                         comparisonFun=NULL, limitConstraints=NULL,
                         rowCap=NULL, keepResults=FALSE) {
    isFactor <- is.factor(v)
    CombinatoricsRcpp(v, m, repetition, constraintFun, 
                      comparisonFun, limitConstraints,
                      rowCap, TRUE, isFactor, keepResults, NULL)
}

permuteGeneral <- function(v, m=NULL, repetition=FALSE, constraintFun=NULL,
                           comparisonFun=NULL, limitConstraints=NULL,
                           rowCap=NULL, keepResults=FALSE, freqs = NULL) {
    isFactor <- is.factor(v)
    CombinatoricsRcpp(v, m, repetition, constraintFun, 
                      comparisonFun, limitConstraints,
                      rowCap, FALSE, isFactor, keepResults, freqs)
}