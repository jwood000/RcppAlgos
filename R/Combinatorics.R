comboGeneral <- function(v,m,repetition=FALSE,constraintFun=NULL,
                         comparisonFun=NULL,limitConstraints=NULL,rowCap=NULL) {
    isFactor <- is.factor(v)
    CombinatoricsRcpp(v, m, repetition, constraintFun, 
                      comparisonFun, limitConstraints, rowCap, TRUE, isFactor)
}

permuteGeneral <- function(v,m,repetition=FALSE,constraintFun=NULL,
                           comparisonFun=NULL,limitConstraints=NULL,rowCap=NULL) {
    isFactor <- is.factor(v)
    CombinatoricsRcpp(v, m, repetition, constraintFun, 
                      comparisonFun, limitConstraints, rowCap, FALSE, isFactor)
}