comboGeneral <- function(v,m,repetition=FALSE,constraintFun=NULL,
                         comparisonFun=NULL,limitConstraints=NULL,rowCap=NULL) {
    CombinatoricsRcpp(v, m, repetition, constraintFun, 
                      comparisonFun, limitConstraints, rowCap, TRUE)
}

permuteGeneral <- function(v,m,repetition=FALSE,constraintFun=NULL,
                           comparisonFun=NULL,limitConstraints=NULL,rowCap=NULL) {
    CombinatoricsRcpp(v, m, repetition, constraintFun, 
                      comparisonFun, limitConstraints, rowCap, FALSE)
}