comboGeneral <- function(v,m,repetition=FALSE,constraintFun=NULL,
                         comparisonFun=NULL,limitConstraints=NULL,
                         rowCap=NULL,keepResults=FALSE) {
    isFactor <- is.factor(v)
    CombinatoricsRcpp(v, m, repetition, constraintFun, 
                      comparisonFun, limitConstraints,
                      rowCap, TRUE, isFactor, keepResults)
}

permuteGeneral <- function(v,m,repetition=FALSE,constraintFun=NULL,
                           comparisonFun=NULL,limitConstraints=NULL,
                           rowCap=NULL,keepResults=FALSE) {
    isFactor <- is.factor(v)
    CombinatoricsRcpp(v, m, repetition, constraintFun, 
                      comparisonFun, limitConstraints,
                      rowCap, FALSE, isFactor, keepResults)
}

permuteRepLen <- function(v, repLengths) {
    isFactor <- is.factor(v)
    PermsRepLenRcpp(v, repLengths, isFactor)
}