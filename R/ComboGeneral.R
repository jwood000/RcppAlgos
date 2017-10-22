comboGeneral <- function(v,m,repetition=FALSE,constraintFun=NULL,comparisonFun=NULL,limitConstraints=NULL,rowCap=NULL) {
    ComboRcpp(v, m, repetition, constraintFun, comparisonFun, limitConstraints, rowCap)
}
