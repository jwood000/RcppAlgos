comboGeneral <- function(n,r,v=1:n,repetition=FALSE,constraintFun=NULL,comparisonFun=NULL,limitConstraints=NULL,rowCap=NULL) {
    ComboRcpp(n, r, v, repetition, constraintFun, comparisonFun, limitConstraints, rowCap)
}
