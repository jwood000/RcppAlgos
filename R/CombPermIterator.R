comboIter <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                      constraintFun = NULL, comparisonFun = NULL,
                      limitConstraints = NULL, FUN = NULL, tolerance = NULL) {
    
    IsFactor <- is.factor(v)
    IsStdRet <- (CheckReturn(v, constraintFun, comparisonFun,
                             limitConstraints, IsFactor, FALSE, FUN) > 0)
    InitVals <- GetClassVals(IsStdRet, v, m, repetition, freqs, TRUE, IsFactor, FUN);
    
    if (IsStdRet) {
        if (InitVals$applyFun) {
            new(ComboFUN, InitVals$Rv, InitVals$Rm, InitVals$nRows,
                InitVals$bVec, InitVals$freqInfo, list(FUN, new.env()))
        } else {
            new(Combo, InitVals$Rv, InitVals$Rm, InitVals$nRows, InitVals$bVec, InitVals$freqInfo)
        }
    } else {
        stop("This feature will be available in future releases")
    }
}

permuteIter <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                        constraintFun = NULL, comparisonFun = NULL,
                        limitConstraints = NULL, FUN = NULL, tolerance = NULL) {
    
    IsFactor <- is.factor(v)
    IsStdRet <- (CheckReturn(v, constraintFun, comparisonFun,
                             limitConstraints, IsFactor, FALSE, FUN) > 0)
    InitVals <- GetClassVals(IsStdRet, v, m, repetition, freqs, FALSE, IsFactor, FUN);
    
    if (IsStdRet) {
        if (InitVals$applyFun) {
            new(ComboFUN, InitVals$Rv, InitVals$Rm, InitVals$nRows,
                InitVals$bVec, InitVals$freqInfo, list(FUN, new.env()))
        } else {
            new(Combo, InitVals$Rv, InitVals$Rm, InitVals$nRows, InitVals$bVec, InitVals$freqInfo)
        }
    } else {
        stop("This feature will be available in future releases")
    }
}

loadModule("Combo", TRUE)
loadModule("ComboFUN", TRUE)
