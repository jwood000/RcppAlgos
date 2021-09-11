setClass(
    "Combo",
    slots = c(
        ptr           = "externalptr",
        startOver     = "function",
        nextIter      = "function",
        nextNIter     = "function",
        nextRemaining = "function",
        prevIter      = "function",
        prevNIter     = "function",
        prevRemaining = "function",
        currIter      = "function",
        randomAccess  = "function",
        sourceVector  = "function",
        front         = "function",
        back          = "function",
        summary       = "function"
    )
)

setMethod(
    "initialize",
    "Combo",
    function(.Object, init) {
        .Object@ptr <- .Call(ComboNew, init$RVals, init$bVec,
                             init$FreqsInfo, PACKAGE = "RcppAlgos")
        .Object@startOver <- function() {
            .Call(StartOverGlue, .Object@ptr, PACKAGE = "RcppAlgos")
            invisible(NULL)
        }
        .Object@nextIter <- function() {
            .Call(NextCombGlue, .Object@ptr, PACKAGE = "RcppAlgos")
        }
        .Object@nextNIter <- function(n = 1) {
            .Call(NextNumCombGlue, .Object@ptr, n, PACKAGE = "RcppAlgos")
        }
        .Object@nextRemaining <- function() {
            .Call(NextGatherGlue, .Object@ptr, PACKAGE = "RcppAlgos")
        }
        .Object@prevIter <- function() {
            .Call(PrevCombGlue, .Object@ptr, PACKAGE = "RcppAlgos")
        }
        .Object@prevNIter <- function(n = 1) {
            .Call(PrevNumCombGlue, .Object@ptr, n, PACKAGE = "RcppAlgos")
        }
        .Object@prevRemaining <- function() {
            .Call(PrevGatherGlue, .Object@ptr, PACKAGE = "RcppAlgos")
        }
        .Object@currIter <- function() {
            .Call(CurrCombGlue, .Object@ptr, PACKAGE = "RcppAlgos")
        }
        .Object@front <- function() {
            .Call(FrontGlue, .Object@ptr, PACKAGE = "RcppAlgos")
        }
        .Object@back <- function() {
            .Call(BackGlue, .Object@ptr, PACKAGE = "RcppAlgos")
        }
        .Object@randomAccess <- function(samp) {
            .Call(RandomAccessGlue, .Object@ptr, samp, PACKAGE = "RcppAlgos")
        }
        .Object@sourceVector <- function() {
            .Call(SourceVectorGlue, .Object@ptr, PACKAGE = "RcppAlgos")
        }
        .Object@summary <- function() {
            .Call(SummaryGlue, .Object@ptr, PACKAGE = "RcppAlgos")
        }
        .Object
    }
)

setClass(
    "ComboApply",
    contains = "Combo"
)

setMethod(
    "initialize",
    "ComboApply",
    function(.Object, init, stdFun, rho, RFunVal) {
        .Object@ptr <- .Call(ComboApplyNew, init$RVals, init$bVec,
                             init$FreqsInfo, stdFun, rho, RFunVal,
                             PACKAGE = "RcppAlgos")
        .Object@startOver <- function() {
            .Call(StartOverGlue, .Object@ptr, PACKAGE = "RcppAlgos")
            invisible(NULL)
        }
        .Object@nextIter <- function() {
            .Call(NextCombGlue, .Object@ptr, PACKAGE = "RcppAlgos")
        }
        .Object@nextNIter <- function(n = 1) {
            .Call(NextNumCombGlue, .Object@ptr, n, PACKAGE = "RcppAlgos")
        }
        .Object@nextRemaining <- function() {
            .Call(NextGatherGlue, .Object@ptr, PACKAGE = "RcppAlgos")
        }
        .Object@prevIter <- function() {
            .Call(PrevCombGlue, .Object@ptr, PACKAGE = "RcppAlgos")
        }
        .Object@prevNIter <- function(n = 1) {
            .Call(PrevNumCombGlue, .Object@ptr, n, PACKAGE = "RcppAlgos")
        }
        .Object@prevRemaining <- function() {
            .Call(PrevGatherGlue, .Object@ptr, PACKAGE = "RcppAlgos")
        }
        .Object@currIter <- function() {
            .Call(CurrCombGlue, .Object@ptr, PACKAGE = "RcppAlgos")
        }
        .Object@front <- function() {
            .Call(FrontGlue, .Object@ptr, PACKAGE = "RcppAlgos")
        }
        .Object@back <- function() {
            .Call(BackGlue, .Object@ptr, PACKAGE = "RcppAlgos")
        }
        .Object@randomAccess <- function(samp) {
            .Call(RandomAccessGlue, .Object@ptr, samp, PACKAGE = "RcppAlgos")
        }
        .Object@sourceVector <- function() {
            .Call(SourceVectorGlue, .Object@ptr, PACKAGE = "RcppAlgos")
        }
        .Object@summary <- function() {
            .Call(SummaryGlue, .Object@ptr, PACKAGE = "RcppAlgos")
        }
        .Object
    }
)

"[.Combo" <- function(x, ...) {
    x@randomAccess(...)
}

"[.ComboApply" <- function(x, ...) {
    x@randomAccess(...)
}

comboIter <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                      constraintFun = NULL, comparisonFun = NULL,
                      limitConstraints = NULL, FUN = NULL,
                      tolerance = NULL, nThreads = NULL,
                      FUN.VALUE = NULL) {

    RetValue <- .Call(CheckReturn, v, constraintFun,
                      comparisonFun, limitConstraints,
                      FALSE, FUN, PACKAGE = "RcppAlgos")
    
    InitVals <- .Call(GetClassVals, v, m, repetition, freqs,
                      TRUE, FUN, nThreads, pkgEnv$nThreads,
                      PACKAGE = "RcppAlgos")

    if (RetValue) {
        if (InitVals$applyFun) {
            new("ComboApply", InitVals, FUN, new.env(), FUN.VALUE)
        } else {
            new("Combo", InitVals)
        }
    } else {
        stop("This feature will be available in future releases")
    }
}

permuteIter <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                        constraintFun = NULL, comparisonFun = NULL,
                        limitConstraints = NULL, FUN = NULL,
                        tolerance = NULL, nThreads = NULL,
                        FUN.VALUE = NULL) {
    
    RetValue <- .Call(CheckReturn, v, constraintFun,
                      comparisonFun, limitConstraints,
                      FALSE, FUN, PACKAGE = "RcppAlgos")
    
    InitVals <- .Call(GetClassVals, v, m, repetition, freqs,
                      FALSE, FUN, nThreads, pkgEnv$nThreads,
                      PACKAGE = "RcppAlgos")
    
    if (RetValue) {
        if (InitVals$applyFun) {
            new("ComboApply", InitVals, FUN, new.env(), FUN.VALUE)
        } else {
            new("Combo", InitVals)
        }
    } else {
        stop("This feature will be available in future releases")
    }
}
