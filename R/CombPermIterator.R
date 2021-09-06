setClass(
    "Combo",
    slots = c(
        ptr = "externalptr",
        startOver = "function",
        nextIter = "function",
        nextNIter = "function",
        nextRemaining = "function",
        prevIter = "function",
        prevNIter = "function",
        prevRemaining = "function",
        currIter = "function",
        sourceVector = "function",
        summary = "function",
        sample = "function"
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
        .Object@sample <- function(samp) {
            .Call(RandomAccessGlue, .Object@ptr, samp, PACKAGE = "RcppAlgos")
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
            .Call(PrevNumCombGlue, .Object@ptr, PACKAGE = "RcppAlgos")
        }
        .Object@prevNIter <- function(n = 1) {
            .Call(PrevCombGlue, .Object@ptr, n, PACKAGE = "RcppAlgos")
        }
        .Object@prevRemaining <- function() {
            .Call(PrevGatherGlue, .Object@ptr, PACKAGE = "RcppAlgos")
        }
        .Object@currIter <- function() {
            .Call(CurrCombGlue, .Object@ptr, PACKAGE = "RcppAlgos")
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

comboIter <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                      constraintFun = NULL, comparisonFun = NULL,
                      limitConstraints = NULL, FUN = NULL,
                      tolerance = NULL, nThreads = NULL) {

    RetValue <- .Call(CheckReturn, v, constraintFun,
                      comparisonFun, limitConstraints,
                      FALSE, FUN, PACKAGE = "RcppAlgos")
    
    InitVals <- .Call(GetClassVals, v, m, repetition, freqs,
                      TRUE, FUN, nThreads, pkgEnv$nThreads,
                      PACKAGE = "RcppAlgos")

    if (RetValue) {
        # if (InitVals$applyFun) {
        #     # new(ComboFUN, InitVals$Rv, InitVals$Rm, InitVals$nRows,
        #     #     InitVals$bVec, InitVals$freqInfo, list(FUN, new.env()))
        # } else {
            # new(Combo, InitVals$Rv, InitVals$Rm, InitVals$nRows, InitVals$bVec, InitVals$freqInfo)
        new("Combo", InitVals)
        # }
    } else {
        stop("This feature will be available in future releases")
    }
}
# 
# permuteIter <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
#                         constraintFun = NULL, comparisonFun = NULL,
#                         limitConstraints = NULL, FUN = NULL, tolerance = NULL) {
#     
#     IsFactor <- is.factor(v)
#     IsStdRet <- (CheckReturn(v, constraintFun, comparisonFun,
#                              limitConstraints, IsFactor, FALSE, FUN) > 0)
#     InitVals <- GetClassVals(IsStdRet, v, m, repetition, freqs, FALSE, IsFactor, FUN);
#     
#     if (IsStdRet) {
#         if (InitVals$applyFun) {
#             new(ComboFUN, InitVals$Rv, InitVals$Rm, InitVals$nRows,
#                 InitVals$bVec, InitVals$freqInfo, list(FUN, new.env()))
#         } else {
#             new(Combo, InitVals$Rv, InitVals$Rm, InitVals$nRows, InitVals$bVec, InitVals$freqInfo)
#         }
#     } else {
#         stop("This feature will be available in future releases")
#     }
# }
