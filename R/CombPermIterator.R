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
    function(.Object, init, Parallel) {
        .Object@ptr <- .Call(Algos_CombClassNew, init$RVals, init$bVec,
                             init$FreqsInfo, Parallel, NULL, NULL, NULL,
                             NULL, NULL, NULL, NULL, NULL, NULL, 1)
        .Object@startOver <- function() {
            .Call(Algos_StartOverGlue, .Object@ptr)
            invisible(NULL)
        }
        .Object@nextIter <- function() {
            .Call(Algos_NextCombGlue, .Object@ptr)
        }
        .Object@nextNIter <- function(n = 1) {
            .Call(Algos_NextNumCombGlue, .Object@ptr, n)
        }
        .Object@nextRemaining <- function() {
            .Call(Algos_NextGatherGlue, .Object@ptr)
        }
        .Object@prevIter <- function() {
            .Call(Algos_PrevCombGlue, .Object@ptr)
        }
        .Object@prevNIter <- function(n = 1) {
            .Call(Algos_PrevNumCombGlue, .Object@ptr, n)
        }
        .Object@prevRemaining <- function() {
            .Call(Algos_PrevGatherGlue, .Object@ptr)
        }
        .Object@currIter <- function() {
            .Call(Algos_CurrCombGlue, .Object@ptr)
        }
        .Object@front <- function() {
            .Call(Algos_FrontGlue, .Object@ptr)
        }
        .Object@back <- function() {
            .Call(Algos_BackGlue, .Object@ptr)
        }
        .Object@randomAccess <- function(samp) {
            .Call(Algos_RandomAccessGlue, .Object@ptr, samp)
        }
        .Object@sourceVector <- function() {
            .Call(Algos_SourceVectorGlue, .Object@ptr)
        }
        .Object@summary <- function() {
            .Call(Algos_SummaryGlue, .Object@ptr)
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
        .Object@ptr <- .Call(Algos_CombClassNew, init$RVals, init$bVec,
                             init$FreqsInfo, FALSE, stdFun, rho, RFunVal,
                             NULL, NULL, NULL, NULL, NULL, NULL, 2)
        .Object@startOver <- function() {
            .Call(Algos_StartOverGlue, .Object@ptr)
            invisible(NULL)
        }
        .Object@nextIter <- function() {
            .Call(Algos_NextCombGlue, .Object@ptr)
        }
        .Object@nextNIter <- function(n = 1) {
            .Call(Algos_NextNumCombGlue, .Object@ptr, n)
        }
        .Object@nextRemaining <- function() {
            .Call(Algos_NextGatherGlue, .Object@ptr)
        }
        .Object@prevIter <- function() {
            .Call(Algos_PrevCombGlue, .Object@ptr)
        }
        .Object@prevNIter <- function(n = 1) {
            .Call(Algos_PrevNumCombGlue, .Object@ptr, n)
        }
        .Object@prevRemaining <- function() {
            .Call(Algos_PrevGatherGlue, .Object@ptr)
        }
        .Object@currIter <- function() {
            .Call(Algos_CurrCombGlue, .Object@ptr)
        }
        .Object@front <- function() {
            .Call(Algos_FrontGlue, .Object@ptr)
        }
        .Object@back <- function() {
            .Call(Algos_BackGlue, .Object@ptr)
        }
        .Object@randomAccess <- function(samp) {
            .Call(Algos_RandomAccessGlue, .Object@ptr, samp)
        }
        .Object@sourceVector <- function() {
            .Call(Algos_SourceVectorGlue, .Object@ptr)
        }
        .Object@summary <- function() {
            .Call(Algos_SummaryGlue, .Object@ptr)
        }
        .Object
    }
)

setClass(
    "Constraints",
    contains = "Combo"
)

# SEXP RmainFun, SEXP RcompFun,
# SEXP Rlimits, SEXP RKeepRes, SEXP Rtarget,
# SEXP Rtolerance, SEXP RmIsNull

setMethod(
    "initialize",
    "Constraints",
    function(.Object, init, Parallel, constraintFun, comparisonFun,
             limitConstraints, keepResults, tolerance, mIsNull) {
        .Object@ptr <- .Call(Algos_CombClassNew, init$RVals, init$bVec,
                             init$FreqsInfo, Parallel, NULL, NULL, NULL,
                             constraintFun, comparisonFun, limitConstraints,
                             keepResults, tolerance, mIsNull, 3)
        .Object@startOver <- function() {
            .Call(Algos_StartOverGlue, .Object@ptr)
            invisible(NULL)
        }
        .Object@nextIter <- function() {
            .Call(Algos_NextCombGlue, .Object@ptr)
        }
        .Object@nextNIter <- function(n = 1) {
            .Call(Algos_NextNumCombGlue, .Object@ptr, n)
        }
        .Object@nextRemaining <- function() {
            .Call(Algos_NextGatherGlue, .Object@ptr)
        }
        .Object@prevIter <- function() {
            .Call(Algos_PrevCombGlue, .Object@ptr)
        }
        .Object@prevNIter <- function(n = 1) {
            .Call(Algos_PrevNumCombGlue, .Object@ptr, n)
        }
        .Object@prevRemaining <- function() {
            .Call(Algos_PrevGatherGlue, .Object@ptr)
        }
        .Object@currIter <- function() {
            .Call(Algos_CurrCombGlue, .Object@ptr)
        }
        .Object@front <- function() {
            .Call(Algos_FrontGlue, .Object@ptr)
        }
        .Object@back <- function() {
            .Call(Algos_BackGlue, .Object@ptr)
        }
        .Object@randomAccess <- function(samp) {
            .Call(Algos_RandomAccessGlue, .Object@ptr, samp)
        }
        .Object@sourceVector <- function() {
            .Call(Algos_SourceVectorGlue, .Object@ptr)
        }
        .Object@summary <- function() {
            .Call(Algos_SummaryGlue, .Object@ptr)
        }
        .Object
    }
)

# setMethod(
#     "initialize",
#     "Constraints",
#     function(.Object, init, constraint, compFun, limits, keepRes, tol) {
#         .Object@ptr <- .Call(Algos_ConstraintsNew, init$RVals, init$bVec,
#                              init$FreqsInfo, constraint, compFun, limits,
#                              keepRes, tol)
#         .Object@startOver <- function() {
#             .Call(Algos_StartOverGlue, .Object@ptr)
#             invisible(NULL)
#         }
#         .Object@nextIter <- function() {
#             .Call(Algos_NextCombGlue, .Object@ptr)
#         }
#         .Object@nextNIter <- function(n = 1) {
#             .Call(Algos_NextNumCombGlue, .Object@ptr, n)
#         }
#         .Object@nextRemaining <- function() {
#             .Call(Algos_NextGatherGlue, .Object@ptr)
#         }
#         .Object@currIter <- function() {
#             .Call(Algos_CurrCombGlue, .Object@ptr)
#         }
#         .Object@sourceVector <- function() {
#             .Call(Algos_SourceVectorGlue, .Object@ptr)
#         }
#         .Object@summary <- function() {
#             .Call(Algos_SummaryGlue, .Object@ptr)
#         }
#         .Object
#     }
# )
# 
# setClass(
#     "Partitions",
#     contains = "Constraints",
#     slots = c(
#         randomAccess  = "function",
#         front         = "function",
#         back          = "function"
#     )
# )
# 
# setMethod(
#     "initialize",
#     "Partitions",
#     function(.Object, init, constraint, compFun, limits, keepRes, tol) {
#         .Object@ptr <- .Call(Algos_ConstraintsNew, init$RVals, init$bVec,
#                              init$FreqsInfo, constraint, compFun, limits,
#                              keepRes, tol)
#         .Object@startOver <- function() {
#             .Call(Algos_StartOverGlue, .Object@ptr)
#             invisible(NULL)
#         }
#         .Object@nextIter <- function() {
#             .Call(Algos_NextCombGlue, .Object@ptr)
#         }
#         .Object@nextNIter <- function(n = 1) {
#             .Call(Algos_NextNumCombGlue, .Object@ptr, n)
#         }
#         .Object@nextRemaining <- function() {
#             .Call(Algos_NextGatherGlue, .Object@ptr)
#         }
#         .Object@currIter <- function() {
#             .Call(Algos_CurrCombGlue, .Object@ptr)
#         }
#         .Object@sourceVector <- function() {
#             .Call(Algos_SourceVectorGlue, .Object@ptr)
#         }
#         .Object@summary <- function() {
#             .Call(Algos_SummaryGlue, .Object@ptr)
#         }
#         .Object
#     }
# )

"[[.Combo" <- function(x, ...) {
    x@randomAccess(...)
}

"[[.ComboApply" <- function(x, ...) {
    x@randomAccess(...)
}

## The '$' accessor is for backwards compatibility only. Moving forward, one
## should prefer using the '@' accessor to avoid calling the C function
## 'duplicate'. For more information see:
## https://cran.r-project.org/doc/manuals/R-exts.html#Profiling-R-code-for-memory-use
setMethod("$", "Combo", function(x, name) {
    function(...) slot(x, name)(...)
})

setMethod("$", "ComboApply", function(x, name) {
    function(...) slot(x, name)(...)
})

comboIter <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                      constraintFun = NULL, comparisonFun = NULL,
                      limitConstraints = NULL, keepResults = NULL,
                      FUN = NULL, Parallel = FALSE, nThreads = NULL,
                      tolerance = NULL, FUN.VALUE = NULL) {

    RetValue <- .Call(Algos_CheckReturn, v, constraintFun,
                      comparisonFun, limitConstraints,
                      keepResults, FUN)
    IsCnstrd <- .Call(Algos_CheckConstrndCpp, constraintFun,
                      comparisonFun, limitConstraints)
    InitVals <- .Call(Algos_GetClassVals, v, m, repetition, freqs,
                      TRUE, FUN, nThreads, pkgEnv$nThreads, IsCnstrd)

    if (RetValue == 1) {
        new("Combo", InitVals, Parallel)
    } else if (RetValue == 2) {
        new("ComboApply", InitVals, FUN, new.env(), FUN.VALUE)
    } else {
        new("Constraints", InitVals, Parallel, constraintFun, comparisonFun,
            limitConstraints, keepResults, tolerance, is.null(m))
    }
}

permuteIter <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                        constraintFun = NULL, comparisonFun = NULL,
                        limitConstraints = NULL, keepResults = NULL,
                        FUN = NULL, Parallel = FALSE, nThreads = NULL,
                        tolerance = NULL, FUN.VALUE = NULL) {
    
    RetValue <- .Call(Algos_CheckReturn, v, constraintFun,
                      comparisonFun, limitConstraints,
                      keepResults, FUN)
    IsCnstrd <- .Call(Algos_CheckConstrndCpp, constraintFun,
                      comparisonFun, limitConstraints)
    InitVals <- .Call(Algos_GetClassVals, v, m, repetition, freqs,
                      FALSE, FUN, nThreads, pkgEnv$nThreads, IsCnstrd)
    
    if (RetValue == 1) {
        new("Combo", InitVals, Parallel)
    } else if (RetValue == 2) {
        new("ComboApply", InitVals, FUN, new.env(), FUN.VALUE)
    } else {
        new("Constraints", InitVals, Parallel, constraintFun, comparisonFun,
            limitConstraints, keepResults, tolerance, is.null(m))
    }
}
