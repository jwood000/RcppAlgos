ALGOS_METHODS <- c(
    startOver     = ".Object@startOver <- function() {
        .Call(`_RcppAlgos_StartOverGlue`, .Object@ptr)
        invisible(NULL)}",
    nextIter      = ".Object@nextIter <- function() {
        .Call(`_RcppAlgos_NextIterGlue`, .Object@ptr)}",
    nextNIter     = ".Object@nextNIter <- function(n = 1) {
        .Call(`_RcppAlgos_NextNumIterGlue`, .Object@ptr, n)}",
    nextRemaining = ".Object@nextRemaining <- function() {
        .Call(`_RcppAlgos_NextGatherGlue`, .Object@ptr)}",
    prevIter      = ".Object@prevIter <- function() {
        .Call(`_RcppAlgos_PrevIterGlue`, .Object@ptr)}",
    prevNIter     = ".Object@prevNIter <- function(n = 1) {
        .Call(`_RcppAlgos_PrevNumIterGlue`, .Object@ptr, n)}",
    prevRemaining = ".Object@prevRemaining <- function() {
        .Call(`_RcppAlgos_PrevGatherGlue`, .Object@ptr)}",
    currIter      = ".Object@currIter <- function() {
        .Call(`_RcppAlgos_CurrIterGlue`, .Object@ptr)}",
    randomAccess  = ".Object@randomAccess <- function(samp) {
        .Call(`_RcppAlgos_RandomAccessGlue`, .Object@ptr, samp)}",
    sourceVector  = ".Object@sourceVector <- function() {
        .Call(`_RcppAlgos_SourceVectorGlue`, .Object@ptr)}",
    front         = ".Object@front <- function() {
        .Call(`_RcppAlgos_FrontGlue`, .Object@ptr)}",
    back          = ".Object@back <- function() {
        .Call(`_RcppAlgos_BackGlue`, .Object@ptr)}",
    summary       = ".Object@summary <- function() {
        .Call(`_RcppAlgos_SummaryGlue`, .Object@ptr)}",
    object        = ".Object"
)

ALGOS_METHODS_SAN_PREV <- ALGOS_METHODS[c(
    "startOver", "nextIter", "nextNIter", "nextRemaining",
    "currIter", "randomAccess", "sourceVector", "front",
    "back", "summary", "object")
]

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

setClass(
    "ComboApply",
    contains = "Combo"
)

setClass(
    "ComboRes",
    contains = "Combo"
)

setClass(
    "Cartesian",
    slots = c(
        ptr           = "externalptr",
        startOver     = "function",
        nextIter      = "function",
        nextNIter     = "function",
        nextRemaining = "function",
        currIter      = "function",
        randomAccess  = "function",
        sourceVector  = "function",
        front         = "function",
        back          = "function",
        summary       = "function"
    )
)

setClass(
    "Constraints",
    slots = c(
        ptr           = "externalptr",
        startOver     = "function",
        nextIter      = "function",
        nextNIter     = "function",
        nextRemaining = "function",
        currIter      = "function",
        sourceVector  = "function",
        summary       = "function"
    )
)

setClass(
    "Partitions",
    slots = c(
        ptr           = "externalptr",
        startOver     = "function",
        nextIter      = "function",
        nextNIter     = "function",
        nextRemaining = "function",
        currIter      = "function",
        randomAccess  = "function",
        sourceVector  = "function",
        front         = "function",
        back          = "function",
        summary       = "function"
    )
)

setClass(
    "ComboGroups",
    slots = c(
        ptr           = "externalptr",
        startOver     = "function",
        nextIter      = "function",
        nextNIter     = "function",
        nextRemaining = "function",
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
        .Object@ptr <- .Call(`_RcppAlgos_CombClassNew`, init$RVals, init$bVec,
                             init$FreqsInfo, Parallel, NULL, NULL, NULL,
                             NULL, NULL, NULL, NULL, NULL, NULL, 1)
        eval(str2expression(text = ALGOS_METHODS))
    }
)

setMethod(
    "initialize",
    "ComboApply",
    function(.Object, init, stdFun, rho, RFunVal) {
        .Object@ptr <- .Call(`_RcppAlgos_CombClassNew`, init$RVals, init$bVec,
                             init$FreqsInfo, FALSE, stdFun, rho, RFunVal,
                             NULL, NULL, NULL, NULL, NULL, NULL, 2)
        eval(str2expression(text = ALGOS_METHODS))
    }
)

setMethod(
    "initialize",
    "ComboRes",
    function(.Object, init, Parallel, constraintFun, comparisonFun,
             limitConstraints, keepResults, tolerance, mIsNull) {
        .Object@ptr <- .Call(`_RcppAlgos_CombClassNew`, init$RVals, init$bVec,
                             init$FreqsInfo, Parallel, NULL, NULL, NULL,
                             constraintFun, comparisonFun, limitConstraints,
                             keepResults, tolerance, mIsNull, 3)
        eval(str2expression(text = ALGOS_METHODS))
    }
)

setMethod(
    "initialize",
    "Constraints",
    function(.Object, init, Parallel, constraintFun, comparisonFun,
             limitConstraints, keepResults, tolerance, mIsNull) {
        .Object@ptr <- .Call(`_RcppAlgos_CombClassNew`, init$RVals, init$bVec,
                             init$FreqsInfo, Parallel, NULL, NULL, NULL,
                             constraintFun, comparisonFun, limitConstraints,
                             keepResults, tolerance, mIsNull, 3)
        eval(
            str2expression(
                text = ALGOS_METHODS[c(
                    "startOver", "nextIter", "nextNIter", "nextRemaining",
                    "currIter", "sourceVector", "summary", "object")]
            )
        )
    }
)

setMethod(
    "initialize",
    "Partitions",
    function(.Object, init, Parallel, constraintFun, comparisonFun,
             limitConstraints, keepResults, tolerance, mIsNull) {
        .Object@ptr <- .Call(`_RcppAlgos_CombClassNew`, init$RVals, init$bVec,
                             init$FreqsInfo, Parallel, NULL, NULL, NULL,
                             constraintFun, comparisonFun, limitConstraints,
                             keepResults, tolerance, mIsNull, 3)
        eval(str2expression(text = ALGOS_METHODS_SAN_PREV))
    }
)

setMethod(
    "initialize",
    "ComboGroups",
    function(.Object, init, Parallel) {
        .Object@ptr <- .Call(`_RcppAlgos_CombClassNew`, init$RVals, init$bVec,
                             init$FreqsInfo, Parallel, NULL, NULL, NULL,
                             NULL, NULL, NULL, NULL, NULL, NULL, 4)
        eval(str2expression(text = ALGOS_METHODS_SAN_PREV))
    }
)

setMethod(
    "initialize",
    "Cartesian",
    function(.Object, Rlist, nThreads, Return_DF) {
        .Object@ptr <- .Call(
            `_RcppAlgos_CartClassNew`, Rlist, nThreads,
            pkgEnv$nThreads, Return_DF
        )
        eval(str2expression(text = ALGOS_METHODS_SAN_PREV))
    }
)

"[[.Combo" <- function(x, ...) {
    x@randomAccess(...)
}

"[[.ComboApply" <- function(x, ...) {
    x@randomAccess(...)
}

"[[.ComboRes" <- function(x, ...) {
    x@randomAccess(...)
}

"[[.Partitions" <- function(x, ...) {
    x@randomAccess(...)
}

"[[.ComboGroups" <- function(x, ...) {
    x@randomAccess(...)
}

"[[.Cartesian" <- function(x, ...) {
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

setMethod("$", "ComboRes", function(x, name) {
    function(...) slot(x, name)(...)
})

setMethod("$", "Constraints", function(x, name) {
    function(...) slot(x, name)(...)
})

setMethod("$", "Partitions", function(x, name) {
    function(...) slot(x, name)(...)
})

setMethod("$", "ComboGroups", function(x, name) {
    function(...) slot(x, name)(...)
})

setMethod("$", "Cartesian", function(x, name) {
    function(...) slot(x, name)(...)
})
