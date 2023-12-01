ComboPermuteGen <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, lower = NULL, upper = NULL,
    constraintFun = NULL, comparisonFun = NULL, limitConstraints = NULL,
    keepResults = NULL, FUN = NULL, Parallel = FALSE, nThreads = NULL,
    tolerance = NULL, FUN.VALUE = NULL, IsComb = TRUE, ...
) {

    ## Credit to user2554330:
    ##      https://stackoverflow.com/a/77573995/4408538
    fun_pass <- if (!is.null(FUN)) {
        FUN <- match.fun(FUN)

        if (as.character(quote(FUN)) != "function" && length(list(...))) {
            force(FUN)
            function(x) FUN(x, ...)
        } else {
            FUN
        }
    }

    RetValue <- .Call(`_RcppAlgos_CheckReturn`, v, constraintFun,
                      comparisonFun, limitConstraints,
                      keepResults, fun_pass)

    if (RetValue == 1) {
        return(.Call(`_RcppAlgos_CombinatoricsStndrd`, v, m, repetition,
                     freqs, lower, upper, Parallel, nThreads,
                     pkgEnv$nThreads, IsComb))
    } else if (RetValue == 2) {
        return(.Call(`_RcppAlgos_CombinatoricsApply`, v, m,
                     repetition, freqs, lower, upper,
                     fun_pass, new.env(), FUN.VALUE, IsComb))
    } else {
        return(.Call(`_RcppAlgos_CombinatoricsCnstrt`, v, m, repetition,
                     freqs, lower, upper, constraintFun, comparisonFun,
                     limitConstraints, IsComb, keepResults, Parallel,
                     nThreads, pkgEnv$nThreads, tolerance, FALSE, FALSE))
    }
}

ComboPermuteSample <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, n = NULL, sampleVec = NULL,
    seed = NULL, FUN = NULL, Parallel = FALSE, nThreads = NULL,
    namedSample = FALSE, FUN.VALUE = NULL, IsComb = TRUE
) {

    if (!is.null(seed)) {
        set.seed(seed)
    }

    if (!is.null(FUN)) {
        FUN <- match.fun(FUN)
    }

    return(.Call(
        `_RcppAlgos_SampleCombPerm`, v, m, repetition, freqs, sampleVec,
        IsComb, seed, n, sample, FUN, new.env(), Parallel, nThreads,
        pkgEnv$nThreads, namedSample, FUN.VALUE
    ))
}

ComboPermuteCount <-  function(v, m = NULL, repetition = FALSE,
                            freqs = NULL, IsComb = TRUE) {
    .Call(`_RcppAlgos_CombinatoricsCount`, v, m,
          repetition, freqs, IsComb);
}

ComboPermuteIter <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, constraintFun = NULL,
    comparisonFun = NULL, limitConstraints = NULL, keepResults = NULL,
    FUN = NULL, Parallel = FALSE, nThreads = NULL, tolerance = NULL,
    FUN.VALUE = NULL, IsComb = TRUE
) {

    if (!is.null(FUN)) {
        FUN <- match.fun(FUN)
    }

    RetValue <- .Call(`_RcppAlgos_CheckReturn`, v, constraintFun,
                      comparisonFun, limitConstraints,
                      keepResults, FUN)
    IsCnstrd <- .Call(`_RcppAlgos_CheckConstrndCpp`, constraintFun,
                      comparisonFun, limitConstraints)
    InitVals <- .Call(`_RcppAlgos_GetClassVals`, v, m, repetition,
                      freqs, IsComb, FUN, nThreads, pkgEnv$nThreads,
                      IsCnstrd, FALSE, FALSE, NULL, NULL, NULL)

    if (RetValue == 1) {
        new("Combo", InitVals, Parallel)
    } else if (RetValue == 2) {
        new("ComboApply", InitVals, FUN, new.env(), FUN.VALUE)
    } else if (IsCnstrd) {
        new("Constraints", InitVals, Parallel, constraintFun, comparisonFun,
            limitConstraints, keepResults, tolerance, is.null(m))
    } else {
        new("ComboRes", InitVals, Parallel, constraintFun, comparisonFun,
            limitConstraints, keepResults, tolerance, is.null(m))
    }
}
