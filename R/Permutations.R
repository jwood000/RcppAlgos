permuteGeneral <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                           lower = NULL, upper = NULL, constraintFun = NULL,
                           comparisonFun = NULL, limitConstraints = NULL,
                           keepResults = NULL, FUN = NULL, Parallel = FALSE,
                           nThreads = NULL, tolerance = NULL,
                           FUN.VALUE = NULL) {

    RetValue <- .Call(Algos_CheckReturn, v, constraintFun,
                      comparisonFun, limitConstraints,
                      keepResults, FUN)

    if (RetValue == 1) {
        return(.Call(Algos_CombinatoricsStndrd, v, m, repetition,
                     freqs, lower, upper, Parallel, nThreads,
                     pkgEnv$nThreads, FALSE))
    } else if (RetValue == 2) {
        return(.Call(Algos_CombinatoricsApply, v, m,
                     repetition, freqs, lower, upper,
                     FUN, new.env(), FUN.VALUE, FALSE))
    } else {
        return(.Call(Algos_CombinatoricsCnstrt, v, m, repetition,
                     freqs, lower, upper, constraintFun, comparisonFun,
                     limitConstraints, FALSE, keepResults, Parallel,
                     nThreads, pkgEnv$nThreads, tolerance))
    }
}

permuteSample <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                          n = NULL, sampleVec = NULL, seed = NULL,
                          FUN = NULL, Parallel = FALSE, nThreads = NULL,
                          namedSample = FALSE, FUN.VALUE = NULL) {

    if (!is.null(seed)) {
        set.seed(seed)
    }

    return(.Call(Algos_SampleCombPerm, v, m, repetition, freqs, sampleVec,
                 FALSE, seed, n, sample, FUN, new.env(), Parallel, nThreads,
                 pkgEnv$nThreads, namedSample, FUN.VALUE))
}

permuteCount <- function(v, m = NULL, repetition = FALSE, freqs = NULL) {
    .Call(Algos_CombinatoricsCount, v, m,
          repetition, freqs, FALSE);
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
    } else if (IsCnstrd) {
        new("Constraints", InitVals, Parallel, constraintFun, comparisonFun,
            limitConstraints, keepResults, tolerance, is.null(m))
    } else {
        new("ComboRes", InitVals, Parallel, constraintFun, comparisonFun,
            limitConstraints, keepResults, tolerance, is.null(m))
    }
}