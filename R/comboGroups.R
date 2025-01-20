comboGroups <- function(v, numGroups = NULL, grpSizes = NULL,
                        retType = "matrix", lower = NULL, upper = NULL,
                        Parallel = FALSE, nThreads = NULL) {

    return(.Call(`_RcppAlgos_ComboGroupsCpp`, v, numGroups, grpSizes, retType,
                 lower, upper, Parallel, nThreads, pkgEnv$nThreads, FALSE,
                 NULL, NULL, NULL, sample, FALSE, NULL))
}

comboGroupsSample <- function(v, numGroups = NULL, grpSizes = NULL,
                              retType = "matrix", n = NULL, sampleVec = NULL,
                              seed = NULL, Parallel = FALSE, nThreads = NULL,
                              namedSample = FALSE) {

    if (!is.null(seed)) {
        set.seed(seed)
    }

    return(.Call(`_RcppAlgos_ComboGroupsCpp`, v, numGroups, grpSizes, retType,
                 NULL, NULL, Parallel, nThreads, pkgEnv$nThreads, TRUE,
                 sampleVec, seed, n, sample, namedSample, new.env()))
}

comboGroupsCount <- function(v, numGroups = NULL, grpSizes = NULL) {
    return(.Call(`_RcppAlgos_ComboGroupsCountCpp`, v, numGroups, grpSizes))
}

comboGroupsIter <- function(v, numGroups = NULL, grpSizes = NULL,
                            retType = "matrix", Parallel = FALSE,
                            nThreads = NULL) {

    InitVals <- .Call(`_RcppAlgos_GetClassVals`, v, length(GetV(v)),
                      FALSE, NULL, TRUE, NULL, nThreads, pkgEnv$nThreads,
                      FALSE, FALSE, FALSE, numGroups, grpSizes, retType)

    new("ComboGroups", InitVals, Parallel)
}
