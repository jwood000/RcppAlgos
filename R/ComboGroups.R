comboGroups <- function(v, numGroups, retType = "matrix", lower = NULL,
                        upper = NULL, Parallel = FALSE, nThreads = NULL) {

    return(.Call(Algos_ComboGroupsCpp, v, numGroups, retType, lower,
                 upper, Parallel, nThreads, pkgEnv$nThreads, FALSE, NULL,
                 NULL, NULL, sample, FALSE, NULL))
}

comboGroupsSample <- function(v, numGroups, retType = "matrix",
                              n = NULL, sampleVec = NULL, seed = NULL,
                              Parallel = FALSE, nThreads = NULL,
                              namedSample = FALSE) {

    if (!is.null(seed)) {
        set.seed(seed)
    }

    return(.Call(Algos_ComboGroupsCpp, v, numGroups, retType, NULL, NULL,
                 Parallel, nThreads, pkgEnv$nThreads, TRUE, sampleVec,
                 seed, n, sample, namedSample, new.env()))
}

comboGroupsCount <- function(v, numGroups) {
    return(.Call(Algos_ComboGroupsCountCpp, v, numGroups))
}