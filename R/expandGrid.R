expandGrid <- function(..., lower = NULL, upper = NULL, nThreads = NULL) {

    lst <- GridInputs(...)

    if (lst$early_return) {
        return(lst$res)
    }

    res <- .Call(
        `_RcppAlgos_ExpandGridCpp`, lst$p, lower, upper, nThreads,
        pkgEnv$nThreads, FALSE, NULL, NULL, NULL, sample, FALSE, NULL
    )

    if (is.matrix(res)) {
        colnames(res) <- lst$nmc
    }

    return(res)
}

expandGridSample <- function(..., n = NULL, sampleVec = NULL, seed = NULL,
                             nThreads = NULL, namedSample = FALSE) {

    if (!is.null(seed)) {
        set.seed(seed)
    }

    lst <- GridInputs(...)

    if (lst$early_return) {
        return(lst$res)
    }

    res <- .Call(
        `_RcppAlgos_ExpandGridCpp`, lst$p, NULL, NULL, nThreads,
        pkgEnv$nThreads, TRUE, sampleVec, seed, n, sample,
        namedSample, new.env()
    )

    if (is.matrix(res)) {
        colnames(res) <- lst$nmc
    }

    return(res)
}

expandGridCount <- function(...) {

    lst <- GridInputs(...)

    if (lst$early_return) {
        return(nrow(lst$res))
    }

    return(.Call(`_RcppAlgos_ExpandGridCountCpp`, lst$p))
}

expandGridIter <- function(..., nThreads = NULL) {

    lst <- GridInputs(...)

    if (lst$early_return) {
        return(nrow(lst$res))
    }

    new("Cartesian", lst$p, nThreads)
}
