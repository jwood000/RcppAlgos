expandGrid <- function(..., lower = NULL, upper = NULL,
                       nThreads = NULL, return_df = FALSE) {

    lst <- GridInputs(...)

    if (lst$early_return) {
        return(lst$res)
    }

    res <- .Call(
        `_RcppAlgos_ExpandGridCpp`, lst$p, lower, upper,
        nThreads, pkgEnv$nThreads, FALSE, NULL, NULL, NULL,
        sample, FALSE, NULL, return_df
    )

    return(res)
}

expandGridSample <- function(
    ..., n = NULL, sampleVec = NULL, seed = NULL,
    nThreads = NULL, namedSample = FALSE, return_df = FALSE
) {

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
        namedSample, new.env(), return_df
    )

    return(res)
}

expandGridCount <- function(...) {

    lst <- GridInputs(...)

    if (lst$early_return) {
        return(nrow(lst$res))
    }

    return(.Call(`_RcppAlgos_ExpandGridCountCpp`, lst$p))
}

expandGridIter <- function(..., nThreads = NULL, return_df = FALSE) {

    lst <- GridInputs(...)

    if (lst$early_return) {
        return(nrow(lst$res))
    }

    new("Cartesian", lst$p, nThreads, return_df)
}
