expandGrid <- function(..., lower = NULL, upper = NULL, nThreads = NULL) {

    lst <- GridInputs(...)

    if (length(lst$i) == lst$n) {
        return(as.data.frame(lst$p))
    }

    res <- .Call(
        `_RcppAlgos_ExpandGridCpp`, lst$p, lower,
        upper, nThreads, pkgEnv$nThreads
    )

    if (is.matrix(res)) {
        colnames(res) <- nmc
    }

    return(res)
}
