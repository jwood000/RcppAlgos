stdThreadMax <- function() {
    nThreads <- .Call(cpp11GetNumThreads, PACKAGE = "RcppAlgos")
    if (nThreads < 1L || is.na(nThreads) || is.null(nThreads)) nThreads = 1L
    return(nThreads)
}