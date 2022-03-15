stdThreadMax <- function() {
    nThreads <- .Call(Algos_cpp11GetNumThreads)
    if (nThreads < 1L || is.na(nThreads) || is.null(nThreads)) nThreads = 1L
    return(nThreads)
}
