primeFactorizeSieve <- function(bound1, bound2 = NULL,
                                namedList = FALSE, nThreads = NULL) {
    return(.Call(`_RcppAlgos_MotleyContainer`, bound1, bound2, FALSE, namedList,
                 nThreads, pkgEnv$maxThreads))
}

eulerPhiSieve <- function(bound1, bound2 = NULL,
                          namedVector = FALSE, nThreads = NULL) {
    return(.Call(`_RcppAlgos_MotleyContainer`, bound1, bound2, TRUE, namedVector,
                 nThreads, pkgEnv$maxThreads))
}

primeSieve <- function(bound1, bound2 = NULL, nThreads = NULL) {
    return(.Call(`_RcppAlgos_PrimeSieveCpp`, bound1, bound2, nThreads,
                 pkgEnv$nCores, pkgEnv$maxThreads))
}

divisorsSieve <- function(bound1, bound2 = NULL,
                          namedList = FALSE, nThreads = NULL) {
    return(.Call(`_RcppAlgos_DivNumSieveCpp`, bound1, bound2, TRUE, namedList,
                 nThreads, pkgEnv$maxThreads))
}

numDivisorSieve <- function(bound1, bound2 = NULL,
                            namedVector = FALSE, nThreads = NULL) {
    return(.Call(`_RcppAlgos_DivNumSieveCpp`, bound1, bound2, FALSE, namedVector,
                 nThreads, pkgEnv$maxThreads))
}

primeFactorize <- function(v, namedList = FALSE, nThreads = NULL) {
    return(.Call(`_RcppAlgos_PollardRhoContainer`, v, namedList, TRUE,
                 FALSE, nThreads, pkgEnv$maxThreads))
}

divisorsRcpp <- function(v, namedList = FALSE, nThreads = NULL) {
    return(.Call(`_RcppAlgos_PollardRhoContainer`, v, namedList, FALSE,
                 TRUE, nThreads, pkgEnv$maxThreads))
}

isPrimeRcpp <- function(v, namedVector = FALSE, nThreads = NULL) {
    return(.Call(`_RcppAlgos_PollardRhoContainer`, v, namedVector, FALSE,
                 FALSE, nThreads, pkgEnv$maxThreads))
}

primeCount <- function(n, nThreads = NULL) {
    return(.Call(`_RcppAlgos_PrimeCountCpp`, n, nThreads, pkgEnv$maxThreads))
}
