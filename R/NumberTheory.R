primeFactorizeSieve <- function(bound1, bound2 = NULL,
                                namedList = FALSE, nThreads = NULL) {
    return(.Call(MotleyContainer, bound1, bound2, FALSE, namedList,
                 nThreads, pkgEnv$nThreads, PACKAGE = "RcppAlgos"))
}

eulerPhiSieve <- function(bound1, bound2 = NULL,
                          namedVector = FALSE, nThreads = NULL) {
    return(.Call(MotleyContainer, bound1, bound2, TRUE, namedVector,
                 nThreads, pkgEnv$nThreads, PACKAGE = "RcppAlgos"))
}

primeSieve <- function(bound1, bound2 = NULL, nThreads = NULL) {
    return(.Call(PrimeSieveCpp, bound1, bound2, nThreads,
                 pkgEnv$nCores, pkgEnv$nThreads, PACKAGE = "RcppAlgos"))
}

divisorsSieve <- function(bound1, bound2 = NULL, namedList = FALSE, nThreads = NULL) {
    return(.Call(DivNumSieveCpp, bound1, bound2, TRUE, namedList, 
                 nThreads, pkgEnv$nThreads, PACKAGE = "RcppAlgos"))
}

numDivisorSieve <- function(bound1, bound2 = NULL, namedVector = FALSE, nThreads = NULL) {
    return(.Call(DivNumSieveCpp, bound1, bound2, FALSE, namedVector,
                 nThreads, pkgEnv$nThreads, PACKAGE = "RcppAlgos"))
}

# primeFactorize <- function(v, namedList = FALSE, nThreads = NULL) {
#     PollardRhoContainer(v, namedList, TRUE, FALSE, nThreads, pkgEnv$nThreads)
# }
# 
# divisorsRcpp <- function(v, namedList = FALSE, nThreads = NULL) {
#     PollardRhoContainer(v, namedList, FALSE, TRUE, nThreads, pkgEnv$nThreads)
# }
# 
# isPrimeRcpp <- function(v, namedVector = FALSE, nThreads = NULL) {
#     PollardRhoContainer(v, namedVector, FALSE, FALSE, nThreads, pkgEnv$nThreads)
# }

primeCount <- function(n, nThreads = NULL) {
    return(.Call(PrimeCountCpp, n, nThreads, pkgEnv$nThreads))
}
