primeFactorizeSieve <- function(bound1, bound2 = NULL, namedList = FALSE, nThreads = NULL) {
    MotleyContainer(bound1, bound2, FALSE, namedList, nThreads, pkgEnv$nThreads)
}

eulerPhiSieve <- function(bound1, bound2 = NULL, namedVector = FALSE, nThreads = NULL) {
    MotleyContainer(bound1, bound2, TRUE, namedVector, nThreads, pkgEnv$nThreads)
}

primeSieve <- function(bound1, bound2 = NULL, nThreads = NULL) {
    EratosthenesRcpp(bound1, bound2, nThreads, pkgEnv$nCores, pkgEnv$nThreads)
}

divisorsSieve <- function(bound1, bound2 = NULL, namedList = FALSE, nThreads = NULL) {
    DivNumSieve(bound1, bound2, TRUE, namedList, nThreads, pkgEnv$nThreads)
}

numDivisorSieve <- function(bound1, bound2 = NULL, namedVector = FALSE, nThreads = NULL) {
    DivNumSieve(bound1, bound2, FALSE, namedVector, nThreads, pkgEnv$nThreads)
}

primeFactorize <- function(v, namedList = FALSE, nThreads = NULL) {
    PollardRhoContainer(v, namedList, TRUE, FALSE, nThreads, pkgEnv$nThreads)
}

divisorsRcpp <- function(v, namedList = FALSE, nThreads = NULL) {
    PollardRhoContainer(v, namedList, FALSE, TRUE, nThreads, pkgEnv$nThreads)
}

isPrimeRcpp <- function(v, namedVector = FALSE, nThreads = NULL) {
    PollardRhoContainer(v, namedVector, FALSE, FALSE, nThreads, pkgEnv$nThreads)
}

primeCount <- function(n, nThreads = NULL) {
    PrimeCountRcpp(n, nThreads, pkgEnv$nThreads)
}

stdThreadMax <- function() {
    nThreads <- cpp11GetNumThreads()
    if (nThreads < 1L || is.na(nThreads) || is.null(nThreads)) nThreads = 1L
    nThreads
}
