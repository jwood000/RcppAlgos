pkgEnv <- new.env(parent = emptyenv())
pkgEnv$nCores <- NULL
pkgEnv$nThreads <- NULL

setPkgVars <- function() {
    if (is.null(pkgEnv$nCores)) {
        pkgEnv$nCores <- physicalCoreCount()
        tempThreads <- parallel::detectCores()
        if (is.na(tempThreads))
            pkgEnv$nThreads <- 1L
        else
            pkgEnv$nThreads <- tempThreads
    }
}

primeFactorizeSieve <- function(bound1 = 100L, bound2 = NULL, namedList = FALSE, nThreads = NULL) {
    setPkgVars()
    MotleyContainer(bound1, bound2, FALSE, namedList, nThreads, pkgEnv$nThreads)
}

eulerPhiSieve <- function(bound1 = 100L, bound2 = NULL, namedVector = FALSE, nThreads = NULL) {
    setPkgVars()
    MotleyContainer(bound1, bound2, TRUE, namedVector, nThreads, pkgEnv$nThreads)
}

primeSieve <- function(bound1 = 100L, bound2 = NULL, nThreads = NULL) {
    setPkgVars()
    EratosthenesRcpp(bound1, bound2, nThreads, pkgEnv$nCores, pkgEnv$nThreads)
}

divisorsSieve <- function(bound1 = 100L, bound2 = NULL, namedList = FALSE, nThreads = NULL) {
    setPkgVars()
    DivNumSieve(bound1, bound2, TRUE, namedList, nThreads, pkgEnv$nThreads)
}

numDivisorSieve <- function(bound1 = 100L, bound2 = NULL, namedVector = FALSE, nThreads = NULL) {
    setPkgVars()
    DivNumSieve(bound1, bound2, FALSE, namedVector, nThreads, pkgEnv$nThreads)
}

primeFactorize <- function(v, namedList = FALSE, nThreads = NULL) {
    setPkgVars()
    PollardRhoContainer(v, namedList, TRUE, FALSE, nThreads, pkgEnv$nThreads)
}

divisorsRcpp <- function(v, namedList = FALSE, nThreads = NULL) {
    setPkgVars()
    PollardRhoContainer(v, namedList, FALSE, TRUE, nThreads, pkgEnv$nThreads)
}

isPrimeRcpp <- function(v, namedVector = FALSE, nThreads = NULL) {
    setPkgVars()
    PollardRhoContainer(v, namedVector, FALSE, FALSE, nThreads, pkgEnv$nThreads)
}

primeCount <- function(n) {PrimeCountRcpp(n)}
