pkgEnv <- new.env(parent = emptyenv())
pkgEnv$nCores <- NULL
pkgEnv$nThreads <- NULL

setPkgVars <- function() {
    if (is.null(pkgEnv$nCores)) pkgEnv$nCores <- physicalCoreCount()
    if (is.null(pkgEnv$nThreads)) {
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

divisorsSieve <- function(bound1 = 100L, bound2 = NULL, namedList = FALSE) {
    DivisorsGeneral(bound1, bound2, TRUE, namedList)
}

numDivisorSieve <- function(bound1 = 100L, bound2 = NULL, namedVector = FALSE) {
    DivisorsGeneral(bound1, bound2, FALSE, namedVector)
}

primeFactorize <- function(v = 100L, namedList = FALSE) {PrimeFactorsContainer(v, namedList)}
divisorsRcpp <- function(v = 100L, namedList = FALSE) {getAllDivisorsRcpp(v, namedList)}
isPrimeRcpp <- function(v = 100L, namedVector = FALSE) {IsPrimeContainer(v, namedVector)}
primeCount <- function(n = 100L) {PrimeCountRcpp(n)}
