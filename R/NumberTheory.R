RcppAlgosPkgEnv <- new.env(parent=emptyenv())
RcppAlgosPkgEnv$NumCores <- NULL
RcppAlgosPkgEnv$NumThreads <- NULL

primeFactorizeSieve <- function(bound1 = 100L, bound2 = NULL, namedList = FALSE) {
    MotleyPrimes(bound1, bound2, TRUE, namedList, NULL)
}

eulerPhiSieve <- function(bound1 = 100L, bound2 = NULL, namedVector = FALSE) {
    MotleyPrimes(bound1, bound2, FALSE, namedVector, NULL)
}

primeSieve <- function(bound1 = 100L, bound2 = NULL, nThreads = NULL) {
    
    if (is.null(RcppAlgosPkgEnv$NumCores)) RcppAlgosPkgEnv$NumCores <- physicalCoreCount()
    if (is.null(RcppAlgosPkgEnv$NumThreads)) RcppAlgosPkgEnv$NumThreads <- parallel::detectCores()
    maxCores = RcppAlgosPkgEnv$NumCores
    maxThreads = RcppAlgosPkgEnv$NumThreads
    if (is.na(maxThreads)) maxThreads = 1L
    
    EratosthenesRcpp(bound1, bound2, nThreads, maxCores, maxThreads)
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
