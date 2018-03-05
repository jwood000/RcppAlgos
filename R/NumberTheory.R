primeFactorizeSieve <- function(bound1 = 100L, bound2 = NULL, namedList = FALSE) {
    EratosthenesRcpp(bound1, bound2, TRUE, FALSE, namedList)
}

eulerPhiSieve <- function(bound1 = 100L, bound2 = NULL, namedVector = FALSE) {
    EratosthenesRcpp(bound1, bound2, FALSE, TRUE, namedVector)
}

primeSieve <- function(bound1 = 100L, bound2 = NULL) {
    EratosthenesRcpp(bound1, bound2, FALSE, FALSE, FALSE)
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
primeCount <- function(n = 100L) {MasterPrimeCount(n)}