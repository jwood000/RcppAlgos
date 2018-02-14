primeFactorizationList <- function(bound1 = 100L, bound2 = NULL, namedVector = FALSE) {
    EratosthenesRcpp(bound1, bound2, TRUE, FALSE, namedVector)
}

eulerPhiSieve <- function(bound1 = 100L, bound2 = NULL, namedVector = FALSE) {
    EratosthenesRcpp(bound1, bound2, FALSE, TRUE, namedVector)
}

primeSieve <- function(bound1 = 100L, bound2 = NULL) {
    EratosthenesRcpp(bound1, bound2, FALSE, FALSE, FALSE)
}

divisorsList <- function(bound1 = 100L, bound2 = NULL, namedVector = FALSE) {
    DivisorsGeneral(bound1, bound2, TRUE, namedVector)
}

numDivisorSieve <- function(bound1 = 100L, bound2 = NULL, namedVector = FALSE) {
    DivisorsGeneral(bound1, bound2, FALSE, namedVector)
}

primeFactorize <- function(n = 100L) {PrimeFactorsContainer(n)}

factorizeAllRcpp <- function(n = 100L) {getAllDivisorsRcpp(n)}

isPrimeRcpp <- function(n = 100L) {IsPrimeContainer(n)}