primeFactorizationList <- function(n = 100L) {PrimeFactorizationListRcpp(n)}

eulerPhiSieve <- function(n = 100L) {EulerPhiSieveRcpp(n)}

primeSieve <- function(bound1 = 100L, bound2 = NULL) {EratosthenesRcpp(bound1, bound2)}

divisorsList <- function(n) {DivisorListRcpp(n)}

numDivisorSieve <- function(n) {NumDivisorsSieve(n)}