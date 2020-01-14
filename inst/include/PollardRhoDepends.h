#ifndef POLLARD_RHO_DEPENDS_H
#define POLLARD_RHO_DEPENDS_H

#include <Rcpp.h>
#include <cmath>

template <typename typeReturn>
void PrimeFacList(std::size_t m, std::size_t n, std::vector<double> myNums,
                  std::vector<std::vector<typeReturn>> &MyPrimeList);

template <typename typeReturn>
std::vector<typeReturn> Factorize(std::vector<typeReturn> &factors);

template <typename typeReturn>
void FactorList(std::size_t m, std::size_t n, std::vector<double> &myNums,
                std::vector<std::vector<typeReturn>> &MyDivList);

void IsPrimeVec(std::size_t m, std::size_t n, std::vector<double> &myNums,
                Rcpp::LogicalVector &primeTest);

#endif
