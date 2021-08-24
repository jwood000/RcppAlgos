#ifndef POLLARD_RHO_DEPENDS_H
#define POLLARD_RHO_DEPENDS_H

#include <vector>
#include <cmath>

template <typename T>
void PrimeFacList(std::size_t m, std::size_t n,
                  const std::vector<double> &myNums,
                  std::vector<std::vector<T>> &MyPrimeList);

template <typename T>
std::vector<T> Factorize(std::vector<T> &factors);

template <typename T>
void FactorList(std::size_t m, std::size_t n,
                const std::vector<double> &myNums,
                std::vector<std::vector<T>> &MyDivList);

void IsPrimeVec(std::size_t m, std::size_t n,
                const std::vector<double> &myNums,
                int* primeTest);

#endif
