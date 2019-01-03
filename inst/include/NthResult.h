#ifndef NTH_RESULT_H
#define NTH_RESULT_H

#include <vector>
#include <gmp.h>

std::vector<int> nthCombination(const int n, const int r, double myIndex, bool isRep,
                                bool isMult, std::vector<int> Reps);

std::vector<int> nthPermutation(const int n, const int r, double myIndex, bool isRep,
                                bool isMult, std::vector<int> Reps, std::vector<int> freqs,
                                bool isStarter = false);

std::vector<int> nthCombinationGmp(const int n, const int r, mpz_t myIndex, bool isRep,
                                   bool isMult, std::vector<int> Reps);

std::vector<int> nthPermutationGmp(const int n, const int r, mpz_t myIndex, bool isRep,
                                   bool isMult, std::vector<int> Reps, std::vector<int> freqs,
                                   bool isStarter = false);

#endif
