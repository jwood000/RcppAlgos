#ifndef NTH_RESULT_H
#define NTH_RESULT_H

#include <vector>
#include <gmp.h>

std::vector<int> nthCombination(int n, int r, double myIndex, bool isRep,
                                bool isMult, std::vector<int> Reps);

std::vector<int> nthPermutation(int n, int r, double myIndex, bool isRep, bool isMult,
                                std::vector<int> Reps, std::vector<int> freqs, bool isStarter = false);

std::vector<int> nthCombinationGmp(int n, int r, mpz_t myIndex, bool isRep,
                                   bool isMult, std::vector<int> Reps);

std::vector<int> nthPermutationGmp(int n, int r, mpz_t myIndex, bool isRep, bool isMult,
                                   std::vector<int> Reps, std::vector<int> freqs, bool isStarter = false);

#endif
