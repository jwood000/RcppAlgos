#include "StandardCount.h"
#include "Cpp14MakeUnique.h"
#include <gmp.h>

// All functions below are exactly the same as the functions
// in StandardCount.cpp. The only difference is that they
// utilize the gmp library and deal mostly with mpz_t types

void NumPermsWithRepGmp(mpz_t result, const std::vector<int> &v) {
    
    mpz_set_ui(result, 1);
    std::vector<int> myLens = rleCpp(v);
    std::sort(myLens.begin(), myLens.end(), std::greater<int>());
    
    const int myMax = myLens[0];
    const int numUni = myLens.size();
    
    for (int i = v.size(); i > myMax; --i)
        mpz_mul_ui(result, result, i);
    
    if (numUni > 1)
        for (int i = 1; i < numUni; ++i)
            for (int j = 2; j <= myLens[i]; ++j)
                mpz_divexact_ui(result, result, j);
}

void NumPermsNoRepGmp(mpz_t result, int n, int k) {
    
    mpz_set_ui(result, 1);
    
    for (int i = n, m = n - k; i > m; --i)
        mpz_mul_ui(result, result, i);
}

void nChooseKGmp(mpz_t result, int n, int k) {
    
    mpz_set_ui(result, 1);
    
    if (k != n && k != 0) {
        for (int i = (n - k + 1), d = 1; d <= k; ++i, ++d) {
            mpz_mul_ui(result, result, i);
            mpz_divexact_ui(result, result, d);
        }
    }
}

void NumCombsWithRepGmp(mpz_t result, int n, int r) {
    nChooseKGmp(result, n + r - 1, r);
}

bool onlyOneCombo(int n, int r, const std::vector<int> &Reps) {
    if (r < 1 || n <= 1)
        return true;
    
    if (r == n)
        if (std::accumulate(Reps.cbegin(), Reps.cend(), 0) == n)
            return true;
    
    return false;
}

void MultisetCombRowNumGmp(mpz_t result, int n, int r, 
                           const std::vector<int> &Reps) {
    
    if (!onlyOneCombo(n, r, Reps)) {
        const int r1 = r + 1;
        int myMax = r1;
        
        if (myMax > Reps[0] + 1)
            myMax = Reps[0] + 1;

        auto triangleVec = FromCpp14::make_unique<mpz_t[]>(r1);
        auto temp = FromCpp14::make_unique<mpz_t[]>(r1);
        
        for (int i = 0; i < r1; ++i) {
            mpz_init(triangleVec[i]);
            mpz_init(temp[i]);
        }
        
        for (int i = 0; i < myMax; ++i) {
            mpz_set_ui(triangleVec[i], 1);
            mpz_set_ui(temp[i], 1);
        }
        
        --myMax;
        int ind = 1;
        
        for (; myMax < r; ++ind) {
            int myMin = std::min(Reps[ind], r);
            
            for (int i = 1; i <= myMin; ++i)
                mpz_add(triangleVec[i], triangleVec[i], triangleVec[i - 1]);
            
            myMin = std::min(Reps[ind] + myMax, r);
            int j = 0;
            
            for (int i = (Reps[ind] + 1); i <= myMin; ++i, ++j) {
                mpz_add(triangleVec[i], triangleVec[i], triangleVec[i - 1]);
                mpz_sub(triangleVec[i], triangleVec[i], temp[j]);
                mpz_set(temp[j], triangleVec[j]);
            }
            
            for (; j <= myMin; ++j)
                mpz_set(temp[j], triangleVec[j]);
            
            myMax = myMin;
        }
        
        const int n1 = n - 1;
        mpz_t mySum, t;
        mpz_init(mySum);
        mpz_init(t);
        
        for (; ind < n1; ++ind) {
            mpz_set(t, triangleVec[r]);
            const int s = std::min(Reps[ind], r);
            
            for (int i = 1; i <= s; ++i)
                mpz_add(triangleVec[r], triangleVec[r], triangleVec[r - i]);
            
            mpz_set(mySum, triangleVec[r]);
            
            for (int i = r - 1; i >= s; --i) {
                mpz_sub(mySum, mySum, t);
                mpz_set(t, triangleVec[i]);
                mpz_add(mySum, mySum, triangleVec[i - s]);
                mpz_set(triangleVec[i], mySum);
            }
            
            for (int i = s - 1; i > 0; --i) {
                mpz_sub(mySum, mySum, t);
                mpz_set(t, triangleVec[i]);
                mpz_set(triangleVec[i], mySum);
            }
        }
        
        if (ind < n) {
            const int myMin2 = std::min(Reps[n1], r);
            
            for (int i = 1; i <= myMin2; ++i)
                mpz_add(triangleVec[r], triangleVec[r], triangleVec[r - i]);
        }
        
        mpz_set(result, triangleVec[r]);
        
        for (int i = 0; i < r1; ++i) {
            mpz_clear(triangleVec[i]);
            mpz_clear(temp[i]);
        }

        mpz_clear(mySum); mpz_clear(t);
    } else {
        mpz_set_ui(result, 1);
    }
}

void MultisetPermRowNumGmp(mpz_t result, int n, int r, const std::vector<int> &myReps) {
    
    const int sumFreqs = std::accumulate(myReps.cbegin(), myReps.cend(), 0);
    
    if (n < 2 || r < 1) {
        mpz_set_ui(result, 1);
    } else if (r > sumFreqs) {
        mpz_set_ui(result, 0);
    } else {
        const int n1 = n - 1;
        int maxFreq = *std::max_element(myReps.cbegin(), myReps.cend());
        
        std::vector<int> seqR(r);
        std::iota(seqR.begin(), seqR.end(), 1);
        
        mpz_t prodR;
        mpz_init(prodR);
        mpz_set_ui(prodR, 1);
        
        for (int i = 0; i < r; ++i)
            mpz_mul_ui(prodR, prodR, seqR[i]);
        
        const std::size_t uR1 = r + 1;
        const int myMax = (r < maxFreq) ? (r + 2) : (maxFreq + 2);
        auto cumProd = FromCpp14::make_unique<mpz_t[]>(myMax);
        auto resV = FromCpp14::make_unique<mpz_t[]>(uR1);
        
        for (int i = 0; i < myMax; ++i)
            mpz_init(cumProd[i]);
        
        // Equivalent to c(1, 1:myMax)
        mpz_set_ui(cumProd[0], 1);
        for (int i = 1; i < myMax; ++i)
            mpz_set_ui(cumProd[i], i);
        
        for (std::size_t i = 0; i < uR1; ++i) {
            mpz_init(resV[i]);
            mpz_set_ui(resV[i], 0);
        }
        
        for (int i = 1; i < myMax; ++i)
            mpz_mul(cumProd[i], cumProd[i], cumProd[i - 1]);
        
        int myMin = std::min(r, myReps[0]);
        
        for (int i = 0; i <= myMin; ++i)
            mpz_divexact(resV[i], prodR, cumProd[i]);
        
        mpz_t temp;
        mpz_init(temp);
        
        for (int i = 1; i < n1; ++i) {
            for (int j = r; j > 0; --j) {
                myMin = std::min(j, myReps[i]);
                mpz_set_ui(result, 0);
                
                for (int k = 0; k <= myMin; ++k) {
                    mpz_divexact(temp, resV[j - k], cumProd[k]);
                    mpz_add(result, result, temp);
                }
                
                mpz_set(resV[j], result);
            }
        }
        
        myMin = std::min(r, myReps[n1]);
        mpz_set_ui(result, 0);
        
        for (int k = 0; k <= myMin; ++k) {
            mpz_divexact(temp, resV[r - k], cumProd[k]);
            mpz_add(result, result, temp);
        }
        
        mpz_clear(temp); mpz_clear(prodR);
        
        for (int i = 0; i < myMax; ++i)
            mpz_clear(cumProd[i]);
        
        for (std::size_t i = 0; i < uR1; ++i)
            mpz_clear(resV[i]);
    }
}

void GetComputedRowMpz(mpz_t computedRowsMpz, bool IsMultiset, bool IsComb, bool IsRep,
                       int n, int m, const SEXP &Rm, const std::vector<int> &freqs, 
                       const std::vector<int> &myReps) {
    
    mpz_init(computedRowsMpz);
    
    if (IsMultiset) {
        if (IsComb) {
            MultisetCombRowNumGmp(computedRowsMpz, n, m, myReps);
        } else {
            if (Rf_isNull(Rm) || m == static_cast<int>(freqs.size()))
                NumPermsWithRepGmp(computedRowsMpz, freqs);
            else
                MultisetPermRowNumGmp(computedRowsMpz, n, m, myReps);
        }
    } else {
        if (IsRep) {
            if (IsComb)
                NumCombsWithRepGmp(computedRowsMpz, n, m);
            else
                mpz_ui_pow_ui(computedRowsMpz, n, m);
        } else {
            if (IsComb)
                nChooseKGmp(computedRowsMpz, n, m);
            else
                NumPermsNoRepGmp(computedRowsMpz, n, m);
        }
    }
}
