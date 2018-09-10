#include <CombPermUtils.h>
#include <gmp.h>

/*
 * All functions below are exactly the same as the functions
 * in CombPermUtils.cpp. The only difference is that they
 * utilize the gmp library and deal mostly with mpz_t types
 */

void NumPermsWithRepGmp(mpz_t result, std::vector<int> &v) {
    mpz_set_ui(result, 1);
    std::vector<std::vector<int> > myRle = rleCpp(v);
    int n = v.size(), myMax;
    std::vector<int> myLens = myRle[0], myUnis = myRle[2];
    std::sort(myLens.begin(), myLens.end(), std::greater<int>());
    
    myMax = myLens[0];
    int numUni = myUnis[0];
    
    for (int i = n; i > myMax; --i)
        mpz_mul_ui(result, result, i);
    
    if (numUni > 0)
        for (int i = 1; i <= numUni; ++i)
            for (int j = 2; j <= myLens[i]; ++j)
                mpz_divexact_ui(result, result, j);
}

void NumPermsNoRepGmp(mpz_t result, int n, int k) {
    mpz_set_ui(result, 1);
    int i, m = n - k;
    for (i = n; i > m; --i)
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

void MultisetCombRowNumGmp(mpz_t result, int n, int r, std::vector<int> &Reps) {
    
    if (r >= 1 && n > 1) {
        int i, k, j, myMax, r1 = r + 1;
        mpz_t* triangleVec;
        mpz_t* temp;
        
        triangleVec = (mpz_t *) malloc(r1 * sizeof(mpz_t));
        temp = (mpz_t *) malloc(r1 * sizeof(mpz_t));
        
        for (i = 0; i < r1; ++i) {
            mpz_init(triangleVec[i]);
            mpz_init(temp[i]);
        }
    
        mpz_t tempSum;
        mpz_init(tempSum);
        
        myMax = r1;
        if (myMax > Reps[0] + 1)
            myMax = Reps[0] + 1;
        
        for (i = 0; i < myMax; ++i) {
            mpz_set_ui(triangleVec[i], 1);
            mpz_set_ui(temp[i], 1);
        }
        
        for (k = 1; k < n; ++k) {
            for (i = r; i > 0; --i) {
                myMax = i - Reps[k];
                if (myMax < 0)
                    myMax = 0;
                
                mpz_set_ui(tempSum, 0);
                
                for (j = myMax; j <= i; ++j)
                    mpz_add(tempSum, tempSum, triangleVec[j]);
                
                mpz_set(temp[i], tempSum);
            }
            
            for (int i = 0; i < r1; ++i)
                mpz_set(triangleVec[i], temp[i]);
        }
        
        mpz_set(result, triangleVec[r]);
        
        for (i = 0; i < r1; ++i) {
            mpz_clear(triangleVec[i]);
            mpz_clear(temp[i]);
        }
        
        free(triangleVec);
        free(temp);
        mpz_clear(tempSum);
    } else {
        mpz_set_ui(result, 1);
    }
}

void MultisetPermRowNumGmp(mpz_t result, int n, int r, std::vector<int> &myReps) {
    
    int sumFreqs = std::accumulate(myReps.begin(), myReps.end(), 0);
    
    if (n < 2 || r < 1) {
        mpz_set_ui(result, 1);
    } else if (r > sumFreqs) {
        mpz_set_ui(result, 0);
    } else {
        int maxFreq, n1 = n - 1;
        maxFreq = *std::max_element(myReps.begin(), myReps.end());
        
        std::vector<int> seqR(r);
        std::iota(seqR.begin(), seqR.end(), 1);
        
        mpz_t prodR;
        mpz_init(prodR);
        mpz_set_ui(prodR, 1);
        unsigned long int uR = r;
        
        for (std::size_t i = 0; i < uR; ++i)
            mpz_mul_ui(prodR, prodR, seqR[i]);
        
        unsigned long int uR1 = uR + 1;
        int myMax = (r < maxFreq) ? r : maxFreq;
        myMax += 2;
        
        mpz_t *cumProd, *resV;
        cumProd = (mpz_t *) malloc(sizeof(mpz_t) * myMax);
        resV = (mpz_t *) malloc(sizeof(mpz_t) * uR1);
        
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
        
        free(cumProd);
        free(resV);
    }
}
