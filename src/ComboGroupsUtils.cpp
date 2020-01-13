#include "GmpDependUtils.h"

// ******* Overview of the Crucial Part of the Algorithm *******
// -------------------------------------------------------------
// last1 is one plus the upper bound in the previous section, so to obtain the current
// current upper bound, we must first add the size of a section (i.e. grpSize) and su-
// bstract one. We can now compute the length we need to reset v by subtracting idx1. E.g.
//
// Given a portion of v w/ s1 = 9, gSize = 4, idx1 = 9, 6 groups (24 subjects) and base 0:
//       
//              prev sections   bound (index = 8)
//                  /  \        |
//            ............. 8 | 9 12 23 24 | 10 20 21 22 | 11 ... 
//                                 |
//                               idx1 (equal to last1, in this case)
//
// Sort v past idx1:
//                      ... 8 | 9 12 10 11 | 13 14 15 16 | 17... 
//
// Determine the index, idx3, such that v[idx3] > v[idx1]
//
//                      ... 8 | 9 12 10 11 | 13 14 15 16 | 17 ... 
//                                 |          |
//                               idx1       idx3
// 
// Swap idx1 and idx3:
//                      ... 8 | 9 13 10 11 | 12 14 15 16 | 17... 
//
// Move enough indices after idx1 to fill that specific group:
//
//                      ... 8 | 9 13 __ __ | 10 11 12 14 | 15 16 ... 
//
// Identify and move indices that are successively incrementing values of v past idx1:
//
//                      ... 8 | 9 13 14 15 | 10 11 12 16 | 17 ... 
//
// The last two steps are accomplished with std::rotate. This completes the algorithm.

bool nextComboGroup(std::vector<int> &z, int nGrps, 
                    int grpSize, int idx1, int idx2, int last1) {
    
    while (idx2 > idx1 && z[idx2] > z[idx1])
        --idx2;
    
    if ((idx2 + 1) < static_cast<int>(z.size())) {
        if (z[idx2 + 1] > z[idx1])
            std::swap(z[idx1], z[idx2 + 1]);
        
        return true;
    } else {
        const auto zbeg = z.begin();
        
        while (idx1 > 0) {
            const int tipPnt = z[idx2];
            
            while (idx1 > last1 && tipPnt < z[idx1])
                --idx1;
            
            if (tipPnt > z[idx1]) { // **Crucial Part**
                int idx3 = idx1 + 1;
                std::sort(zbeg + idx3, z.end());
                const int xtr = last1 + grpSize - idx3;
                
                while (z[idx3] < z[idx1])
                    ++idx3;
                
                std::swap(z[idx3], z[idx1]);
                std::rotate(zbeg + idx1 + 1, 
                            zbeg + idx3 + 1, zbeg + idx3 + xtr);
                return true;
            } else {
                idx1 -= 2;
                idx2 -= grpSize;
                last1 -= grpSize;
            }
        }
    }
    
    return false;
}

double numGroupCombs(int n, int numGroups, int grpSize) {
    
    double result = 1;
    
    for (double i = n; i > numGroups; --i)
        result *= i;
    
    if (result < std::numeric_limits<double>::max()) {
        double myDiv = 1;
        
        for (double i = 2; i <= grpSize; ++i)
            myDiv *= i;
        
        result /= std::pow(myDiv, numGroups);
        return std::round(result);
    } else {
        return std::numeric_limits<double>::infinity();
    }
}

void numGroupCombsGmp(mpz_t result, int n, 
                      int numGroups, int grpSize) {
    
    for (int i = n; i > numGroups; --i)
        mpz_mul_ui(result, result, i);
    
    mpz_t myDiv;
    mpz_init(myDiv);
    mpz_set_ui(myDiv, 1);
    
    for (int i = 2; i <= grpSize; ++i)
        mpz_mul_ui(myDiv, myDiv, i);
    
    mpz_pow_ui(myDiv, myDiv, numGroups);
    mpz_divexact(result, result, myDiv);
    mpz_clear(myDiv);
}

std::vector<int> nthComboGroup(int n, int gSize, int r,
                               double myIndex, double total) {
    
    double ind1 = myIndex, ind2 = myIndex;
    int s = n - 1;
    const int g = gSize - 1;
    int temp = static_cast<int>(nChooseK(s, g));
    int secLen = total / temp;
    
    std::vector<int> res(n, 0);
    std::vector<int> v(s);
    std::iota(v.begin(), v.end(), 1);
    
    int myMin = 0;
    mpz_t mpzDefault;
    mpz_init(mpzDefault);
    
    for (int j = 0; j < (r - 1); ++j) {
        ind2 = std::floor(ind2 / secLen);
        res[j * gSize] = myMin;
        const std::vector<int> comb = nthComb(s, g, ind2, mpzDefault, v);
        
        for (int k = j * gSize + 1, i = 0; k < ((j + 1) * gSize); ++k, ++i)
            res[k] = v[comb[i]];
        
        v.clear();
        
        for (int i = 1; i <= n; ++i) {
            const auto it = std::find(res.begin(), res.end(), i);
            
            if (it == res.end())
                v.push_back(i);
        }
        
        myMin = v.front();
        v.erase(v.begin());
        ind1 -= ind2 * secLen;
        ind2 = ind1;
        
        s -= gSize;
        temp = static_cast<int>(nChooseK(s, g));
        secLen /= temp;
    }
    
    res[(r - 1) * gSize] = myMin;
    
    for (int k = (r - 1) * gSize + 1, i = 0; k < (r * gSize); ++k, ++i)
        res[k] = v[i];
    
    return res;
}

std::vector<int> nthComboGroupGmp(int n, int gSize, int r,
                                  mpz_t lowerMpz, mpz_t computedRowMpz) {
    mpz_t ind1, ind2;
    mpz_init(ind1); mpz_init(ind2);
    mpz_set(ind1, lowerMpz); mpz_set(ind2, lowerMpz);
    
    int s = n - 1;
    const int g = gSize - 1;
    
    mpz_t temp, secLen;
    mpz_init(temp); mpz_init(secLen);
    
    nChooseKGmp(temp, s, g);
    mpz_divexact(secLen, computedRowMpz, temp);
    
    std::vector<int> res(n, 0);
    std::vector<int> v(s);
    std::iota(v.begin(), v.end(), 1);
    
    int myMin = 0;
    constexpr double dblDefault = 0;
    
    for (int j = 0; j < (r - 1); ++j) {
        mpz_tdiv_q(ind2, ind2, secLen);
        res[j * gSize] = myMin;
        const std::vector<int> comb = nthCombGmp(s, g, dblDefault, ind2, v);
        
        for (int k = j * gSize + 1, i = 0; k < ((j + 1) * gSize); ++k, ++i)
            res[k] = v[comb[i]];
        
        v.clear();
        
        for (int i = 1; i <= n; ++i) {
            const auto it = std::find(res.begin(), res.end(), i);
            
            if (it == res.end())
                v.push_back(i);
        }
        
        myMin = v.front();
        v.erase(v.begin());
        mpz_mul(temp, ind2, secLen);
        mpz_sub(ind1, ind1, temp);
        mpz_set(ind2, ind1);
        
        s -= gSize;
        nChooseKGmp(temp, s, g);
        mpz_divexact(secLen, secLen, temp);
    }
    
    res[(r - 1) * gSize] = myMin;
    
    for (int k = (r - 1) * gSize + 1, i = 0; k < (r * gSize); ++k, ++i)
        res[k] = v[i];
    
    return res;
}
