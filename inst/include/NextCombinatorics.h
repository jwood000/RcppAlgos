#ifndef NEXT_COMBINATORICS_H
#define NEXT_COMBINATORICS_H

#include <Rcpp.h>

inline void nextComb(std::vector<int> &z, int m1, int nMinusM) {
    
    for (int i = m1 - 1; i >= 0; --i) {
        if (z[i] != (nMinusM + i)) {
            ++z[i];
            
            for (int j = i; j < m1; ++j) 
                z[j + 1] = z[j] + 1;
            
            break;
        }
    }
}

inline void nextCombRep(std::vector<int> &z, int m1, int n1) {
    
    for (int i = m1 - 1; i >= 0; --i) {
        if (z[i] != n1) {
            ++z[i];
            
            for (int j = i; j < m1; ++j)
                z[j + 1] = z[j];
            
            break;
        }
    }
}

inline void nextCombMulti(const std::vector<int> &freqs, const std::vector<int> &zIndex,
                          std::vector<int> &zGroup, std::vector<int> &z, int m, int pentExtreme) {
    
    for (int i = m - 2; i >= 0; --i) {
        if (freqs[zIndex[z[i]]] != freqs[pentExtreme + i]) {
            ++z[i];
            zGroup[i] = zIndex[z[i]];
            
            for (int j = (i + 1); j < m; ++j) {
                zGroup[j] = zGroup[j - 1] + 1;
                z[j] = freqs[zGroup[j]];
            }
            
            break;
        }
    }
}

// This algorithm is nearly identical to the
// one found in the standard algorithm library
inline void nextFullPerm(int *const myArray, int maxInd) {
    
    int p1 = maxInd - 1;
    int p2 = maxInd;
    
    while (myArray[p1 + 1] <= myArray[p1])
        --p1;
    
    while (myArray[p2] <= myArray[p1])
        --p2;
    
    int temp = myArray[p1];
    myArray[p1] = myArray[p2];
    myArray[p2] = temp;
    
    for (int k = p1 + 1, q = maxInd; k < q; ++k, --q) {
        int temp = myArray[k];
        myArray[k] = myArray[q];
        myArray[q] = temp;
    }
}

// This algorithm is the same as above except that
// since we are not using the entire vector, we have
// to first check that the rth element is the largest.
// If it is, we have to reverse all of the elements
// to the right of the rth position before finding
// the next permutation. This is so because if we
// didn't, all of the next perms of the entire vector
// would produce many duplicate r-length perms. If it
// isn't the largest, we find the element to the right
// and swap them. We can then proceed to the next perm.
// We can do this because the standard algo would end
// up performing two unnecessary reversings.
inline void nextPartialPerm(int *const myArray, int r1, int maxInd) {
    
    int p1 = r1;
    
    while (p1 <= maxInd && myArray[r1] >= myArray[p1])
        ++p1;
    
    if (p1 <= maxInd) {
        int temp = myArray[p1];
        myArray[p1] = myArray[r1];
        myArray[r1] = temp;
    } else {
        for (int k = r1 + 1, q = maxInd; k < q; ++k, --q) {
            int temp = myArray[k];
            myArray[k] = myArray[q];
            myArray[q] = temp;
        }
        
        p1 = r1;
        while (myArray[p1 + 1] <= myArray[p1])
            --p1;
        
        int p2 = maxInd;
        
        while (myArray[p2] <= myArray[p1])
            --p2;
        
        int temp = myArray[p1];
        myArray[p1] = myArray[p2];
        myArray[p2] = temp;
        
        for (int k = p1 + 1, q = maxInd; k < q; ++k, --q) {
            int temp = myArray[k];
            myArray[k] = myArray[q];
            myArray[q] = temp;
        }
    }
}


#endif
