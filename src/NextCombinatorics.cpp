#include "Combinations/NextComboSection.h"

using nextIterPtr = void (*const)(const std::vector<int> &freqs,
                          std::vector<int> &z, int n1, int m1);

void nextCombDistinct(const std::vector<int> &freqs,
                      std::vector<int> &z, int n1, int m1) {
    
    if (z[m1] != n1) {
        ++z[m1];
    } else {
        nextCombSec(z, m1, n1 - m1);
    }
}

void nextCombRep(const std::vector<int> &freqs,
                 std::vector<int> &z, int n1, int m1) {
    
    if (z[m1] != n1) {
        ++z[m1];
    } else {
        nextCombSecRep(z, m1, n1);
    }
}

void nextCombMulti(const std::vector<int> &freqs,
                   std::vector<int> &z, int n1, int m1) {
    
    if (z[m1] != n1) {
        ++z[m1];
    } else {
        std::vector<int> zIndex(n1 + 1);
        
        for (int i = 0; i <= n1; ++i) {
            zIndex[i] = std::find(freqs.cbegin(),
                                  freqs.cend(), i) - freqs.cbegin();
        }

        nextCombSecMulti(freqs, zIndex, z, m1, freqs.size() - m1 - 1);
    }
}

// These two algorithms are essentially the same as the two defined in
// NextPermutation.cpp, however they are defined here to have the same
// signature as the other functions in this file so that they can be
// used by nextIterPtr
void nextPermFull(const std::vector<int> &freqs,
                  std::vector<int> &z, int n1, int m1) {
    
    int p1 = n1 - 1;
    int p2 = n1;
    
    while (z[p1 + 1] <= z[p1])
        --p1;
    
    while (z[p2] <= z[p1])
        --p2;
    
    std::swap(z[p1], z[p2]);
    std::reverse(z.begin() + p1 + 1, z.end());
}

void nextPermPartial(const std::vector<int> &freqs,
                     std::vector<int> &z, int n1, int m1) {
    
    int p1 = m1 + 1;
    
    while (p1 <= n1 && z[m1] >= z[p1])
        ++p1;
    
    if (p1 <= n1) {
        std::swap(z[p1], z[m1]);
    } else {
        std::reverse(z.begin() + m1 + 1, z.end());
        p1 = m1;
        
        while (z[p1 + 1] <= z[p1]) {
            --p1;
        }

        int p2 = n1;
        
        while (z[p2] <= z[p1]) {
            --p2;
        }

        std::swap(z[p1], z[p2]);
        std::reverse(z.begin() + p1 + 1, z.end());
    }
}

void nextPermRep(const std::vector<int> &freqs,
                 std::vector<int> &z, int n1, int m1) {
    
    for (int i = m1; i >= 0; --i) {
        if (z[i] != n1) {
            ++z[i];
            break;
        } else {
            z[i] = 0;
        }
    }
}

nextIterPtr GetNextIterPtr(bool IsComb, bool IsMult,
                           bool IsRep, bool IsFull) {
    
    if (IsComb) {
        if (IsMult) {
            return(nextIterPtr(nextCombMulti));
        } else if (IsRep) {
            return(nextIterPtr(nextCombRep));
        } else {
            return(nextIterPtr(nextCombDistinct));
        }
    } else {
        if (IsRep) {
            return(nextIterPtr(nextPermRep));
        } else if (IsFull) {
            return(nextIterPtr(nextPermFull));
        } else {
            return(nextIterPtr(nextPermPartial));
        }
    }
}
