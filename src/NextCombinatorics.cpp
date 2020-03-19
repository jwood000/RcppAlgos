#include "NextStandard.h"

using nextIterPtr = void (*const)(const std::vector<int> &freqs,
                          std::vector<int> &z,int n1, int m1);

void nextCombCpp(const std::vector<int> &freqs, std::vector<int> &z,int n1, int m1) {
    
    if (z[m1] != n1) {
        ++z[m1];
    } else {
        nextCombSec(z, m1, n1 - m1);
    }
}

void nextCombRepCpp(const std::vector<int> &freqs, std::vector<int> &z, int n1, int m1) {
    
    if (z[m1] != n1) {
        ++z[m1];
    } else {
        nextCombSecRep(z, m1, n1);
    }
}

void nextCombMultiCpp(const std::vector<int> &freqs, std::vector<int> &z, int n1, int m1) {
    
    if (z[m1] != n1) {
        ++z[m1];
    } else {
        std::vector<int> zIndex(n1 + 1);
        
        for (int i = 0; i <= n1; ++i)
            zIndex[i] = std::find(freqs.cbegin(), freqs.cend(), i) - freqs.cbegin();
        
        nextCombSecMulti(freqs, zIndex, z, m1, freqs.size() - m1 - 1);
    }
}

// These two algorithms are essentially the same as the two above, however they are
// defined to have the same signature as the other next*Cpp algorithtms so that
// they can be used for the class Combo
void nextFullPermCpp(const std::vector<int> &freqs, std::vector<int> &z, int n1, int m1) {
    
    int p1 = n1 - 1;
    int p2 = n1;
    
    while (z[p1 + 1] <= z[p1])
        --p1;
    
    while (z[p2] <= z[p1])
        --p2;
    
    std::swap(z[p1], z[p2]);
    std::reverse(z.begin() + p1 + 1, z.end());
}

void nextPartialPermCpp(const std::vector<int> &freqs, std::vector<int> &z, int n1, int m1) {
    
    int p1 = m1 + 1;
    
    while (p1 <= n1 && z[m1] >= z[p1])
        ++p1;
    
    if (p1 <= n1) {
        std::swap(z[p1], z[m1]);
    } else {
        std::reverse(z.begin() + m1 + 1, z.end());
        p1 = m1;
        
        while (z[p1 + 1] <= z[p1])
            --p1;
        
        int p2 = n1;
        
        while (z[p2] <= z[p1])
            --p2;
        
        std::swap(z[p1], z[p2]);
        std::reverse(z.begin() + p1 + 1, z.end());
    }
}

void nextRepPermCpp(const std::vector<int> &freqs, std::vector<int> &z, int n1, int m1) {
    
    for (int i = m1; i >= 0; --i) {
        if (z[i] != n1) {
            ++z[i];
            break;
        } else {
            z[i] = 0;
        }
    }
}

Rcpp::XPtr<nextIterPtr> putNextIterPtrInXPtr(bool IsComb, bool IsMult, bool IsRep, bool IsFull) {
    
    if (IsComb) {
        if (IsMult) {
            return(Rcpp::XPtr<nextIterPtr>(new nextIterPtr(&nextCombMultiCpp)));
        } else if (IsRep) {
            return(Rcpp::XPtr<nextIterPtr>(new nextIterPtr(&nextCombRepCpp)));
        } else {
            return(Rcpp::XPtr<nextIterPtr>(new nextIterPtr(&nextCombCpp)));
        }
    } else {
        if (IsRep) {
            return(Rcpp::XPtr<nextIterPtr>(new nextIterPtr(&nextRepPermCpp)));
        } else if (IsFull) {
            return(Rcpp::XPtr<nextIterPtr>(new nextIterPtr(&nextFullPermCpp)));
        } else {
            return(Rcpp::XPtr<nextIterPtr>(new nextIterPtr(&nextPartialPermCpp)));
        }
    }
}
