#include <Rcpp.h>

using prevIterPtr = void (*const)(const std::vector<int> &freqs,
                          std::vector<int> &z,int n1, int m1);

void prevCombCpp(const std::vector<int> &freqs, std::vector<int> &z, int n1, int m1) {
    
    for (int i = 0; i <= m1; ++i) {
        if (z[m1] - z[i] == m1 - i) {
            --z[i];
            
            for (int j = i + 1, k = n1 - m1 + i + 1; j <= m1 && z[j] != k; ++j, ++k) 
                z[j] = k;
            
            break;
        }
    }
}

void prevCombRepCpp(const std::vector<int> &freqs, std::vector<int> &z, int n1, int m1) {
    
    for (int i = 0; i <= m1; ++i) {
        if (z[i] == z[m1]) {
            --z[i];
            
            for (int j = i + 1; j <= m1; ++j)
                z[j] = n1;
            
            break;
        }
    }
}

void prevCombMultiCpp(const std::vector<int> &freqs, std::vector<int> &z, int n1, int m1) {
    
    const int pentExtreme = freqs.size() - m1 - 1;
    std::vector<int> myReps(n1 + 1);
    myReps.back() = freqs.size();
    
    for (int i = n1; i > 0; --i) {
        myReps[i - 1] = std::find(freqs.cbegin(), freqs.cend(), i) - freqs.cbegin();
        myReps[i] -= myReps[i - 1];
    };
    
    for (int i = m1; i >= 0; --i)
        --myReps[z[i]];
    
    bool keepGoing = true;
    
    for (int i = m1; i > 0; --i) {
        while (i && z[i] == z[i - 1])
            --i;
        
        if (myReps[z[i] - 1]) {
            --z[i];
            
            for (int j = i + 1, k = pentExtreme + i + 1; j <= m1; ++j, ++k)
                z[j] = freqs[k];
            
            keepGoing = false;
            break;
        }
    }
    
    if (keepGoing) {
        --z[0];
        
        for (int j = 1, k = pentExtreme + 1; j <= m1; ++j, ++k)
            z[j] = freqs[k];
    }
}

void prevFullPermCpp(const std::vector<int> &freqs, std::vector<int> &z, int n1, int m1) {
    
    int p1 = n1 - 1;
    int p2 = n1;
    
    while (z[p1 + 1] >= z[p1])
        --p1;
    
    while (z[p2] >= z[p1])
        --p2;
    
    std::swap(z[p1], z[p2]);
    std::reverse(z.begin() + p1 + 1, z.end());
}

void prevPartialPermCpp(const std::vector<int> &freqs, std::vector<int> &z, int n1, int m1) {
    
    int p1 = n1;
    
    while (p1 > m1 && z[p1] >= z[m1])
        --p1;
    
    if (p1 > m1) {
        std::swap(z[p1], z[m1]);
    } else {
        while (z[p1 + 1] >= z[p1])
            --p1;
        
        std::reverse(z.begin() + p1 + 1, z.end());
        
        int p2 = p1 + 1;
        
        while (z[p2] >= z[p1])
            ++p2;
        
        std::swap(z[p1], z[p2]);
        std::reverse(z.begin() + m1 + 1, z.end());
    }
}

void prevRepPermCpp(const std::vector<int> &freqs, std::vector<int> &z, int n1, int m1) {
    
    for (int i = m1; i >= 0; --i) {
        if (z[i] != 0) {
            --z[i];
            break;
        } else {
            z[i] = n1;
        }
    }
}

Rcpp::XPtr<prevIterPtr> putPrevIterPtrInXPtr(bool IsComb, bool IsMult, bool IsRep, bool IsFull) {
    
    if (IsComb) {
        if (IsMult) {
            return(Rcpp::XPtr<prevIterPtr>(new prevIterPtr(&prevCombMultiCpp)));
        } else if (IsRep) {
            return(Rcpp::XPtr<prevIterPtr>(new prevIterPtr(&prevCombRepCpp)));
        } else {
            return(Rcpp::XPtr<prevIterPtr>(new prevIterPtr(&prevCombCpp)));
        }
    } else {
        if (IsRep) {
            return(Rcpp::XPtr<prevIterPtr>(new prevIterPtr(&prevRepPermCpp)));
        } else if (IsFull) {
            return(Rcpp::XPtr<prevIterPtr>(new prevIterPtr(&prevFullPermCpp)));
        } else {
            return(Rcpp::XPtr<prevIterPtr>(new prevIterPtr(&prevPartialPermCpp)));
        }
    }
}
