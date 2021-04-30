#include "Partitions/PartitionsCountSection.h"
#include <vector>
#include <cmath>

int width_dist;
int bSize_dist;
static std::vector<double> mem_dist;

double pDistCap(int n, int m, int myMax) {
    
    if (m > n || myMax < m) return 0;
    
    if (m == n) {
        if (n == 1 && myMax >= 1) {
            return 1;
        } else {
            return 0;
        }
    }
    
    if (m == 1) {
        if (myMax >= n) {
            return 1;
        } else {
            return 0;
        }
    }
    
    const int block = myMax * bSize_dist + (n - m) * width_dist + m - 2;
    
    if (mem_dist[block])
        return mem_dist[block];
    
    int test = (m * (m + 1)) / 2;
    
    if (test > n) return 0;
    if (test == n) return 1;
    
    const int low = myMax - m;
    test = (myMax * (myMax + 1) - low * (low + 1)) / 2;
    
    if (test < n) return 0;
    if (test == n) return 1;
    
    double count = (mem_dist[block] = pDistCap(n - m, m - 1, myMax - 1) + pDistCap(n - m, m, myMax - 1));
    return count;
}

double CountPartDistinctLenCap(int n, int m, int myMax) {
    
    if (m > n || myMax < m) return 0;
    
    if (m == n) {
        if (n == 1 && myMax >= 1) {
            return 1;
        } else {
            return 0;
        }
    }
    
    if (m == 1) {
        if (myMax >= n) {
            return 1;
        } else {
            return 0;
        }
    }
    
    // Ensure max is large enough given the width
    //
    // Below, we have an expression that represents the
    // absolute maximum value we could obtain with a 
    // given max and width:
    //
    // max + (max - 1) + (max - 2) + ... + (max - (m - 1))
    //
    // (max * m) - (0 + 1 + 2 + ... + (m - 1))
    //
    // (max * m) - ((m - 1) * m) / 2
    
    const int limit = (myMax * m) - ((m - 1) * m) / 2;
    
    if (limit <= n) {
        if (limit == n) {
            return 1.0;
        } else {
            return 0.0;
        }
    }
    
    width_dist = m;
    bSize_dist = m * (n - m + 1);
    mem_dist = std::vector<double>((myMax + 1) * bSize_dist, 0.0);
    return pDistCap(n, m, myMax);
}

double pDist(int n, int m) {
    
    if (mem_dist[n * width_dist + m]) 
        return mem_dist[n * width_dist + m];
    
    if (m == 3)
        return SumSection(n - 3);
    
    int myLim = n / m;
    double count = 0;
    
    for (; --myLim; n -= m)
        count += (mem_dist[(n - m) * width_dist + m - 1] = pDist(n - m, m - 1));
    
    return count;
}

double CountPartDistinctLen(int n, int m) {
    
    if (m < 2)
        return 1;
    
    if (m == 2)
        return std::floor((n - 1) / 2);
    
    width_dist = m + 1;
    mem_dist.assign((n + 1) * (m + 1), 0.0);
    return pDist(n, m);
}

// Credit to Robin K. S. Hankin, author of the excellent partitions package.
// From the partitions.c, here are Hankin's comments for c_numbdiffparts:
//      "the recursion on p826 of Abramowitz and Stegun"
double CountPartDistinct(int n) {
    
    const int n1 = n + 1;
    std::vector<double> qq(n1, 1);
    
    for(int i = 2 ; i < n1; ++i) {
        qq[i] = 0;
        
        for (int s = 1, f = 5, r = 2; i >= r; r += f, f += 3, s *= -1) {
            qq[i] += s * qq[i - r];
            if(i == r * 2) {qq[i] -= s;}
        }
        
        for (int s = 1, f = 4, r = 1; i >= r; r += f, f += 3, s *= -1) {
            qq[i] += s * qq[i - r];
            if(i == r * 2) {qq[i] -= s;}
        }
    }
    
    return qq.back();
}
