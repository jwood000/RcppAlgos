#include "Partitions/PartitionsCountSection.h"
#include <vector>
#include <cmath>

int width;
int blockSize;
static std::vector<double> memoize;

double pStdCap(int n, int m, int myMax) {
    
    if (myMax * m < n || n < m) return 0;
    if (myMax * m == n || n <= m + 1) return 1;
    if (m < 2) return m;
    
    const int block = myMax * blockSize + (n - m) * width + m - 2;
    if (memoize[block]) return memoize[block];
    
    int niter = n / m;
    
    if (m == 2) {
        if (myMax * 2 >= n) {
            myMax = std::min(myMax, n - 1);
            return niter - (n - 1 - myMax);
        } else {
            return 0;
        }
    }
    
    double count = 0;
    
    for (; niter--; n -= m, --myMax) {
        count += (memoize[myMax * blockSize + (n - m) * width + m - 3] = pStdCap(n - 1, m - 1, myMax));
    }
    
    return count;
}

double CountPartRepLenCap(int n, int m, int myMax) {
    
    if (myMax * m < n || n < m) return 0;
    if (myMax * m == n || n <= m + 1) return 1;
    if (m < 2) return m;
    
    if (m == 2) {
        if (myMax * 2 >= n) {
            myMax = std::min(myMax, n - 1);
            return n / m - (n - 1 - myMax);
        } else {
            return 0;
        }
    }
    
    width = m;
    blockSize = m * (n - m + 1);
    memoize = std::vector<double>((myMax + 1) * blockSize, 0.0);
    return pStdCap(n, m, myMax);
}

// This algorithm can be derived as follows:
double CountPartRepLen(int n, int m) {
    
    if (m == 0)
        return (n == 0) ? 1.0 : 0.0;
    
    if (n < m)
        return 0.0;
    
    if (n == m)
        return 1.0;
    
    if (m < 2)
        return 1.0;
    
    if (n - m == 1)
        return 1.0;
    
    // If n > 3, we have the following:
    // 1st part: n - 3 1's followed by a 3
    // 2nd part: n - 4 1's followed by two 2's
    //
    // N.B. We have taken care of every case where n <= 2 above
    // i.e. n = 2, m = 2; n = 2, m = 1; n = 1, m = 1; n = 0, m = 0
    if (n - m == 2)
        return 2.0;
    
    if (m == 2)
        return std::floor(n / 2);
    
    if (m == 3) {
        const double res = SumSection(n);
        return(res);
    } else {
        const int limit = std::min(n - m, m);
        n = (n < 2 * m) ? 2 * limit : n;
        
        std::vector<double> p1(n + 1);
        std::vector<double> p2(n + 1);
        
        for (int i = 3; i <= n; ++i) {
            p1[i] = SumSection(i);
        }
        
        for (int i = 4; i <= limit; ++i) {
            const int m2 = i * 2;
            
            if (i % 2) { 
                p1[i] = 1;
                
                for (int j = i + 1; j < m2; ++j) {
                    p1[j] = p2[j - 1];
                }
                
                for (int j = m2; j <= n; ++j) {
                    p1[j] = p2[j - 1] + p1[j - i];
                }
            } else {
                p2[i] = 1;
                
                for (int j = i + 1; j < m2; ++j) {
                    p2[j] = p1[j - 1];
                }
                
                for (int j = m2; j <= n; ++j) {
                    p2[j] = p1[j - 1] + p2[j - i];
                }
            }
        }
        
        return (limit % 2) ? p1.back() : p2.back();
    }
}

// Similar to CountPartDistinct
double CountPartRep(int n) {
    
    if (n < 2)
        return 1.0;
    
    const int n1 = n + 1;
    std::vector<double> qq(n1, 1);
    qq[0] = qq[1] = 1;
    
    for(int i = 2; i <= n; ++i) {
        qq[i] = 0;
        
        for (int s = 1, f = 1, r = 1; i >= r; f += 3, r += f, s *= -1) {
            qq[i] += s * qq[i - r];
        }
        
        for (int s = 1, f = 2, r = 2; i >= r; f += 3, r += f, s *= -1) {
            qq[i] += s * qq[i - r];
        }
    }
    
    return qq.back();
}
