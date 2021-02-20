#include "StandardCount.h"
#include "PartitionEnums.h"

// Credit to Robin K. S. Hankin, author of the excellent partitions package.
// From the partitions.c, here are Hankin's comments for c_numbdiffparts:
//      "the recursion on p826 of Abramowitz and Stegun"
double c_numbdiffparts(int n) {
    
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

// Consider finding partitions of 50 of length 5 with no zeros.
// The first few results are:
//                  1 1 1 1 46
//                  1 1 1 2 45
//                  1 1 1 3 44
//
// If we cont. until the num in the 3rd col changes, we have:
//              [22,] 1 1 1 22 25
//              [23,] 1 1 1 23 24
//              [24,] 1 1 2  2 44
//
// That's 23 results, which is expected as we are decrementing
// the right most number and incrementing the number in the
// 4th column. Thus we have:
//                   (46 - n) > (1 + n)
//              ==>        45 > 2n
//              ==>      22.5 > n
//              ==>         n = 22
//
// Including the first result (i.e. 1 46), we have n + 1 = 23.
// In general, in can be show that when we are presented with
// a situaiton of the form: 1 1 ... 1 m, then the number of
// results for the "1st" group (i.e. while the 2nd to the last
// column stays the same), is given by floor((1 + m) / 2).
//
// Continuing on with the next group, we have:
//                   1 1 2 2 44
//
// We need to find the number of result of (2 44). This is
// equivalent to (1 43). The next few groups are (3 42), (4, 40),
// (5, 38), and (6, 36). These all get converted to (1, 40),
// (1, 37), (1, 34), and (1, 31). The number of results for
// each are given by (assume integer division and don't 
// forget to add 1):
//          44 / 2, 41 / 2, 38 / 2, 35 / 2, 32 / 2
//
// The total number of result for all groups is given by:
//          47 / 2 + 44 / 2 + ... + 5 / 2 + 2 / 2
//
// This can be converted to:
// {(47 - 0) + (47 - 3) + ... + (47 - 42) + (47 - 45)} / 2
//
// The 8 below comes from adding 1 to every other term
// that needs to be removed (leave as a proof to the reader)
//   ==>>  (47 * 16 - 3 * (0 + 1 + 2 + .... + 14 + 15) - 8) / 2
//   ==>>    (47 * 16 - 3 * (16 * 15 / 2) - 8) / 2

double SumSection(int n) {
    const int numIter = n / 3;
    const double sumOne = numIter * (n - 2);
    const double sumTwo = 3 * (numIter * (numIter - 1)) / 2;
    return std::floor((sumOne - sumTwo - numIter / 2) / 2);
}

double CountPartLen(int n, int m) {
    
    if (m < 2)
        return 1;
    
    if (m == 2)
        return std::floor((n - 1) / 2);
    
    const int limit = std::min(n - m, m);
    n = (n < 2 * m) ? 2 * limit : n;
    
    std::vector<double> p1(n + 1);
    std::vector<double> p2(n + 1);
    
    for (int i = 3; i <= n; ++i) {
        p1[i] = SumSection(i + 1);
    }
    
    for (int i = 4; i <= limit; ++i) {
        const int m2 = i * 2;
        int j = i + 1;
        
        if (i % 2) { 
            p1[i] = 1;
            
            for (; j < m2; ++j) {
                p1[j] = p2[j - 1];
            }
            
            for (; j <= n; ++j) {
                p1[j] = p2[j - 1] + p1[j - i];
            }
        } else {
            p2[i] = 1;
            
            for (; j < m2; ++j) {
                p2[j] = p1[j - 1];
            }
            
            for (; j <= n; ++j) {
                p2[j] = p1[j - 1] + p2[j - i];
            }
        }
    }
    
    return (m % 2) ? p1.back() : p2.back();
}

// *********************************************************************** //

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

double CountPartLenCap(int n, int m, int myMax) {
    
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

double pDist(int n, int m) {
    
    if (memoize[n * width + m]) 
        return memoize[n * width + m];
    
    if (m == 3)
        return SumSection(n - 2);
    
    int myLim = n / m;
    double count = 0;
    
    for (; --myLim; n -= m)
        count += (memoize[(n - m) * width + m - 1] = pDist(n - m, m - 1));
    
    return count;
}

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
    
    const int block = myMax * blockSize + (n - m) * width + m - 2;
    
    if (memoize[block])
        return memoize[block];
    
    int test = (m * (m + 1)) / 2;
    
    if (test > n) return 0;
    if (test == n) return 1;
    
    const int low = myMax - m;
    test = (myMax * (myMax + 1) - low * (low + 1)) / 2;
    
    if (test < n) return 0;
    if (test == n) return 1;
    
    double count = (memoize[block] = pDistCap(n - m, m - 1, myMax - 1) + pDistCap(n - m, m, myMax - 1));
    return count;
}

double CountDistPartLen(int n, int m) {
    
    if (m < 2)
        return 1;
    
    if (m == 2)
        return std::floor((n - 1) / 2);
    
    width = m + 1;
    memoize.assign((n + 1) * (m + 1), 0.0);
    return pDist(n, m);
}

double CountDistPartLenCap(int n, int m, int myMax) {
    
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
    
    width = m;
    blockSize = m * (n - m + 1);
    memoize = std::vector<double>((myMax + 1) * blockSize, 0.0);
    
    return pDistCap(n, m, myMax);
}

double GetComputedPartsComps(const std::vector<int> &z, PartitionType PartType, 
                             int target, int m, bool IsComb, bool IncludeZero, bool mIsNull) {
    
    double partCountTest = 0;
    const bool IsDistinct = PartType > PartitionType::TradNoZero;
    const int startLen = std::count_if(z.cbegin(), z.cend(), 
                                       [](int i){return i > 0;});
    
    if (IsDistinct) {
        if (IsComb) {
            if (PartType == PartitionType::DstctStdAll) {
                partCountTest = c_numbdiffparts(target);
            } else if (PartType == PartitionType::DistCapped) {
                partCountTest = CountDistPartLenCap(target, m, z.back());
            } else if (PartType == PartitionType::TradCapped) {
                partCountTest = CountPartLenCap(target, m, z.back());
            } else {
                for (int i = startLen; i <= m; ++i)
                    partCountTest += CountDistPartLen(target, i);
            }
        } else {
            if (IncludeZero) {
                // Given z = c(0, 0, 3, 4)
                // When m is Null, here is the output (only permute last
                // 2 indices):
                //                  c(0, 0, 3, 4)
                //                  c(0, 0, 4, 3)
                //
                // And when m is give, we have (permute all indices):
                //                   [,1] [,2] [,3] [,4]
                //              [1,]    0    0    3    4
                //              [2,]    0    0    4    3
                //              [3,]    0    3    0    4
                //                 .    .    .    .    .
                //             [10,]    4    0    0    3
                //             [11,]    4    0    3    0
                //             [12,]    4    3    0    0
                
                if (mIsNull) {
                    
                    for (int i = startLen; i <= m; ++i)
                        partCountTest += (CountDistPartLen(target, i) * NumPermsNoRep(i, i));
                    
                } else {
                    std::vector<int> permCountVec(m);
                    std::iota(permCountVec.begin(), permCountVec.begin() + startLen, 1);
                    
                    for (int i = startLen; i <= m; ++i) {
                        permCountVec[i - 1] = i;
                        partCountTest += (CountDistPartLen(target, i)
                                              * NumPermsWithRep(permCountVec));
                    }
                }
            } else {
                partCountTest = CountDistPartLen(target, m) * NumPermsNoRep(m, m);
            }
        }
    } else {
        if (IsComb) {
            for (int i = startLen; i <= m; ++i)
                partCountTest += CountPartLen(target, i);
            
        } else if (mIsNull && IncludeZero) {
            partCountTest = std::pow(2.0, static_cast<double>(target - 1));
        } else {
            partCountTest = (IncludeZero) ? nChooseK(target + m - 1, m - 1)
                                          : nChooseK(target - 1, m - 1);
        }
    }
    
    return partCountTest;
}
