#include "CombPermUtils.h"

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

// Credit to user @m69 on stackoverflow:
//       https://stackoverflow.com/a/32918426/4408538
// *********************************************************************** //
static Rcpp::NumericMatrix memoize;

double p(int n, int m) {
    
    if (n <= m + 1)
        return 1;
    
    if (memoize(n - m, m - 2))
        return memoize(n - m, m - 2);
    
    int myLim = n / m;
    
    if (m == 2)
        return myLim;
    
    double count = 0;
    
    for (; myLim--; n -= m)
        count += (memoize(n - m, m - 3) = p(n - 1, m - 1));
    
    return count;
}

double partitionCount(int n, int m, bool IncludeZero) {
    double count = 0;
    
    if (IncludeZero) {
        Rcpp::NumericMatrix bigRefill(n, m);
        memoize = bigRefill;
        
        for (int i = m; i > 1; --i) {
            for (int j = 0; j < n - i + 1; ++j) {
                for (int k = 0; k < i; ++k) {
                    memoize(j, k) = 0;
                }
            }
            
            count += p(n, i);
        }
        
        ++count;  // Add 1 for the case p(n, 1)
    } else {
        Rcpp::NumericMatrix refill(n - m + 1, m); // Initialize matrix to zero
        memoize = refill;
        count = p(n, m);
    }
    
    memoize = Rcpp::no_init_matrix(0, 0);
    return count;
}

// *********************************************************************** //

double SumSection(int high, int low) {
    
    const int numIter = (high - low + 2) / 3;
    const bool oddFlag = (high % 2);
    int n1 = oddFlag ? high - 1 : high;
    
    n1 >>= 1;
    double sumOne = (n1 * (n1 + 1)) / 2;
    const int n2 = n1 - numIter;
    
    sumOne -= (n2 * (n2 + 1)) / 2;
    sumOne *= 2;
    if (oddFlag) sumOne += numIter;
    
    const int n3 = low - 1;
    const int n4 = n3 + numIter;
    const double sumTwo = (n4 * (n4 + 1) - n3 * (n3 + 1)) / 2;
    
    return std::floor((1 + sumOne - sumTwo + numIter / 2) / 2);
}

double pDist(int n, int m) {
    
    if (memoize(n, m)) 
        return memoize(n, m);
    
    if (m == 3)
        return SumSection(n - 3, 2);
    
    int myLim = n / m;
    double count = 0;
    
    for (; --myLim; n -= m)
        count += (memoize(n - m, m - 1) = pDist(n - m, m - 1));
    
    return count;
}

double CountDistinctPartLen(int n, int m) {
    
    if (m < 2)
        return 1;
    
    if (m == 2)
        return std::floor((n - 1) / 2);
    
    Rcpp::NumericMatrix refill(n + 1, m + 1); // Initialize matrix to zero
    memoize = refill;
    double count = pDist(n, m);
    memoize = Rcpp::no_init_matrix(0, 0);
    
    return count;
}

double GetComputedPartsComps(const std::vector<int> &z, PartitionType PartType, 
                             int target, int m, bool IsComb, bool IncludeZero, bool mIsNull) {
    
    double partCountTest = 0;
    const bool IsDistinct = PartType > PartitionType::PartTradNoZero;
    
    if (IsDistinct) {
        const int startLen = std::count_if(z.cbegin(), z.cend(), 
                                           [](int i){return i > 0;});
        if (IsComb) {
            if (PartType == PartitionType::PartDstctStdAll) {
                partCountTest = c_numbdiffparts(target);
            } else {
                for (int i = startLen; i <= m; ++i)
                    partCountTest += CountDistinctPartLen(target, i);
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
                        partCountTest += (CountDistinctPartLen(target, i) * NumPermsNoRep(i, i));
                    
                } else {
                    std::vector<int> permCountVec(m);
                    std::iota(permCountVec.begin(), permCountVec.begin() + startLen, 1);
                    
                    for (int i = startLen; i <= m; ++i) {
                        permCountVec[i - 1] = i;
                        partCountTest += (CountDistinctPartLen(target, i)
                                              * NumPermsWithRep(permCountVec));
                    }
                }
            } else {
                partCountTest = CountDistinctPartLen(target, m) * NumPermsNoRep(m, m);
            }
        }
    } else {
        if (IsComb) {
            partCountTest = partitionCount(target, m, IncludeZero);
        } else if (mIsNull && IncludeZero) {
            partCountTest = std::pow(2.0, static_cast<double>(target - 1));
        } else {
            partCountTest = (IncludeZero) ? nChooseK(target + m - 1, m - 1)
                                          : nChooseK(target - 1, m - 1);
        }
    }
    
    return partCountTest;
}
