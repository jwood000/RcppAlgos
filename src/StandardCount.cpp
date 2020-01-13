#include "GmpDependUtils.h"

// Most of the code for rleCpp was obtained from Hadley Wickham's
// article titled "High Performance functions with Rcpp" found:
//             http://adv-r.had.co.nz/Rcpp.html
std::vector<int> rleCpp(const std::vector<int> &x) {
    std::vector<int> lengths;
    int prev = x[0];
    std::size_t i = 0;
    lengths.push_back(1);
    
    for(auto it = x.cbegin() + 1; it != x.cend(); ++it) {
        if (prev == *it) {
            ++lengths[i];
        } else {
            lengths.push_back(1);
            prev = *it;
            ++i;
        }
    }
    
    return lengths;
}

double NumPermsWithRep(const std::vector<int> &v) {
    std::vector<int> myLens = rleCpp(v);
    std::sort(myLens.begin(), myLens.end(), std::greater<int>());
    
    const int myMax = myLens[0];
    const int numUni = myLens.size();
    double result = 1;
    
    for (int i = v.size(); i > myMax; --i)
        result *= i;
    
    if (numUni > 1)
        for (int i = 1; i < numUni; ++i)
            for (int j = 2; j <= myLens[i]; ++j)
                result /= j;
    
    return result;
}

double NumPermsNoRep(int n, int k) {
    const double dblN = static_cast<double>(n);
    const double nMinusK = dblN - static_cast<double>(k);
    double result = 1.0;
    
    for (double i = n; i > nMinusK; --i)
        result *= i;
    
    return result;
}

// Returns number of k-combinations from n elements.
// Mathematically speaking, we have: n!/(k!*(n-k)!)
double nChooseK(int n, int k) {
    
    if (k == n || k == 0)
        return 1.0;
    
    double nCk = 1;
    
    for (double i = (n - k + 1), d = 1; d <= k; ++i, ++d) {
        nCk *= i;
        nCk /= d;
    }
    
    return std::round(nCk);
}

double NumCombsWithRep(int n, int r) {
    return nChooseK(n + r - 1, r);
}

// The resulting vector, "triangleVec" resembles triangle numbers. In
// fact, this vector is obtained in a very similar method as generating
// triangle numbers, albeit in a repeating fashion. Two things to keep
// in mind is that we can't guarantee the following:
//      1) the repetition of each element is greater than or equal to n
//      2) that the repetition of each element isn't the same
double MultisetCombRowNumFast(int n, int r, const std::vector<int> &Reps) {
    
    if (r < 1 || n <= 1)
        return 1.0;
    
    if (r == n)
        if (std::accumulate(Reps.cbegin(), Reps.cend(), 0) == n)
            return 1.0;
        
    const int r1 = r + 1;
    std::vector<double> triangleVec(r1);
    std::vector<double> temp(r1);
    
    int myMax = r1;
    
    if (myMax > Reps[0] + 1)
        myMax = Reps[0] + 1;
    
    for (int i = 0; i < myMax; ++i)
        triangleVec[i] = temp[i] = 1;
    
    --myMax;
    int ind = 1;
    
    for (; myMax < r; ++ind) {
        int myMin = std::min(Reps[ind], r);
        
        for (int i = 1; i <= myMin; ++i)
            triangleVec[i] += triangleVec[i - 1];
        
        myMin = std::min(Reps[ind] + myMax, r);
        int j = 0;
        
        for (int i = (Reps[ind] + 1); i <= myMin; ++i, ++j) {
            triangleVec[i] += triangleVec[i - 1];
            triangleVec[i] -= temp[j];
            temp[j] = triangleVec[j];
        }
        
        for (; j <= myMin; ++j)
            temp[j] = triangleVec[j];
        
        myMax = myMin;
    }
    
    const int n1 = n - 1;
    
    for (; ind < n1; ++ind) {
        double t = triangleVec[r];
        const int s = std::min(Reps[ind], r);
        
        for (int i = 1; i <= s; ++i)
            triangleVec[r] += triangleVec[r - i];
        
        double mySum = triangleVec[r];
        
        for (int i = r - 1; i >= s; --i) {
            mySum -= t;
            t = triangleVec[i];
            mySum += triangleVec[i - s];
            triangleVec[i] = mySum;
        }
        
        for (int i = s - 1; i > 0; --i) {
            mySum -= t;
            t = triangleVec[i];
            triangleVec[i] = mySum;
        }
    }
    
    if (ind < n) {
        const int myMin2 = std::min(Reps[n1], r);
        
        for (int i = 1; i <= myMin2; ++i)
            triangleVec[r] += triangleVec[r - i];
    }
    
    return triangleVec[r];
}

// The algorithm below is credited to Randy Lai, author of arrangements
// and iterpc. It is much faster than the original naive approach whereby
// we create all combinations of the multiset, thensubsequently count the
// number of permutations of each of those combinations.
double MultisetPermRowNum(int n, int r, const std::vector<int> &Reps) {
    
    if (n < 2 || r < 1)
        return 1.0;
    
    int sumFreqs = std::accumulate(Reps.cbegin(), Reps.cend(), 0);
    
    if (r > sumFreqs)
        return 0.0;
    
    const int n1 = n - 1;
    const int maxFreq = *std::max_element(Reps.cbegin(), Reps.cend());
    const int myMax = (r < maxFreq) ? (r + 1) : (maxFreq + 1);
    
    // factorial(171)
    // [1] Inf 
    // factorial(170)
    // [1] 7.257416e+306
    if (myMax > 170 || r > 170) {
        mpz_t result;
        mpz_init(result);
        MultisetPermRowNumGmp(result, n, r, Reps);
        
        if (mpz_cmp_d(result, Significand53) > 0)
            return std::numeric_limits<double>::infinity();
        else
            return mpz_get_d(result);
    }
    
    std::vector<int> seqR(r);
    std::iota(seqR.begin(), seqR.end(), 1);
    const double prodR = std::accumulate(seqR.cbegin(), seqR.cend(), 
                                         1.0, std::multiplies<double>());
    
    // Create seqeunce from 1 to myMax, then add another
    // 1 at the front... equivalent to c(1, 1:myMax)
    std::vector<double> cumProd(myMax), resV(r + 1, 0.0);
    std::iota(cumProd.begin(), cumProd.end(), 1);
    
    cumProd.insert(cumProd.begin(), 1);
    std::partial_sum(cumProd.begin(), cumProd.end(), 
                     cumProd.begin(), std::multiplies<double>());
    
    double numPerms = 0.0;
    int myMin = std::min(r, Reps[0]);
    
    for (int i = 0; i <= myMin; ++i)
        resV[i] = prodR / cumProd[i];
    
    for (int i = 1; i < n1; ++i) {
        for (int j = r; j > 0; --j) {
            myMin = std::min(j, Reps[i]);
            numPerms = 0;
            
            for (int k = 0; k <= myMin; ++k)
                numPerms += resV[j - k] / cumProd[k];
            
            resV[j] = numPerms;
        }
    }
    
    myMin = std::min(r, Reps[n1]);
    numPerms = 0;
    
    for (int i = 0; i <= myMin; ++i)
        numPerms += resV[r - i] / cumProd[i];
    
    return numPerms;
}

// This function will be used in the main function to determine whether
// gmp analogs are needed as the fast algorithm above could potentially
// produce negative results because of issues with double precision
double MultisetCombRowNum(int n, int r, const std::vector<int> &Reps) {
    
    if (r < 1 || n <= 1)
        return 1;
    
    int i, k, j, myMax, r1 = r + 1;
    std::vector<double> triangleVec(r1);
    std::vector<double> temp(r1);
    
    myMax = r1;
    
    if (myMax > Reps[0] + 1)
        myMax = Reps[0] + 1;
    
    for (i = 0; i < myMax; ++i)
        triangleVec[i] = 1;
    
    temp = triangleVec;
    
    for (k = 1; k < n; ++k) {
        for (i = r; i > 0; --i) {
            myMax = i - Reps[k];
            
            if (myMax < 0) {myMax = 0;}
            
            double tempSum = 0;
            
            for (j = myMax; j <= i; ++j)
                tempSum += triangleVec[j];
            
            temp[i] = tempSum;
        }
        
        triangleVec = temp;
    }
    
    return triangleVec[r];
}

double GetComputedRows(bool IsMultiset, bool IsComb, bool IsRep, int n, int &m, SEXP Rm,
                       int lenFreqs, std::vector<int> &freqs, std::vector<int> &Reps) {
    
    double computedRows = 0;
    
    if (IsMultiset) {
        if (n != lenFreqs)
            Rcpp::stop("the length of freqs must equal the length of v");
        
        if (m > static_cast<int>(freqs.size()))
            m = freqs.size();
        
        if (IsComb) {
            computedRows = MultisetCombRowNum(n, m, Reps);
        } else {
            if (Rf_isNull(Rm) || m == static_cast<int>(freqs.size()))
                computedRows = NumPermsWithRep(freqs);
            else
                computedRows = MultisetPermRowNum(n, m, Reps);
        }
    } else {
        if (IsRep) {
            if (IsComb)
                computedRows = NumCombsWithRep(n, m);
            else
                computedRows = std::pow(static_cast<double>(n), static_cast<double>(m));
        } else {
            if (m > n)
                Rcpp::stop("m must be less than or equal to the length of v");
            
            if (IsComb)
                computedRows = nChooseK(n, m);
            else
                computedRows = NumPermsNoRep(n, m);
        }
    }
    
    return computedRows;
}
