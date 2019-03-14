#include <Rcpp.h>

// Most of the code for rleCpp was obtained
// from Hadley Wickham's article titled 
// "High Performance functions with Rcpp"
// found: http://adv-r.had.co.nz/Rcpp.html
std::vector<int> rleCpp(const std::vector<int> &x) {
    std::vector<int> lengths;
    int prev = x[0];
    unsigned long int i = 0;
    lengths.push_back(1);
    
    for(auto it = x.cbegin() + 1; it < x.cend(); ++it) {
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

double NumPermsNoRep(const int n, const int k) {
    double dblN = static_cast<double>(n), result = 1;
    double m = dblN - static_cast<double>(k);
    for (double i = n; i > m; --i) {result *= i;}
    return result;
}

// Returns number of k-combinations from n elements.
// Mathematically speaking, we have: n!/(k!*(n-k)!)
double nChooseK(const int n, const int k) {
    
    if (k == n || k == 0)
        return 1.0;
    
    double nCk = 1;
    
    for (double i = (n - k + 1), d = 1; d <= k; ++i, ++d) {
        nCk *= i;
        nCk /= d;
    }
    
    return round(nCk);
}

double NumCombsWithRep(const int n, const int r) {
    return nChooseK(n + r - 1, r);
}

// The resulting vector, "triangleVec" resembles triangle
// numbers. In fact, this vector is obtained in a very
// similar method as generating triangle numbers, albeit
// in a repeating fashion. Two things to keep in mind is
// that we can't guarantee the following:
//      1) the repetition of each element is greater
//         than or equal to n
//      2) that the repetition of the each element 
//         isn't the same
double MultisetCombRowNumFast(const int n, const int r, 
                              const std::vector<int> &Reps) {
    
    if (r < 1 || n <= 1)
        return 1.0;
    
    if (r == n)
        if (std::accumulate(Reps.begin(), Reps.end(), 0) == n)
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

// The algorithm below is credited to Randy Lai,
// author of arrangements and iterpc. It is much
// faster than the original naive approach whereby
// we create all combinations of the multiset, then
// subsequently count the number of permutations
// of each of those combinations.
double MultisetPermRowNum(const int n, const int r, 
                          const std::vector<int> &myReps) {
    
    if (n < 2 || r < 1)
        return 1.0;
    
    int sumFreqs = std::accumulate(myReps.begin(), myReps.end(), 0);
    
    if (r > sumFreqs)
        return 0.0;
    
     const int n1 = n - 1;
     const int maxFreq = *std::max_element(myReps.begin(), myReps.end());
    
    std::vector<int> seqR(r);
    std::iota(seqR.begin(), seqR.end(), 1);
    const double prodR = std::accumulate(seqR.cbegin(), seqR.cend(), 
                                         1.0, std::multiplies<double>());
    
    const int myMax = (r < maxFreq) ? (r + 1) : (maxFreq + 1);
    std::vector<double> cumProd(myMax), resV(r + 1, 0.0);
    
    // Create seqeunce from 1 to myMax, then add another
    // 1 at the front... equivalent to c(1, 1:myMax)
    std::iota(cumProd.begin(), cumProd.end(), 1);
    cumProd.insert(cumProd.begin(), 1);
    
    std::partial_sum(cumProd.begin(), cumProd.end(), 
                     cumProd.begin(), std::multiplies<double>());
    
    double numPerms = 0.0;
    int myMin = std::min(r, myReps[0]);
    
    for (int i = 0; i <= myMin; ++i)
        resV[i] = prodR / cumProd[i];
    
    for (int i = 1; i < n1; ++i) {
        for (int j = r; j > 0; --j) {
            myMin = std::min(j, myReps[i]);
            numPerms = 0;
            for (int k = 0; k <= myMin; ++k)
                numPerms += resV[j - k] / cumProd[k];
            
            resV[j] = numPerms;
        }
    }
    
    myMin = std::min(r, myReps[n1]);
    numPerms = 0;
    for (int i = 0; i <= myMin; ++i)
        numPerms += resV[r - i] / cumProd[i];
    
    return numPerms;
}

// This function will be used in the main function to
// determine whether gmp analogs are needed as the fast
// algorithm above could potentionally produce negative
// results because of issues with double precision
double MultisetCombRowNum(const int n, const int r, 
                          const std::vector<int> &Reps) {
    
    if (r < 1 || n <= 1)
        return 1;
    
    int i, k, j, myMax, r1 = r + 1;
    std::vector<double> triangleVec(r1);
    std::vector<double> temp(r1);
    double tempSum;
    
    myMax = r1;
    if (myMax > Reps[0] + 1)
        myMax = Reps[0] + 1;
    
    for (i = 0; i < myMax; ++i)
        triangleVec[i] = 1;
    
    temp = triangleVec;
    
    for (k = 1; k < n; ++k) {
        for (i = r; i > 0; --i) {
            myMax = i - Reps[k];
            if (myMax < 0)
                myMax = 0;
            
            tempSum = 0;
            for (j = myMax; j <= i; ++j)
                tempSum += triangleVec[j];
            
            temp[i] = tempSum;
        }
        triangleVec = temp;
    }
    
    return triangleVec[r];
}

// This algorithm is nearly identical to the
// one found in the standard algorithm library
void nextFullPerm(int *myArray, const unsigned long int n1,
                  const unsigned long int n2) {
    
    unsigned long int p1 = n2, p2 = n1;
    int temp;
    
    while (myArray[p1 + 1] <= myArray[p1])
        --p1;
    
    while (myArray[p2] <= myArray[p1])
        --p2;
    
    temp = myArray[p1];
    myArray[p1] = myArray[p2];
    myArray[p2] = temp;
    
    for (std::size_t k = p1 + 1, q = n1; k < q; ++k, --q) {
        temp = myArray[k];
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
// didn't, all of the next perms. of the entire vector
// would produce many duplicate r-length perms. If it
// isn't the largest, we find the element to the right
// and swap them. We can then proceed to the next perm.
// We can do this because the standard algo would end
// up performing two unnecessary reversings.
void nextPartialPerm(int *myArray, const unsigned long int r, 
                     const unsigned long int r1, const unsigned long int n,
                     const unsigned long int lastElem) {
    
    int temp;
    unsigned long int p1 = r1;
    
    while (p1 < n && myArray[r1] >= myArray[p1])
        ++p1;
    
    if (p1 < n) {
        temp = myArray[p1];
        myArray[p1] = myArray[r1];
        myArray[r1] = temp;
    } else {
        for (std::size_t k = r, q = lastElem; k < q; ++k, --q) {
            temp = myArray[k];
            myArray[k] = myArray[q];
            myArray[q] = temp;
        }
        
        p1 = r1;
        while (myArray[p1 + 1] <= myArray[p1])
            --p1;
        
        unsigned long int p2 = lastElem;
        
        while (myArray[p2] <= myArray[p1])
            --p2;
        
        temp = myArray[p1];
        myArray[p1] = myArray[p2];
        myArray[p2] = temp;
        
        for (std::size_t k = p1 + 1, q = lastElem; k < q; ++k, --q) {
            temp = myArray[k];
            myArray[k] = myArray[q];
            myArray[q] = temp;
        }
    }
}
