#include <Rcpp.h>

// Most of the code for rleCpp was obtained
// from Hadley Wickham's article titled 
// "High Performance functions with Rcpp"
// found: http://adv-r.had.co.nz/Rcpp.html
std::vector<std::vector<int> > rleCpp(std::vector<int> x) {
    std::vector<int> lengths, numUni, values;
    std::vector<int>::iterator it, xBeg, xEnd;
    xBeg = x.begin() + 1; xEnd = x.end();
    int prev = x[0];
    unsigned long int i = 0;
    
    values.push_back(prev);
    lengths.push_back(1);
    
    for(it = xBeg; it < xEnd; ++it) {
        if (prev == *it) {
            ++lengths[i];
        } else {
            values.push_back(*it);
            lengths.push_back(1);
            ++i;
            prev = *it;
        }
    }
    
    numUni.push_back((int) i);
    std::vector<std::vector<int> > myList {lengths, values, numUni};
    
    return myList;
}

double NumPermsWithRep(std::vector<int> v) {
    std::vector<std::vector<int> > myRle = rleCpp(v);
    int n = v.size(), myMax;
    std::vector<int> myLens = myRle[0], myUnis = myRle[2];
    std::sort(myLens.begin(), myLens.end(), std::greater<int>());
    
    myMax = myLens[0];
    int numUni = myUnis[0];
    double result = 1;
    
    for (int i = n; i > myMax; --i)
        result *= i;
    
    if (numUni > 0)
        for (int i = 1; i <= numUni; ++i)
            for (int j = 2; j <= myLens[i]; ++j)
                result /= j;
    
    return result;
}

double NumPermsNoRep(int n, int k) {
    double dblN = (double) n, result = 1;
    double i, m = dblN - (double) k;
    for (i = n; i > m; --i) {result *= i;}
    return result;
}

// Returns number of k-combinations from n elements.
// Mathematically speaking, we have: n!/(k!*(n-k)!)
double nChooseK(double n, double k) {
    
    if (k == n || k == 0)
        return 1;
    
    double nCk;
    double temp = 1;
    for(int i = 1; i <= k; ++i)
        temp *= (n - k + i)/i;
    
    nCk = round(temp);
    return nCk;
}

// For combinations where repetition is allowed, this
// function returns the number of combinations for
// a given n and r. The resulting vector, "triangleVec"
// resembles triangle numbers. In fact, this vector
// is obtained in a very similar method as generating
// triangle numbers, albeit in a repeating fashion.
double NumCombsWithRep(int n, int r) {

    if (r == 0)
        return 1;
    
    int i, k;
    std::vector<double> temp(n), triangleVec(n);
    std::iota(triangleVec.begin(), triangleVec.end(), 1.0);
    
    for (i = 1; i < r; ++i) {
        for (k = 1; k <= n; ++k)
            temp[k-1] = std::accumulate(triangleVec.begin(), triangleVec.begin() + k, 0.0);

        triangleVec = temp;
    }
    
    return triangleVec[n-1];
}

// Slightly different than CombsWithRep above as we can't
// guarantee 1) the repetition of each element is
// greater than or equal to n, and 2) that the
// repetition of the each element isn't the same
double MultisetCombRowNum(int n, int r, std::vector<int> Reps) {
    
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

// The algorithm below is credited to Randy Lai,
// author of arrangements and iterpc. It is much
// faster than the original naive approach whereby
// we create all combinations of the multiset, then
// subsequently count the number of permutations
// of each of those combinations.
double MultisetPermRowNum(int n, int r, std::vector<int> myReps) {
    
    if (n < 2 || r < 1)
        return 1.0;
    
    int sumFreqs = std::accumulate(myReps.begin(), myReps.end(), 0);
    
    if (r > sumFreqs)
        return 0.0;
    
    int maxFreq, n1 = n - 1;
    maxFreq = *std::max_element(myReps.begin(), myReps.end());
    
    std::vector<int> seqR(r);
    std::iota(seqR.begin(), seqR.end(), 1);
    
    double prodR, numPerms = 0.0;
    prodR = std::accumulate(seqR.begin(), seqR.end(), 
                            1.0, std::multiplies<double>());
    
    int myMax = (r < maxFreq) ? r : maxFreq;
    ++myMax;
    
    std::vector<double> cumProd(myMax), resV(r + 1, 0.0);
    
    // Create seqeunce from 1 to myMax, then add another
    // 1 at the front... equivalent to c(1, 1:myMax)
    std::iota(cumProd.begin(), cumProd.end(), 1);
    cumProd.insert(cumProd.begin(), 1);
    
    std::partial_sum(cumProd.begin(), cumProd.end(), 
                     cumProd.begin(), std::multiplies<double>());
    
    int myMin = std::min(r, myReps[0]);
    
    for (int i = 0; i <= myMin; ++i)
        resV[i] = prodR / cumProd[i];
    
    numPerms = resV[r];
    
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

// This algorithm is nearly identical to the
// one found in the standard algorithm library
void nextFullPerm(uint16_t *myArray, unsigned long int n1) {
    
    unsigned long int p1 = n1, p2 = n1;
    uint16_t temp;
    
    while (myArray[p1] <= myArray[p1 - 1])
        --p1;
    
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
void nextPartialPerm(uint16_t *myArray, unsigned long int nCols, 
                     unsigned long int r1, unsigned long int r,
                     unsigned long int n1, unsigned long int n) {
    
    uint16_t temp;
    unsigned long int p1 = nCols;
    
    while (p1 < n && myArray[r1] >= myArray[p1])
        ++p1;
    
    if (p1 < n) {
        temp = myArray[p1];
        myArray[p1] = myArray[r1];
        myArray[r1] = temp;
    } else {
        for (std::size_t k = r, q = n1; k < q; ++k, --q) {
            temp = myArray[k];
            myArray[k] = myArray[q];
            myArray[q] = temp;
        }
        
        p1 = n1;
        while (myArray[p1] <= myArray[p1 - 1])
            --p1;
        
        --p1;
        unsigned long int p2 = n1;
        
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
}
