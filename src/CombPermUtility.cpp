#include <Rcpp.h>
#include "CombPermUtility.h"

// Most of the code for rleCpp was obtained
// from Hadley Wickham's article titled 
// "High Performance functions with Rcpp"
// found: http://adv-r.had.co.nz/Rcpp.html
Rcpp::List rleCpp(std::vector<int> x) {
    std::vector<unsigned long int> lengths, numUni;
    std::vector<int> values;
    std::vector<int>::iterator it, xBeg, xEnd;
    xBeg = x.begin() + 1; xEnd = x.end();
    int prev = x[0];
    unsigned long int n = x.size(), i = 0;
    lengths.reserve(n);
    values.reserve(n);
    values.push_back(prev);
    lengths.push_back(1);
    for(it = xBeg; it < xEnd; it++) {
        if (prev == *it) {
            lengths[i]++;
        } else {
            values.push_back(*it);
            lengths.push_back(1);
            i++;
            prev = *it;
        }
    }
    
    numUni.push_back(i);
    return Rcpp::List::create(lengths,values,numUni);
}

double NumPermsWithRep(std::vector<int> v) {
    Rcpp::List myRle = rleCpp(v);
    unsigned long int n = v.size(), myMax;
    std::vector<unsigned long int> myLens = myRle[0], myUnis = myRle[2];
    std::sort(myLens.begin(), myLens.end(),
              std::greater<unsigned long int>());
    
    myMax = myLens[0];
    unsigned long int i, j, numUni = myUnis[0];
    double result = 1;
    
    for (i = n; i > myMax; i--) {result *= i;}
    
    if (numUni > 0) {
        for (i = 1; i <= numUni; i++) {
            // No need to divide by 1.
            // Start j at 2 instead.
            for (j = 2; j <= myLens[i]; j++) {
                result /= j;
            }
        }
    }
    
    return result;
}

double NumPermsNoRep(int n, int k) {
    double dblN = (double)n, result = 1;
    double i, m = dblN - (double)k;
    for (i = n; i > m; i--) {result *= i;}
    return result;
}

double nChooseK(double n, double k) {
    // returns the number of k-combinations from a set
    // of n elements. Mathematically speaking, 
    //  we have: n!/(k!*(n-k)!)
    double nCk;
    double temp = 1;
    for(int i = 1; i <= k; i++) {temp *= (n - k + i)/i;}
    nCk = round(temp);
    return nCk;
}

double GetRowNum(int n, int r) {
    // for combinations where repetition is allowed, this
    // function returns the number of combinations for
    // a given n and r. The resulting vector, "triangleVec"
    // resembles triangle numbers. In fact, this vector
    // is obtained in a very similar method as generating
    // triangle numbers, albeit in a repeating fashion.
    int i, k;
    std::vector<double> triangleVec(n);
    std::vector<double> temp(n);
    for (i = 0; i < n; i++) {triangleVec[i] = i+1;}
    
    for (i = 1; i < r; i++) {
        for (k = 1; k <= n; k++) {
            temp[k-1] = std::accumulate(triangleVec.begin(), triangleVec.begin() + k, 0.0);
        }
        triangleVec = temp;
    }
    
    return triangleVec[n-1];
}

// Slightly different than GetRowNum above as we can't
// guarantee 1) the repetition of each element is
// greater than or equal to n, and 2) that the
// repetition of the each element isn't the same
double MultisetCombRowNum(int n, int r, std::vector<int> Reps) {
    int i, k, j, myMax, r1 = r+1;
    std::vector<double> triangleVec(r1);
    std::vector<double> temp(r1);
    double tempSum;
    
    myMax = r1;
    if (myMax > Reps[0] + 1) {myMax = Reps[0] + 1;}
    
    for (i = 0; i < myMax; i++) {triangleVec[i] = 1;}
    temp = triangleVec;
    
    for (k = 1; k < n; k++) {
        for (i = r; i > 0; i--) {
            myMax = i - Reps[k];
            if (myMax < 0) {myMax = 0;}
            
            tempSum = 0;
            for (j = myMax; j <= i; j++) {tempSum += triangleVec[j];}
            temp[i] = tempSum;
        }
        triangleVec = temp;
    }
    
    return triangleVec[r];
}


Rcpp::IntegerMatrix MakeIndexHeaps(unsigned long int indRows, unsigned long int r) {
    unsigned long int j, i = 0, count = 0;
    Rcpp::IntegerMatrix indexMatrix(indRows, r);
    
    std::vector<unsigned long int> vecInd(r, 0);
    Rcpp::IntegerVector mySeq = Rcpp::seq(0, r-1);
    for (j = 0; j < r; j++) {indexMatrix(count, j) = j;}
    
    while (i < r) {
        if (vecInd[i] < i) {
            if (i % 2 == 0)
                std::swap(mySeq[0], mySeq[i]);
            else
                std::swap(mySeq[vecInd[i]], mySeq[i]);
            count++;
            for (j = 0; j < r; j++) {indexMatrix(count,j) = mySeq[j];}
            vecInd[i]++;
            i = 0;
        } else {
            vecInd[i] = 0;
            i++;
        }
    }
    
    return indexMatrix;
}

// Below, we define five functions that will be utilized
// as constraint functions. We also define five comparison
// operations (<, <=, >, >=, ==). The framework is based on
// the information posted by Dirk Eddelbuettel from
// this Rcpp Gallery (http://gallery.rcpp.org/articles/passing-cpp-function-pointers/)

double prodCpp(std::vector<double>& v) {
    std::vector<double>::iterator it, vEnd = v.end();
    double myProduct = 1.0;
    for (it = v.begin(); it < vEnd; it++) {myProduct *= *it;}
    return(myProduct);
}

double sumCpp(std::vector<double>& v) {
    std::vector<double>::iterator it, vEnd = v.end();
    double mySum = 0.0;
    for (it = v.begin(); it < vEnd; it++) {mySum += *it;}
    return(mySum);
}

double meanCpp(std::vector<double>& v){
    double s = v.size();
    double mySum = sumCpp(v);
    return (mySum/s);
}

double maxCpp(std::vector<double>& v) {
    std::vector<double>::iterator y = std::max_element(v.begin(), v.end());
    return v[std::distance(v.begin(), y)];
}

double minCpp(std::vector<double>& v) {
    std::vector<double>::iterator y = std::min_element(v.begin(), v.end());
    return v[std::distance(v.begin(), y)];
}
