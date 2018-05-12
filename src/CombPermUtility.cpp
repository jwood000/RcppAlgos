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
    
    for (i = n; i > myMax; i--)
        result *= i;
    
    if (numUni > 0)
        for (i = 1; i <= numUni; i++)
            for (j = 2; j <= myLens[i]; j++)
                result /= j;
    
    return result;
}

double NumPermsNoRep(int n, int k) {
    double dblN = (double)n, result = 1;
    double i, m = dblN - (double)k;
    for (i = n; i > m; i--) {result *= i;}
    return result;
}

double nChooseK(double n, double k) {
    // Returns number of k-combinations from n elements.
    // Mathematically speaking, we have: n!/(k!*(n-k)!)
    
    if (k == n || k == 0)
        return 1;
    
    double nCk;
    double temp = 1;
    for(int i = 1; i <= k; i++)
        temp *= (n - k + i)/i;
    
    nCk = round(temp);
    return nCk;
}

double NumCombsWithRep(int n, int r) {
    // For combinations where repetition is allowed, this
    // function returns the number of combinations for
    // a given n and r. The resulting vector, "triangleVec"
    // resembles triangle numbers. In fact, this vector
    // is obtained in a very similar method as generating
    // triangle numbers, albeit in a repeating fashion.
    
    if (r == 0)
        return 1;
    
    int i, k;
    std::vector<double> triangleVec(n);
    std::vector<double> temp(n);
    
    for (i = 0; i < n; i++)
        triangleVec[i] = i+1;
    
    for (i = 1; i < r; i++) {
        for (k = 1; k <= n; k++)
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
            if (myMax < 0)
                myMax = 0;
            
            tempSum = 0;
            for (j = myMax; j <= i; j++)
                tempSum += triangleVec[j];
            
            temp[i] = tempSum;
        }
        triangleVec = temp;
    }
    
    return triangleVec[r];
}

double MultisetPermRowNum(int n, int r, 
                          std::vector<int> Reps, 
                          Rcpp::IntegerMatrix myCombs) {
    
    int combRows = (int) MultisetCombRowNum(n, r, Reps);
    std::vector<int> rowVec(r);
    Rcpp::IntegerVector temp(1);
    bool KeepGoing = true;
    double numRows = 0, computedRows;
    int rowCount;
    
    for (std::size_t i = 0; i < combRows; i++) {
        int j = 0, k = 0;
        while (j < r) {
            rowVec[j] = k;
            temp[0] = myCombs(i, j);
            KeepGoing = true;
            while (KeepGoing) {
                rowVec[j] = k;
                j++;
                if (j >= r)
                    KeepGoing = false;
                else if (myCombs(i, j) != temp[0])
                    KeepGoing = false;
            }
            k++;
        }

        computedRows = NumPermsWithRep(rowVec);
        rowCount = (int) computedRows;
        numRows += rowCount;
    }
    
    return numRows;
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

// Comparison functions
bool lessCpp(double& x, double& y) {return x < y;}
bool greaterCpp(double& x, double& y) {return x > y;}
bool lessEqualCpp(double& x, double& y) {return x <= y;}
bool greaterEqualCpp(double& x, double& y) {return x >= y;}
bool equalCpp(double& x, double& y) {return x == y;}

Rcpp::XPtr<funcPtr> putFunPtrInXPtr(std::string fstr) {
    if (fstr == "prod")
        return(Rcpp::XPtr<funcPtr>(new funcPtr(&prodCpp)));
    else if (fstr == "sum")
        return(Rcpp::XPtr<funcPtr>(new funcPtr(&sumCpp)));
    else if (fstr == "mean")
        return(Rcpp::XPtr<funcPtr>(new funcPtr(&meanCpp)));
    else if (fstr == "max")
        return(Rcpp::XPtr<funcPtr>(new funcPtr(&maxCpp)));
    else if (fstr == "min")
        return(Rcpp::XPtr<funcPtr>(new funcPtr(&minCpp)));
    else
        return Rcpp::XPtr<funcPtr>(R_NilValue); // runtime error as NULL no XPtr
}

Rcpp::XPtr<compPtr> putCompPtrInXPtr(std::string fstr) {
    if (fstr == "<")
        return(Rcpp::XPtr<compPtr>(new compPtr(&lessCpp)));
    else if (fstr == ">")
        return(Rcpp::XPtr<compPtr>(new compPtr(&greaterCpp)));
    else if (fstr == "<=")
        return(Rcpp::XPtr<compPtr>(new compPtr(&lessEqualCpp)));
    else if (fstr == ">=")
        return(Rcpp::XPtr<compPtr>(new compPtr(&greaterEqualCpp)));
    else if (fstr == "==")
        return(Rcpp::XPtr<compPtr>(new compPtr(&equalCpp)));
    else
        return Rcpp::XPtr<compPtr>(R_NilValue); // runtime error as NULL no XPtr
}
