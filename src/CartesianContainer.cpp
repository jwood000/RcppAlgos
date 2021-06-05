#include "CleanConvert.h"
#include "ComboCartesian.h"
#include "Eratosthenes.h"

enum rcppType {
    tInt = 0,
    tDbl = 1,
    tStr = 2,
    tLog = 3,
    tFac = 4
};

void convertToString(std::vector<std::string> &tempVec,
                     SEXP ListElement, rcppType &typePass, bool bFac) {
    
    switch(TYPEOF(ListElement)) {
        case INTSXP: {
            if (bFac) {
                Rcpp::IntegerVector facVec = Rcpp::as<Rcpp::IntegerVector>(ListElement);
                std::vector<std::string> strVec = Rcpp::as<std::vector<std::string>>(facVec.attr("levels"));
                typePass = tFac;
                
                for (auto v: strVec) {
                    bool isNum = (!v.empty() && v.find_first_not_of("0123456789.") == std::string::npos);
                    
                    if (isNum) {
                        double n = std::atof(v.c_str());
                        if (std::floor(n) != n) v += "dbl";
                    } else {
                        v += "str";
                    }
                    
                    tempVec.push_back(v);
                }
            } else {
                std::vector<int> intVec = Rcpp::as<std::vector<int>>(ListElement);
                typePass = tInt;
                
                for (auto v: intVec)
                    tempVec.push_back(std::to_string(v));
            }
            
            break;
        }
        case LGLSXP: {
            std::vector<int> intVec = Rcpp::as<std::vector<int>>(ListElement);
            typePass = tLog;
            
            for (auto v: intVec) {
                std::string r = std::to_string(v);
                tempVec.push_back(r + "log");
            }
            
            break;
        }
        case REALSXP: {
            std::vector<double> dblVec = Rcpp::as<std::vector<double>>(ListElement);
            typePass = tDbl;
            
            for (auto v: dblVec) {
                std::string r = std::to_string(v);
                
                while (r.back() == '0')
                    r.pop_back();
                
                if (std::floor(v) != v) r += "dbl";
                tempVec.push_back(r);
            }
            
            break;
        }
        case STRSXP: {
            std::vector<std::string> strVec = Rcpp::as<std::vector<std::string>>(ListElement);
            typePass = tStr;
            
            for (auto v: strVec) {
                bool isNum = (!v.empty() && v.find_first_not_of("0123456789.") == std::string::npos);
                
                if (isNum) {
                    double n = std::atof(v.c_str());
                    if (std::floor(n) != n) v += "dbl";
                } else {
                    v += "str";
                }
                
                tempVec.push_back(v);
            }
            
            break;
        }
    }
}

template <typename typeVec, typename typeRcpp>
void GetPureOutput(const std::vector<int> &cartCombs,
                   const Rcpp::List &RList, const typeVec &standardVec, 
                   typeRcpp &result, std::size_t nCols, std::size_t nRows) {
    
    for (std::size_t i = 0, row = 0; i < nRows; ++i, row += nCols)
        for (std::size_t j = 0; j < nCols; ++j)
            result(i, j) = standardVec[cartCombs[row + j]];
    
    if (!Rf_isNull(RList.attr("names"))) {
        Rcpp::CharacterVector myNames = RList.attr("names");
        colnames(result) = myNames;
    }
}

SEXP GlueComboCart(const std::vector<int> &cartCombs,
                   Rcpp::List &RList, const Rcpp::IntegerVector &intVec,
                   const Rcpp::LogicalVector &boolVec, const Rcpp::NumericVector &dblVec,
                   const Rcpp::CharacterVector &charVec, 
                   const std::vector<std::vector<int>> &facList,
                   std::vector<int> &typeCheck, bool IsDF, std::size_t nRows,
                   std::size_t nCols, std::vector<int> IsFactor) {
    
    if (IsDF) {
        Rcpp::List resList(nCols);
        
        for (std::size_t i = 0, facInd = 0; i < nCols; ++i) {
            switch (TYPEOF(RList[i])) {
                case INTSXP: {
                    Rcpp::IntegerVector rcppVec = Rcpp::no_init_vector(nRows);

                    if (IsFactor[i]) {
                        Rcpp::IntegerVector facVec(facList[facInd].cbegin(), facList[facInd].cend());
                        
                        for (std::size_t j = 0, row = i; j < nRows; ++j, row += nCols)
                            rcppVec[j] = facVec[cartCombs[row]];

                        Rcpp::IntegerVector testFactor = Rcpp::as<Rcpp::IntegerVector>(RList[i]);
                        Rcpp::CharacterVector myClass = testFactor.attr("class");
                        Rcpp::CharacterVector myLevels = testFactor.attr("levels");
                        rcppVec.attr("class") = myClass;
                        rcppVec.attr("levels") = myLevels;
                        ++facInd;
                    } else {
                        for (std::size_t j = 0, row = i; j < nRows; ++j, row += nCols)
                            rcppVec[j] = intVec[cartCombs[row]];
                    }

                    resList[i] = rcppVec;
                    break;
                }
                case LGLSXP: {
                    Rcpp::LogicalVector rcppVec = Rcpp::no_init_vector(nRows);

                    for (std::size_t j = 0, row = i; j < nRows; ++j, row += nCols)
                        rcppVec[j] = boolVec[cartCombs[row]];

                    resList[i] = rcppVec;
                    break;
                }
                case REALSXP: {
                    Rcpp::NumericVector rcppVec = Rcpp::no_init_vector(nRows);

                    for (std::size_t j = 0, row = i; j < nRows; ++j, row += nCols)
                        rcppVec[j] = dblVec[cartCombs[row]];

                    resList[i] = rcppVec;
                    break;
                }
                case STRSXP: {
                    Rcpp::CharacterVector rcppVec = Rcpp::no_init_vector(nRows);

                    for (std::size_t j = 0, row = i; j < nRows; ++j, row += nCols)
                        rcppVec[j] = charVec[cartCombs[row]];

                    resList[i] = rcppVec;
                    break;
                }
            }
        }

        Rcpp::DataFrame myDF(resList);
        myDF.attr("names") = RList.attr("names");
        return myDF;

    } else {
        if (typeCheck[tInt]) {
            Rcpp::IntegerMatrix intMat = Rcpp::no_init_matrix(nRows, nCols);
            GetPureOutput(cartCombs, RList, intVec, intMat, nCols, nRows);
            return intMat;
        } else if (typeCheck[tFac]) {
            Rcpp::IntegerMatrix intMat = Rcpp::no_init_matrix(nRows, nCols);
            GetPureOutput(cartCombs, RList, facList.front(), intMat, nCols, nRows);
            Rcpp::IntegerVector testFactor = Rcpp::as<Rcpp::IntegerVector>(RList[0]);
            Rcpp::CharacterVector myClass = testFactor.attr("class");
            Rcpp::CharacterVector myLevels = testFactor.attr("levels");
            intMat.attr("class") = myClass;
            intMat.attr("levels") = myLevels;
            return intMat;
        } else if (typeCheck[tLog]) {
            Rcpp::LogicalMatrix boolMat = Rcpp::no_init_matrix(nRows, nCols);
            GetPureOutput(cartCombs, RList, boolVec, boolMat, nCols, nRows);
            return boolMat;
        } else if (typeCheck[tDbl]) {
            Rcpp::NumericMatrix dblMat = Rcpp::no_init_matrix(nRows, nCols);
            GetPureOutput(cartCombs, RList, dblVec, dblMat, nCols, nRows);
            return dblMat;
        } else {
            Rcpp::CharacterMatrix charMat = Rcpp::no_init_matrix(nRows, nCols);
            GetPureOutput(cartCombs, RList, charVec, charMat, nCols, nRows);
            return charMat;
        }
    }
}

void getAtLeastNPrimes(std::vector<int> &primes, std::size_t sumLength) {

    double limit = 100;
    std::size_t guess = PrimeSieve::EstimatePiPrime(1.0, limit);

    while (guess < (1.1 * sumLength)) {
        limit *= 2;
        guess = PrimeSieve::EstimatePiPrime(1.0, limit);
    }
    
    bool tempPar = false;
    int_fast64_t intMin = static_cast<int_fast64_t>(1);
    int_fast64_t intMax = static_cast<int_fast64_t>(limit);
    std::vector<std::vector<int>> tempList;
    
    PrimeSieve::PrimeSieveMaster(intMin, intMax, primes, tempList, tempPar);
}

// [[Rcpp::export]]
SEXP comboGridRcpp(Rcpp::List RList, std::vector<int> IsFactor,
                   bool IsRep, std::size_t sumLength) {
    
    std::vector<int> primes;
    getAtLeastNPrimes(primes, sumLength);
    std::size_t numFactorVec = std::accumulate(IsFactor.cbegin(), IsFactor.cend(), 0);
    
    // All duplicates have been removed from RList via lapply(RList, function(x) sort(unique(x)))
    std::size_t nCols = RList.size();
    std::vector<std::vector<int>> myVec(nCols);
    std::unordered_map<std::string, int> mapIndex;
    std::vector<int> typeCheck(5, 0);
    
    Rcpp::CharacterVector charVec(sumLength);
    Rcpp::NumericVector dblVec(sumLength);
    Rcpp::IntegerVector intVec(sumLength);
    Rcpp::LogicalVector boolVec(sumLength);
    
    std::vector<std::vector<int>> facList(numFactorVec,
                                          std::vector<int>(sumLength, 0));
    
    for (std::size_t i = 0, total = 0, myIndex = 0, facInd = 0; i < nCols; ++i) {
        std::vector<std::string> tempVec;
        rcppType myType;
        
        convertToString(tempVec, RList[i], myType, IsFactor[i]);
        std::size_t j = 0;
        
        for (const auto &v: tempVec) {
            if (mapIndex.empty()) {
                mapIndex.insert({v, total});
            } else {
                auto mapCheck = mapIndex.find(v);
                
                if (mapCheck == mapIndex.end()) {
                    ++total;
                    mapIndex.insert({v, total});
                    myIndex = total;
                } else {
                    myIndex = mapCheck->second;
                }
            }
            
            myVec[i].push_back(myIndex);
            
            switch(myType) {
                case tInt: {
                    intVec[myIndex] = Rcpp::as<Rcpp::IntegerVector>(RList[i])[j];
                    typeCheck[tInt] = 1;
                    break;
                }
                case tFac: {
                    facList[facInd][myIndex] = Rcpp::as<std::vector<int>>(RList[i])[j];
                    typeCheck[tFac] = 1;
                    break;
                }
                case tDbl: {
                    dblVec[myIndex] = Rcpp::as<Rcpp::NumericVector>(RList[i])[j];
                    typeCheck[tDbl] = 1;
                    break;
                }
                case tStr: {
                    charVec[myIndex] = Rcpp::as<Rcpp::CharacterVector>(RList[i])[j];
                    typeCheck[tStr] = 1;
                    break;
                }
                case tLog: {
                    boolVec[myIndex] = Rcpp::as<Rcpp::LogicalVector>(RList[i])[j];
                    typeCheck[tLog] = 1;
                    break;
                }
            }
            
            ++j;
        }
        
        facInd += IsFactor[i];
    }
    
    int mySum = std::accumulate(typeCheck.cbegin(), typeCheck.cend(), 0);
    std::vector<std::string> testLevels;
    Rcpp::IntegerVector facVec(sumLength);
    
    if (typeCheck[tFac] && mySum == 1) {
        for (std::size_t i = 0; i < nCols; ++i) {
            if (IsFactor[i]) {
                facVec = Rcpp::as<Rcpp::IntegerVector>(RList[i]);
                std::vector<std::string> strVec = Rcpp::as<std::vector<std::string>>(facVec.attr("levels"));
                
                if (testLevels.size()) {
                    if (strVec != testLevels) {
                        ++mySum;
                        break;
                    }
                } else {
                    testLevels = strVec;
                }
            }
        }
    }
    
    bool IsDF = (mySum > 1) ? true : false;
    std::vector<int> cartCombs;
    comboGrid(cartCombs, IsRep, myVec, primes);
    
    const std::size_t nRows = cartCombs.size() / nCols;
    return GlueComboCart(cartCombs, RList, intVec, boolVec, dblVec,
                         charVec, facList, typeCheck, IsDF, nRows, nCols, IsFactor);
}
