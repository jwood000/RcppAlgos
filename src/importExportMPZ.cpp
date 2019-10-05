// This file contains utility functions that
// are used for converting to and from type mpz_t,
// as well as sorting an array of type mpz_t.
// 
// createMPZArray and myRaw are slightly modified versions
// of "bigvec create_vector(const SEXP & param)" and 
// "int biginteger::as_raw(char* raw) const", respectively,
// from the source files bigintegerR.cc/ biginteger.cc from
// the R gmp package.

#include <gmp.h>
#include <Rcpp.h>
#include "Cpp14MakeUnique.h"

constexpr std::size_t intSize = sizeof(int);
constexpr std::size_t numb = 8u * intSize;

void createMPZArray(SEXP input, mpz_t *myVec, std::size_t vecSize, 
                    const std::string &nameOfObject, bool negPoss) {
    
    const std::string suffix = (vecSize > 1) ? 
                               "Each element in " + nameOfObject : nameOfObject;
    
    switch (TYPEOF(input)) {
        case RAWSXP: {
            // deserialise the vector. first int is the size.
            const char* raw = (char*)RAW(input);
            const std::size_t numb = 8 * intSize;
            int pos = intSize; // position in raw[]. Starting after header.
            
            for (std::size_t i = 0; i < vecSize; ++i) {
                const int* r = (int*)(&raw[pos]);
                
                if (r[0] > 0) {
                    mpz_import(myVec[i], r[0], 1, intSize, 0, 0, (void*)&(r[2]));
                    
                    if(r[1] == -1) {
                        if (negPoss)
                            mpz_neg(myVec[i], myVec[i]);
                        else
                            Rcpp::stop(suffix + " must be a positive number");
                    }
                } else {
                    Rcpp::stop(suffix + " cannot be NA or NaN");
                }
                
                pos += intSize * (2 + (mpz_sizeinbase(myVec[i], 2) + numb - 1) / numb);
            }
            
            break;
        }
        case REALSXP: {
            std::vector<double> dblVec = Rcpp::as<std::vector<double>>(input);
            constexpr double Sig53 = 9007199254740991.0;
            
            for (std::size_t j = 0; j < vecSize; ++j) {
                if (Rcpp::NumericVector::is_na(dblVec[j]) || std::isnan(dblVec[j]))
                    Rcpp::stop(suffix + " cannot be NA or NaN");
                
                if (negPoss) {
                    if (std::abs(dblVec[j]) > Sig53) {
                        Rcpp::stop("Number is too large for double precision. Consider "
                                       "using gmp::as.bigz or as.character for " + nameOfObject);
                    }
                } else {
                    if (dblVec[j] < 1)
                        Rcpp::stop(suffix + " must be a positive number");
                    
                    if (dblVec[j] > Sig53) {
                        Rcpp::stop("Number is too large for double precision. Consider "
                                       "using gmp::as.bigz or as.character for " + nameOfObject);
                    }
                }
                        
                if (static_cast<int64_t>(dblVec[j]) != dblVec[j])
                    Rcpp::stop(suffix + " must be a whole number.");
                
                mpz_set_d(myVec[j], dblVec[j]);
            }
            
            break;
        }
        case INTSXP:
        case LGLSXP: {
            std::vector<int> intVec = Rcpp::as<std::vector<int>>(input);
            std::vector<double> dblVec = Rcpp::as<std::vector<double>>(input);
            
            for (std::size_t j = 0; j < vecSize; ++j) {
                if (Rcpp::NumericVector::is_na(dblVec[j]) || std::isnan(dblVec[j]))
                    Rcpp::stop(suffix + " cannot be NA or NaN");
                
                if (!negPoss)
                    if (intVec[j] < 1)
                        Rcpp::stop(suffix + " must be a positive number");
                
                mpz_set_si(myVec[j], intVec[j]);
            }
            
            break;
        }
        case STRSXP: {
            for (std::size_t i = 0; i < vecSize; ++i) {
                if (STRING_ELT(input, i) == NA_STRING) {
                    Rcpp::stop(suffix + " cannot be NA or NaN");
                } else {
                    mpz_set_str(myVec[i], CHAR(STRING_ELT(input, i)), 10);
                    
                    if (!negPoss)
                        if (mpz_sgn(myVec[i]) < 1)
                            Rcpp::stop(suffix + " must be a positive whole number");
                }
            }
            
            break;
        }
        default:
            Rcpp::stop("This type is not supported! No conversion possible for " + nameOfObject);
    }
}

int myRaw(char* raw, mpz_t value, std::size_t totals) {
    memset(raw, 0, totals);
    
    int* r = (int*)raw;
    r[0] = totals / intSize - 2;
    
    r[1] = (int) mpz_sgn(value);
    mpz_export(&r[2], 0, 1, intSize, 0, 0, value);
    
    return totals;
}

SEXP BigMatrix(Rcpp::IntegerMatrix &indexMat, int nRows, int nCol, SEXP Rv) {
    
    Rcpp::RawVector raw(Rv);
    const int myVecSize = static_cast<int>(raw[0]);
    auto myVec = FromCpp14::make_unique<mpz_t[]>(myVecSize);
    createMPZArray(Rv, myVec.get(), myVecSize, "v", true);
    
    std::vector<std::size_t> mySizes(myVecSize);
    std::vector<std::size_t> vecIndFreq(myVecSize, 0u);
    
    for (std::size_t j = 0; j < myVecSize; j++) // adding each bigint's needed size
        mySizes[j] = intSize * (2 + (mpz_sizeinbase(myVec[j], 2) + numb - 1) / numb);
    
    for (int j = 0; j < nCol; ++j)
        for (int i = 0; i < nRows; ++i)
            ++vecIndFreq[indexMat(i, j)];
    
    // Add the size for vector size header
    std::size_t bigMatCount = intSize;
    
    for (std::size_t j = 0; j < myVecSize; ++j)
        bigMatCount += (mySizes[j] * vecIndFreq[j]);
    
    Rcpp::RawVector bigMat = Rcpp::no_init_vector(bigMatCount);
    char* rPos = (char*)(RAW(bigMat));
    ((int*)(rPos))[0] = nRows * nCol; // first int is vector-size-header
    
    // current position in rPos[] (starting after vector-
    // size-header) N.B. Only need for the first row.
    std::size_t posPos = intSize;
    
    for (std::size_t j = 0; j < nCol; ++j)
        for (std::size_t i = 0; i < nRows; ++i)
            posPos += myRaw(&rPos[posPos], myVec[indexMat(i, j)], mySizes[indexMat(i, j)]);
    
    bigMat.attr("class") = Rcpp::CharacterVector::create("bigz");
    bigMat.attr("nrow") = nRows;
    return(bigMat);
}
