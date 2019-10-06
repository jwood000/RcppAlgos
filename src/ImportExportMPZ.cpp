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

constexpr std::size_t intSize = sizeof(int);

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
    
    r[1] = static_cast<int>(mpz_sgn(value));
    mpz_export(&r[2], 0, 1, intSize, 0, 0, value);
    
    return totals;
}
