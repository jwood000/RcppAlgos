#ifndef CLEAN_CONVERT_H
#define CLEAN_CONVERT_H

#include "importExportMPZ.h"

const double Significand53 = 9007199254740991.0;

namespace CleanConvert {
    
    inline bool convertLogical(SEXP boolInput, std::string nameOfBool) {
        bool result = false;
        
        if (!Rf_isNull(boolInput)) {
            if (TYPEOF(boolInput) == LGLSXP) {
                double dblInp = Rcpp::as<double>(boolInput);
                int intTest = static_cast<int>(dblInp);
                
                if (dblInp != intTest)
                    Rcpp::stop("NA values are not allowed for " + nameOfBool);
                
                result = Rcpp::as<bool>(boolInput);
            } else {
                Rcpp::stop("Only logical values are supported for " + nameOfBool);
            }
        }
        
        return result;
    }

    template <typename stdType>
    inline void convertPrimitive(SEXP input, stdType &result, std::string nameOfObject,
                          bool strPoss = false, bool checkWhole = true) {
        
        switch(TYPEOF(input)) {
            case REALSXP:
            case INTSXP: {
                double dblInp = Rcpp::as<double>(input);
                int64_t intTest = static_cast<int64_t>(dblInp);
                if (checkWhole && intTest != dblInp)
                    Rcpp::stop(nameOfObject + " must be a whole number");
                
                result = Rcpp::as<stdType>(input);
                break;
            }
            case RAWSXP:
            case STRSXP: {
                if (!strPoss)
                    Rcpp::stop(nameOfObject + " must be of type numeric or integer");
                
                mpz_t temp[1];
                mpz_init(temp[0]);
                createMPZArray(input, temp, 1);
                double dblTemp = mpz_get_d(temp[0]);
                result = dblTemp;
                mpz_clear(temp[0]);
                break;
            }
            default:
                Rcpp::stop("This type is not supported! No conversion possible for " + nameOfObject);
        }
    }
    
    template <typename stdType>
    inline void convertVector(SEXP input, std::vector<stdType> &result, std::string nameOfObject,
                       bool checkWhole = true, bool numOnly = true) {
        
        int total = Rf_length(input);
        
        switch(TYPEOF(input)) {
            case REALSXP:
            case INTSXP: {
                if (checkWhole) {
                    std::vector<double> vecCheck = Rcpp::as<std::vector<double>>(input);
                    result.resize(vecCheck.size());
                    
                    for (std::size_t i = 0; i < vecCheck.size(); ++i) {
                        if (static_cast<int64_t>(vecCheck[i]) != vecCheck[i])
                            Rcpp::stop("Each element in " + nameOfObject + " must be a whole number");
                        
                        result[i] = static_cast<stdType>(vecCheck[i]);
                    }
                } else {
                    result = Rcpp::as<std::vector<stdType>>(input);
                }
                
                break;
            }
            case RAWSXP: {
                if (numOnly)
                    Rcpp::stop(nameOfObject + " must be of type numeric or integer");
                
                const char* raw = (char*)RAW(input);
                total = ((int*)raw)[0];
                // do not put a break here. Fall to
                // the next case for complete conversion
            }
            case STRSXP: {
                if (numOnly)
                    Rcpp::stop(nameOfObject + " must be of type numeric or integer");
                
                mpz_t *temp;
                temp = (mpz_t *) malloc(sizeof(mpz_t) * total);
                
                for (int i = 0; i < total; ++i)
                    mpz_init(temp[i]);
                
                createMPZArray(input, temp, total);
                std::vector<double> dblTemp(total);
                
                for (int i = 0; i < total; ++i) {
                    dblTemp[i] = mpz_get_d(temp[i]);
                    result.push_back(static_cast<stdType>(dblTemp[i]));
                }
                
                for (int i = 0; i < total; ++i)
                    mpz_clear(temp[i]);
                
                free(temp);
                break;
            }
            default:
                Rcpp::stop("This type is not supported! No conversion possible for " + nameOfObject);
        }
    }
}

#endif
