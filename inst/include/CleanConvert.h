#ifndef CLEAN_CONVERT_H
#define CLEAN_CONVERT_H

#include "importExportMPZ.h"
#include <memory>

constexpr double Significand53 = 9007199254740991.0;

namespace CleanConvert {
    
    inline bool convertLogical(SEXP boolInput, std::string nameOfBool) {
        bool result = false;
        
        if (!Rf_isNull(boolInput)) {
            if (TYPEOF(boolInput) == LGLSXP) {
                double dblInp = Rcpp::as<double>(boolInput);
                
                if (Rcpp::NumericVector::is_na(dblInp) || std::isnan(dblInp))
                    Rcpp::stop(nameOfBool + " cannot be NA or NaN");
                
                if (std::abs(dblInp) > Significand53)
                    Rcpp::stop("Only logical values are allowed for " + nameOfBool);
                
                result = Rcpp::as<bool>(boolInput);
            } else {
                Rcpp::stop("Only logical values are supported for " + nameOfBool);
            }
        }
        
        return result;
    }

    template <typename stdType>
    inline void convertPrimitive(SEXP input, stdType &result, std::string nameOfObject,
                                 bool numOnly = true, bool checkWhole = true, bool negPoss = false) {
        
        const stdType maxType = std::numeric_limits<stdType>::max();
        
        switch(TYPEOF(input)) {
            case REALSXP:
            case INTSXP: {
                const double dblInp = Rcpp::as<double>(input);
                
                if (Rcpp::NumericVector::is_na(dblInp) || std::isnan(dblInp))
                    Rcpp::stop(nameOfObject + " cannot be NA or NaN");
                
                if (negPoss) {
                    if (std::abs(dblInp) > Significand53)
                        Rcpp::stop("The abs value of " + nameOfObject + " must be less than 2^53");
                } else {
                    if (dblInp < 1)
                        Rcpp::stop(nameOfObject + " must be a positive number");
                    
                    if (dblInp > maxType)
                        Rcpp::stop(nameOfObject + " must be less than or equal to " + std::to_string(maxType));
                    
                    if (dblInp > Significand53)
                        Rcpp::stop(nameOfObject + " must be less than 2^53");
                }

                if (checkWhole)
                    if (static_cast<int64_t>(dblInp) != dblInp)
                        Rcpp::stop(nameOfObject + " must be a whole number");
                
                result = Rcpp::as<stdType>(input);
                break;
            }
            case RAWSXP:
            case STRSXP: {
                if (numOnly)
                    Rcpp::stop(nameOfObject + " must be of type numeric or integer");
                
                mpz_t temp[1];
                mpz_init(temp[0]);
                createMPZArray(input, temp, 1, nameOfObject, negPoss);
                double dblTemp = mpz_get_d(temp[0]);
                
                if (Rcpp::NumericVector::is_na(dblTemp) || std::isnan(dblTemp))
                    Rcpp::stop(nameOfObject + " cannot be NA or NaN");
                
                if (negPoss) {
                    if (std::abs(dblTemp) > Significand53)
                        Rcpp::stop("The abs value of " + nameOfObject + " must be less than 2^53");
                } else {
                    if (dblTemp < 1)
                        Rcpp::stop(nameOfObject + " must be a positive number");
                    
                    if (dblTemp > maxType)
                        Rcpp::stop(nameOfObject + " must be less than or equal to " + std::to_string(maxType));
                    
                    if (dblTemp > Significand53)
                        Rcpp::stop(nameOfObject + " must be less than 2^53");
                }
                
                if (checkWhole)
                    if (static_cast<int64_t>(dblTemp) != dblTemp)
                        Rcpp::stop(nameOfObject + " must be a whole number");
                
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
                              bool numOnly = true, bool checkWhole = true, bool negPoss = false) {
        
        int total = Rf_length(input);
        const stdType maxType = std::numeric_limits<stdType>::max();
        
        switch(TYPEOF(input)) {
            case REALSXP:
            case INTSXP: {
                std::vector<double> vecCheck = Rcpp::as<std::vector<double>>(input);
                result.resize(vecCheck.size());
                
                for (std::size_t i = 0; i < vecCheck.size(); ++i) {
                    if (Rcpp::NumericVector::is_na(vecCheck[i]) || std::isnan(vecCheck[i]))
                        Rcpp::stop(nameOfObject + " cannot be NA or NaN");
                    
                    if (negPoss) {
                        if (std::abs(vecCheck[i]) > Significand53)
                            Rcpp::stop("The abs value of each element in " + nameOfObject + " must be less than 2^53");
                    } else {
                        if (vecCheck[i] < 1)
                            Rcpp::stop("Each element in " + nameOfObject + " must be a positive number");
                        
                        if (vecCheck[i] > maxType)
                            Rcpp::stop("Each element in " + nameOfObject + " must be less than or equal to " + std::to_string(maxType));
                        
                        if (vecCheck[i] > Significand53)
                            Rcpp::stop("Each element in " + nameOfObject + " must be less than 2^53");
                    }
                    
                    if (checkWhole)
                        if (static_cast<int64_t>(vecCheck[i]) != vecCheck[i])
                            Rcpp::stop("Each element in " + nameOfObject + " must be a whole number");
                    
                    result[i] = static_cast<stdType>(vecCheck[i]);
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
                
                auto temp = std::make_unique<mpz_t[]>(total);
                
                for (int i = 0; i < total; ++i)
                    mpz_init(temp[i]);
                
                createMPZArray(input, temp.get(), total, nameOfObject, negPoss);
                std::vector<double> dblTemp(total);
                result.resize(total);
                
                for (int i = 0; i < total; ++i) {
                    dblTemp[i] = mpz_get_d(temp[i]);
                    
                    if (Rcpp::NumericVector::is_na(dblTemp[i]) || std::isnan(dblTemp[i]))
                        Rcpp::stop(nameOfObject + " cannot be NA or NaN");
                    
                    if (negPoss) {
                        if (std::abs(dblTemp[i]) > Significand53)
                            Rcpp::stop("The abs value of each element in " + nameOfObject + " must be less than 2^53");
                    } else {
                        if (dblTemp[i] < 1)
                            Rcpp::stop("Each element in " + nameOfObject + " must be a positive number");
                        
                        if (dblTemp[i] > maxType)
                            Rcpp::stop("Each element in " + nameOfObject + " must be less than or equal to " + std::to_string(maxType));
                        
                        if (dblTemp[i] > Significand53)
                            Rcpp::stop("Each element in " + nameOfObject + " must be less than 2^53");
                    }
                    
                    if (checkWhole)
                        if (static_cast<int64_t>(dblTemp[i]) != dblTemp[i])
                            Rcpp::stop("Each element in " + nameOfObject + " must be a whole number");
                    
                    result[i] = static_cast<stdType>(dblTemp[i]);
                }
                
                for (int i = 0; i < total; ++i)
                    mpz_clear(temp[i]);
                
                break;
            }
            default:
                Rcpp::stop("This type is not supported! No conversion possible for " + nameOfObject);
        }
    }
}

#endif
