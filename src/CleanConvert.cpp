#include "CleanConvert.h"
#include "Cpp14MakeUnique.h"
#include "ImportExportMPZ.h"
#include <string>
#include <cmath>

namespace CleanConvert {
    
    template <typename T>
    std::vector<T> GetNumVec(SEXP Rv) {
        std::vector<T> v;
        int len = Rf_length(Rv);
        
        if (TYPEOF(Rv) == REALSXP) {
            double* dblRv = REAL(Rv);
            v.assign(dblRv, dblRv + len);
        } else {
            int* intRv = INTEGER(Rv);
            v.assign(intRv, intRv + len);
        }
        
        return v;
    }

    bool convertLogical(SEXP boolInput, const std::string &nameOfBool) {
        bool result = false;
        std::string msg;
        
        if (!Rf_isNull(boolInput)) {
            if (TYPEOF(boolInput) == LGLSXP) {
                if (Rf_length(boolInput) > 1) {
                    Rf_error("Expecting a single value for %s", nameOfBool.c_str());
                }
                
                double dblInp = Rf_asReal(boolInput);
                
                if (ISNAN(dblInp) || std::isnan(dblInp))
                    Rf_error("%s cannot be NA or NaN", nameOfBool.c_str());
                
                if (std::abs(dblInp) > Significand53)
                    Rf_error("Only logical values are allowed for %s", nameOfBool.c_str());
                
                result = Rf_asLogical(boolInput);
            } else {
                Rf_error("Only logical values are supported for %s", nameOfBool.c_str());
            }
        }
        
        return result;
    }
    
    template <typename stdType>
    void convertPrimitive(SEXP input, stdType &result,
                          const std::string &nameOfObject,
                          bool numOnly, bool checkWhole, 
                          bool negPoss, bool decimalFraction) {
        
        const stdType maxType = std::numeric_limits<stdType>::max();
        
        switch(TYPEOF(input)) {
            case REALSXP:
            case INTSXP: {
                const double dblInp = Rf_asReal(input);
                const double posDblInp = std::abs(dblInp);
                
                if (ISNAN(dblInp) || std::isnan(dblInp))
                    Rf_error("%s cannot be NA or NaN", nameOfObject.c_str());
                
                if (!negPoss) {
                    if (decimalFraction) {
                        if (dblInp < 0) Rf_error("%s must be a positive number", nameOfObject.c_str());
                    } else if (dblInp < 1) {
                        Rf_error("%s must be a positive whole number", nameOfObject.c_str());
                    }
                }
                
                if (posDblInp > maxType) {
                    Rf_error("The abs value of %s must be less than or equal to %s",
                             nameOfObject.c_str(), std::to_string(maxType).c_str());
                }
                
                if (posDblInp > Significand53)
                    Rf_error("The abs value of %s must be less than 2^53", nameOfObject.c_str());
                
                if (checkWhole && static_cast<int64_t>(dblInp) != dblInp)
                    Rf_error("%s must be a whole number", nameOfObject.c_str());
                    
                result = static_cast<stdType>(Rf_asReal(input));
                break;
            }
            case RAWSXP:
            case STRSXP: {
                if (numOnly) {
                    Rf_error("%s must be of type numeric or integer", nameOfObject.c_str());
                }
                
                mpz_t temp[1];
                mpz_init(temp[0]);
                createMPZArray(input, temp, 1, nameOfObject, negPoss);
                const double dblTemp = mpz_get_d(temp[0]);
                const double posDblTemp = std::abs(dblTemp);
                
                if (ISNAN(dblTemp) || std::isnan(dblTemp))
                    Rf_error("%s cannot be NA or NaN", nameOfObject.c_str());
                
                if (!negPoss) {
                    if (decimalFraction) {
                        if (dblTemp < 0) Rf_error("%s must be a positive number", nameOfObject.c_str());
                    } else if (dblTemp < 1) {
                        Rf_error("%s must be a positive whole number", nameOfObject.c_str());
                    }
                }
                
                if (posDblTemp > maxType) {
                    Rf_error("The abs value of %s must be less than or equal to %s",
                             std::to_string(maxType).c_str(), nameOfObject.c_str());
                }
                
                if (posDblTemp > Significand53)
                    Rf_error("The abs value of %s must be less than 2^53", nameOfObject.c_str());
                
                if (checkWhole && static_cast<int64_t>(dblTemp) != dblTemp)
                    Rf_error("%s must be a whole number", nameOfObject.c_str());
                    
                result = dblTemp;
                mpz_clear(temp[0]);
                break;
            }
            default: {
                Rf_error("This type is not supported! No conversion possible for %s", nameOfObject.c_str());
            }
        }
    }
    
    template <typename stdType>
    void convertVector(SEXP input, std::vector<stdType> &result,
                       const std::string &nameOfObject,
                       bool numOnly, bool checkWhole,
                       bool negPoss) {
        
        int total = Rf_length(input);
        const stdType maxType = std::numeric_limits<stdType>::max();
        
        switch(TYPEOF(input)) {
            case REALSXP:
            case INTSXP: {
                std::vector<double> vecCheck = GetNumVec<double>(input);
                result.resize(vecCheck.size());
                
                for (std::size_t i = 0; i < vecCheck.size(); ++i) {
                    const double posDblInp = std::abs(vecCheck[i]);
                    
                    if (ISNAN(vecCheck[i]) || std::isnan(vecCheck[i]))
                        Rf_error("%s cannot be NA or NaN", nameOfObject.c_str());
                    
                    if (!negPoss && vecCheck[i] < 1)
                        Rf_error("Each element in %s must be a positive number", nameOfObject.c_str());
                    
                    if (posDblInp > maxType) {
                        Rf_error("The abs value of each element in %s must be less than or equal to %s",
                                 std::to_string(maxType).c_str(), nameOfObject.c_str());
                    }
                    
                    if (posDblInp > Significand53)
                        Rf_error("The abs value of each element in %s must be less than 2^53", nameOfObject.c_str());
                    
                    if (checkWhole && static_cast<int64_t>(vecCheck[i]) != vecCheck[i])
                        Rf_error("Each element in %s must be a whole number", nameOfObject.c_str());
                        
                    result[i] = static_cast<stdType>(vecCheck[i]);
                }
                
                break;
            }
            case RAWSXP: {
                if (numOnly)
                    Rf_error("%s must be of type numeric or integer", nameOfObject.c_str());
                
                const char* raw = (char*) RAW(input);
                total = ((int*) raw)[0];
                // do not put a break here. Fall to
                // the next case for complete conversion
            }
            case STRSXP: {
                if (numOnly)
                    Rf_error("%s must be of type numeric or integer", nameOfObject.c_str());
                
                auto temp = FromCpp14::make_unique<mpz_t[]>(total);
                
                for (int i = 0; i < total; ++i)
                    mpz_init(temp[i]);
                
                createMPZArray(input, temp.get(), total, nameOfObject, negPoss);
                std::vector<double> dblTemp(total);
                result.resize(total);
                
                for (int i = 0; i < total; ++i) {
                    dblTemp[i] = mpz_get_d(temp[i]);
                    const double posDblInp = std::abs(dblTemp[i]);
                    
                    if (ISNAN(dblTemp[i]) || std::isnan(dblTemp[i]))
                        Rf_error("%s cannot be NA or NaN", nameOfObject.c_str());
                    
                    if (!negPoss && dblTemp[i] < 1)
                        Rf_error("Each element in %s must be a positive number", nameOfObject.c_str());
                    
                    if (posDblInp > maxType) {
                        Rf_error("The abs value of each element in %s must be less than or equal to %s",
                                 std::to_string(maxType).c_str(), nameOfObject.c_str());
                    }
                    
                    if (posDblInp > Significand53)
                        Rf_error("The abs value of each element in %s must be less than 2^53", nameOfObject.c_str());
                    
                    if (checkWhole && static_cast<int64_t>(dblTemp[i]) != dblTemp[i])
                        Rf_error("Each element in %s must be a whole number", nameOfObject.c_str());
                    
                    result[i] = static_cast<stdType>(dblTemp[i]);
                }
                
                for (int i = 0; i < total; ++i)
                    mpz_clear(temp[i]);
                
                break;
            }
            default: {
                Rf_error("This type is not supported! No conversion possible for %s", nameOfObject.c_str());
            }
        }
    }
}

template void CleanConvert::convertPrimitive(SEXP, int&, const std::string&, bool, bool, bool, bool);
template void CleanConvert::convertPrimitive(SEXP, double&, const std::string&, bool, bool, bool, bool);

template void CleanConvert::convertVector(SEXP, std::vector<int>&, const std::string&, bool, bool, bool);
template void CleanConvert::convertVector(SEXP, std::vector<double>&, const std::string&, bool, bool, bool);

template std::vector<int> CleanConvert::GetNumVec(SEXP);
template std::vector<double> CleanConvert::GetNumVec(SEXP);
