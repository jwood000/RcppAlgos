#include "CleanConvert.h"
#include "Cpp14MakeUnique.h"
#include "ImportExportMPZ.h"
#include <string>
#include <cmath>

namespace CleanConvert {

    bool CheckNA(double val, VecType myType) {
        if (myType == VecType::Integer) {
            return (ISNAN(val) || val == NA_INTEGER);
        } else {
            return ISNAN(val);
        }
    }

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

    SEXP GetCount(bool IsGmp, const mpz_t computedRowsMpz, double computedRows) {

        if (IsGmp) {
            constexpr std::size_t numb = 8 * intSize;
            const std::size_t sizeNum = intSize *
                (2 + (mpz_sizeinbase(computedRowsMpz, 2) + numb - 1) / numb);
            const std::size_t size = intSize + sizeNum;

            SEXP ans = PROTECT(Rf_allocVector(RAWSXP, size));
            char* rPos = (char*) RAW(ans);
            ((int*) rPos)[0] = 1; // first int is vector-size-header

            // current position in rPos[] (starting after vector-size-header)
            myRaw(&rPos[intSize], computedRowsMpz, sizeNum);
            Rf_setAttrib(ans, R_ClassSymbol, Rf_mkString("bigz"));
            UNPROTECT(1);
            return(ans);
        } else {
            if (computedRows > std::numeric_limits<int>::max()) {
                return Rf_ScalarReal(computedRows);
            } else {
                return Rf_ScalarInteger(static_cast<int>(computedRows));
            }
        }
    }

    bool convertFlag(SEXP boolInput, const std::string &nameOfBool) {
        bool result = false;
        std::string msg;

        if (!Rf_isNull(boolInput)) {
            if (TYPEOF(boolInput) == LGLSXP) {
                if (Rf_length(boolInput) > 1) {
                    Rf_error("Expecting a single value for %s", nameOfBool.c_str());
                }

                double dblInp = Rf_asReal(boolInput);

                if (CheckNA(dblInp, VecType::Integer)) {
                    Rf_error("%s cannot be NA or NaN", nameOfBool.c_str());
                }

                if (std::abs(dblInp) > Significand53) {
                    Rf_error("Only logical values are allowed for %s", nameOfBool.c_str());
                }

                result = Rf_asLogical(boolInput);
            } else {
                Rf_error("Only logical values are supported for %s", nameOfBool.c_str());
            }
        }

        return result;
    }

    template <typename T>
    void convertPrimitive(SEXP input, T &result, VecType myType,
                          const std::string &nameOfObject, bool numOnly,
                          bool checkWhole, bool negPoss, bool decimalFraction) {

        const T maxType = std::numeric_limits<T>::max();

        switch(TYPEOF(input)) {
            case REALSXP:
            case LGLSXP:
            case INTSXP: {
                const double dblInp = Rf_asReal(input);
                const double posDblInp = std::abs(dblInp);

                if (CheckNA(dblInp, myType)) {
                    Rf_error("%s cannot be NA or NaN", nameOfObject.c_str());
                }

                if (!negPoss) {
                    if (decimalFraction) {
                        if (dblInp < 0) Rf_error("%s must be a positive number", nameOfObject.c_str());
                    } else if (dblInp < 1) {
                        Rf_error("%s must be a positive whole number", nameOfObject.c_str());
                    }
                }

                if (checkWhole && static_cast<int64_t>(dblInp) != dblInp) {
                    Rf_error("%s must be a whole number", nameOfObject.c_str());
                }

                if (posDblInp > maxType) {
                    Rf_error("The abs value of %s must be less than or equal to %s",
                             nameOfObject.c_str(), std::to_string(maxType).c_str());
                }

                if (posDblInp > Significand53) {
                    Rf_error("The abs value of %s must be less than 2^53", nameOfObject.c_str());
                }

                result = static_cast<T>(Rf_asReal(input));
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

                if (CheckNA(dblTemp, myType)) {
                    Rf_error("%s cannot be NA or NaN", nameOfObject.c_str());
                }

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

                if (posDblTemp > Significand53) {
                    Rf_error("The abs value of %s must be less than 2^53", nameOfObject.c_str());
                }

                if (checkWhole && static_cast<int64_t>(dblTemp) != dblTemp) {
                    Rf_error("%s must be a whole number", nameOfObject.c_str());
                }

                result = dblTemp;
                mpz_clear(temp[0]);
                break;
            } default: {
                Rf_error("This type is not supported! No conversion possible for %s", nameOfObject.c_str());
            }
        }
    }

    template <typename T>
    void convertVector(SEXP input, std::vector<T> &result,
                       VecType myType, const std::string &nameOfObject,
                       bool numOnly, bool checkWhole, bool negPoss) {

        int total = Rf_length(input);
        const T maxType = std::numeric_limits<T>::max();

        switch(TYPEOF(input)) {
            case REALSXP:
            case LGLSXP:
            case INTSXP: {
                std::vector<double> vecCheck = GetNumVec<double>(input);
                result.resize(vecCheck.size());

                for (std::size_t i = 0; i < vecCheck.size(); ++i) {
                    const double posDblInp = std::abs(vecCheck[i]);

                    if (CheckNA(vecCheck[i], myType)) {
                        Rf_error("%s cannot be NA or NaN", nameOfObject.c_str());
                    }

                    if (!negPoss && vecCheck[i] < 1) {
                        Rf_error("Each element in %s must be a positive number", nameOfObject.c_str());
                    }

                    if (posDblInp > maxType) {
                        Rf_error("The abs value of each element in %s must be less than or equal to %s",
                                 std::to_string(maxType).c_str(), nameOfObject.c_str());
                    }

                    if (posDblInp > Significand53) {
                        Rf_error("The abs value of each element in %s must be less than 2^53", nameOfObject.c_str());
                    }

                    if (checkWhole && static_cast<int64_t>(vecCheck[i]) != vecCheck[i]) {
                        Rf_error("Each element in %s must be a whole number", nameOfObject.c_str());
                    }

                    result[i] = static_cast<T>(vecCheck[i]);
                }

                break;
            } case RAWSXP: {
                if (numOnly) {
                    Rf_error("%s must be of type numeric or integer", nameOfObject.c_str());
                }

                const char* raw = (char*) RAW(input);
                total = ((int*) raw)[0];
                // do not put a break here. Fall to
                // the next case for complete conversion
            } case STRSXP: {
                if (numOnly) {
                    Rf_error("%s must be of type numeric or integer", nameOfObject.c_str());
                }

                auto temp = FromCpp14::make_unique<mpz_t[]>(total);

                for (int i = 0; i < total; ++i) {
                    mpz_init(temp[i]);
                }

                createMPZArray(input, temp.get(), total, nameOfObject, negPoss);
                std::vector<double> dblTemp(total);
                result.resize(total);

                for (int i = 0; i < total; ++i) {
                    dblTemp[i] = mpz_get_d(temp[i]);
                    const double posDblInp = std::abs(dblTemp[i]);

                    if (CheckNA(dblTemp[i], myType)) {
                        Rf_error("%s cannot be NA or NaN", nameOfObject.c_str());
                    }

                    if (!negPoss && dblTemp[i] < 1) {
                        Rf_error("Each element in %s must be a positive number", nameOfObject.c_str());
                    }

                    if (posDblInp > maxType) {
                        Rf_error("The abs value of each element in %s must be less than or equal to %s",
                                 std::to_string(maxType).c_str(), nameOfObject.c_str());
                    }

                    if (posDblInp > Significand53) {
                        Rf_error("The abs value of each element in %s must be less than 2^53", nameOfObject.c_str());
                    }

                    if (checkWhole && static_cast<int64_t>(dblTemp[i]) != dblTemp[i]) {
                        Rf_error("Each element in %s must be a whole number", nameOfObject.c_str());
                    }

                    result[i] = static_cast<T>(dblTemp[i]);
                }

                for (int i = 0; i < total; ++i) {
                    mpz_clear(temp[i]);
                }

                break;
            } default: {
                Rf_error("This type is not supported! No conversion possible for %s", nameOfObject.c_str());
            }
        }
    }
}

template void CleanConvert::convertPrimitive(SEXP, int&, VecType, const std::string&,
                                             bool, bool, bool, bool);
template void CleanConvert::convertPrimitive(SEXP, double&, VecType, const std::string&,
                                             bool, bool, bool, bool);

template void CleanConvert::convertVector(SEXP, std::vector<int>&, VecType,
                                          const std::string&, bool, bool, bool);
template void CleanConvert::convertVector(SEXP, std::vector<double>&, VecType,
                                          const std::string&, bool, bool, bool);

template std::vector<int> CleanConvert::GetNumVec(SEXP);
template std::vector<double> CleanConvert::GetNumVec(SEXP);
