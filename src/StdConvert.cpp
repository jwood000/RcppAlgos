#include "CppConvert/StdConvert.h"
#include "CppConvert/GmpConvert.h"

namespace CppConvert {

    bool convertFlag(SEXP boolInput, const std::string &nameOfBool) {
        bool result = false;
        std::string msg;

        if (!Rf_isNull(boolInput)) {
            if (TYPEOF(boolInput) == LGLSXP) {
                if (Rf_length(boolInput) > 1) {
                    cpp11::stop("Expecting a single value for %s",
                                nameOfBool.c_str());
                }

                double dblInp = Rf_asReal(boolInput);

                if (CheckNA(dblInp, VecType::Integer)) {
                    cpp11::stop("%s cannot be NA or NaN",
                                nameOfBool.c_str());
                }

                if (std::abs(dblInp) > Significand53) {
                    cpp11::stop("Only logical values are allowed for %s",
                                nameOfBool.c_str());
                }

                result = Rf_asLogical(boolInput);
            } else {
                cpp11::stop("Only logical values are supported for %s",
                            nameOfBool.c_str());
            }
        }

        return result;
    }

    template <typename T>
    void convertPrimitive(SEXP input, T &result, VecType myType,
                          const std::string &nameOfObject, bool numOnly,
                          bool checkWhole, bool negPoss,
                          bool decimalFraction) {

        const T maxType = std::numeric_limits<T>::max();

        switch(TYPEOF(input)) {
            case REALSXP:
            case LGLSXP:
            case INTSXP: {
                const double dblInp = Rf_asReal(input);
                const double posDblInp = std::abs(dblInp);

                if (CheckNA(dblInp, myType)) {
                    cpp11::stop("%s cannot be NA or NaN",
                                nameOfObject.c_str());
                }

                if (!negPoss) {
                    if (decimalFraction && dblInp < 0) {
                        cpp11::stop("%s must be a positive number",
                                    nameOfObject.c_str());
                    } else if (!decimalFraction && dblInp < 1) {
                        cpp11::stop("%s must be a positive whole number",
                                    nameOfObject.c_str());
                    }
                }

                if (checkWhole && static_cast<int64_t>(dblInp) != dblInp) {
                    cpp11::stop("%s must be a whole number",
                                nameOfObject.c_str());
                }

                if (posDblInp > maxType) {
                    std::string msg = "The abs value of " + nameOfObject +
                        " must be less than or equal to " +
                        std::to_string(maxType);
                    cpp11::stop(msg.c_str());
                }

                if (posDblInp > Significand53) {
                    cpp11::stop("The abs value of %s must be less than 2^53",
                                nameOfObject.c_str());
                }

                result = static_cast<T>(Rf_asReal(input));
                break;
            }
            case RAWSXP:
            case STRSXP: {
                if (numOnly) {
                    cpp11::stop("%s must be of type numeric or integer",
                                nameOfObject.c_str());
                }

                mpz_class temp;
                CppConvert::convertMpzClass(input, temp, nameOfObject,
                                            negPoss);
                const double dblTemp = temp.get_d();
                const double posDblTemp = std::abs(dblTemp);

                if (CheckNA(dblTemp, myType)) {
                    cpp11::stop("%s cannot be NA or NaN",
                                nameOfObject.c_str());
                }

                if (!negPoss) {
                    if (decimalFraction && dblTemp < 0) {
                        cpp11::stop("%s must be a positive number",
                                    nameOfObject.c_str());
                    } else if (!decimalFraction && dblTemp < 1) {
                        cpp11::stop("%s must be a positive whole number",
                                    nameOfObject.c_str());
                    }
                }

                if (posDblTemp > maxType) {
                    std::string msg = "The abs value of " + nameOfObject +
                        " must be less than or equal to " +
                        std::to_string(maxType);
                    cpp11::stop(msg.c_str());
                }

                if (posDblTemp > Significand53) {
                    cpp11::stop("The abs value of %s must be less than 2^53",
                                nameOfObject.c_str());
                }

                if (checkWhole && static_cast<int64_t>(dblTemp) != dblTemp) {
                    cpp11::stop("%s must be a whole number",
                                nameOfObject.c_str());
                }

                result = dblTemp;
                break;
            } default: {
                cpp11::stop(
                    "This type is not supported! No conversion",
                    " possible for %s", nameOfObject.c_str()
                );
            }
        }
    }

    template <typename T>
    void convertVector(SEXP input, std::vector<T> &result, VecType myType,
                       const std::string &nameOfObject, bool numOnly,
                       bool checkWhole, bool negPoss) {

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
                        cpp11::stop("%s cannot be NA or NaN",
                                    nameOfObject.c_str());
                    }

                    if (!negPoss && vecCheck[i] < 1) {
                        cpp11::stop(
                            "Each element in %s must be a positive number",
                            nameOfObject.c_str()
                        );
                    }

                    if (posDblInp > maxType) {
                        std::string msg = "The abs value of each element "
                        "in " + nameOfObject + " must be less than " +
                            std::to_string(maxType);
                        cpp11::stop(msg.c_str());
                    }

                    if (posDblInp > Significand53) {
                        std::string msg = "The abs value of each element in " +
                            nameOfObject + " must be less than 2^53";
                        cpp11::stop(msg.c_str());
                    }

                    if (checkWhole &&
                        static_cast<int64_t>(vecCheck[i]) != vecCheck[i]) {
                        cpp11::stop(
                            "Each element in %s must be a whole number",
                            nameOfObject.c_str()
                        );
                    }

                    result[i] = static_cast<T>(vecCheck[i]);
                }

                break;
            } case RAWSXP: {
                if (numOnly) {
                    cpp11::stop("%s must be of type numeric or integer",
                                nameOfObject.c_str());
                }

                const char* raw = (char*) RAW(input);
                total = ((int*) raw)[0];
                // do not put a break here. Fall to
                // the next case for complete conversion
            } case STRSXP: {
                if (numOnly) {
                    cpp11::stop("%s must be of type numeric or integer",
                                nameOfObject.c_str());
                }

                std::vector<mpz_class> temp(total);
                CppConvert::convertMPZVector(input, temp, total,
                                             nameOfObject, negPoss);

                std::vector<double> dblTemp(total);
                result.resize(total);

                for (int i = 0; i < total; ++i) {
                    dblTemp[i] = temp[i].get_d();
                    const double posDblInp = std::abs(dblTemp[i]);

                    if (CheckNA(dblTemp[i], myType)) {
                        cpp11::stop("%s cannot be NA or NaN",
                                    nameOfObject.c_str());
                    }

                    if (!negPoss && dblTemp[i] < 1) {
                        cpp11::stop(
                            "Each element in %s must be a positive number",
                            nameOfObject.c_str()
                        );
                    }

                    if (posDblInp > maxType) {
                        std::string msg = "The abs value of each element "
                        "in " + nameOfObject + " must be less than " +
                            std::to_string(maxType);
                        cpp11::stop(msg.c_str());
                    }

                    if (posDblInp > Significand53) {
                        std::string msg = "The abs value of each element "
                        "in " + nameOfObject + " must be less than 2^53";
                        cpp11::stop(msg.c_str());
                    }

                    if (checkWhole &&
                        static_cast<int64_t>(dblTemp[i]) != dblTemp[i]) {
                        cpp11::stop(
                            "Each element in %s must be a whole number",
                            nameOfObject.c_str()
                        );
                    }

                    result[i] = static_cast<T>(dblTemp[i]);
                }

                break;
            } default: {
                std::string msg = "This type is not supported! No "
                "conversion possible for " + nameOfObject;
                cpp11::stop(msg.c_str());
            }
        }
    }
}
template void CppConvert::convertPrimitive(
    SEXP, int&, VecType, const std::string&, bool, bool, bool, bool
);

template void CppConvert::convertPrimitive(
    SEXP, double&, VecType, const std::string&, bool, bool, bool, bool
);

template void CppConvert::convertVector(
    SEXP, std::vector<int>&, VecType, const std::string&, bool, bool, bool
);

template void CppConvert::convertVector(
    SEXP, std::vector<double>&, VecType, const std::string&, bool, bool, bool
);
