#include "CppConvert/GmpConvert.h"

namespace CppConvert {

    void convertMPZVector(SEXP input, std::vector<mpz_class> &myVec,
                          std::size_t vecSize, const std::string &nameOfObject,
                          bool negPoss) {

        const std::string suffix = (vecSize > 1) ?
        "Each element in " + nameOfObject : nameOfObject;
        std::string myError;
        bool foundError = false;

        switch (TYPEOF(input)) {
            case RAWSXP: {
                // deserialise the vector. first int is the size.
                const char* raw = (char*) RAW(input);
                const std::size_t numb = 8 * intSize;
                int pos = intSize; // position in raw[]. Starting after header.

                for (std::size_t i = 0; i < vecSize; ++i) {
                    const int* r = (int*) (&raw[pos]);

                    if (r[0] > 0) {
                        mpz_import(myVec[i].get_mpz_t(), r[0], 1,
                                   intSize, 0, 0, (void*) & (r[2]));

                        if(r[1] == -1) {
                            if (negPoss) {
                                mpz_neg(myVec[i].get_mpz_t(),
                                        myVec[i].get_mpz_t());
                            } else {
                                myError = suffix + " must be a positive number";
                                foundError = true;
                                break;
                            }
                        }
                    } else {
                        myError = suffix + " cannot be NA or NaN";
                        foundError = true;
                        break;
                    }

                    pos += intSize * (2 + (
                        mpz_sizeinbase(myVec[i].get_mpz_t(), 2) + numb - 1
                    ) / numb);
                }

                break;
            } case REALSXP: {
                double* dblInput = REAL(input);
                std::vector<double> dblVec(dblInput, dblInput + vecSize);

                for (std::size_t j = 0; j < vecSize; ++j) {
                    if (ISNAN(dblVec[j])) {
                        myError = suffix + " cannot be NA or NaN";
                        foundError = true;
                        break;
                    }

                    if (negPoss) {
                        if (std::abs(dblVec[j]) > Significand53) {
                            myError = "Number is too large for double precision."
                            " Consider using gmp::as.bigz or "
                            "as.character for " + nameOfObject;
                            foundError = true;
                            break;
                        }
                    } else {
                        if (dblVec[j] < 1) {
                            myError = suffix + " must be a positive number";
                            foundError = true;
                            break;
                        }

                        if (dblVec[j] > Significand53) {
                            myError = "Number is too large for double precision."
                            " Consider using gmp::as.bigz or "
                            "as.character for " + nameOfObject;
                            foundError = true;
                            break;
                        }
                    }

                    if (static_cast<int64_t>(dblVec[j]) != dblVec[j]) {
                        myError = suffix + " must be a whole number";
                        foundError = true;
                        break;
                    }

                    myVec[j] = dblVec[j];
                }

                break;
            }
            case INTSXP:
            case LGLSXP: {
                int* intInput = INTEGER(input);
                std::vector<double> dblVec(intInput, intInput + vecSize);
                std::vector<int> intVec(intInput, intInput + vecSize);

                for (std::size_t j = 0; j < vecSize; ++j) {
                    if (ISNAN(dblVec[j]) || dblVec[j] == NA_INTEGER) {
                        myError = suffix + " cannot be NA or NaN";
                        foundError = true;
                        break;
                    }

                    if (!negPoss && intVec[j] < 1) {
                        myError = suffix + " must be a positive number";
                        foundError = true;
                        break;
                    }

                    myVec[j] = intVec[j];
                }

                break;
            } case STRSXP: {
                for (std::size_t i = 0; i < vecSize; ++i) {
                    if (STRING_ELT(input, i) == NA_STRING) {
                        myError = suffix + " cannot be NA or NaN";
                        foundError = true;
                        break;
                    } else {
                        mpz_set_str(myVec[i].get_mpz_t(),
                                    CHAR(STRING_ELT(input, i)), 10);

                        if (!negPoss && mpz_sgn(myVec[i].get_mpz_t()) < 1) {
                            myError = suffix + " must be a positive whole number";
                            foundError = true;
                            break;
                        }
                    }
                }

                break;
            } default: {
                myError = "This type is not supported! No conversion"
                " possible for " + nameOfObject;
                foundError = true;
                break;
            }
        }

        if (foundError) cpp11::stop(myError.c_str());
    }

    void convertMpzClass(SEXP input, mpz_class &result,
                         const std::string &nameOfObject, bool negPoss) {

        std::string myError;
        bool foundError = false;

        switch (TYPEOF(input)) {
            case RAWSXP: {
                // deserialise the vector. first int is the size.
                const char* raw = (char*) RAW(input);
                const int* r = (int*) (&raw[intSize]);

                if (r[0] > 0) {
                    mpz_import(result.get_mpz_t(), r[0], 1,
                               intSize, 0, 0, (void*) & (r[2]));

                    if(r[1] == -1) {
                        if (negPoss) {
                            mpz_neg(result.get_mpz_t(), result.get_mpz_t());
                        } else {
                            myError = nameOfObject + " must be a positive number";
                            foundError = true;
                            break;
                        }
                    }
                } else {
                    myError = nameOfObject + " cannot be NA or NaN";
                    foundError = true;
                    break;
                }

                break;
            } case REALSXP: {
                const double dblInput = Rf_asReal(input);

                if (ISNAN(dblInput)) {
                    myError = nameOfObject + " cannot be NA or NaN";
                    foundError = true;
                    break;
                }

                if (negPoss) {
                    if (std::abs(dblInput) > Significand53) {
                        myError = "Number is too large for double precision."
                        " Consider using gmp::as.bigz or "
                        "as.character for " + nameOfObject;
                        foundError = true;
                        break;
                    }
                } else {
                    if (dblInput < 1) {
                        myError = nameOfObject + " must be a positive number";
                        foundError = true;
                        break;
                    }

                    if (dblInput > Significand53) {
                        myError = "Number is too large for double precision."
                        " Consider using gmp::as.bigz or "
                        "as.character for " + nameOfObject;
                        foundError = true;
                        break;
                    }
                }

                if (static_cast<int64_t>(dblInput) != dblInput) {
                    myError = nameOfObject + " must be a whole number";
                    foundError = true;
                    break;
                }

                result = dblInput;
                break;
            }
            case INTSXP:
            case LGLSXP: {
                const int intInput = Rf_asInteger(input);
                const int dblInput = Rf_asReal(input);

                if (ISNAN(dblInput)) {
                    myError = nameOfObject + " cannot be NA or NaN";
                    foundError = true;
                    break;
                }

                if (!negPoss && intInput < 1) {
                    myError = nameOfObject + " must be a positive number";
                    foundError = true;
                    break;
                }

                result = intInput;
                break;
            } case STRSXP: {
                if (STRING_ELT(input, 0) == NA_STRING) {
                    myError = nameOfObject + " cannot be NA or NaN";
                    foundError = true;
                    break;
                } else {
                    mpz_set_str(result.get_mpz_t(),
                                CHAR(STRING_ELT(input, 0)), 10);

                    if (!negPoss && mpz_sgn(result.get_mpz_t()) < 1) {
                        myError = nameOfObject + " must be a positive whole number";
                        foundError = true;
                        break;
                    }
                }

                break;
            } default: {
                myError = "This type is not supported! No conversion"
                " possible for " + nameOfObject;
                foundError = true;
                break;
            }
        }

        if (foundError) cpp11::stop(myError.c_str());
    }
}
