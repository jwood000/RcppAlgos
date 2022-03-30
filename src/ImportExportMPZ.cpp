// This file contains utility functions that
// are used for converting to and from type mpz_t,
// as well as sorting an array of type mpz_t.
//
// createMPZArray and myRaw are slightly modified versions
// of "bigvec create_vector(const SEXP & param)" and
// "int biginteger::as_raw(char* raw) const", respectively,
// from the source files bigintegerR.cc/ biginteger.cc from
// the R gmp package.

#include "ImportExportMPZ.h"
#include "SetUpUtils.h"
#include <cstring>

void createMPZArray(SEXP input, mpz_t *myVec, std::size_t vecSize,
                    const std::string &nameOfObject, bool negPoss) {

    const std::string suffix = (vecSize > 1) ?
                               "Each element in " + nameOfObject : nameOfObject;
    std::string myError;
    bool foundError = false;

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
                        if (negPoss) {
                            mpz_neg(myVec[i], myVec[i]);
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

                pos += intSize * (2 + (mpz_sizeinbase(myVec[i], 2) + numb - 1) / numb);
            }

            break;
        } case REALSXP: {
            double* dblInput = REAL(input);
            std::vector<double> dblVec(dblInput, dblInput + vecSize);
            constexpr double Sig53 = 9007199254740991.0;

            for (std::size_t j = 0; j < vecSize; ++j) {
                if (ISNAN(dblVec[j])) {
                    myError = suffix + " cannot be NA or NaN";
                    foundError = true;
                    break;
                }

                if (negPoss) {
                    if (std::abs(dblVec[j]) > Sig53) {
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

                    if (dblVec[j] > Sig53) {
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

                mpz_set_d(myVec[j], dblVec[j]);
            }

            break;
        }
        case INTSXP:
        case LGLSXP: {
            int* intInput = INTEGER(input);
            std::vector<double> dblVec(intInput, intInput + vecSize);
            std::vector<int> intVec(intInput, intInput + vecSize);

            for (std::size_t j = 0; j < vecSize; ++j) {
                if (ISNAN(dblVec[j])) {
                    myError = suffix + " cannot be NA or NaN";
                    foundError = true;
                    break;
                }

                if (!negPoss && intVec[j] < 1) {
                    myError = suffix + " must be a positive number";
                    foundError = true;
                    break;
                }

                mpz_set_si(myVec[j], intVec[j]);
            }

            break;
        } case STRSXP: {
            for (std::size_t i = 0; i < vecSize; ++i) {
                if (STRING_ELT(input, i) == NA_STRING) {
                    myError = suffix + " cannot be NA or NaN";
                    foundError = true;
                    break;
                } else {
                    mpz_set_str(myVec[i], CHAR(STRING_ELT(input, i)), 10);

                    if (!negPoss && mpz_sgn(myVec[i]) < 1) {
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

    if (foundError) {
        MpzClearVec(myVec, vecSize);
        cpp11::stop(myError.c_str());
    }
}

int myRaw(char* raw, const mpz_t value, std::size_t totals) {
    std::memset(raw, 0, totals);

    int* r = (int*)raw;
    r[0] = totals / intSize - 2;

    r[1] = static_cast<int>(mpz_sgn(value));
    mpz_export(&r[2], 0, 1, intSize, 0, 0, value);

    return totals;
}
