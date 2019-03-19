/* This file contains utility functions that
 * are used for converting to and from type mpz_t,
 * as well as sorting an array of type mpz_t.
 * 
 * createMPZArray and myRaw are slightly modified versions
 * of "bigvec create_vector(const SEXP & param)" and 
 * "int biginteger::as_raw(char* raw) const", respectively,
 * from the source files bigintegerR.cc/ biginteger.cc from
 * the R gmp package.
 */

#include <gmp.h>
#include <Rcpp.h>

void createMPZArray(SEXP v, mpz_t myVec[], unsigned long int sizevec) {
    switch (TYPEOF(v)) {
        case RAWSXP: {
            // deserialise the vector. first int is the size.
            const char* raw = (char*)RAW(v);
            const unsigned long int intSize = sizeof(int);
            const unsigned long int numb = 8 * intSize;
            int pos = intSize; // position in raw[]. Starting after header.
            
            for (std::size_t i = 0; i < sizevec; ++i) {
                const int* r = (int*)(&raw[pos]);
                
                if (r[0] > 0) {
                    mpz_import(myVec[i], r[0], 1, intSize, 0, 0, (void*)&(r[2]));
                    if(r[1] == -1) mpz_neg(myVec[i], myVec[i]);
                } else {
                    mpz_set_si(myVec[i], 0);
                }
                
                pos += intSize * (2 + (mpz_sizeinbase(myVec[i], 2) + numb - 1) / numb);
            }
            
            break;
        }
        case REALSXP: {
            double* myDbl = REAL(v);
            
            for (std::size_t j = 0; j < sizevec; ++j) {
                /// New : numeric '+- Inf' give +- "Large" instead of NA
                double dj = myDbl[j];
                
                if (std::isnan(dj) || static_cast<int64_t>(dj) != dj) {
                    Rcpp::stop("Must be a whole number. If number is too large for "
                               "double precision, consider using gmp::as.bigz or "
                               "passing sampleVec as a character vector.");
                } else {
                    mpz_set_d(myVec[j], dj);
                }
            }
            
            break;
        }
        case INTSXP:
        case LGLSXP: {
            int* myInt = INTEGER(v);
            for (std::size_t j = 0; j < sizevec; ++j)
                mpz_set_si(myVec[j], myInt[j]);
            
            break;
        }
        case STRSXP: {
            for (std::size_t i = 0; i < sizevec; ++i) {
                if (STRING_ELT(v, i) == NA_STRING)
                    mpz_set_si(myVec[i], 0);
                else
                    mpz_set_str(myVec[i], CHAR(STRING_ELT(v, i)), 10);
            }
            
            break;
        }
        default:
            Rcpp::stop("only logical, numeric or character (atomic) vectors can be coerced to 'bigz'");
    }
}

int myRaw(char* raw, mpz_t value, unsigned long int totals) {
    memset(raw, 0, totals);
    
    const unsigned long int intSize = sizeof(int);
    int* r = (int*)raw;
    r[0] = totals / intSize - 2;
    
    r[1] = (int) mpz_sgn(value);
    mpz_export(&r[2], 0, 1, intSize, 0, 0, value);
    
    return totals;
}
