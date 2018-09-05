#ifndef PERMUTATIONS_H
#define PERMUTATIONS_H

#include <CombPermUtils.h>

namespace Permutations {
    
    template <typename typeMatrix, typename typeVector>
    typeMatrix PermuteGeneral(int n, int r, typeVector &v, bool repetition,
                              int numRows, bool xtraCol, std::vector<int> &z, bool nonTrivial) {
        
        unsigned long int uN = n, uR = r, uRowN = numRows;
        unsigned long int numCols, lastElem = n - 1;
        unsigned long int lastCol = r - 1, pentultimate = n - 2;
        
        numCols = xtraCol ? (uR + 1) : uR;
        typeMatrix permuteMatrix = Rcpp::no_init_matrix(uRowN, numCols);
        
        if (repetition) {
            int lastElemInt = lastElem;
            
            for (std::size_t i = 0; i < uRowN; ++i) {
                for (std::size_t j = 0; j < uR; ++j)
                    permuteMatrix(i, j) = v[z[j]];
                
                for (int k = lastCol; k >= 0; --k) {
                    if (z[k] != lastElemInt) {
                        ++z[k];
                        break;
                    } else {
                        z[k] = 0;
                    }
                }
            }
            
        } else if (nonTrivial) {
            
            unsigned long int numR1 = numRows - 1;
            uint16_t *arrPerm = new uint16_t[uN];

            for (std::size_t i = 0; i < uN; ++i)
                arrPerm[i] = (uint16_t) z[i];
            
            if (r == n) {
                for (std::size_t i = 0; i < numR1; ++i) {
                    for (std::size_t j = 0; j < uR; ++j)
                        permuteMatrix(i, j) = v[arrPerm[j]];
                    
                    nextFullPerm(arrPerm, lastElem, pentultimate);
                }
            } else {
                for (std::size_t i = 0; i < numR1; ++i) {
                    for (std::size_t j = 0; j < uR; ++j)
                        permuteMatrix(i, j) = v[arrPerm[j]];
                    
                    nextPartialPerm(arrPerm, uR, lastCol, uN, lastElem);
                }
            }
            
            // Get last permutation
            for (std::size_t j = 0; j < uR; ++j)
                permuteMatrix(numR1, j) = v[arrPerm[j]];
            
            delete[] arrPerm;
            
        } else {
            
            if (n > 1) {
                unsigned long int phaseOne, maxN = NumPermsNoRep(n, r);
                unsigned long int segment = maxN / uN;
                phaseOne = (uRowN < segment) ? uRowN : segment;
    
                uint16_t *indexMat = new uint16_t[phaseOne * uR];
                uint16_t *arrPerm = new uint16_t[uN];
                
                for (std::size_t i = 0; i < uN; ++i)
                    arrPerm[i] = (uint16_t) i;
    
                if (r == n) {
                    for (std::size_t i = 0; i < phaseOne; ++i) {
                        for (std::size_t j = 0, k = i; j < uR; ++j, k += phaseOne) {
                            permuteMatrix(i, j) = v[arrPerm[j]];
                            indexMat[k] = arrPerm[j];
                        }
                        
                        nextFullPerm(arrPerm, lastElem, pentultimate);
                    }
                } else {
                    for (std::size_t i = 0; i < phaseOne; ++i) {
                        for (std::size_t j = 0, k = i; j < uR; ++j, k += phaseOne) {
                            permuteMatrix(i, j) = v[arrPerm[j]];
                            indexMat[k] = arrPerm[j];
                        }
                        
                        nextPartialPerm(arrPerm, uR, lastCol, uN, lastElem);
                    }
                }
                
                unsigned long int start = segment, last = 2 * segment;
                unsigned long int ind = 1;
                typeVector vTemp(1);
                
                for (; start < uRowN; start += segment, last += segment, ++ind) {
                    vTemp[0] = v[0];
                    v[0] = v[ind];
                    v[ind] = vTemp[0];
        
                    if (last > uRowN) {
                        unsigned long int skip = last - uRowN;
                        last = uRowN;
                        
                        for (std::size_t j = 0, k = 0; j < uR; ++j, k += skip)
                            for (std::size_t i = start; i < last; ++i, ++k)
                                permuteMatrix(i, j) = v[indexMat[k]];
                    } else {
                        for (std::size_t j = 0, k = 0; j < uR; ++j)
                            for (std::size_t i = start; i < last; ++i, ++k)
                                permuteMatrix(i, j) = v[indexMat[k]];
                    }
                }
        
                delete[] indexMat;
                delete[] arrPerm;
            } else {
                permuteMatrix(0, 0) = v[0];
            }
        }
        
        return permuteMatrix;
    }
    
    template <typename typeMatrix, typename typeVector>
    typeMatrix MultisetPermutation(int n, int r, typeVector &v, int numRows,
                                   bool xtraCol, std::vector<int> &z) {
        
        unsigned long int numCols = r;
        unsigned long int lenFreqs = z.size();
        
        if (xtraCol)
            ++numCols;
        
        typeMatrix permuteMatrix = Rcpp::no_init_matrix(numRows, numCols);
        uint16_t *arrPerm = new uint16_t[lenFreqs];
        
        unsigned long int uN = n, numR1 = numRows - 1;
        unsigned long int uR = r, lastCol = r - 1;
        unsigned long int lastElem = lenFreqs - 1;
        
        for (std::size_t j = 0; j < lenFreqs; ++j)
            arrPerm[j] = (uint16_t) z[j];
        
        if (uR == lenFreqs) {
            unsigned long int pentultimate = lenFreqs - 2;
            
            for (std::size_t i = 0; i < numR1; ++i) {
                for (std::size_t j = 0; j < uR; ++j)
                    permuteMatrix(i, j) = v[arrPerm[j]];
                
                nextFullPerm(arrPerm, lastElem, pentultimate);
            }
        } else {
            for (std::size_t i = 0; i < numR1; ++i) {
                for (std::size_t j = 0; j < uR; ++j)
                    permuteMatrix(i, j) = v[arrPerm[j]];
                
                nextPartialPerm(arrPerm, uR, lastCol, uN, lastElem);
            }
        }
        
        // Get last permutation
        for (std::size_t j = 0; j < uR; ++j)
            permuteMatrix(numR1, j) = v[arrPerm[j]];
        
        delete[] arrPerm;
        return permuteMatrix;
    }
    
    template <typename typeVector>
    SEXP PermutationApplyFun(int n, int r, typeVector &v, bool repetition,
                              int numRows, bool Multi, std::vector<int> &z,
                              SEXP func, SEXP rho) {
        
        unsigned long int lenFreqs = 0, uR = r, uN = n;
        
        if (Multi)
            lenFreqs = z.size();
        
        SEXP ans = PROTECT(Rf_allocVector(VECSXP, numRows));
        SEXP R_fcall = PROTECT(Rf_lang2(func, R_NilValue));
        typeVector vectorPass(r);
        
        unsigned long int numR1 = numRows - 1;
        unsigned long int uRowN = numRows, lastCol = uR - 1;
        unsigned long int lastElem = (Multi) ? (lenFreqs - 1) : (n - 1);
        
        if (repetition) {
            int lastElemInt = lastElem;
            
            for (std::size_t i = 0; i < uRowN; ++i) {
                for (std::size_t j = 0; j < uR; ++j)
                    vectorPass[j] = v[z[j]];

                SETCADR(R_fcall, vectorPass);
                SET_VECTOR_ELT(ans, i, Rf_eval(R_fcall, rho));

                for (int k = lastCol; k >= 0; --k) {
                    if (z[k] != lastElemInt) {
                        ++z[k];
                        break;
                    } else {
                        z[k] = 0;
                    }
                }
            }
        } else {
            unsigned long int arrLength = lastElem + 1;
            uint16_t *arrPerm = new uint16_t[arrLength];
            
            for (std::size_t i = 0; i < arrLength; ++i)
                arrPerm[i] = (uint16_t) z[i];
            
            if (uR == uN || uR == lenFreqs) {
                unsigned long int pentultimate = lastElem - 1;
                
                for (std::size_t i = 0; i < numR1; ++i) {
                    for (std::size_t j = 0; j < uR; ++j)
                        vectorPass[j] = v[arrPerm[j]];
                    
                    SETCADR(R_fcall, vectorPass);
                    SET_VECTOR_ELT(ans, i, Rf_eval(R_fcall, rho));
                    nextFullPerm(arrPerm, lastElem, pentultimate);
                }
            } else {
                for (std::size_t i = 0; i < numR1; ++i) {
                    for (std::size_t j = 0; j < uR; ++j)
                        vectorPass[j] = v[arrPerm[j]];
                        
                    SETCADR(R_fcall, vectorPass);
                    SET_VECTOR_ELT(ans, i, Rf_eval(R_fcall, rho));
                    nextPartialPerm(arrPerm, uR, lastCol, uN, lastElem);
                }
            }
            
            // Get last permutation
            for (std::size_t j = 0; j < uR; ++j)
                vectorPass[j] = v[arrPerm[j]];

            SETCADR(R_fcall, vectorPass);
            SET_VECTOR_ELT(ans, numR1, Rf_eval(R_fcall, rho));
            delete[] arrPerm;
        }
        
        UNPROTECT(2);
        return ans;
    }
}

#endif
