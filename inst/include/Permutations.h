#ifndef PERMUTATIONS_H
#define PERMUTATIONS_H

#include <CombPermUtils.h>

namespace Permutations {
    
    template <typename typeMatrix, typename typeVector>
    typeMatrix PermuteGeneral(int n, int r, typeVector v, bool repetition,
                              int numRows, bool xtraCol, std::vector<int> z, bool nonTrivial) {
        
        unsigned long int uN = n, uR = r, uRowN = numRows;
        unsigned long int chunk, numCols, segment;
        numCols = xtraCol ? (uR + 1) : uR;
        typeMatrix permuteMatrix = Rcpp::no_init_matrix(uRowN, numCols);
        
        if (repetition) {
            int r1 = r - 1, n1 = n - 1;
            
            for (std::size_t i = 0; i < uRowN; ++i) {
                for (std::size_t j = 0; j < uR; ++j)
                    permuteMatrix(i, j) = v[z[j]];
                
                for (int k = r1; k >= 0; --k) {
                    if (z[k] != n1) {
                        ++z[k];
                        break;
                    } else {
                        z[k] = 0;
                    }
                }
            }
            
        } else if (nonTrivial) {
            
            uint16_t *arrPerm = new uint16_t[uN];

            for (std::size_t i = 0; i < uN; ++i)
                arrPerm[i] = (uint16_t) z[i];
            
            unsigned long int numR1 = numRows - 1;
            unsigned long int r1 = uR - 1, n1 = uN - 1;
            
            if (r == n) {
                for (std::size_t i = 0; i < numR1; ++i) {
                    for (std::size_t j = 0; j < uR; ++j)
                        permuteMatrix(i, j) = v[arrPerm[j]];
                    
                    nextFullPerm(arrPerm, n1);
                }
            } else {
                for (std::size_t i = 0; i < numR1; ++i) {
                    for (std::size_t j = 0; j < uR; ++j)
                        permuteMatrix(i, j) = v[arrPerm[j]];
                    
                    nextPartialPerm(arrPerm, uR, r1, uR, n1, uN);
                }
            }
            
            // Get last permutation
            for (std::size_t j = 0; j < uR; ++j)
                permuteMatrix(numR1, j) = v[arrPerm[j]];
            
            delete[] arrPerm;
            return permuteMatrix;
            
        } else {
            
            if (r == 1) {
                for (std::size_t i = 0; i < uN; i++)
                    permuteMatrix(i, 0) = v[i];
                return permuteMatrix;
            }
            
            unsigned long int phaseOne, maxN = NumPermsNoRep(n, r);
            chunk = segment = maxN / uN;
            phaseOne = (uRowN < segment) ? uRowN : segment;
            uint16_t *indexMat = new uint16_t[phaseOne * uR];
            unsigned long int start, last, colInd = 0;
            
            if (r < n) {
                std::vector<uint16_t> indexVec(uN), origSeqeunce(uN);
                std::iota(origSeqeunce.begin(), origSeqeunce.end(), 0);
                indexVec = origSeqeunce;
                
                std::vector<char> whichPerm(uN, 1);
                int resetInd = 0, resetSize = uN;
                unsigned long int colLim = 1;
                
                for (std::size_t i = 0; i < uR; ++i) {
                    colInd = i;
                    start = last = resetInd = 0;
                    
                    while (last < phaseOne) {
                        ++resetInd;
                        
                        for (std::size_t j = 0; (j < colLim) && (start < uRowN); ++j, start += chunk) {
                            last += chunk;
                            if (last > uRowN)
                                last = uRowN;
                            
                            for (std::size_t k = start; k < last; ++k, colInd += uR) {
                                permuteMatrix(k, i) = v[indexVec[j]];
                                indexMat[colInd] = indexVec[j];
                            }
                        }
                        
                        if (start < phaseOne) {
                            if (i == 2) {
                                --indexVec[resetInd - 1];
                            } else if (i > 2) {
                                if (resetInd <= resetSize) {
                                    indexVec[resetInd - 1] = indexMat[(i - 1) + (start - 1) * uR];
                                } else {
                                    resetInd = 0;
                                    std::fill(whichPerm.begin() + 1, whichPerm.end(), 1);
                                    for (std::size_t j = 1, k = 1 + (start * uR); j < i; ++j, ++k)
                                        whichPerm[indexMat[k]] = 0;
                                    
                                    for (std::size_t j = 1, k = 0; j < uN; ++j)
                                        if (whichPerm[j])
                                            indexVec[k++] = j;
                                }
                            }
                        }
                    }
                    
                    --resetSize;
                    colLim = resetSize;
                    chunk /= resetSize;
                    
                    origSeqeunce.erase(origSeqeunce.begin());
                    indexVec = origSeqeunce;
                }
    
            } else {
                
                uint16_t *arrPerm = new uint16_t[uN];
                
                for (std::size_t i = 0; i < uN; ++i)
                    arrPerm[i] = (uint16_t) i;
                
                unsigned long int n1 = n - 1;
                
                for (std::size_t i = 0; i < phaseOne; ++i) {
                    for (std::size_t j = 0; j < uN; ++j, ++colInd) {
                        permuteMatrix(i, j) = v[arrPerm[j]];
                        indexMat[colInd] = arrPerm[j];
                    }
                    
                    nextFullPerm(arrPerm, n1);
                }
                
                delete[] arrPerm;
            }
            
            start = last = segment;
            last += segment;
            unsigned long int ind = 0;
            typeVector vTemp(1);
            
            for (; start < uRowN; start += segment, last += segment) {
                ++ind;
                vTemp[0] = v[0];
                v[0] = v[ind];
                v[ind] = vTemp[0];
    
                if (last > uRowN)
                    last = uRowN;

                for (std::size_t i = start, k = 0; i < last; ++i)
                    for (std::size_t j = 0; j < uR; ++j, ++k)
                        permuteMatrix(i, j) = v[indexMat[k]];
            }
    
            delete[] indexMat;
        }
        
        return permuteMatrix;
    }
    
    template <typename typeMatrix, typename typeVector>
    typeMatrix MultisetPermutation(int n, int r, typeVector v,
                                   std::vector<int> Reps,
                                   int numRows, bool xtraCol,
                                   std::vector<int> z) {
        
        unsigned long int numCols, sumReps;
        sumReps = std::accumulate(Reps.begin(), Reps.end(), 0);
        
        bool retAllPerms = true;
        
        if (r < (int) sumReps) {
            numCols = r;
            retAllPerms = false;
        } else {
            numCols = sumReps;
        }
        
        unsigned long int arrLength = numCols;
        
        if (xtraCol)
            ++numCols;
        
        typeMatrix permuteMatrix = Rcpp::no_init_matrix(numRows, numCols);
        uint16_t *arrPerm = new uint16_t[sumReps];
        
        for (std::size_t j = 0; j < sumReps; ++j)
            arrPerm[j] = (uint16_t) z[j];
        
        unsigned long int numR1 = numRows - 1;
        unsigned long int r1 = r - 1;
        unsigned long int n1 = sumReps - 1;
        
        if (retAllPerms) {
            for (std::size_t i = 0; i < numR1; ++i) {
                for (std::size_t j = 0; j < arrLength; ++j)
                    permuteMatrix(i, j) = v[arrPerm[j]];
                
                nextFullPerm(arrPerm, n1);
            }
            
        } else {
            for (std::size_t i = 0; i < numR1; ++i) {
                for (std::size_t j = 0; j < arrLength; ++j)
                    permuteMatrix(i, j) = v[arrPerm[j]];
                
                nextPartialPerm(arrPerm, arrLength, r1, r, n1, n);
            }
        }
        
        // Get last permutation
        for (std::size_t j = 0; j < arrLength; ++j)
            permuteMatrix(numR1, j) = v[arrPerm[j]];
        
        delete[] arrPerm;
        return permuteMatrix;
    }
}

#endif
