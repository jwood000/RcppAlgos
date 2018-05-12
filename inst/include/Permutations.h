#ifndef RcppAlgos_Permutations_h
#define RcppAlgos_Permutations_h

#include <CombPermUtility.h>

template <typename typeMatrix, typename typeVector>
typeMatrix PermuteGeneral(int n, int r, typeVector v,
                          bool repetition, int numRows, bool xtraCol);

template <typename typeMatrix, typename typeVector>
typeMatrix MultisetPermutation(int n, int r, typeVector v,
                               std::vector<int> Reps,
                               int numRows, bool xtraCol);

SEXP PermutationsRcpp(int n, int m, bool repetition, Rcpp::CharacterVector vStr,
                      int nRows, std::vector<int> vInt, std::vector<double> vNum,
                      bool isMult, bool isFac, bool keepRes, bool isChar, SEXP Rv,
                      bool isInt, std::vector<int> myReps, SEXP f1, SEXP f2);

template <typename typeMatrix, typename typeVector>
typeMatrix PermuteGeneral(int n, int r, typeVector v,
                          bool repetition, int numRows, bool xtraCol) {
    
    unsigned long int uN = n, uR = r, uRowN = numRows;
    unsigned long int chunk, numCols, segment;
    numCols = xtraCol ? (uR + 1) : uR;
    typeMatrix permuteMatrix(uRowN, numCols);
    
    if (repetition) {
        
        unsigned long int groupLen = 1, repLen = 1, myCol;
        typename typeVector::iterator it, vBeg, vEnd;
        vBeg = v.begin(); vEnd = v.end();
        
        for (std::size_t i = uR; i > 0; i--) {
            myCol = i - 1;
            groupLen *= uN;
            chunk = 0;
            bool KeepGoing = true;
            segment = repLen;
            for (std::size_t k = 0; k < uRowN && KeepGoing; k += groupLen) {
                for (it = vBeg; it < vEnd && KeepGoing; 
                        it++, segment += repLen, chunk += repLen) {

                    if (segment > uRowN) {
                        segment = uRowN;
                        KeepGoing = false;
                    }

                    for (std::size_t j = chunk; j < segment; j++)
                        permuteMatrix(j, myCol) = *it;
                }
            }
            repLen *= uN;
        }
        
    } else {
        
        if (n < 3) {
            if (n == 1) {
                permuteMatrix(0, 0) = v[0];
            } else if (r == 2) {
                permuteMatrix(1, 1) = permuteMatrix(0, 0) = v[0];
                permuteMatrix(1, 0) = permuteMatrix(0, 1) = v[1];
            } else {
                permuteMatrix(0, 0) = v[0];
                permuteMatrix(1, 0) = v[1];
            }
            return permuteMatrix;
        }
        
        unsigned long int firstCol = (r == n) ? (r - 1) : r;
        unsigned long int phaseOne, maxN = NumPermsNoRep(n, r);
        
        chunk = segment = maxN / uN;
        
        phaseOne = (uRowN < segment) ? uRowN : segment;
        int *indexMat = new int[phaseOne * r];
        
        std::vector<int> indexVec(n), origSeqeunce(n);
        std::iota(origSeqeunce.begin(), origSeqeunce.end(), 0);
        indexVec = origSeqeunce;
        
        std::vector<char> whichPerm(n, 1);
        unsigned long int resetInd = 0, resetSize = n;
        unsigned long int start, last, colInd, colLim = 1;
        
        for (std::size_t i = 0; i < firstCol; i++) {
            colInd = phaseOne * i;
            start = last = resetInd = 0;

            while (last < phaseOne) {
                resetInd++;

                for (std::size_t j = 0; j < colLim && start < uRowN; j++, start += chunk) {
                    last += chunk;
                    if (last > uRowN)
                        last = uRowN;

                    for (std::size_t k = start; k < last; k++, colInd++) {
                        permuteMatrix(k, i) = v[indexVec[j]];
                        indexMat[colInd] = indexVec[j];
                    }
                }

                if (start < phaseOne) {
                    if (i == 2) {
                        indexVec[resetInd - 1]--;
                    } else if (i > 2) {
                        if (resetInd <= resetSize) {
                            indexVec[resetInd - 1] = indexMat[(i - 1) * phaseOne + start - 1];
                        } else {
                            resetInd = 0;
                            std::fill(whichPerm.begin() + 1, whichPerm.end(), 1);
                            for (std::size_t j = 1, k = phaseOne; j < i; j++, k += phaseOne)
                                whichPerm[indexMat[k + start]] = 0;

                            for (std::size_t j = 1, k = 0; j < uN; j++)
                                if (whichPerm[j])
                                    indexVec[k++] = j;
                        }
                    }
                }
            }

            resetSize--;
            colLim = resetSize;
            chunk /= resetSize;

            origSeqeunce.erase(origSeqeunce.begin());
            indexVec = origSeqeunce;
        }

        if (r == n) {
            unsigned long int colNew, r1 = r - 1, r2 = r - 2;
            colInd = r2 * phaseOne;
            colNew = r1 * phaseOne;

            for (std::size_t i = 0, j = 1; i < phaseOne; i += 2, j += 2) {
                permuteMatrix(i, r1) = v[indexMat[colNew + i] = indexMat[colInd + j]];
                permuteMatrix(j, r1) = v[indexMat[colNew + j] = indexMat[colInd + i]];
            }

            if (uRowN < segment) {
                std::fill(whichPerm.begin() + 1, whichPerm.end(), 1);
                for (std::size_t j = 2, k = 2 * phaseOne; j <= r1; j++, k += phaseOne)
                    whichPerm[indexMat[k - 1]] = 0;

                for (std::size_t j = 1; j < uN; j++)
                    if (whichPerm[j])
                        permuteMatrix(uRowN - 1, r1) = v[j];
            }
        }

        if (phaseOne < uRowN) {
            start = last = segment;
            last += segment;
            unsigned long int ind = 0;
            typeVector vTemp(1);

            for (; start < uRowN; start += segment, last += segment) {
                ind++;
                vTemp[0] = v[0];
                v[0] = v[ind];
                v[ind] = vTemp[0];

                if (last > uRowN)
                    last = uRowN;

                for (std::size_t i = 0, q = 0; i < uR; i++, q += phaseOne)
                    for (std::size_t j = start, k = q; j < last; j++, k++)
                        permuteMatrix(j, i) = v[indexMat[k]];
            }
        }

        delete[] indexMat;
    }
    
    return permuteMatrix;
}

template <typename typeMatrix, typename typeVector>
typeMatrix MultisetPermutation(int n, int r, typeVector v,
                               std::vector<int> Reps,
                               int numRows, bool xtraCol) {
    
    unsigned long int uN = n, count = 0, numCols = 100;
    typeVector vTemp(1);
    // sort v and order Reps by the ordering of v.
    for (std::size_t i = 0; i < (uN - 1); i++) {
        for (std::size_t j = (i+1); j < uN; j++) {
            if (v[i] > v[j]) {
                vTemp[0] = v[i];
                v[i] = v[j];
                v[j] = vTemp[0];
                std::swap(Reps[i], Reps[j]);
            }
        }
    }
    
    unsigned long int sumReps = std::accumulate(Reps.begin(), Reps.end(), 0);
    bool retAllPerms = true;
    typeMatrix myCombs;
    
    if (r < sumReps) {
        numCols = r;
        retAllPerms = false;
    } else {
        numCols = sumReps;
    }
    
    if (xtraCol)
        numCols++;
    
    typeMatrix permuteMatrix(numRows, numCols);
    uint16_t *arrPerm = new uint16_t[sumReps];
    
    for (std::size_t i = 0; i < n; i++)
        for (std::size_t j = 0; j < Reps[i]; j++, count++)
            arrPerm[count] = (uint16_t) i;
    
    int p1, p2, numR1 = numRows - 1;
    int temp, n1 = sumReps - 1;
    
    if (retAllPerms) {
        for (std::size_t i = 0; i < numR1; i++) {
            for (std::size_t j = 0; j < numCols; j++)
                permuteMatrix(i, j) = v[arrPerm[j]];
            
            // This algorithm is nearly identical to the
            // one found in the standard algorithm library
            
            p1 = n1;
            while (arrPerm[p1] <= arrPerm[p1 - 1])
                p1--;
            
            p1--;
            p2 = n1;
            
            while (arrPerm[p2] <= arrPerm[p1])
                p2--;
            
            // swap
            temp = arrPerm[p1];
            arrPerm[p1] = arrPerm[p2];
            arrPerm[p2] = temp;
            
            // reverse
            for (std::size_t k = p1 + 1, q = n1; k < q; k++, q--) {
                temp = arrPerm[k];
                arrPerm[k] = arrPerm[q];
                arrPerm[q] = temp;
            }
        }
        
    } else {
        
        int r1 = r - 1;
        
        for (std::size_t i = 0; i < numR1; i++) {
            for (std::size_t j = 0; j < numCols; j++)
                permuteMatrix(i, j) = v[arrPerm[j]];
            
            // This algorithm is the same as above except that
            // since we are not using the entire vector, we have
            // to check first that the rth element is the largest.
            // If it is, we have to reverse all of the elements
            // to the right of the rth position before finding
            // the next permutation. This is so because if we
            // didn't, all of the next perms. of the entire vector
            // would produce many duplicate r-length perms. If it
            // isn't the largest, we find the element to the right
            // and swap them. We can then proceed to the next perm.
            // We can do this because the standard algo would end
            // up performing two unnecessary reversings.
            
            p1 = numCols;
            while(p1 < n && arrPerm[r1] >= arrPerm[p1])
                p1++;
            
            if (p1 < n) {
                temp = arrPerm[p1];
                arrPerm[p1] = arrPerm[r1];
                arrPerm[r1] = temp;
            } else {
                for (std::size_t k = r, q = n1; k < q; k++, q--) {
                    temp = arrPerm[k];
                    arrPerm[k] = arrPerm[q];
                    arrPerm[q] = temp;
                }
                
                p1 = n1;
                while (arrPerm[p1] <= arrPerm[p1 - 1])
                    p1--;
                
                p1--;
                p2 = n1;
                
                while (arrPerm[p2] <= arrPerm[p1])
                    p2--;
                
                temp = arrPerm[p1];
                arrPerm[p1] = arrPerm[p2];
                arrPerm[p2] = temp;
                
                for (std::size_t k = p1 + 1, q = n1; k < q; k++, q--) {
                    temp = arrPerm[k];
                    arrPerm[k] = arrPerm[q];
                    arrPerm[q] = temp;
                }
            }
        }
    }
    
    // Get last permutation
    for (std::size_t j = 0; j < numCols; j++)
        permuteMatrix(numR1, j) = v[arrPerm[j]];
    
    delete[] arrPerm;
    return permuteMatrix;
}

#endif
