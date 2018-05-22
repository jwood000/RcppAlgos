#ifndef RcppAlgos_Combinations_h
#define RcppAlgos_Combinations_h

#include <CombPermUtils.h>

namespace Combinations {

    template <typename typeRcpp>
    typeRcpp SubMat(typeRcpp m, int n) {
        int k = m.ncol();
        typeRcpp subMatrix(n,k);
        
        for (int i = 0; i < n; i++)
            subMatrix(i, Rcpp::_) = m(i, Rcpp::_);
        
        return subMatrix;
    }
    
    template <typename typeMatrix, typename typeVector>
    typeMatrix ComboGeneral(int n, int r, typeVector v, bool repetition, 
                            int numRows, bool xtraCol, std::vector<int> z) {
        
        int r1 = r - 1, r2 = r - 2;
        int k, i, numIter, vSize;
        int numCols, maxZ, count = 0;
        bool needsSubsetting = false;
        numCols = xtraCol ? (r + 1) : r;
        typeMatrix combinationMatrix(numRows, numCols);
        
        if (repetition) {
            v.erase(std::unique(v.begin(), v.end()), v.end());
            vSize = v.size();
            int testRows = NumCombsWithRep(vSize, r);
            
            if (testRows < numRows) {
                needsSubsetting = true;
                numRows = testRows;
            }
            
            maxZ = vSize - 1;
            
            while (count < numRows) {
                numIter = vSize - z[r1];
                if (numIter + count > numRows)
                    numIter = numRows - count;
                
                for (i = 0; i < numIter; i++, count++, z[r1]++)
                    for (k = 0; k < r; k++)
                        combinationMatrix(count, k) = v[z[k]];
                
                for (i = r2; i >= 0; i--) {
                    if (z[i] != maxZ) {
                        z[i]++;
                        for (k = (i + 1); k < r; k++)
                            z[k] = z[k - 1];
                        
                        break;
                    }
                }
            }
        } else {
            while (count < numRows) {
                numIter = n - z[r1];
                if ((numIter + count) > numRows)
                    numIter = numRows - count;
                
                for (i = 0; i < numIter; i++, count++, z[r1]++)
                    for (k = 0; k < r; k++)
                        combinationMatrix(count, k) = v[z[k]];
                
                for (i = r2; i >= 0; i--) {
                    if (z[i] != (n - r + i)) {
                        z[i]++;
                        for (k = (i + 1); k < r; k++) 
                            z[k] = z[k - 1] + 1;
                        
                        break;
                    }
                }
            }
        }
        
        if (needsSubsetting) {
            return SubMat(combinationMatrix, count);
        } else {
            return combinationMatrix;
        }
    }
    
    template <typename typeMatrix, typename typeVector>
    typeMatrix MultisetCombination(int n, int r, typeVector v,
                                   std::vector<int> Reps,
                                   int numRows, bool xtraCol,
                                   std::vector<int> z) {
        
        int i, j, numIter, numCols;
        int count = 0, zExpSize = 0;
        
        std::vector<int> zExpand, zIndex, zGroup(r);
        int r1 = r - 1, r2 = r - 2, k = 0;
        
        for (i = 0; i < n; i++)
            zExpSize += Reps[i];
        
        zIndex.reserve(n);
        zExpand.reserve(zExpSize);
        for (i = 0; i < n; i++) {
            zIndex.push_back(k);
            for (j = 0; j < Reps[i]; j++, k++)
                zExpand.push_back(i);
        }
        
        numCols = xtraCol ? (r + 1) : r;
        typeMatrix combinationMatrix(numRows, numCols);
        
        while (count < numRows) {
            numIter = n - z[r1];
            
            if (numIter + count > numRows)
                numIter = numRows - count;
            
            for (i = 0; i < numIter; i++, count++, z[r1]++)
                for (k = 0; k < r; k++)
                    combinationMatrix(count, k) = v[zExpand[zIndex[z[k]]]];
            
            for (i = r2; i >= 0; i--) {
                if (zExpand[zIndex[z[i]]] != zExpand[zExpSize - r + i]) {
                    z[i]++;
                    zGroup[i] = zIndex[z[i]];
                    for (k = (i + 1); k < r; k++) {
                        zGroup[k] = zGroup[k - 1] + 1;
                        z[k] = zExpand[zGroup[k]];
                    }
                    break;
                }
            }
        }
        
        return combinationMatrix;
    }
}

#endif
