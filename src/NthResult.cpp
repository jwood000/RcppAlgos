#include <CombPermUtils.h>

std::vector<int> nonZeroVec(std::vector<int> &v) {
    std::vector<int> nonZero;
    
    for (std::size_t i = 0; i < v.size(); i++)
        if (v[i] > 0)
            nonZero.push_back(v[i]);
        
    return nonZero;
}

std::vector<int> nthPermutation(int n, int r, double myIndex, bool isRep,
                                bool isMult, std::vector<int> Reps) {
    
    int j = 0, n1 = n;
    double temp, index1 = myIndex;
    std::vector<int> res(r);
    
    if (isMult) {
        
        ++index1;
        std::vector<int> Counts;
        int r1 = r - 1;
        double test, index2 = index1;
        
        for (int k = 0; k < r; ++k, --r1) {
            
            j = 0;
            while (Reps[j] == 0)
                ++j;
            
            --Reps[j];
            Counts = nonZeroVec(Reps);
            test = temp = MultisetPermRowNum(Counts.size(), r1, Counts);
            
            while (test < index1) {
                index2 -= temp;
                ++Reps[j];
                ++j;
                
                while (Reps[j] == 0)
                    ++j;
                
                --Reps[j];
                
                Counts = nonZeroVec(Reps);
                temp = MultisetPermRowNum(Counts.size(), r1, Counts);
                test += temp;
            }
            
            res[k] = j;
            index1 = index2;
        }
    } else if (isRep) {
        temp = std::pow((double) n, (double) r);
        
        for (int k = 0; k < r; ++k) {
            temp /= n;
            j = (int) std::trunc(index1 / temp);
            res[k] = j;
            index1 -= (temp * (double) j);
        }
    } else {
        temp = NumPermsNoRep(n, r);
        std::vector<int> indexVec(n);
        std::iota(indexVec.begin(), indexVec.end(), 0);
        
        for (int k = 0; k < r; ++k, --n1) {
            temp /= n1;
            j = (int) std::trunc(index1 / temp);
            res[k] = indexVec[j];
            index1 -= (temp * (double) j);
            indexVec.erase(indexVec.begin() + j);
        }
    }
    
    return res;
}

std::vector<int> nthCombination(int n, int r, double myIndex, bool isRep,
                                bool isMult, std::vector<int> Reps) {
    
    int j = 0, n1 = n, r1 = r - 1;
    double test, temp, index1 = myIndex, index2 = myIndex;
    std::vector<int> res(r);
    unsigned long int uR = r;
    
    if (isMult) {
        
        std::vector<int> Counts = Reps;
        
        for (int k = 0; k < r; ++k, --r1) {
            
            --Counts[0];
            if (Counts[0] == 0 && Counts.size() > 1) {
                --n1;
                Counts.erase(Counts.begin());
            }
            
            test = temp = MultisetCombRowNum(n1, r1, Counts);
            
            while (test <= index1) {
                index2 -= temp;
                Reps[j] = 0;
                
                if ((int) Counts.size() == (n - j)) {
                    --n1;
                    Counts.erase(Counts.begin());
                }
                
                --Counts[0];
                if (Counts[0] == 0 && Counts.size() > 1) {
                    --n1;
                    Counts.erase(Counts.begin());
                }
                
                ++j;
                temp = MultisetCombRowNum(n1, r1, Counts);
                test += temp;
            }
            
            res[k] = j;
            index1 = index2;
            
            --Reps[j];
            if (Reps[j] <= 0)
                ++j;
        }
    } else if (isRep) {
        
        for (std::size_t k = 0; k < uR; ++k, --r1) {
            temp = test = NumCombsWithRep(n1, r1);
            while (test <= index1) {
                index2 -= temp;
                --n1;
                ++j;
                temp = NumCombsWithRep(n1, r1);
                test += temp;
            }
            res[k] = j;
            index1 = index2;
        }
        
    } else {
        
        --n1;
        
        for (std::size_t k = 0; k < uR; ++k, --n1, --r1, ++j) {
            temp = test = nChooseK(n1, r1);
            while (test <= index1) {
                index2 -= temp;
                --n1;
                ++j;
                temp = nChooseK(n1, r1);
                test += temp;
            }
            res[k] = j;
            index1 = index2;
        }
    }
    
    return res;
}
