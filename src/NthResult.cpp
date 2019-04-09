#include "CombPermUtils.h"
#include "CountGmp.h"

std::vector<int> nonZeroVec(std::vector<int> v) {
    std::vector<int> nonZero;
    
    for (std::size_t i = 0; i < v.size(); i++)
        if (v[i] > 0)
            nonZero.push_back(v[i]);
        
    return nonZero;
}

std::vector<int> nthPermutation(int n, int r, double myIndex, bool isRep, 
                                bool isMult, std::vector<int> Reps, std::vector<int> freqs,
                                bool isStarter = false) {
    
    double temp, index1 = myIndex;
    std::vector<int> res(r); 
    
    if (isMult) {
        
        ++index1;
        std::vector<int> Counts;
        double test, index2 = index1;
        
        for (int k = 0, r1 = r - 1; k < r; ++k, --r1) {
            int j = 0;
            while (Reps[j] == 0)
                ++j;
            
            --Reps[j];
            Counts = nonZeroVec(Reps);
            test = temp = MultisetPermRowNum(Counts.size(), r1, Counts);
            
            for (; test < index1; test += temp) {
                index2 -= temp;
                ++Reps[j];
                ++j;
                
                while (Reps[j] == 0)
                    ++j;
                
                --Reps[j];
                Counts = nonZeroVec(Reps);
                temp = MultisetPermRowNum(Counts.size(), r1, Counts);
            }
            
            res[k] = j;
            index1 = index2;
        }
        
        if (isStarter) {
            for (std::size_t j = 0; j < res.size(); ++j) {
                for (std::size_t i = 0; i < freqs.size(); ++i) {
                    if (freqs[i] == res[j]) {
                        freqs.erase(freqs.begin() + i);
                        break;
                    }
                }
            }
            
            for (std::size_t i = 0; i < freqs.size(); ++i)
                res.push_back(freqs[i]);
        }
    } else if (isRep) {
        temp = std::pow(static_cast<double>(n), static_cast<double>(r));
        
        for (int k = 0; k < r; ++k) {
            temp /= n;
            int j = static_cast<int>(index1 / temp);
            res[k] = j;
            index1 -= (temp * j);
        }
        
    } else {
        temp = NumPermsNoRep(n, r);
        std::vector<int> indexVec(n);
        std::iota(indexVec.begin(), indexVec.end(), 0);
        
        for (int k = 0, n1 = n; k < r; ++k, --n1) {
            temp /= n1;
            int j = static_cast<int>(index1 / temp);
            res[k] = indexVec[j];
            index1 -= (temp * j);
            indexVec.erase(indexVec.begin() + j);
        }
        
        if (isStarter) {
            if (r < n) {
                for (int i = 0; i < n; ++i) {
                    bool bExist = false;
                    for (std::size_t j = 0; j < res.size(); ++j) {
                        if (res[j] == i) {
                            bExist = true;
                            break;
                        }
                    }
                    if (!bExist)
                        res.push_back(i);
                }
            }
        }
    }
    
    return res;
}

std::vector<int> nthCombination(int n, int r, double myIndex, bool isRep, 
                                bool isMult, std::vector<int> Reps) {
    
    double test, temp;
    double index1 = myIndex, index2 = myIndex;
    std::vector<int> res(r);
    
    if (isMult) {
        std::vector<int> Counts = Reps;
        
        for (int k = 0, j = 0, n1 = n, r1 = r - 1; k < r; ++k, --r1) {
            
            --Counts[0];
            if (Counts[0] == 0 && Counts.size() > 1) {
                --n1;
                Counts.erase(Counts.begin());
            }
            
            test = temp = MultisetCombRowNumFast(n1, r1, Counts);
            
            for (; test <= index1; ++j, test += temp) {
                index2 -= temp;
                Reps[j] = 0;
                
                if (static_cast<int>(Counts.size()) == (n - j)) {
                    --n1;
                    Counts.erase(Counts.begin());
                }
                
                --Counts[0];
                if (Counts[0] == 0 && Counts.size() > 1) {
                    --n1;
                    Counts.erase(Counts.begin());
                }
                
                temp = MultisetCombRowNumFast(n1, r1, Counts);
            }
            
            res[k] = j;
            index1 = index2;
            
            --Reps[j];
            if (Reps[j] <= 0) ++j;
        }
    } else if (isRep) {
        
        temp = NumCombsWithRep(n, r - 1);
        
        for (int k = 0, j = 0, n1 = n, r1 = r - 1; k < r; ++k, --r1) {
            test = temp;
            
            for (; test <= index1; --n1, ++j, test += temp) {
                index2 -= temp;
                temp *= (n1 - 1);
                temp /= (n1 + r1 - 1);
            }
            
            temp *= r1;
            temp /= (n1 + r1 - 1);
            res[k] = j;
            index1 = index2;
        }
        
    } else {
        
        temp = nChooseK(n - 1, r - 1);
        
        for (int k = 0, j = 0, n1 = n - 1, r1 = r - 1; 
                                k < r; ++k, --n1, --r1, ++j) {
            test = temp;
            
            for (int rTemp = n1 - r1; test <= index1; 
                        --n1, ++j, --rTemp, test += temp) {
                index2 -= temp;
                temp *= rTemp;
                temp /= n1;
            }
            
            temp *= r1;
            temp /= n1;
            res[k] = j;
            index1 = index2;
        }
    }
    
    return res;
}

std::vector<int> nthPermutationGmp(int n, int r, mpz_t myIndex, bool isRep, bool isMult,
                                   std::vector<int> Reps, std::vector<int> freqs,
                                   bool isStarter = false) {
    
    mpz_t temp, temp2, index1;
    mpz_init(temp); mpz_init(temp2); mpz_init(index1);
    mpz_set(index1, myIndex);
    std::vector<int> res(r);
    
    if (isMult) {
        
        mpz_add_ui(index1, index1, 1);
        std::vector<int> Counts;
        mpz_t test, index2;
        mpz_init(test); mpz_init(index2);
        mpz_set(index2, index1);
        
        for (int k = 0, r1 = r - 1; k < r; ++k, --r1) {
            
            int j = 0;
            while (Reps[j] == 0)
                ++j;
            
            --Reps[j];
            Counts = nonZeroVec(Reps);
            MultisetPermRowNumGmp(temp, static_cast<int>(Counts.size()), r1, Counts);
            mpz_set(test, temp);
            
            while (mpz_cmp(test, index1) < 0) {
                mpz_sub(index2, index2, temp);
                ++Reps[j];
                ++j;

                while (Reps[j] == 0)
                    ++j;

                --Reps[j];
                Counts = nonZeroVec(Reps);
                MultisetPermRowNumGmp(temp, static_cast<int>(Counts.size()), r1, Counts);
                mpz_add(test, test, temp);
            }
            
            res[k] = j;
            mpz_set(index1, index2);
        }
        
        mpz_clear(test); mpz_clear(index2);
        
        if (isStarter) {
            for (std::size_t j = 0; j < res.size(); ++j) {
                for (std::size_t i = 0; i < freqs.size(); ++i) {
                    if (freqs[i] == res[j]) {
                        freqs.erase(freqs.begin() + i);
                        break;
                    }
                }
            }
            
            for (std::size_t i = 0; i < freqs.size(); ++i)
                res.push_back(freqs[i]);
        }
        
    } else if (isRep) {
        mpz_ui_pow_ui(temp, n, r);
        
        for (int k = 0; k < r; ++k) {
            mpz_divexact_ui(temp, temp, n);
            mpz_tdiv_q(temp2, index1, temp);
            int j = mpz_get_si(temp2);
            res[k] = j;
            mpz_submul_ui(index1, temp, j);
        }
        
    } else {
        mpz_set_ui(temp, 1u);
        NumPermsNoRepGmp(temp, n, r);
        std::vector<int> indexVec(n);
        std::iota(indexVec.begin(), indexVec.end(), 0);
        
        for (int k = 0, n1 = n; k < r; ++k, --n1) {
            mpz_divexact_ui(temp, temp, n1);
            mpz_tdiv_q(temp2, index1, temp);
            int j = mpz_get_si(temp2);
            res[k] = indexVec[j];
            mpz_submul_ui(index1, temp, j);
            indexVec.erase(indexVec.begin() + j);
        }
        
        if (isStarter) {
            if (r < n) {
                for (int i = 0; i < n; ++i) {
                    bool bExist = false;
                    for (std::size_t j = 0; j < res.size(); ++j) {
                        if (res[j] == i) {
                            bExist = true;
                            break;
                        }
                    }
                    if (!bExist)
                        res.push_back(i);
                }
            }
        }
    }
    
    mpz_clear(temp); mpz_clear(temp2); mpz_clear(index1);
    return res;
}

std::vector<int> nthCombinationGmp(int n, int r, mpz_t myIndex, bool isRep,
                                   bool isMult, std::vector<int> Reps) {
    ;
    mpz_t test, temp, index1, index2;
    mpz_init(test); mpz_init(temp);
    mpz_init(index1); mpz_init(index2);
    mpz_set(index1, myIndex);
    mpz_set(index2, myIndex);
    std::vector<int> res(r);
    
    if (isMult) {
        std::vector<int> Counts = Reps;
        
        for (int k = 0, n1 = n, j = 0, r1 = r - 1; k < r; ++k, --r1) {
            
            --Counts[0];
            if (Counts[0] == 0 && Counts.size() > 1) {
                --n1;
                Counts.erase(Counts.begin());
            }
            
            MultisetCombRowNumGmp(temp, n1, r1, Counts);
            mpz_set(test, temp);
            
            for (; mpz_cmp(test, index1) <= 0; ++j) {
                mpz_sub(index2, index2, temp);
                Reps[j] = 0;
                
                if (static_cast<int>(Counts.size()) == (n - j)) {
                    --n1;
                    Counts.erase(Counts.begin());
                }
                
                --Counts[0];
                if (Counts[0] == 0 && Counts.size() > 1) {
                    --n1;
                    Counts.erase(Counts.begin());
                }
                
                MultisetCombRowNumGmp(temp, n1, r1, Counts);
                mpz_add(test, test, temp);
            }
            
            res[k] = j;
            mpz_set(index1, index2);
            
            --Reps[j];
            if (Reps[j] <= 0) ++j;
        }
        
    } else if (isRep) {
        
        NumCombsWithRepGmp(temp, n, r - 1);
        
        for (int k = 0, j = 0, n1 = n, r1 = r - 1; k < r; ++k, --r1) {
            mpz_set(test, temp);
            
            for (; mpz_cmp(test, index1) <= 0; ++j, --n1) {
                mpz_sub(index2, index2, temp);
                mpz_mul_ui(temp, temp, n1 - 1);
                mpz_divexact_ui(temp, temp, n1 + r1 - 1);
                mpz_add(test, test, temp);
            }
            
            mpz_mul_ui(temp, temp, r1);
            if ((n1 + r1) > 1) mpz_divexact_ui(temp, temp, n1 + r1 - 1);
            res[k] = j;
            mpz_set(index1, index2);
        }
        
    } else {
        
        nChooseKGmp(temp, n - 1, r - 1);
        
        for (int k = 0, j = 0, n1 = n - 1, r1 = r - 1; 
                            k < r; ++k, --n1, --r1, ++j) {
            mpz_set(test, temp);
            
            for (int rTemp = n1 - r1; 
                 mpz_cmp(test, index1) <= 0; --rTemp, ++j, --n1) {
                
                mpz_sub(index2, index2, temp);
                mpz_mul_ui(temp, temp, rTemp);
                mpz_divexact_ui(temp, temp, n1);
                mpz_add(test, test, temp);
            }
            
            mpz_mul_ui(temp, temp, r1);
            if (n1 > 0) mpz_divexact_ui(temp, temp, n1);
            res[k] = j;
            mpz_set(index1, index2);
        }
    }
    
    mpz_clear(index1); mpz_clear(index2);
    mpz_clear(temp); mpz_clear(test);
    return res;
}
