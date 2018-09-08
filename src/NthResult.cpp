#include <CombPermUtils.h>
#include <CountGmp.h>

std::vector<int> nonZeroVec(std::vector<int> v) {
    std::vector<int> nonZero;
    
    for (std::size_t i = 0; i < v.size(); i++)
        if (v[i] > 0)
            nonZero.push_back(v[i]);
        
    return nonZero;
}

std::vector<int> nthPermutation(int n, int r, double myIndex, bool isRep, bool isMult,
                                std::vector<int> Reps, std::vector<int> freqs, bool isStarter = false) {
    
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
        temp = std::pow((double) n, (double) r);
        
        for (int k = 0; k < r; ++k) {
            temp /= n;
            j = (int) std::trunc(index1 / temp);
            res[k] = j;
            index1 -= (temp * j);
        }
    } else {
        temp = NumPermsNoRep(n, r);
        std::vector<int> indexVec(n);
        std::iota(indexVec.begin(), indexVec.end(), 0);
        
        for (int k = 0; k < r; ++k, --n1) {
            temp /= n1;
            j = (int) std::trunc(index1 / temp);
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
    
    int j = 0, n1 = n, r1 = r - 1, rTemp;
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
            rTemp = n1 - 1;
            while (test <= index1) {
                index2 -= temp;
                temp *= rTemp;
                temp /= (n1 + r1 - 1);
                --n1;
                ++j;
                --rTemp;
                test += temp;
            }
            res[k] = j;
            index1 = index2;
        }
        
    } else {
        
        --n1;
        
        for (std::size_t k = 0; k < uR; ++k, --n1, --r1, ++j) {
            temp = test = nChooseK(n1, r1);
            rTemp = n1 - r1;
            while (test <= index1) {
                index2 -= temp;
                temp *= rTemp;
                temp /= n1;
                --n1;
                ++j;
                --rTemp;
                test += temp;
            }
            res[k] = j;
            index1 = index2;
        }
    }
    
    return res;
}

std::vector<int> nthPermutationGmp(int n, int r, mpz_t myIndex, bool isRep, bool isMult,
                                   std::vector<int> Reps, std::vector<int> freqs, bool isStarter = false) {
    int j = 0, n1 = n;
    mpz_t temp, temp2, index1;
    mpz_init(temp); mpz_init(temp2); mpz_init(index1);
    mpz_set(index1, myIndex);
    std::vector<int> res(r);
    
    if (isMult) {
        
        mpz_add_ui(index1, index1, 1);
        std::vector<int> Counts;
        int r1 = r - 1;
        mpz_t test, index2;
        mpz_init(test); mpz_init(index2);
        mpz_set(index2, index1);
        
        for (int k = 0; k < r; ++k, --r1) {
            
            j = 0;
            while (Reps[j] == 0)
                ++j;
            
            --Reps[j];
            Counts = nonZeroVec(Reps);
            MultisetPermRowNumGmp(temp, (int) Counts.size(), r1, Counts);
            mpz_set(test, temp);
            
            while (mpz_cmp(test, index1) < 0) {
                mpz_sub(index2, index2, temp);
                ++Reps[j];
                ++j;

                while (Reps[j] == 0)
                    ++j;

                --Reps[j];
                
                Counts = nonZeroVec(Reps);
                MultisetPermRowNumGmp(temp, (int) Counts.size(), r1, Counts);
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
            j = mpz_get_si(temp2);
            res[k] = j;
            mpz_submul_ui(index1, temp, j);
        }
    } else {
        mpz_set_ui(temp, 1);
        NumPermsNoRepGmp(temp, n, r);
        std::vector<int> indexVec(n);
        std::iota(indexVec.begin(), indexVec.end(), 0);
        
        for (int k = 0; k < r; ++k, --n1) {
            mpz_divexact_ui(temp, temp, n1);
            mpz_tdiv_q(temp2, index1, temp);
            j = mpz_get_si(temp2);
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
    
    int j = 0, n1 = n, r1 = r - 1, rTemp;
    mpz_t test, temp, index1, index2;
    mpz_init(test); mpz_init(temp);
    mpz_init(index1); mpz_init(index2);
    mpz_set(index1, myIndex);
    mpz_set(index2, myIndex);
    
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
            
            MultisetCombRowNumGmp(temp, n1, r1, Counts);
            mpz_set(test, temp);
            
            while (mpz_cmp(test, index1) <= 0) {
                mpz_sub(index2, index2, temp);
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
                
                MultisetCombRowNumGmp(temp, n1, r1, Counts);
                mpz_add(test, test, temp);
            }
            
            res[k] = j;
            mpz_set(index1, index2);
            
            --Reps[j];
            if (Reps[j] <= 0)
                ++j;
        }
    } else if (isRep) {
        
        for (std::size_t k = 0; k < uR; ++k, --r1) {
            NumCombsWithRepGmp(temp, n1, r1);
            mpz_set(test, temp);
            rTemp = n1 - 1;
            
            while (mpz_cmp(test, index1) <= 0) {
                mpz_sub(index2, index2, temp);
                mpz_mul_ui(temp, temp, rTemp);
                mpz_divexact_ui(temp, temp, n1 + r1 - 1);
                --n1;
                ++j;
                --rTemp;
                mpz_add(test, test, temp);
            }
            res[k] = j;
            mpz_set(index1, index2);
        }
        
    } else {
        
        --n1;
        
        for (std::size_t k = 0; k < uR; ++k, --n1, --r1, ++j) {
            nChooseKGmp(temp, n1, r1);
            mpz_set(test, temp);
            rTemp = n1 - r1;
            
            while (mpz_cmp(test, index1) <= 0) {
                mpz_sub(index2, index2, temp);
                mpz_mul_ui(temp, temp, rTemp);
                mpz_divexact_ui(temp, temp, n1);
                --n1;
                ++j;
                --rTemp;
                mpz_add(test, test, temp);
            }
            
            res[k] = j;
            mpz_set(index1, index2);
        }
    }
    
    mpz_clear(index1); mpz_clear(index2);
    mpz_clear(temp); mpz_clear(test);
    
    return res;
}
