#include "GmpDependUtils.h"

using nthResultPtr = std::vector<int> (*const)(int n, int r, double dblIdx,
                                      mpz_t mpzIdx, const std::vector<int> &Reps);

std::vector<int> nonZeroVec(std::vector<int> v) {
    std::vector<int> nonZero;
    
    for (std::size_t i = 0; i < v.size(); i++)
        if (v[i] > 0)
            nonZero.push_back(v[i]);
        
    return nonZero;
}

std::vector<int> nthComb(int n, int r, double dblIdx, 
                         mpz_t mpzIdx, const std::vector<int> &Reps) {
    
    double index1 = dblIdx, index2 = dblIdx;
    std::vector<int> res(r);
    double temp = nChooseK(n - 1, r - 1);
    
    for (int k = 0, j = 0, n1 = n - 1, 
         r1 = r - 1; k < r; ++k, --n1, --r1, ++j) {
        double test = temp;
        
        for (int s = n1 - r1; test <= index1; 
                    --n1, ++j, --s, test += temp) {
            index2 -= temp;
            temp *= s;
            temp /= n1;
        }
        
        temp *= r1;
        temp /= n1;
        res[k] = j;
        index1 = index2;
    }
    
    return res;
}

std::vector<int> nthCombRep(int n, int r, double dblIdx, 
                            mpz_t mpzIdx, const std::vector<int> &Reps) {
    
    double index1 = dblIdx, index2 = dblIdx;
    std::vector<int> res(r);
    double temp = NumCombsWithRep(n, r - 1);
    
    for (int k = 0, j = 0, n1 = n, r1 = r - 1; k < r; ++k, --r1) {
        double test = temp;
        
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
        
    return res;
}

std::vector<int> nthCombMult(int n, int r, double dblIdx, 
                             mpz_t mpzIdx, const std::vector<int> &Reps) {
    
    double index1 = dblIdx, index2 = dblIdx;
    std::vector<int> res(r);
    std::vector<int> Counts = Reps;
    std::vector<int> TempReps = Reps;
    
    for (int k = 0, j = 0, n1 = n, r1 = r - 1; k < r; ++k, --r1) {
        
        --Counts[0];
        if (Counts[0] == 0 && Counts.size() > 1) {
            --n1;
            Counts.erase(Counts.begin());
        }
        
        double test = MultisetCombRowNumFast(n1, r1, Counts);
        double temp = test;
        
        for (; test <= index1; ++j, test += temp) {
            index2 -= temp;
            TempReps[j] = 0;
            
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
        
        --TempReps[j];
        if (TempReps[j] <= 0) ++j;
    }
    
    return res;
}

std::vector<int> nthCombGmp(int n, int r, double dblIdx, 
                            mpz_t mpzIdx, const std::vector<int> &Reps) {
    
    mpz_t test, temp, index1, index2;
    mpz_init(test); mpz_init(temp);
    mpz_init(index1); mpz_init(index2);
    mpz_set(index1, mpzIdx);
    mpz_set(index2, mpzIdx);
    
    std::vector<int> res(r);
    nChooseKGmp(temp, n - 1, r - 1);
    
    for (int k = 0, j = 0, n1 = n - 1, r1 = r - 1; 
         k < r; ++k, --n1, --r1, ++j) {
        mpz_set(test, temp);
        
        for (int s = n1 - r1; mpz_cmp(test, index1) <= 0; --s, ++j, --n1) {
            mpz_sub(index2, index2, temp);
            mpz_mul_ui(temp, temp, s);
            mpz_divexact_ui(temp, temp, n1);
            mpz_add(test, test, temp);
        }
        
        mpz_mul_ui(temp, temp, r1);
        if (n1 > 0) mpz_divexact_ui(temp, temp, n1);
        res[k] = j;
        mpz_set(index1, index2);
    }
    
    mpz_clear(index1); mpz_clear(index2);
    mpz_clear(temp); mpz_clear(test);
    return res;
}

std::vector<int> nthCombRepGmp(int n, int r, double dblIdx, 
                               mpz_t mpzIdx, const std::vector<int> &Reps) {
    
    mpz_t test, temp, index1, index2;
    mpz_init(test); mpz_init(temp);
    mpz_init(index1); mpz_init(index2);
    mpz_set(index1, mpzIdx);
    mpz_set(index2, mpzIdx);
    
    std::vector<int> res(r);
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
    
    mpz_clear(index1); mpz_clear(index2);
    mpz_clear(temp); mpz_clear(test);
    return res;
}

std::vector<int> nthCombMultGmp(int n, int r, double dblIdx, 
                                mpz_t mpzIdx, const std::vector<int> &Reps) {
    
    mpz_t test, temp, index1, index2;
    mpz_init(test); mpz_init(temp);
    mpz_init(index1); mpz_init(index2);
    mpz_set(index1, mpzIdx);
    mpz_set(index2, mpzIdx);
    
    std::vector<int> res(r);
    std::vector<int> Counts = Reps;
    std::vector<int> TempReps = Reps;
    
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
            TempReps[j] = 0;
            
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
        
        --TempReps[j];
        if (TempReps[j] <= 0) ++j;
    }
    
    mpz_clear(index1); mpz_clear(index2);
    mpz_clear(temp); mpz_clear(test);
    return res;
}

std::vector<int> nthPerm(int n, int r, double dblIdx, 
                         mpz_t mpzIdx, const std::vector<int> &Reps) {
    
    double index1 = dblIdx;
    std::vector<int> res(r); 
    double temp = NumPermsNoRep(n, r);
    std::vector<int> indexVec(n);
    std::iota(indexVec.begin(), indexVec.end(), 0);
    
    for (int k = 0, n1 = n; k < r; ++k, --n1) {
        temp /= n1;
        int j = static_cast<int>(index1 / temp);
        res[k] = indexVec[j];
        index1 -= (temp * j);
        indexVec.erase(indexVec.begin() + j);
    }
    
    return res;
}

std::vector<int> nthPermRep(int n, int r, double dblIdx, 
                            mpz_t mpzIdx, const std::vector<int> &Reps) {
    
    double index1 = dblIdx;
    std::vector<int> res(r); 
    double temp = std::pow(static_cast<double>(n), static_cast<double>(r));
    
    for (int k = 0; k < r; ++k) {
        temp /= n;
        int j = static_cast<int>(index1 / temp);
        res[k] = j;
        index1 -= (temp * j);
    }
    
    return res;
}

std::vector<int> nthPermMult(int n, int r, double dblIdx, 
                             mpz_t mpzIdx, const std::vector<int> &Reps) {
    
    double index1 = dblIdx + 1;
    double index2 = index1;
    
    std::vector<int> res(r);
    std::vector<int> Counts;
    std::vector<int> TempReps = Reps;
    
    for (int k = 0, r1 = r - 1; k < r; ++k, --r1) {
        
        int j = 0;
        while (TempReps[j] == 0)
            ++j;
        
        --TempReps[j];
        Counts = nonZeroVec(TempReps);
        double test = MultisetPermRowNum(Counts.size(), r1, Counts);
        double temp = test;
        
        for (; test < index1; test += temp) {
            index2 -= temp;
            ++TempReps[j];
            ++j;
            
            while (TempReps[j] == 0)
                ++j;
            
            --TempReps[j];
            Counts = nonZeroVec(TempReps);
            temp = MultisetPermRowNum(Counts.size(), r1, Counts);
        }
        
        res[k] = j;
        index1 = index2;
    }
    
    return res;
}

std::vector<int> nthPermGmp(int n, int r, double dblIdx, 
                            mpz_t mpzIdx, const std::vector<int> &Reps) {
    
    mpz_t temp, temp2, index1;
    mpz_init(temp); mpz_init(temp2); mpz_init(index1);
    mpz_set(index1, mpzIdx);
    std::vector<int> res(r);
    
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
    
    mpz_clear(temp); mpz_clear(temp2); mpz_clear(index1);
    return res;
}

std::vector<int> nthPermRepGmp(int n, int r, double dblIdx, 
                               mpz_t mpzIdx, const std::vector<int> &Reps) {
    
    mpz_t temp, temp2, index1;
    mpz_init(temp); mpz_init(temp2); mpz_init(index1);
    mpz_set(index1, mpzIdx);
    std::vector<int> res(r);
    mpz_ui_pow_ui(temp, n, r);
    
    for (int k = 0; k < r; ++k) {
        mpz_divexact_ui(temp, temp, n);
        mpz_tdiv_q(temp2, index1, temp);
        int j = mpz_get_si(temp2);
        res[k] = j;
        mpz_submul_ui(index1, temp, j);
    }
    
    mpz_clear(temp); mpz_clear(temp2); mpz_clear(index1);
    return res;
}

std::vector<int> nthPermMultGmp(int n, int r, double dblIdx, 
                                mpz_t mpzIdx, const std::vector<int> &Reps) {
    
    mpz_t temp, index1;
    mpz_init(temp); mpz_init(index1);
    mpz_set(index1, mpzIdx);
    std::vector<int> res(r);
    
    mpz_add_ui(index1, index1, 1);
    std::vector<int> Counts;
    std::vector<int> TempReps = Reps;
    
    mpz_t test, index2;
    mpz_init(test); mpz_init(index2);
    mpz_set(index2, index1);
    
    for (int k = 0, r1 = r - 1; k < r; ++k, --r1) {
        
        int j = 0;
        while (TempReps[j] == 0)
            ++j;
        
        --TempReps[j];
        Counts = nonZeroVec(TempReps);
        MultisetPermRowNumGmp(temp, static_cast<int>(Counts.size()), r1, Counts);
        mpz_set(test, temp);
        
        while (mpz_cmp(test, index1) < 0) {
            mpz_sub(index2, index2, temp);
            ++TempReps[j];
            ++j;
            
            while (TempReps[j] == 0)
                ++j;
            
            --TempReps[j];
            Counts = nonZeroVec(TempReps);
            MultisetPermRowNumGmp(temp, static_cast<int>(Counts.size()), r1, Counts);
            mpz_add(test, test, temp);
        }
        
        res[k] = j;
        mpz_set(index1, index2);
    }
    
    mpz_clear(test); mpz_clear(index2);
    mpz_clear(temp); mpz_clear(index1);
    return res;
}

Rcpp::XPtr<nthResultPtr> putNthResPtrInXPtr(bool IsComb, bool IsMult,
                                            bool IsRep, bool IsGmp) {
    
    if (IsGmp) {
        if (IsComb) {
            if (IsMult) {
                return(Rcpp::XPtr<nthResultPtr>(new nthResultPtr(&nthCombMultGmp)));
            } else if (IsRep) {
                return(Rcpp::XPtr<nthResultPtr>(new nthResultPtr(&nthCombRepGmp)));
            } else {
                return(Rcpp::XPtr<nthResultPtr>(new nthResultPtr(&nthCombGmp)));
            }
        } else {
            if (IsMult) {
                return(Rcpp::XPtr<nthResultPtr>(new nthResultPtr(&nthPermMultGmp)));
            } else if (IsRep) {
                return(Rcpp::XPtr<nthResultPtr>(new nthResultPtr(&nthPermRepGmp)));
            } else {
                return(Rcpp::XPtr<nthResultPtr>(new nthResultPtr(&nthPermGmp)));
            }
        }
    } else {
        if (IsComb) {
            if (IsMult) {
                return(Rcpp::XPtr<nthResultPtr>(new nthResultPtr(&nthCombMult)));
            } else if (IsRep) {
                return(Rcpp::XPtr<nthResultPtr>(new nthResultPtr(&nthCombRep)));
            } else {
                return(Rcpp::XPtr<nthResultPtr>(new nthResultPtr(&nthComb)));
            }
        } else {
            if (IsMult) {
                return(Rcpp::XPtr<nthResultPtr>(new nthResultPtr(&nthPermMult)));
            } else if (IsRep) {
                return(Rcpp::XPtr<nthResultPtr>(new nthResultPtr(&nthPermRep)));
            } else {
                return(Rcpp::XPtr<nthResultPtr>(new nthResultPtr(&nthPerm)));
            }
        }
    }
}
