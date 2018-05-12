#include <NthCombinations.h>
#include <NthPermutations.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP NthResultRcpp(SEXP Rv, SEXP Rm, SEXP Rind, SEXP Rrepetition, 
                   SEXP RIsComb, SEXP RIsFactor, SEXP RFreqs) {
    
    int n, m = 0, i, j, m1, m2, lenFreqs = 0;
    double ind, computedRows, seqEnd;
    bool IsRepetition, IsInteger, IsFactor;
    bool IsComb, IsMultiset, IsCharacter;
    
    std::vector<double> vNum;
    std::vector<int> vInt, myReps, freqsExpanded;
    std::vector<std::string > vStr;
    
    switch(TYPEOF(Rv)) {
        case INTSXP: {
            IsCharacter = false;
            IsInteger = true;
            break;
        }
        case REALSXP: {
            IsInteger = IsCharacter = false;
            break;
        }
        case STRSXP: {
            IsCharacter = true;
            IsInteger = false;
            break;
        }
        default: {
            stop("Only integers, numerical, character, and factor classes are supported for v");   
        }
    }
    
    IsComb = as<bool>(RIsComb);
    
    if (Rf_isNull(RFreqs)) {
        IsMultiset = false;
        myReps.push_back(1);
    } else {
        IsMultiset = true;
        switch(TYPEOF(RFreqs)) {
            case REALSXP: {
                myReps = as<std::vector<int> >(RFreqs);
                break;
            }
            case INTSXP: {
                myReps = as<std::vector<int> >(RFreqs);
                break;
            }
            default: {
                stop("freqs must be of type numeric or integer");
            }
        }
        lenFreqs = myReps.size();
        for (i = 0; i < lenFreqs; i++) {
            if (myReps[i] < 1)
                stop("each element in freqs must be a positive number");
            
            for (j = 0; j < myReps[i]; j++)
                freqsExpanded.push_back(i);
        }
    }
    
    if (Rf_isNull(Rm)) {
        if (IsMultiset)
            m = freqsExpanded.size();
        else
            stop("m and freqs cannot both be NULL");
    } else {
        if (Rf_length(Rm) > 1)
            stop("length of m must be 1");
        
        switch(TYPEOF(Rm)) {
            case REALSXP: {
                m = as<int>(Rm);
                break;
            }
            case INTSXP: {
                m = as<int>(Rm);
                break;
            }
            default: {
                stop("m must be of type numeric or integer");
            }
        }
    }
    
    if (m < 1)
        stop("m must be positive");
    
    switch(TYPEOF(Rind)) {
        case REALSXP: {
            ind = as<double>(Rind);
            break;
        }
        case INTSXP: {
            ind = as<double>(Rind);
            break;
        }
        default: {
            stop("ind must be of type numeric or integer");
        }
    }
    
    ind--;
    if (ind < 0)
        stop("index must be positive");

    if (!Rf_isLogical(Rrepetition)) 
        stop("repetitions must be a logical value");
    
    IsRepetition = as<bool>(Rrepetition);
    IsFactor = as<bool>(RIsFactor);
    
    if (IsCharacter) {
        vStr = as<std::vector<std::string > >(Rv);
        n = vStr.size();
    } else {
        if (Rf_length(Rv) == 1) {
            seqEnd = as<double>(Rv);
            if (NumericVector::is_na(seqEnd)) {seqEnd = 1;}
            if (seqEnd > 1) {m1 = 1; m2 = seqEnd;} else {m1 = seqEnd; m2 = 1;}
            IntegerVector vTemp = seq(m1, m2);
            IsInteger = true;
            vNum = as<std::vector<double> >(vTemp);
        } else {
            vNum = as<std::vector<double> >(Rv);
        }
        
        for (i = (vNum.size() - 1); i >= 0; i--)
            if (NumericVector::is_na(vNum[i]))
                vNum.erase(vNum.begin() + i);

        n = vNum.size();
    }
    
    if (IsInteger)
        vInt.assign(vNum.begin(), vNum.end());
    
    if (IsFactor)
        IsCharacter = IsInteger = false;
    
    if (IsMultiset) {
        if (n != lenFreqs)
            stop("the length of freqs must equal the length of v");
        
        if (IsComb) {
            if (m > (int) freqsExpanded.size())
                m = freqsExpanded.size();
            
            computedRows = MultisetCombRowNum(n, m, myReps);
        } else {
            if (Rf_isNull(Rm)) {
                computedRows = NumPermsWithRep(freqsExpanded);
            } else if (m == (int) freqsExpanded.size()) {
                computedRows = NumPermsWithRep(freqsExpanded);
            } else {
                if (m > (int) freqsExpanded.size())
                    m = freqsExpanded.size();
                
                computedRows = MultisetCombRowNum(n, m, myReps);
            }
        }
    } else {
        if (IsRepetition) {
            if (IsComb)
                computedRows = NumCombsWithRep(n, m);
            else
                computedRows = std::pow((double) n, (double) m);
        } else {
            if (m > n)
                stop("m must be less than or equal to the length of v");
            
            if (IsComb)
                computedRows = nChooseK(n, m);
            else
                computedRows = NumPermsNoRep(n, m);
        }
    }
    
    if (ind >= computedRows)
        stop("The requested index exceeds the number of possible results");
    
    if (IsCharacter) {
        if (IsMultiset) {
            if (IsComb)
                return wrap(nthCombination(n, m, ind, vStr, IsRepetition));
            // else
            //     return wrap(nthPermNoRep(n, m, ind, vStr, IsRepetition));
        } else {
            if (IsComb)
                return wrap(nthCombination(n, m, ind, vStr, IsRepetition));
            else
                return wrap(nthPermutation(n, m, ind, vStr, IsRepetition));
        }
    } else if (IsInteger) {
        if (IsMultiset) {
            if (IsComb)
                return wrap(nthCombination(n, m, ind, vInt, IsRepetition));
            // else
            //     return wrap(nthPermNoRep(n, m, ind, vInt, IsRepetition));
        } else {
            if (IsComb)
                return wrap(nthCombination(n, m, ind, vInt, IsRepetition));
            else
                return wrap(nthPermutation(n, m, ind, vInt, IsRepetition));
        }
    } else if (IsFactor) {
        if (IsMultiset) {
            if (IsComb)
                return wrap(nthCombination(n, m, ind, vStr, IsRepetition));
            // else
            //     return wrap(nthPermNoRep(n, m, ind, vStr, IsRepetition));
        } else {
            if (IsComb)
                return wrap(nthCombination(n, m, ind, vStr, IsRepetition));
            else
                return wrap(nthCombination(n, m, ind, vStr, IsRepetition));
        }
    } else {
        if (IsMultiset) {
            if (IsComb)
                return wrap(nthCombination(n, m, ind, vNum, IsRepetition));
            // else
            //     return wrap(nthPermNoRep(n, m, ind, vNum, IsRepetition));
        } else {
            if (IsComb)
                return wrap(nthCombination(n, m, ind, vNum, IsRepetition));
            else
                return wrap(nthPermutation(n, m, ind, vNum, IsRepetition));
        }
    }
    return wrap(1);
}
