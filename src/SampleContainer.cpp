#include <NthResult.h>
#include <Rcpp.h>

template <typename typeRcpp, typename typeVector>
typeRcpp SampleResults(typeVector v, unsigned long int m, 
                       bool IsRep, std::vector<int> myReps, 
                       unsigned long int n, bool IsComb,
                       std::vector<double> mySample) {
    
    typeRcpp sampleMatrix = Rcpp::no_init_matrix(n, m);
    int lenV = v.size();
    std::vector<int> z(m);
    bool IsMult = false;
    
    if ((int) myReps.size() == lenV)
        IsMult = true;
    
    if (IsComb) {
        for (std::size_t i = 0; i < n; ++i) {
            z = nthCombination(lenV, m, mySample[i] - 1, IsRep, IsMult, myReps);
            for (std::size_t j = 0; j < m; ++j)
                sampleMatrix(i, j) = v[z[j]];
        }
    } else {
        for (std::size_t i = 0; i < n; ++i) {
            z = nthPermutation(lenV, m, mySample[i] - 1, IsRep, IsMult, myReps);
            for (std::size_t j = 0; j < m; ++j)
                sampleMatrix(i, j) = v[z[j]];
        }
    }
    
    return sampleMatrix;
}


// This function only contains a couple of error messages (for RSample 
// & maxNum) because most of the inputs have already been sanitized by
// their respective count functions (i.e. combo/permuteCount).

// [[Rcpp::export]]
SEXP SampleRcpp(SEXP Rv, SEXP Rm, SEXP Rrepetition, 
                SEXP RFreqs, SEXP Rsample, bool IsComb,
                bool IsFactor, double computedRows) {
    
    int m = 0, m1, m2;
    int lenFreqs = 0;
    double seqEnd;
    bool IsRepetition = false;
    bool IsCharacter = false;
    bool IsInteger = false;
    
    std::vector<double> vNum, mySample;
    std::vector<int> vInt, myReps, freqsExpanded;
    Rcpp::CharacterVector vStr;
    
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
    }
    
    if (!Rf_isNull(RFreqs)) {
        myReps = Rcpp::as<std::vector<int> >(RFreqs);

        lenFreqs = (int) myReps.size();
        for (int i = 0; i < lenFreqs; ++i)
            for (int j = 0; j < myReps[i]; ++j)
                freqsExpanded.push_back(i);
    }
    
    if (Rf_isNull(Rm)) {
        m = freqsExpanded.size();
    } else {
        m = Rcpp::as<int>(Rm);
    }
    
    switch(TYPEOF(Rsample)) {
        case REALSXP: {
            mySample = Rcpp::as<std::vector<double> >(Rsample);
            break;
        }
        case INTSXP: {
            mySample = Rcpp::as<std::vector<double> >(Rsample);
            break;
        }
        default: {
            Rcpp::stop("sampleVec must be of type numeric or integer");
        }
    }

    IsRepetition = Rcpp::as<bool>(Rrepetition);
    
    if (IsCharacter) {
        vStr = Rcpp::as<Rcpp::CharacterVector>(Rv);
    } else {
        if (Rf_length(Rv) == 1) {
            seqEnd = Rcpp::as<double>(Rv);
            if (Rcpp::NumericVector::is_na(seqEnd)) {seqEnd = 1;}
            if (seqEnd > 1) {m1 = 1; m2 = seqEnd;} else {m1 = seqEnd; m2 = 1;}
            Rcpp::IntegerVector vTemp = Rcpp::seq(m1, m2);
            IsInteger = true;
            vNum = Rcpp::as<std::vector<double> >(vTemp);
        } else {
            vNum = Rcpp::as<std::vector<double> >(Rv);
        }
    }
    
    if (IsInteger)
        vInt.assign(vNum.begin(), vNum.end());
    
    if (IsFactor)
        IsCharacter = IsInteger = false;
    
    double myMax = *std::max_element(mySample.begin(), mySample.end());
    if (myMax > computedRows)
        Rcpp::stop("One or more of the requested values in sampleVec "
                 "exceeds the maximum number of possible results");
    
    unsigned long int sampSize = mySample.size();
    
    if (IsCharacter) {
        return SampleResults<Rcpp::CharacterMatrix>(vStr, m, IsRepetition, myReps, 
                                              sampSize, IsComb, mySample);
    } else if (IsFactor) {
        Rcpp::IntegerMatrix factorMat;
        Rcpp::IntegerVector testFactor = Rcpp::as<Rcpp::IntegerVector>(Rv);
        Rcpp::CharacterVector myClass = testFactor.attr("class");
        Rcpp::CharacterVector myLevels = testFactor.attr("levels");
        
        factorMat = SampleResults<Rcpp::IntegerMatrix>(vInt, m, IsRepetition, myReps, 
                                                     sampSize, IsComb, mySample);
        
        factorMat.attr("class") = myClass;
        factorMat.attr("levels") = myLevels;
        
        return factorMat;
    } else if (IsInteger) {
        return SampleResults<Rcpp::IntegerMatrix>(vInt, m, IsRepetition, myReps, 
                                            sampSize, IsComb, mySample);
    } else {
        return SampleResults<Rcpp::NumericMatrix>(vNum, m, IsRepetition, myReps, 
                                            sampSize, IsComb, mySample);
    }
}

