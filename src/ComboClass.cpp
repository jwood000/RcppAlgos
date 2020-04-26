#include "ComboClass.h"

SEXP Combo::VecReturn(const SEXP &v) {
    
    switch (TYPEOF(v)) {
        case LGLSXP: {
            Rcpp::LogicalVector boolVec(Rcpp::clone(v));
            return VecPopulate(boolVec);
        } case INTSXP: {
            Rcpp::IntegerVector intVec(Rcpp::clone(v));
            return VecPopulate(intVec);
        } case STRSXP: {
            Rcpp::CharacterVector charVec(Rcpp::clone(v));
            return VecPopulate(charVec);
        } case CPLXSXP: {
            Rcpp::ComplexVector cmplxVec(Rcpp::clone(v));
            return VecPopulate(cmplxVec);
        } case RAWSXP: {
            Rcpp::RawVector rawVec(Rcpp::clone(v));
            return VecPopulate(rawVec);
        } default: {
            Rcpp::NumericVector numVec(Rcpp::clone(v));
            return VecPopulate(numVec);
        }
    }
}

SEXP Combo::MatReturn(const SEXP &v, int nr, bool IsReverse) {
    
    switch (TYPEOF(v)) {
        case LGLSXP: {
            Rcpp::LogicalVector bv(Rcpp::clone(v));
            return (IsReverse) ? MatComboReverse(bv, nr) : MatCombo(bv, nr);
        } case INTSXP: {
            Rcpp::IntegerVector iv(Rcpp::clone(v));
            return (IsReverse) ? MatComboReverse(iv, nr) : MatCombo(iv, nr);
        } case STRSXP: {
            Rcpp::CharacterVector cv(Rcpp::clone(v));
            return (IsReverse) ? MatComboReverse(cv, nr) : MatCombo(cv, nr);
        } case CPLXSXP: {
            Rcpp::ComplexVector xv(Rcpp::clone(v));
            return (IsReverse) ? MatComboReverse(xv, nr) : MatCombo(xv, nr);
        } case RAWSXP: {
            Rcpp::RawVector rv(Rcpp::clone(v));
            return (IsReverse) ? MatComboReverse(rv, nr) : MatCombo(rv, nr);
        } default: {
            Rcpp::NumericVector nv(Rcpp::clone(v));
            return (IsReverse) ? MatComboReverse(nv, nr) : MatCombo(nv, nr);
        }
    }
}

SEXP Combo::SampReturn(const SEXP &v, const std::vector<double> &mySample,
                       mpz_t *const myBigSamp, std::size_t sampSize) {
    
    switch (TYPEOF(v)) {
        case LGLSXP: {
            Rcpp::LogicalVector boolVec(Rcpp::clone(v));
            return SampMatrix(boolVec, mySample, myBigSamp, sampSize);
        } case INTSXP: {
            Rcpp::IntegerVector intVec(Rcpp::clone(v));
            return SampMatrix(intVec, mySample, myBigSamp, sampSize);
        } case STRSXP: {
            Rcpp::CharacterVector charVec(Rcpp::clone(v));
            return SampMatrix(charVec, mySample, myBigSamp, sampSize);
        } case CPLXSXP: {
            Rcpp::ComplexVector cmplxVec(Rcpp::clone(v));
            return SampMatrix(cmplxVec, mySample, myBigSamp, sampSize);
        } case RAWSXP: {
            Rcpp::RawVector rawVec(Rcpp::clone(v));
            return SampMatrix(rawVec, mySample, myBigSamp, sampSize);
        } default: {
            Rcpp::NumericVector numVec(Rcpp::clone(v));
            return SampMatrix(numVec, mySample, myBigSamp, sampSize);
        }
    }
}

template <int RTYPE>
Rcpp::Vector<RTYPE> Combo::VecPopulate(const Rcpp::Vector<RTYPE> &v) {
    
    Rcpp::Vector<RTYPE> vOut(m);
    
    for (int j = 0; j < m; ++j)
        vOut[j] = v[z[j]];
    
    if (IsFactor) {
        Rf_setAttrib(vOut, R_ClassSymbol, myClass);
        Rf_setAttrib(vOut, R_LevelsSymbol, myLevels);
    }
    
    return vOut;
}

template <int RTYPE>
Rcpp::Matrix<RTYPE> Combo::MatCombo(const Rcpp::Vector<RTYPE> &v, int nRows) {
    
    Rcpp::Matrix<RTYPE> matRcpp = Rcpp::no_init_matrix(nRows, m);
    
    Rcpp::XPtr<combPermPtr<Rcpp::Matrix<RTYPE>, Rcpp::Vector<RTYPE>>> xpFunCoPePtr =
        putCombPtrInXPtr<Rcpp::Matrix<RTYPE>, Rcpp::Vector<RTYPE>>(IsComb, IsMult, IsRep, true);

    const combPermPtr<Rcpp::Matrix<RTYPE>, Rcpp::Vector<RTYPE>> myFunCombPerm = *xpFunCoPePtr;
    myFunCombPerm(matRcpp, v, z, n, m, 0, nRows, freqs);

    if (IsFactor) {
        Rf_setAttrib(matRcpp, R_ClassSymbol, myClass);
        Rf_setAttrib(matRcpp, R_LevelsSymbol, myLevels);
    }

    // Update z. We cannot pass by reference since this would cause
    // issues with the parallel version in the general function
    Rcpp::Vector<RTYPE> myRow = matRcpp.row(nRows - 1);
    z = zUpdateIndex(v, myRow, m);
    TopOffPartialPerm(z, myReps, n, m, IsComb, IsRep, IsMult);
    return matRcpp;
}

template <int RTYPE>
Rcpp::Matrix<RTYPE> Combo::MatComboReverse(const Rcpp::Vector<RTYPE> &v, int nRows) {
    
    Rcpp::Matrix<RTYPE> matRcpp = Rcpp::no_init_matrix(nRows, m);
    
    const int loc_n1 = IsComb ? n - 1 : (IsMult ? freqs.size() - 1 : n - 1);
    const int lastRow = nRows - 1;
    const int loc_m = m;
    
    // We iterater to the pentultimate row to avoid iterating z one too many times
    for (int count = 0, loc_m1 = m - 1; count < lastRow; ++count) {
        for (int j = 0; j < loc_m; ++j)
            matRcpp(count, j) = v[z[j]];
        
        prevIter(freqs, z, loc_n1, loc_m1);
    }
    
    // Get the last result
    for (int j = 0; j < loc_m; ++j)
        matRcpp(lastRow, j) = v[z[j]];

    if (IsFactor) {
        Rf_setAttrib(matRcpp, R_ClassSymbol, myClass);
        Rf_setAttrib(matRcpp, R_LevelsSymbol, myLevels);
    }

    return matRcpp;
}

template <int RTYPE>
Rcpp::Matrix<RTYPE> Combo::SampMatrix(const Rcpp::Vector<RTYPE> &v,
                                      const std::vector<double> &mySample,
                                      mpz_t *const myBigSamp, std::size_t sampSize) {

    Rcpp::Matrix<RTYPE> matRcpp = Rcpp::no_init_matrix(sampSize, m);
    
    if (IsGmp) {
        constexpr double dblDefault = 0;
        
        for (std::size_t i = 0; i < sampSize; ++i) {
            const std::vector<int> tempZ = nthResFun(n, m, dblDefault, myBigSamp[i], myReps);
            
            for (int j = 0; j < m; ++j)
                matRcpp(i, j) = v[tempZ[j]];
        }
    } else {
        for (std::size_t i = 0; i < sampSize; ++i) {
            const std::vector<int> tempZ = nthResFun(n, m, mySample[i], mpzTemp, myReps);
            
            for (int j = 0; j < m; ++j)
                matRcpp(i, j) = v[tempZ[j]];
        }
    }

    if (IsFactor) {
        Rf_setAttrib(matRcpp, R_ClassSymbol, myClass);
        Rf_setAttrib(matRcpp, R_LevelsSymbol, myLevels);
    }

    return matRcpp;
}

// The bVec Vector represents IsFac, IsComb, IsMult, IsRep, IsGmp, & IsFull
Combo::Combo(SEXP Rv, int Rm, SEXP RcompRows,
             Rcpp::LogicalVector bVec, std::vector<std::vector<int>> freqInfo)
            : n(Rf_length(Rv)), m(Rm), m1(Rm - 1), sexpVec(Rv),
              IsFactor(bVec[0]), IsComb(bVec[1]), IsMult(bVec[2]),
              IsRep(bVec[3]), IsGmp(bVec[4]), 
              computedRows(IsGmp ? 0 : Rcpp::as<double>(RcompRows)),
              freqs(freqInfo[0]), myReps(freqInfo[1]),
              testFactor(bVec[0] ? Rcpp::as<Rcpp::IntegerVector>(Rv) :
                             Rcpp::IntegerVector::create(1)),
              myClass(bVec[0] ? testFactor.attr("class") :
                          Rcpp::CharacterVector::create("empty")),
              myLevels(bVec[0] ? testFactor.attr("levels") :
                           Rcpp::CharacterVector::create("empty")),
              nthResFun(*putNthResPtrInXPtr(bVec[1], bVec[2], bVec[3], bVec[4])),
              nextIter(*putNextIterPtrInXPtr(bVec[1], bVec[2], bVec[3], bVec[5])),
              prevIter(*putPrevIterPtrInXPtr(bVec[1], bVec[2], bVec[3], bVec[5])),
              n1(IsComb ? n - 1 : (IsMult ? freqs.size() - 1 : n - 1)) {

    z.resize(Rm);

    // Initialize trivial mpz_t value for functions which require mpz_t
    mpz_init(mpzTemp);
    mpz_init(mpzIndex[0]);
    mpz_init(computedRowsMpz[0]);
    
    if (IsGmp)
        createMPZArray(RcompRows, computedRowsMpz, 1, "index vector");
    
    dblIndex = 0;

    SetStartZ(n, m, dblIndex, 0, mpzIndex[0], IsRep, 
              IsComb, IsMult, IsGmp, myReps, freqs, z, nthResFun);
}

void Combo::startOver() {
    if (IsGmp) {
        mpz_set_ui(mpzIndex[0], 0u);
    } else {
        dblIndex = 0;
    }
    
    SetStartZ(n, m, dblIndex, 0, mpzIndex[0], IsRep, 
              IsComb, IsMult, IsGmp, myReps, freqs, z, nthResFun);
}

SEXP Combo::nextComb() {
    
    if (CheckEqSi(IsGmp, mpzIndex[0], dblIndex, 0)) {
        increment(IsGmp, mpzIndex[0], dblIndex);
        return VecReturn(sexpVec);
    } else if (CheckIndLT(IsGmp, mpzIndex[0], dblIndex, computedRowsMpz[0], computedRows)) {
        increment(IsGmp, mpzIndex[0], dblIndex);
        nextIter(freqs, z, n1, m1);
        return VecReturn(sexpVec);
    } else if (CheckEqInd(IsGmp, mpzIndex[0], dblIndex, computedRowsMpz[0], computedRows)) {
        Rcpp::Rcout << "No more results. To see the last result,"
                       " use the prevIter method(s)\n" << std::endl;
        increment(IsGmp, mpzIndex[0], dblIndex);
        return Rcpp::wrap(false);
    } else {
        return Rcpp::wrap(false);
    }
}

SEXP Combo::prevComb() {
    
    if (CheckIndGrT(IsGmp, mpzIndex[0], dblIndex, computedRowsMpz[0], computedRows)) {
        decrement(IsGmp, mpzIndex[0], dblIndex);
        return VecReturn(sexpVec);
    } else if (CheckGrTSi(IsGmp, mpzIndex[0], dblIndex, 1)) {
        decrement(IsGmp, mpzIndex[0], dblIndex);
        prevIter(freqs, z, n1, m1);
        return VecReturn(sexpVec);
    } else if (CheckEqSi(IsGmp, mpzIndex[0], dblIndex, 1)) {
        Rcpp::Rcout << "Iterator Initialized. To see the first result,"
                       " use the nextIter method(s)\n" << std::endl;
        decrement(IsGmp, mpzIndex[0], dblIndex);
        return Rcpp::wrap(false);
    } else {
        return Rcpp::wrap(false);
    }
}

SEXP Combo::nextNumCombs(SEXP RNum) {
    
    int num;
    CleanConvert::convertPrimitive(RNum, num, "The number of results");

    if (CheckIndLT(IsGmp, mpzIndex[0], dblIndex, computedRowsMpz[0], computedRows)) {
        int nRows = 0;
        int numIncrement = 0;
        
        if (IsGmp) {
            mpz_sub(mpzTemp, computedRowsMpz[0], mpzIndex[0]);
            nRows = mpz_cmp_si(mpzTemp, num) < 0 ? mpz_get_si(mpzTemp) : num;
            numIncrement = mpz_cmp_si(mpzTemp, num) < 0 ? (nRows + 1) : nRows;
        } else {
            dblTemp = computedRows - dblIndex;
            nRows = num > dblTemp ? dblTemp : num;
            numIncrement = num > dblTemp ? (nRows + 1) : nRows;
        }
        
        if (CheckGrTSi(IsGmp, mpzIndex[0], dblIndex, 0))
            nextIter(freqs, z, n1, m1);
        
        increment(IsGmp, mpzIndex[0], dblIndex, numIncrement);
        return MatReturn(sexpVec, nRows, false);
    } else if (CheckEqInd(IsGmp, mpzIndex[0], dblIndex, computedRowsMpz[0], computedRows)) {
        Rcpp::Rcout << "No more results. To see the last result,"
                       " use the prevIter method(s)\n" << std::endl;
        increment(IsGmp, mpzIndex[0], dblIndex);
        return Rcpp::wrap(false);
    } else {
        return Rcpp::wrap(false);
    }
}

SEXP Combo::prevNumCombs(SEXP RNum) {
    
    int num;
    CleanConvert::convertPrimitive(RNum, num, "The number of results");

    if (CheckGrTSi(IsGmp, mpzIndex[0], dblIndex, 2)) {
        int nRows = 0;
        int numDecrement = 0;
        
        if (IsGmp) {
            mpz_sub_ui(mpzTemp, mpzIndex[0], 1u);
            nRows = mpz_cmp_si(mpzTemp, num) < 0 ? mpz_get_si(mpzTemp) : num;
            numDecrement = mpz_cmp_si(mpzTemp, num) < 0 ? (nRows + 1) : nRows;
        } else {
            dblTemp = dblIndex - 1;
            nRows = num > dblTemp ? dblTemp : num;
            numDecrement = num > dblTemp ? (nRows + 1) : nRows;
        }
        
        if (CheckIndLT(IsGmp, mpzIndex[0], dblIndex, computedRowsMpz[0], computedRows, true))
            prevIter(freqs, z, n1, m1);
        
        decrement(IsGmp, mpzIndex[0], dblIndex, numDecrement);
        return MatReturn(sexpVec, nRows, true);
    } else if (CheckEqSi(IsGmp, mpzIndex[0], dblIndex, 2)) {
        Rcpp::Rcout << "No more results. To see the last result,"
                       " use the prevIter method(s)\n" << std::endl;
        decrement(IsGmp, mpzIndex[0], dblIndex);
        return Rcpp::wrap(false);
    } else {
        return Rcpp::wrap(false);
    }
}

SEXP Combo::nextGather() {
    
    if (IsGmp) {
        mpz_sub(mpzTemp, computedRowsMpz[0], mpzIndex[0]);
        
        if (mpz_cmp_si(mpzTemp, std::numeric_limits<int>::max()) > 0) {
            Rcpp::stop("The number of requested rows is greater than " +
                std::to_string(std::numeric_limits<int>::max()));
        }
    } else {
        dblTemp = computedRows - dblIndex;
        
        if (dblTemp > std::numeric_limits<int>::max()) {
            Rcpp::stop("The number of requested rows is greater than " +
                std::to_string(std::numeric_limits<int>::max()));
        }
    }

    const int nRows = (IsGmp) ? mpz_get_si(mpzTemp) : dblTemp;

    if (nRows > 0) {
        if (CheckGrTSi(IsGmp, mpzIndex[0], dblIndex, 0))
            nextIter(freqs, z, n1, m1);
        
        if (IsGmp) {
            mpz_add_ui(mpzIndex[0], computedRowsMpz[0], 1u);
        } else {
            dblIndex = computedRows + 1;
        }
        
        return MatReturn(sexpVec, nRows, false);
    } else {
        return Rcpp::wrap(false);
    }
}

SEXP Combo::prevGather() {
    
    if (IsGmp) {
        mpz_sub_ui(mpzTemp, mpzIndex[0], 1);
        
        if (mpz_cmp_si(mpzTemp, std::numeric_limits<int>::max()) > 0) {
            Rcpp::stop("The number of requested rows is greater than " +
                std::to_string(std::numeric_limits<int>::max()));
        }
    } else {
        dblTemp = dblIndex - 1;
        
        if (dblTemp > std::numeric_limits<int>::max()) {
            Rcpp::stop("The number of requested rows is greater than " +
                std::to_string(std::numeric_limits<int>::max()));
        }
    }

    const int nRows = (IsGmp) ? mpz_get_si(mpzTemp) : dblTemp;

    if (nRows) {
        if (CheckIndLT(IsGmp, mpzIndex[0], dblIndex, computedRowsMpz[0], computedRows, true))
            prevIter(freqs, z, n1, m1);
        
        if (IsGmp) {
            mpz_set_si(mpzIndex[0], 0u);
        } else {
            dblIndex = 0;
        }
        
        return MatReturn(sexpVec, nRows, true);
    } else {
        return Rcpp::wrap(false);
    }
}

SEXP Combo::currComb() {
    
    if (CheckIndGrT(IsGmp, mpzIndex[0], dblIndex, computedRowsMpz[0], computedRows)) {
        Rcpp::Rcout << "No more results. To see the last result,"
                        " use the prevIter method(s)\n" << std::endl;
        return Rcpp::wrap(false);
    } else if (CheckGrTSi(IsGmp, mpzIndex[0], dblIndex, 0)) {
        return VecReturn(sexpVec);
    } else {
        Rcpp::Rcout << "Iterator Initialized. To see the first result,"
                       " use the nextIter method(s)\n" << std::endl;
        return Rcpp::wrap(false);
    }
}

SEXP Combo::combIndex(SEXP RindexVec) {
    
    std::size_t sampSize;
    std::vector<double> mySample;
    SetIndexVec(RindexVec, mySample, sampSize, IsGmp, computedRows);
    
    const std::size_t bigSampSize = (IsGmp) ? sampSize : 1;
    auto mpzVec = FromCpp14::make_unique<mpz_t[]>(bigSampSize);
    
    for (std::size_t i = 0; i < bigSampSize; ++i)
        mpz_init(mpzVec[i]);
    
    if (IsGmp)
        SetIndexVecMpz(RindexVec, mpzVec.get(), sampSize, computedRowsMpz[0]);

    if (sampSize > 1) {
        return SampReturn(sexpVec, mySample, mpzVec.get(), sampSize);
    } else {
        if (IsGmp) {
            mpz_add_ui(mpzIndex[0], mpzVec[0], 1u);
            mpz_set(mpzTemp, mpzVec[0]);
        } else {
            dblIndex = mySample.front() + 1;
            dblTemp = mySample.front();
        }
        
        z = nthResFun(n, m, dblTemp, mpzTemp, myReps);
        TopOffPartialPerm(z, myReps, n, m, IsComb, IsRep, IsMult);
        return VecReturn(sexpVec);
    }
}

SEXP Combo::front() {
    
    if (IsGmp) {
        mpz_set_ui(mpzIndex[0], 1u);
        mpz_set_ui(mpzTemp, 0u);
    } else {
        dblIndex = 1;
        dblTemp = 0;
    }
    
    z = nthResFun(n, m, dblTemp, mpzTemp, myReps);
    TopOffPartialPerm(z, myReps, n, m, IsComb, IsRep, IsMult);
    return VecReturn(sexpVec);
}

SEXP Combo::back() {
    
    if (IsGmp) {
        mpz_set(mpzIndex[0], computedRowsMpz[0]);
        mpz_sub_ui(mpzTemp, computedRowsMpz[0], 1u);
    } else {
        dblIndex = computedRows;
        dblTemp = computedRows - 1;
    }
    
    z = nthResFun(n, m, dblTemp, mpzTemp, myReps);
    TopOffPartialPerm(z, myReps, n, m, IsComb, IsRep, IsMult);
    return VecReturn(sexpVec);
}

SEXP Combo::sourceVector() const {
    return sexpVec;
}

Rcpp::List Combo::summary() {
    const std::string CoPerm = (IsComb) ? "Combinations " : "Permutations ";
    const std::string RepStr = (IsRep) ? "with repetition " : "";
    const std::string MultiStr = (IsMult) ? "of a multiset " : "";
    const std::string strDesc = CoPerm + RepStr + MultiStr + "of "
                                + std::to_string(n) + " choose " + std::to_string(m);
    const double dblDiff = (IsGmp) ? 0 : computedRows - dblIndex;
    
    if (IsGmp)
        mpz_sub(mpzTemp, computedRowsMpz[0], mpzIndex[0]);

    return Rcpp::List::create(Rcpp::Named("description") = Rcpp::wrap(strDesc),
                              Rcpp::Named("currentIndex") = GetCount(IsGmp, mpzIndex[0], dblIndex),
                              Rcpp::Named("totalResults") = GetCount(IsGmp, computedRowsMpz[0], computedRows),
                              Rcpp::Named("totalRemaining") = GetCount(IsGmp, mpzTemp, dblDiff));
}

RCPP_MODULE(Combo) {
    Rcpp::class_<Combo>("Combo")
    
        .constructor<SEXP, int, SEXP, 
                     Rcpp::LogicalVector,
                     std::vector<std::vector<int>>>("constructor")
    
        .method("nextIter", &Combo::nextComb, "Get the next combination/permutation")
        .method("nextNIter", &Combo::nextNumCombs, "Get the next n combinations/permutations")
        .method("nextRemaining", &Combo::nextGather, "Get all of the remaining combinations/permutations")
        .method("currIter", &Combo::currComb, "Get the current combination/permutation")
        .method("prevIter", &Combo::prevComb, "Get the previous combination/permutation")
        .method("prevNIter", &Combo::prevNumCombs, "Get the previous n combinations/permutations")
        .method("prevRemaining", &Combo::prevGather, "Get all of the previous combinations/permutations")
        .method("startOver", &Combo::startOver, "Resets the iterator")
        .method("sourceVector", &Combo::sourceVector, "View the source vector")
        .method("summary", &Combo::summary, "See a list of summary information about the iterator")
        .method("front", &Combo::front, "Get the first lexicographical combination/permutation")
        .method("back", &Combo::back, "Get the last lexicographical combination/permutation")
        .method("[[", &Combo::combIndex, "Random access method. Pass a single index or vector of indices")
        ;
}
