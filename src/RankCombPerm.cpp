#include "ImportExportMPZ.h"
#include "Cpp14MakeUnique.h"
#include "ComputedCount.h"
#include "SetUpUtils.h"
#include "RankResult.h"
#include <unordered_map>
#include <thread>

constexpr std::size_t numb = 8 * intSize;

void RankResults(mpz_t* bigRes, int* intRes, double* dblRes,
                 std::vector<int> &idx, const std::vector<int> &myReps,
                 rankResultPtr rankFun, int m, int lenV, int numResults,
                 bool IsGmp, bool IsInteger) {

    mpz_t mpzIdx;
    mpz_init_set_ui(mpzIdx, 0u);

    if (IsGmp) {
        for (int i = 0, j = 0; i < numResults; ++i, j += m) {
            double dblIdx = 0;
            rankFun(idx.begin() + j, lenV, m, dblIdx, mpzIdx, myReps);
            mpz_add_ui(mpzIdx, mpzIdx, 1u);
            mpz_set(bigRes[i], mpzIdx);
        }
    } else {
        for (int i = 0, j = 0; i < numResults; ++i, j += m) {
            double dblIdx = 0;
            rankFun(idx.begin() + j, lenV, m, dblIdx, mpzIdx, myReps);
            ++dblIdx;
            if (IsInteger) intRes[i] = dblIdx; else dblRes[i] = dblIdx;
        }
    }

    mpz_clear(mpzIdx);
}

[[cpp11::register]]
SEXP RankCombPerm(SEXP RIdx, SEXP Rv, SEXP RisRep,
                  SEXP RFreqs, SEXP Rm, SEXP RIsComb) {

    int m = 0;
    VecType myType = VecType::Integer;

    SetType(myType, Rv);
    CleanConvert::convertPrimitive(Rm, m, VecType::Integer, "m");

    bool IsRep = CleanConvert::convertFlag(RisRep, "repetition");
    const bool IsComb = CleanConvert::convertFlag(RIsComb, "IsComb");
    bool IsMult = false;

    std::vector<int> freqs;
    std::vector<int> myReps;
    std::vector<int> idx;
    CleanConvert::convertVector(RIdx, idx, VecType::Integer, "idx");

    for (auto &i: idx) {
        --i;
    }

    if (IsComb) {
        for (auto it = idx.begin(); it != idx.end(); it += m) {
            std::sort(it, it + m);
        }
    }

    int n = GetLength(Rv, myType);
    SetFreqsAndM(myReps, freqs, RFreqs, Rm, n, m, IsMult, IsRep);

    if (IsMult) {
        // See comment in SetUpUtils.cpp in the SetFinalValues function
        if (n != static_cast<int>(myReps.size())) {
            cpp11::stop("the length of freqs must equal the length of v");
        }

        if (m > static_cast<int>(freqs.size())) {
            m = freqs.size();
        }
    } else if (!IsRep && m > n) {
        cpp11::stop("m must be less than or equal to the length of v");
    }

    if (IsMult || !IsRep) {
        for (auto it = idx.begin(); it != idx.end();) {
            std::unordered_map<int, int> table;

            for (int i = 0; i < m; ++i, ++it) {
                ++table[*it];
            }

            if (IsMult) {
                for (auto&& tbl: table) {
                    if (myReps[tbl.first] < tbl.second) {
                        cpp11::stop("Input frequencies do not "
                                        "match supplied freqs");
                    }
                }
            } else {
                for (auto&& tbl: table) {
                    if (tbl.second > 1) {
                        cpp11::stop("No duplicates allowed when "
                                        "repetition = FALSE and freqs = NULL");
                    }
                }
            }
        }
    }

    const double computedRows = GetComputedRows(IsMult, IsComb, IsRep,
                                                n, m, Rm, freqs, myReps);
    const bool IsGmp     = (computedRows > Significand53);
    const bool IsInteger = (computedRows <= std::numeric_limits<int>::max());

    const int numResults = Rf_length(RIdx) / m;
    const rankResultPtr rankFun = GetRankResultFunc(
        IsComb, IsMult, IsRep, IsGmp
    );
    cpp11::sexp res_std_int = Rf_allocVector(INTSXP, (IsInteger) ?
                                                 numResults : 0);
    int* res_int = INTEGER(res_std_int);

    cpp11::sexp res_std_dbl = Rf_allocVector(REALSXP, (!IsInteger && !IsGmp)
                                                 ? numResults : 0);
    double* res_dbl = REAL(res_std_dbl);

    const int bigNumResults = (IsGmp) ? numResults : 0;
    auto myVec = FromCpp14::make_unique<mpz_t[]>(bigNumResults);

    for (int i = 0; i < bigNumResults; ++i) {
        mpz_init(myVec[i]);
    }

    RankResults(myVec.get(), res_int, res_dbl, idx, myReps,
                rankFun, m, n, numResults, IsGmp, IsInteger);

    if (IsInteger) {
        return res_std_int;
    } else if (IsGmp) {
        std::size_t size = intSize;
        std::vector<std::size_t> mySizes(numResults);

        for (std::size_t i = 0; i < numResults; ++i) { // adding each bigint's needed size
            const std::size_t tempSize = intSize * (
                2 + (mpz_sizeinbase(myVec[i], 2) + numb - 1) / numb
            );
            size += tempSize;
            mySizes[i] = tempSize;
        }

        cpp11::sexp res_bigz = Rf_allocVector(RAWSXP, size);

        char* rPos = (char*) (RAW(res_bigz));
        ((int*) (rPos))[0] = numResults; // first int is vector-size-header

        // current position in rPos[] (starting after vector-size-header)
        std::size_t posPos = intSize;

        for (std::size_t i = 0; i < numResults; ++i) {
            posPos += myRaw(&rPos[posPos], myVec[i], mySizes[i]);
        }

        res_bigz.attr("class") = "bigz";
        MpzClearVec(myVec.get(), numResults);
        return res_bigz;
    } else {
        return res_std_dbl;
    }
}
