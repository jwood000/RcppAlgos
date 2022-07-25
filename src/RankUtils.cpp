#include "ImportExportMPZ.h"
#include "SetUpUtils.h"
#include <unordered_map>
#include <algorithm>

constexpr std::size_t numb = 8 * intSize;

void SetUpRank(SEXP RIdx, SEXP Rv, SEXP RisRep, SEXP RFreqs, SEXP Rm,
               std::vector<int> &idx, std::vector<int> &freqs,
               std::vector<int> &myReps, VecType &myType, int &n, int &m,
               bool IsComb, bool &IsMult, bool &IsRep) {

    SetType(myType, Rv);
    CleanConvert::convertPrimitive(Rm, m, VecType::Integer, "m");
    CleanConvert::convertVector(RIdx, idx, VecType::Integer, "idx");

    for (auto &i: idx) {
        --i;
    }

    if (IsComb) {
        for (auto it = idx.begin(); it != idx.end(); it += m) {
            std::sort(it, it + m);
        }
    }

    n = GetLength(Rv, myType);
    SetFreqsAndM(myReps, freqs, RFreqs, Rm, n, m, IsMult, IsRep);

    if (IsMult) {
        // See comment in SetUpUtils.cpp in the SetFinalValues function
        if (n != static_cast<int>(myReps.size())) {
            cpp11::stop("The length of freqs must equal the length of v");
        }

        if (m > static_cast<int>(freqs.size())) {
            cpp11::stop("The input width is too large for the given freqs");
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
}

SEXP MpzReturn(mpz_t *myVec, int numResults) {
    std::size_t size = intSize;
    std::vector<std::size_t> mySizes(numResults);

    for (int i = 0; i < numResults; ++i) { // adding each bigint's needed size
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

    for (int i = 0; i < numResults; ++i) {
        posPos += myRaw(&rPos[posPos], myVec[i], mySizes[i]);
    }

    res_bigz.attr("class") = "bigz";
    MpzClearVec(myVec, numResults);
    return res_bigz;
}
