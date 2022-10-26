#include "Partitions/PartitionsManager.h"
#include "Partitions/NthPartition.h"
#include <thread>

void StandardPartitions(int* mat, std::vector<int> &z, PartitionType ptype,
                        double lower, mpz_class &lowerMpz, int nCols,
                        int width, int nRows, int nThreads, int lastCol,
                        int lastElem, int tar, int strtLen, int cap,
                        bool IsRep, bool IsMult, bool IsGmp, bool IsComb,
                        bool includeZero, bool IsComp, bool zero_spesh) {

    if (nThreads > 1 && (IsComb || IsComp)) {
        RcppParallel::RMatrix<int> parMat(mat, nRows, nCols);
        std::vector<std::thread> threads;

        const int stepSize = nRows / nThreads;
        int nextStep = stepSize;
        int step = 0;

        const nthPartsPtr nthPartFun = GetNthPartsFunc(ptype, IsGmp, IsComp);
        std::vector<std::vector<int>> zs(nThreads, z);

        for (int j = 0; j < (nThreads - 1);
             ++j, step += stepSize, nextStep += stepSize) {

            threads.emplace_back(
                std::cref(PartsStdParallel), std::ref(parMat),
                std::ref(zs[j]), step, width, lastElem, lastCol,
                nextStep, IsRep, IsComp, zero_spesh
            );

            if (IsGmp) {
                lowerMpz += stepSize;
            } else {
                lower    += stepSize;
            }

            zs[j + 1] = nthPartFun(tar, width, cap,
                                   strtLen, lower, lowerMpz);

            if (!includeZero) {
                for (auto &z_i: zs[j + 1]) {
                    ++z_i;
                }
            }
        }

        threads.emplace_back(
            std::cref(PartsStdParallel), std::ref(parMat),
            std::ref(zs.back()), step, width, lastElem,
            lastCol, nRows, IsRep, IsComp, zero_spesh
        );

        for (auto& thr: threads) {
            thr.join();
        }
    } else {
        PartsStdManager(mat, z, width, lastElem, lastCol, nRows,
                        IsComb, IsRep, IsComp, zero_spesh);
    }
}

template <typename T>
void GeneralPartitions(T* mat, const std::vector<T> &v, std::vector<int> &z,
                       const PartDesign &part, double lower,
                       mpz_class &lowerMpz, int nCols, int nRows,
                       int nThreads, int lastCol, int lastElem, int strtLen,
                       int cap, bool IsComb) {

    const bool zero_spesh = part.includeZero && !part.isWeak;

    if (nThreads > 1 && (IsComb || part.isComp)) {
        RcppParallel::RMatrix<T> parMat(mat, nRows, nCols);
        std::vector<std::thread> threads;

        const int stepSize = nRows / nThreads;
        int nextStep = stepSize;
        int step = 0;

        const nthPartsPtr nthPartFun = GetNthPartsFunc(
            part.ptype, part.isGmp, part.isComp
        );
        std::vector<std::vector<int>> zs(nThreads, z);

        for (int j = 0; j < (nThreads - 1);
             ++j, step += stepSize, nextStep += stepSize) {

            threads.emplace_back(
                std::cref(PartsGenParallel<T>), std::ref(parMat),
                std::cref(v), std::ref(zs[j]), step, part.width, lastElem,
                lastCol, nextStep, part.isRep, part.isComp, zero_spesh
            );

            if (part.isGmp) {
                lowerMpz += stepSize;
            } else {
                lower    += stepSize;
            }

            zs[j + 1] = nthPartFun(part.mapTar, part.width,
                                   cap, strtLen, lower, lowerMpz);
        }

        threads.emplace_back(
            std::cref(PartsGenParallel<T>), std::ref(parMat), std::cref(v),
            std::ref(zs.back()), step, part.width, lastElem, lastCol, nRows,
            part.isRep, part.isComp, zero_spesh
        );

        for (auto& thr: threads) {
            thr.join();
        }
    } else {
        PartsGenManager(mat, v, z, part.width, lastElem, lastCol, nRows,
                        IsComb, part.isRep, part.isComp, zero_spesh);
    }
}

template void GeneralPartitions(int*, const std::vector<int>&,
                                std::vector<int>&, const PartDesign&,
                                double, mpz_class&, int, int, int, int,
                                int, int, int, bool);

template void GeneralPartitions(double*, const std::vector<double>&,
                                std::vector<int>&, const PartDesign&,
                                double, mpz_class&, int, int, int, int,
                                int, int, int, bool);
