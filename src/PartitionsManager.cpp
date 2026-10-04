#include "cpp11/protect.hpp"

#include "Partitions/CompositionsDistinct.h"
#include "Partitions/PartitionsMultiset.h"
#include "Partitions/PartitionsDistinct.h"
#include "Partitions/PartitionsTypes.h"
#include "Partitions/CompositionsRep.h"
#include "Partitions/PartitionsRep.h"
#include <algorithm> // std::find
#include <numeric>
#include "RMatrix.h"

int PartsStdManager(
    int* mat, std::vector<int> &z, int width, int lastElem,
    int lastCol, int nRows, PartitionType ptype, int zeroBudget
) {

    switch (ptype) {
        case PartitionType::NoSolution:
            // Do Nothing
            return 1;

        case PartitionType::LengthOne:
            if (nRows) mat[0] = z.front();
            return 1;

            // ----- PartsRep -----
        case PartitionType::RepStdAll:
        case PartitionType::RepNoZero:
        case PartitionType::RepShort:
            return PartsRep(mat, z, width, lastElem, lastCol, nRows);

            // ----- PartsDistinct -----
        case PartitionType::DistinctStdAll:
        case PartitionType::DistinctMZ:
        case PartitionType::DistinctOneZero:
        case PartitionType::DistinctNoZero:
            return PartsDistinct(mat, z, width, lastElem, lastCol, nRows);

            // ----- CompsRep -----
        case PartitionType::CompRepNoZero:
        case PartitionType::CompRepZero:
            return CompsRep<1>(mat, z, width, nRows);

        case PartitionType::CompRepWeak:
            return CompsRep<0>(mat, z, width, nRows);

            // ----- PartsPermRep -----
        case PartitionType::PrmRepPartNoZ:
        case PartitionType::PrmRepPart:
            return PartsPermRep(mat, z, width, lastElem, lastCol, nRows);

            // ----- PartsPermDistinct -----
        case PartitionType::PrmDstPartNoZ:
        case PartitionType::PrmDstPrtOneZ:
            return PartsPermDistinct(mat, z, width, lastElem, lastCol, nRows);

        case PartitionType::PrmDstPartMZ:
            return PartsPermZeroDistinct(
                mat, z, width, lastElem, lastCol, nRows
            );

            // ----- CompsDistinct -----
        case PartitionType::CompDistinctNoZero:
        case PartitionType::CompDistinctZero:
            return CompsDistinct(mat, z, width, nRows, false, zeroBudget);

        case PartitionType::CompDistinctWeak:
        case PartitionType::CompDistinctMZWeak:
            return CompsDistinct(mat, z, width, nRows, true, zeroBudget);

        default:
            cpp11::stop("Case not supported");
    }
}

template <typename T>
int PartsGenManager(T* mat, const std::vector<T> &v,
                    const std::vector<int> &Reps, std::vector<int> &z,
                    int width, int lastElem, int lastCol, int nRows,
                    PartitionType ptype, int zeroBudget) {

    // PrmRepPart, PrmDstPrtOneZ, DistinctStdAll, DistinctOneZero, RepStdAll,
    // RepShort and CompRepWeak shouldn't happen. In any mapped case that
    // would trigger PartsGenManager, the zero would be mapped to some
    // non-zero integer.

    switch (ptype) {
        case PartitionType::NoSolution:
            // Do Nothing
            return 1;

        case PartitionType::LengthOne:
            if (nRows) mat[0] = v[z.front()];
            return 1;

            // ----- PartsGenRep -----
        // case PartitionType::RepShort:
        // case PartitionType::RepStdAll:
        case PartitionType::RepNoZero:
        case PartitionType::RepCapped:
            return PartsGenRep(mat, v, z, width, lastElem, lastCol, nRows);

            // ----- PartsGenDistinct -----
        // case PartitionType::DistinctStdAll:
        // case PartitionType::DistinctOneZero:
        case PartitionType::DistinctNoZero:
        case PartitionType::DistinctCapped:
        case PartitionType::DistinctCappedMZ:
        case PartitionType::DistinctMZ:
            return PartsGenDistinct(mat, v, z, width, lastElem, lastCol, nRows);

            // ----- PartsGenMultiset -----
        case PartitionType::Multiset:
            return PartsGenMultiset(
                mat, v, Reps, z, width, lastElem, lastCol, nRows
            );

            // ----- CompsGenRep -----
        case PartitionType::CompRepZero:
            return CompsGenRep<1>(mat, v, z, width, nRows);

        // case PartitionType::CompRepWeak:
        case PartitionType::CompRepNoZero:
            return CompsGenRep<0>(mat, v, z, width, nRows);

        case PartitionType::CompRepCapped:
        case PartitionType::CompRepWeakCap:
            return CompsGenRepCapped<0>(mat, v, z, width, nRows);

        case PartitionType::CompRepCapZero:
            return CompsGenRepCapped<1>(mat, v, z, width, nRows);

            // ----- PartsGenPerm -----
        // case PartitionType::PrmRepPart:
        case PartitionType::PrmRepPartNoZ:
            return PartsStdPermRep(mat, v, z, width, lastElem, lastCol, nRows);

        case PartitionType::PrmRepCapped:
            return PartsGenPermRep(mat, v, z, width, nRows);

        // case PartitionType::PrmDstPrtOneZ:
        case PartitionType::PrmDstPartNoZ:
        case PartitionType::PrmDstPrtCap:
            return PartsGenPermDistinct(
                mat, v, z, width, lastElem, lastCol, nRows
            );

        case PartitionType::PrmDstPartMZ:
        case PartitionType::PrmDstPrtCapMZ:
            return PartsGenPermZeroDistinct(
                mat, v, z, width, lastElem, lastCol, nRows
            );

        case PartitionType::PrmMultiset:
            return PartsGenPermMultiset(
                mat, v, Reps, z, width, lastElem, lastCol, nRows
            );

            // ----- CompsGenDistinct -----
        case PartitionType::CompDistinctNoZero:
        case PartitionType::CompDistinctZero:
        case PartitionType::CompDistinctCapped:
        case PartitionType::CompDistinctCapMZ:
            return CompsGenDistinct(mat, v, z, width, nRows, false, zeroBudget);

        case PartitionType::CompDistinctWeak:
        case PartitionType::CompDistinctCapWeak:
        case PartitionType::CompDistinctMZWeak:
        case PartitionType::CompDistinctCapMZWeak:
            return CompsGenDistinct(mat, v, z, width, nRows, true, zeroBudget);

        default:
            cpp11::stop("Case not supported");
    }
}

template <typename T>
int PartsGenManager(std::vector<T> &partsVec, const std::vector<T> &v,
                    const std::vector<int> &Reps, std::vector<int> &z,
                    int width, int nRows, PartitionType ptype) {

    // LengthOne, RepCapped, DistinctCapped, and PrmDstPrtCap should not happen.
    // In the past, we used a different algorithm for counting these types
    // of partitions that relied on allocating a large array that represented
    // three dimensions. We moved away from this to a more memory-friendly
    // data structure that should make these cases impossible.

    switch (ptype) {
        case PartitionType::NoSolution: {
            // Do Nothing
            return 1;
        } case PartitionType::LengthOne: {
            if (nRows) partsVec.push_back(v[z.front()]);
            return 1;
        } case PartitionType::DistinctCapped: {
            return PartsGenDistinct(partsVec, v, z, width, nRows, true);
        } case PartitionType::PrmRepCapped : {
            return PartsGenRep(partsVec, v, z, width, nRows, false);
        } case PartitionType::RepCapped: {
            return PartsGenRep(partsVec, v, z, width, nRows, true);
        } case PartitionType::PrmDstPrtCap: {
            return PartsGenDistinct(partsVec, v, z, width, nRows, false);
        } default: {
            cpp11::stop("Case not supported");
        }
    }
}

int PartsStdParallel(RcppParallel::RMatrix<int> &mat, std::vector<int> &z,
                     int strt, int width, int lastElem, int lastCol, int nRows,
                     PartitionType ptype, int zeroBudget) {

    switch (ptype) {
        case PartitionType::NoSolution:
            // Do Nothing
            return 1;

            // ----- PartsRep -----
        case PartitionType::RepStdAll:
        case PartitionType::RepNoZero:
        case PartitionType::RepShort:
            return PartsRep(mat, z, strt, width, lastElem, lastCol, nRows);

            // ----- PartsDistinct -----
        case PartitionType::DistinctStdAll:
        case PartitionType::DistinctMZ:
        case PartitionType::DistinctOneZero:
        case PartitionType::DistinctNoZero:
            return PartsDistinct(mat, z, strt, width, lastElem, lastCol, nRows);

            // ----- CompsRep -----
        case PartitionType::CompRepNoZero:
        case PartitionType::CompRepZero:
            return CompsRep<1>(mat, z, strt, width, nRows);

        case PartitionType::CompRepWeak:
            return CompsRep<0>(mat, z, strt, width, nRows);

            // ----- CompsDistinct -----
        case PartitionType::CompDistinctNoZero:
        case PartitionType::CompDistinctZero:
            return CompsDistinct(
                mat, z, strt, width, nRows, false, zeroBudget
            );

        case PartitionType::CompDistinctWeak:
        case PartitionType::CompDistinctMZWeak:
            return CompsDistinct(
                mat, z, strt, width, nRows, true, zeroBudget
            );

        default:
            cpp11::stop("Case not supported");
    }
}

template <typename T>
int PartsGenParallel(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                     const std::vector<int> &Reps, std::vector<int> &z,
                     int strt, int width, int lastElem, int lastCol,
                     int nRows, PartitionType ptype, int zeroBudget) {

    // DistinctStdAll, DistinctOneZero, RepStdAll, RepShort and CompRepWeak
    // shouldn't happen just as in the non-Parallel case. In any mapped case
    // that would trigger PartsGenManager, the zero would be mapped to some
    // non-zero integer.

    switch (ptype) {
        case PartitionType::NoSolution:
            // Do Nothing
            return 1;

            // ----- PartsGenRep -----
        // case PartitionType::RepStdAll:
        // case PartitionType::RepShort:
        case PartitionType::RepNoZero:
        case PartitionType::RepCapped:
            return PartsGenRep(
                mat, v, z, strt, width, lastElem, lastCol, nRows
            );

            // ----- PartsGenDistinct -----
        // case PartitionType::DistinctStdAll:
        // case PartitionType::DistinctOneZero:

        case PartitionType::DistinctNoZero:
        case PartitionType::DistinctCapped:
        case PartitionType::DistinctCappedMZ:
        case PartitionType::DistinctMZ:
            return PartsGenDistinct(
                mat, v, z, strt, width, lastElem, lastCol, nRows
            );

            // ----- PartsGenMultiset ----- Coming Soon!
        // case PartitionType::Multiset:
        //     return PartsGenMultiset(
        //         mat, v, Reps, z, strt, width, lastElem, lastCol, nRows
        //     );

            // ----- CompsGenRep -----
        case PartitionType::CompRepZero:
            return CompsGenRep<1>(mat, v, z, strt, width, nRows);

        // case PartitionType::CompRepWeak:
        case PartitionType::CompRepNoZero:
            return CompsGenRep<0>(mat, v, z, strt, width, nRows);

        case PartitionType::CompRepCapped:
        case PartitionType::CompRepWeakCap:
            return CompsGenRepCapped<0>(mat, v, z, strt, width, nRows);

        case PartitionType::CompRepCapZero:
            return CompsGenRepCapped<1>(mat, v, z, strt, width, nRows);

            // ----- CompsGenDistinct -----
        case PartitionType::CompDistinctNoZero:
        case PartitionType::CompDistinctZero:
        case PartitionType::CompDistinctCapped:
        case PartitionType::CompDistinctCapMZ:
            return CompsGenDistinct(
                mat, v, z, strt, width, nRows, false, zeroBudget
            );

        case PartitionType::CompDistinctWeak:
        case PartitionType::CompDistinctCapWeak:
        case PartitionType::CompDistinctMZWeak:
        case PartitionType::CompDistinctCapMZWeak:
            return CompsGenDistinct(
                mat, v, z, strt, width, nRows, true, zeroBudget
            );

        default:
            cpp11::stop("Case not supported");
    }
}

template int PartsGenManager(int*, const std::vector<int>&,
                             const std::vector<int>&, std::vector<int>&,
                             int, int, int, int, PartitionType, int);
template int PartsGenManager(double*, const std::vector<double>&,
                             const std::vector<int>&, std::vector<int>&,
                             int, int, int, int, PartitionType, int);

template int PartsGenManager(std::vector<int>&, const std::vector<int>&,
                              const std::vector<int>&, std::vector<int>&,
                              int, int, PartitionType);
template int PartsGenManager(std::vector<double>&, const std::vector<double>&,
                              const std::vector<int>&, std::vector<int>&,
                              int, int, PartitionType);

template int PartsGenParallel(
    RcppParallel::RMatrix<int>&, const std::vector<int>&,
    const std::vector<int>&, std::vector<int>&, int, int,
    int, int, int, PartitionType, int
);
template int PartsGenParallel(
    RcppParallel::RMatrix<double>&, const std::vector<double>&,
    const std::vector<int>&, std::vector<int>&, int, int,
    int, int, int, PartitionType, int
);
