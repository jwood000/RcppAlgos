#include "cpp11/protect.hpp"

#include "Partitions/CompositionsDistinctUtils.h"
#include "Partitions/NextComposition.h"
#include "Partitions/NextPartition.h"
#include <algorithm>  // std::find
#include <numeric>

void NextDistinctPart(std::vector<int> &z, int &boundary,
                      int &edge, int &tarDiff, int lastCol) {

    if (z[boundary] - z[edge] != tarDiff) {
        boundary = edge + 1;
    }

    ++z[edge];
    --z[boundary];

    for (int edge2 = z[edge] + boundary - edge;
         boundary < lastCol; ++boundary, ++edge2) {

        z[lastCol] += (z[boundary] - edge2);
        z[boundary] = edge2;
    }

    while (boundary > 1 && (z[boundary] - z[boundary - 1]) < 2) {
        --boundary;
    }

    edge = boundary - 1;
    tarDiff = 3;

    while (edge > 0 && (z[boundary] - z[edge]) < tarDiff) {
        --edge;
        ++tarDiff;
    }
}

void NextRepPart(std::vector<int> &z, int &boundary,
                 int &edge, int lastCol) {

    if (z[boundary] - z[edge] != 2) {
        boundary = edge + 1;
    }

    ++z[edge];
    --z[boundary];

    for (int edge2 = z[edge]; boundary < lastCol; ++boundary) {
        z[lastCol] += (z[boundary] - edge2);
        z[boundary] = edge2;
    }

    const int currMax = z[boundary];

    while (boundary > 1 && z[boundary - 1] == currMax) {
        --boundary;
    }

    edge = boundary - 1;
    const int edgeTest = z[boundary] - 2;

    while (edge > 0 && edgeTest < z[edge]) {
        --edge;
    }
}

// BndDecrementPossible, VtxDecrementPossible, EdgeIncrementPossible, GetPivotExtr
// and PivotDecrementPossible are all helper functions for PartitionsMultiset
// ************************** Start Helper Functions ***************************
inline bool BndDecrementPossible(const std::vector<int> &rpsCnt,
                                 const std::vector<int> &z,
                                 int boundary) {

    if (boundary > 1) {
        if ((z[boundary] - z[boundary - 1]) < 2) {
            return ((z[boundary] != z[boundary - 1]) ?
                        !rpsCnt[z[boundary] - 1] : true);
        } else {
            return false;
        }
    } else {
        return false;
    }
}

inline bool VtxDecrementPossible(const std::vector<int> &rpsCnt,
                                 const std::vector<int> &z,
                                 int lastCol, int vertex, int edge) {

    if (vertex < lastCol) {
        const int myDiff = z[vertex] - z[edge];

        if (myDiff < 2) {
            return true;
        } else  {
            return ((myDiff != 2 || rpsCnt[z[edge] + 1] <= 1) &&
                        (
                            myDiff <= 2 ||
                            !rpsCnt[z[edge] + 1] ||
                            !rpsCnt[z[vertex] - 1]
                        )
                    );
        }
    } else {
        return false;
    }
}

inline bool EdgeIncrementPossible(const std::vector<int> &rpsCnt,
                                  const std::vector<int> &z,
                                  int edge, int boundary) {

    if (edge) {
        const int myDiff = z[boundary] - z[edge];

        if (myDiff < 2) {
            return true;
        } else {
            return (myDiff != 2 || rpsCnt[z[edge] + 1] < 2) &&
                   (myDiff <= 2 || !rpsCnt[z[edge] + 1]);
        }
    } else {
        return false;
    }
}

inline int GetPivotExtr(const std::vector<int> &rpsCnt,
                        const std::vector<int> &z,
                        int lastCol, int lastElem) {

    int res = lastCol - 1;

    while (res > 0 && z[res] == lastElem) {
        --res;
    }

    while (res > 0 && !rpsCnt[z[res] + 1]) {
        --res;
    }

    return res;
}

inline bool PivotDecrementPossible(const std::vector<int> &rpsCnt,
                                   const std::vector<int> &z,
                                   int lastElem, int pivot, int vertex) {
    if (pivot > vertex) {
        if (z[pivot] == lastElem) {
            return true;
        } else {
            return !rpsCnt[z[pivot] + 1];
        }
    } else {
        return false;
    }
}

// ********** End PartitionsMultiset Helper Functions *****************

void GetLastPart(int* mat, std::vector<int> &z, int m, int nRows) {
    for (std::size_t k = 0; k < static_cast<std::size_t>(m); ++k) {
        mat[static_cast<std::size_t>(nRows - 1) +
            static_cast<std::size_t>(nRows) * k] = z[k];
    }
}

void PrepareMultisetPart(std::vector<int> &rpsCnt,
                         const std::vector<int> &z, int &b,
                         int &p, int &e, int lastCol, int lastElem) {

    b = lastCol;

    for (const auto ind: z) {
        --rpsCnt[ind];
    }

    while (BndDecrementPossible(rpsCnt, z, b)) {
        --b;
    }

    p = (z[lastCol] < lastElem) ? lastCol :
        GetPivotExtr(rpsCnt, z, lastCol, lastElem);

    e = b - 1;
    e = b > 0 ? b - 1 : 0;

    while (EdgeIncrementPossible(rpsCnt, z, e, b)) {
        --e;
    }
}

void PrepareRepPart(const std::vector<int> &z, int &boundary,
                    int &pivot, int &edge, int lastElem, int lastCol) {

    // smallest index such that z[boundary] == currMax
    boundary = lastCol;

    while (boundary > 1 && z[boundary - 1] == z[lastCol]) {
        --boundary;
    }

    // pivot is the greatest index that can be incremented.
    // We know that if z[boundary] < lastElem ==>> pivot = lastCol
    pivot = (z[boundary] < lastElem) ? lastCol : boundary - 1;

    // edge is the greatest index such that z[boundary] - z[edge] >= 2
    // This is the index that will be be used as a starting point
    // to determine the next combination that meets the criteria
    edge = boundary > 0 ? boundary - 1 : 0;
    int edgeTest = z[boundary] - 2;

    while (edge > 0 && edgeTest < z[edge]) {
        --edge;
    }
}

void PrepareDistinctPart(const std::vector<int> &z, int &boundary,
                         int &pivot, int &edge, int &tarDiff,
                         int lastElem, int lastCol) {

    // Largest index such that z[boundary] - z[boundary - 1] > 1
    // This is the index that will be decremented at the same
    // time edge is incremented. See below.
    boundary = lastCol;

    while (boundary > 1 && (z[boundary] - z[boundary - 1]) < 2) {
        --boundary;
    }

    // pivot is the greatest index that can be incremented...
    // Either z[pivot + 1] - z[pivot] > 1 or if z[lastCol] < lastElem
    // pivot = lastCol since incrementing z[lastCol] is possible
    pivot = (z[lastCol] < lastElem) ? lastCol : boundary - 1;

    // edge is the greatest index such that when incremented
    // the result will be at least one less than its neighbor
    // even if its neighbor is decremented
    edge = boundary > 0 ? boundary - 1 : 0;
    tarDiff = 3;

    while (edge > 0 && (z[boundary] - z[edge]) < tarDiff) {
        --edge;
        ++tarDiff;
    }
}

void NextMultisetGenPart(std::vector<int> &rpsCnt,
                         std::vector<int> &z, int &e, int &b,
                         int &p, int lastCol, int lastElem) {

    // vertex is the smallest index greater than edge that can be decremented
    int v = e + 1;

    // Find next vertex v that can't be decremented
    while (VtxDecrementPossible(rpsCnt, z, lastCol, v, e)) {
        ++v;
    }

    // Perform transfer between edge and vertex
    ++rpsCnt[z[e]];
    ++z[e];
    --rpsCnt[z[e]];

    ++rpsCnt[z[v]];
    --z[v];
    --rpsCnt[z[v]];

    // Update boundary and pivot if vertex hits boundary
    if (v == b) {
        if (b < lastCol) {
            ++b;
        }

        while (BndDecrementPossible(rpsCnt, z, b)) {
            --b;
        }

        p = (z[lastCol] < lastElem) ? lastCol :
            GetPivotExtr(rpsCnt, z, lastCol, lastElem);
    }

    // Skip vertex forward over zones that are equal to previous
    while (
        (v < lastCol) &&
        (
            z[v] == z[v - 1] ||
            z[v] == z[e] ||
            (
                z[v] - z[v - 1] == 1 &&
                !rpsCnt[z[v - 1]]
            )
        )) {

        ++v;
    }

    // Balance values between vertex and pivot where possible
    while (v < p && rpsCnt[z[v] - 1] && rpsCnt[z[p] + 1]) {
        ++rpsCnt[z[v]];
        --z[v];
        --rpsCnt[z[v]];

        ++rpsCnt[z[p]];
        ++z[p];
        --rpsCnt[z[p]];

        while (z[v] == z[v - 1] ||
               (z[v] - z[v - 1] == 1 && !rpsCnt[z[v - 1]])) {
            ++v;
        }

        while (PivotDecrementPossible(rpsCnt, z, lastElem, p, v)) {
            --p;
        }
    }

    b = p;

    // Move boundary right over zones that are equal or "close and movable"
    while (
        (b < lastCol) &&
        (
            (z[b] == z[b + 1]) ||
            (
                z[b + 1] > z[b] &&
                (
                    rpsCnt[z[b + 1] - 1] ||
                    rpsCnt[z[b + 1]]
                )
            )
        )) {

        ++b;
    }

    // Refine boundary
    while (BndDecrementPossible(rpsCnt, z, b)) {
        --b;
    }

    e = b - 1;

    // Move edge left as long as it's incrementable
    while (EdgeIncrementPossible(rpsCnt, z, e, b)) {
        --e;
    }
}

void NextRepGenPart(std::vector<int> &z, int &boundary, int &edge,
                    int &pivot, int lastCol, int lastElem) {

    int vertex = (z[boundary] - z[edge] == 2) ? boundary : edge + 1;

    ++z[edge];
    --z[vertex];

    // If vertex == boundary, handle boundary update and pivot assignment
    if (vertex == boundary) {
        if (boundary < lastCol) {
            ++boundary;
        }

        int currMax = z[boundary];

        // Move boundary left as long as values equal currMax
        while (boundary > 1 && z[boundary - 1] == currMax) {
            --boundary;
        }

        pivot = (z[boundary] < lastElem) ? lastCol : boundary - 1;
    } else if (z[vertex] == z[edge]) {
        ++vertex;
    }

    // Main loop: balance values between vertex and pivot
    while (vertex < pivot) {
        int distVert = z[vertex] - z[edge];
        int distPivot = lastElem - z[pivot];

        // Adjust indices
        if (distVert == distPivot) {
            z[vertex] -= distVert;
            z[pivot] += distVert;

            ++vertex;
            --pivot;
        } else if (distVert < distPivot) {
            z[vertex] -= distVert;
            z[pivot] += distVert;

            ++vertex;
        } else {
            z[vertex] -= distPivot;
            z[pivot] += distPivot;

            --pivot;
        }
    }

    boundary = pivot;

    // Recheck and adjust boundary
    if (boundary < lastCol && z[boundary] < z[boundary + 1]) {
        ++boundary;
    }

    int currMax = z[boundary];

    while (boundary > 1 && z[boundary - 1] == currMax) {
        --boundary;
    }

    // Compute edge from boundary
    edge = boundary - 1;
    int edgeTest = z[boundary] - 2;

    // Move edge left while condition holds
    while (edge > 0 && edgeTest < z[edge]) {
        --edge;
    }
}

void NextDistinctGenPart(std::vector<int> &z, int &boundary,
                         int &edge, int &pivot, int &tarDiff,
                         int lastCol, int lastElem) {

    int vertex = edge + 1;
    tarDiff = 3;

    // Find the vertex where the difference from z[edge] is >= tarDiff
    while (vertex < lastCol && (z[vertex] - z[edge]) < tarDiff) {
        ++vertex;
        ++tarDiff;
    }

    ++z[edge];
    --z[vertex];

    if (vertex == boundary) {
        if (boundary < lastCol) {
            ++boundary;
        }

        // Move boundary left until the difference between adjacent values
        // is at least 2
        while (boundary > 1 && (z[boundary] - z[boundary - 1]) < 2) {
            --boundary;
        }

        pivot = (z[lastCol] < lastElem) ? lastCol : boundary - 1;
    }

    // Redistribution logic if vertex is before boundary
    if (vertex < boundary) {
        if (z[vertex] - z[vertex - 1] == 1) {
            ++vertex;
        }

        while (vertex < pivot) {
            --z[vertex];
            ++z[pivot];

            if (z[vertex] - z[vertex - 1] == 1) {
                ++vertex;
            }

            // Move pivot left if maxed out or next is too close
            if (z[pivot] == lastElem ||
                (pivot < lastCol && z[pivot + 1] - z[pivot] == 1)) {

                --pivot;
            }
        }

        boundary = pivot;

        if (boundary < lastCol && z[boundary + 1] - z[boundary] > 1) {
            ++boundary;
        }
    }

    // Update edge based on new boundary
    edge = boundary - 1;
    tarDiff = 3;

    // Find new edge where the difference with boundary is >= tarDiff
    while (edge > 0 && (z[boundary] - z[edge]) < tarDiff) {
        --edge;
        ++tarDiff;
    }
}

void NextDistinct(std::vector<int> &rpsCnt, std::vector<int> &z, int &e,
                  int &b, int &p, int &tarDiff, int lastCol, int lastElem) {
    NextDistinctPart(z, b, e, tarDiff, lastCol);
}

void NextRep(std::vector<int> &rpsCnt, std::vector<int> &z, int &e,
             int &b, int &p, int &tarDiff, int lastCol, int lastElem) {
    NextRepPart(z, b, e, lastCol);
}

void NextMultisetGen(std::vector<int> &rpsCnt, std::vector<int> &z,
                     int &e, int &b, int &p, int &tarDiff, int lastCol,
                     int lastElem) {
    NextMultisetGenPart(rpsCnt, z, e, b, p, lastCol, lastElem);
}

void NextRepGen(std::vector<int> &rpsCnt, std::vector<int> &z, int &e,
                int &b, int &p, int &tarDiff, int lastCol, int lastElem) {
    NextRepGenPart(z, b, e, p, lastCol, lastElem);
}

void NextDistinctGen(std::vector<int> &rpsCnt, std::vector<int> &z,
                     int &e, int &b, int &p, int &tarDiff, int lastCol,
                     int lastElem) {
    NextDistinctGenPart(z, b, e, p, tarDiff, lastCol, lastElem);
}

void NextRepCompZero(std::vector<int> &rpsCnt,
                     std::vector<int> &z, int &e, int &b, int &p,
                     int &tarDiff, int lastCol, int lastElem) {
    NextCompositionRep<0>(z, lastCol);
}

void NextRepCompOne(std::vector<int> &rpsCnt,
                    std::vector<int> &z, int &e, int &b, int &p,
                    int &tarDiff, int lastCol, int lastElem) {
    NextCompositionRep<1>(z, lastCol);
}

void NextDistinctComp(std::vector<int> &complement,
                      std::vector<int> &z, int &i1, int &i2, int &myMax,
                      int &target, int lastCol, int lastIdx) {

    if (complement.size() < 2) {
        std::next_permutation(z.begin(), z.end());
        return;
    }

    std::vector<int> idx;
    std::vector<int> tailSum;
    lastIdx = complement.size() - 1;

    NextCompositionDistinct(
        z, complement, idx, tailSum, i1, i2, myMax, lastCol, lastIdx, target
    );
}

void NextDistMZNotWeakComp(std::vector<int> &complement,
                           std::vector<int> &z, int &i1, int &i2, int &myMax,
                           int &target, int lastCol, int lastIdx) {

    if (complement.size() < 2) {
        std::next_permutation(z.begin(), z.end());
        return;
    }

    std::vector<int> idx;
    std::vector<int> tailSum;
    lastIdx = complement.size() - 1;

    int nz = 0;
    const int z_size = z.size();

    while (nz < z_size && z[nz] == 0) {
        ++nz;
    }

    if (nz == z_size) {
        cpp11::stop("This should not happen!");
    }

    int k = nz + 1;

    for (; k < z_size; ++k) {
        if (z[k] > z[k - 1]) {
            break;
        }
    }

    // means descending
    if (k == z_size && nz) {
        const int non_zero_size = z_size - nz;

        // theoretical maximum if there are no constraints
        int max_possval = target - static_cast<int>(
            (non_zero_size * (non_zero_size - 1)) / 2
        );

        int z_max = *std::max_element(z.begin(), z.end());
        int cmp_max = *std::max_element(complement.begin(), complement.end());
        int cap = std::max(cmp_max, z_max);
        cap = std::min(max_possval, cap);

        if (z[nz] == cap && IsMaximizedGreedySuffix(z, target, nz)) {
            if (nz > 1) z.erase(z.begin(), z.begin() + (nz - 1));

            if (cap == max_possval) {
                std::iota(z.begin(), z.end(), 1);
                z.back() = target - std::accumulate(z.begin(), z.end() - 1, 0);
            } else {
                std::vector<int> allowed(cap);
                std::iota(allowed.begin(), allowed.end(), 1);
                GetFirstPartitionDistinct(allowed, z, target, z.size(), cap);
                for (auto &z_i: z) ++z_i;
            }

            if (nz > 1) z.insert(z.begin(), nz - 1, 0);
            int zeroBudget = std::count(z.cbegin(), z.cend(), 0);

            CompsDistinctSetup(
                z, complement, target, i1, i2, myMax, cap, false, zeroBudget
            );

            return;
        }
    }

    NextCompositionDistinct(
        z, complement, idx, tailSum, i1, i2, myMax, lastCol, lastIdx, target
    );
}

void EmptyReturn(std::vector<int> &rpsCnt, std::vector<int> &z, int &e,
             int &b, int &p, int &tarDiff, int lastCol, int lastElem) {
    // Do nothing
}

nextPartsPtr GetNextPartsPtr(PartitionType ptype, ConstraintType ctype) {

    if (ctype == ConstraintType::PartStandard) {
        switch (ptype) {
            case PartitionType::LengthOne:
            case PartitionType::DstctStdAll:
            case PartitionType::DstctMultiZero:
            case PartitionType::DstctOneZero:
            case PartitionType::DstctNoZero:
                return(nextPartsPtr(NextDistinct));

            case PartitionType::RepStdAll:
            case PartitionType::RepNoZero:
            case PartitionType::RepShort:
                return(nextPartsPtr(NextRep));

            case PartitionType::CompRepNoZero:
            case PartitionType::CmpRpZroNotWk:
                return(nextPartsPtr(NextRepCompOne));

            case PartitionType::CompRepWeak:
                return(nextPartsPtr(NextRepCompZero));

            case PartitionType::CmpDstctWeak:
            case PartitionType::CmpDstCapWeak:
            case PartitionType::CmpDstctMZWeak:
            case PartitionType::CmpDstctNoZero:
            case PartitionType::CmpDstctCapped:
            case PartitionType::CmpDstCapMZWeak:
                return(nextPartsPtr(NextDistinctComp));

            case PartitionType::CmpDstctZNotWk:
            case PartitionType::CmpDstCapMZNotWk:
                return(nextPartsPtr(NextDistMZNotWeakComp));

            case PartitionType::NoSolution:
                return(nextPartsPtr(EmptyReturn));

            default:
                // This should not happen
                cpp11::stop(
                    "This should not happen. Please open an issue here:\n\t"
                    "https://github.com/jwood000/RcppAlgos/issues"
                );
        }
    } else {
        switch (ptype) {
            case PartitionType::LengthOne:
            case PartitionType::DstctMultiZero:
            case PartitionType::DstctNoZero:
            case PartitionType::DstctCapped:
            case PartitionType::DstctCappedMZ:
                return(nextPartsPtr(NextDistinctGen));

            case PartitionType::RepNoZero:
            case PartitionType::RepCapped:
                return(nextPartsPtr(NextRepGen));

            case PartitionType::Multiset:
                return(nextPartsPtr(NextMultisetGen));

            case PartitionType::CompRepNoZero:
                return(nextPartsPtr(NextRepCompZero));

            case PartitionType::CmpRpZroNotWk:
                return(nextPartsPtr(NextRepCompOne));

            case PartitionType::CmpDstctWeak:
            case PartitionType::CmpDstCapWeak:
            case PartitionType::CmpDstctMZWeak:
            case PartitionType::CmpDstctNoZero:
            case PartitionType::CmpDstctCapped:
            case PartitionType::CmpDstCapMZWeak:
                return(nextPartsPtr(NextDistinctComp));

            case PartitionType::CmpDstctZNotWk:
            case PartitionType::CmpDstCapMZNotWk:
                return(nextPartsPtr(NextDistMZNotWeakComp));

            case PartitionType::NoSolution:
                return(nextPartsPtr(EmptyReturn));

            default:
                // This should not happen
                cpp11::stop(
                    "This should not happen. Please open an issue here:\n\t"
                    "https://github.com/jwood000/RcppAlgos/issues"
                );
        }
    }
}
