#include "Partitions/PartitionsTypes.h"
#include <algorithm>  // std::find

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

    while (edge && (z[boundary] - z[edge]) < tarDiff) {
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

    while (edge && edgeTest < z[edge]) {
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
    for (int k = 0; k < m; ++k) {
        mat[nRows - 1 + nRows * k] = z[k];
    }
}

void PrepareMultisetPart(std::vector<int> &rpsCnt,
                         const std::vector<int> &z, int &b,
                         int &p, int &e, int lastCol, int lastElem) {

    b = lastCol;

    for (const auto ind: z)
        --rpsCnt[ind];

    while (BndDecrementPossible(rpsCnt, z, b)) {
        --b;
    }

    p = (z[lastCol] < lastElem) ? lastCol :
        GetPivotExtr(rpsCnt, z, lastCol, lastElem);

    e = b - 1;

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
    edge = boundary - 1;
    int edgeTest = z[boundary] - 2;

    while (edge && edgeTest < z[edge]) {
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
    edge = boundary - 1;
    tarDiff = 3;

    while (edge && (z[boundary] - z[edge]) < tarDiff) {
        --edge;
        ++tarDiff;
    }
}

void NextMultisetGenPart(std::vector<int> &rpsCnt,
                         std::vector<int> &z, int &e, int &b,
                         int &p, int lastCol, int lastElem) {

    // vertex is the smallest index greater than edge that can be decremented
    int v = e + 1;

    while (VtxDecrementPossible(rpsCnt, z, lastCol, v, e)) {
        ++v;
    }

    ++rpsCnt[z[e]];
    ++z[e];
    --rpsCnt[z[e]];

    ++rpsCnt[z[v]];
    --z[v];
    --rpsCnt[z[v]];

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

    while (BndDecrementPossible(rpsCnt, z, b)) {
        --b;
    }

    e = b - 1;

    while (EdgeIncrementPossible(rpsCnt, z, e, b)) {
        --e;
    }
}

void NextRepGenPart(std::vector<int> &z, int &boundary, int &edge,
                    int &pivot, int lastCol, int lastElem) {

    int vertex = (z[boundary] - z[edge] == 2) ? boundary : edge + 1;

    ++z[edge];
    --z[vertex];

    if (vertex == boundary) {
        if (boundary < lastCol) {
            ++boundary;
        }

        const int currMax = z[boundary];

        while (boundary > 1 && z[boundary - 1] == currMax) {
            --boundary;
        }

        pivot = (z[boundary] < lastElem) ? lastCol : boundary - 1;

    } else if (z[vertex] == z[edge]) {
        ++vertex;
    }

    while (vertex < pivot) {
        const int distVert = z[vertex] - z[edge];
        const int distPivot = lastElem - z[pivot];

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

    if (boundary < lastCol &&
        z[boundary] < z[boundary + 1]) {

        ++boundary;
    }

    const int currMax = z[boundary];

    while (boundary > 1 && z[boundary - 1] == currMax)
        --boundary;

    edge = boundary - 1;
    int edgeTest = z[boundary] - 2;

    while (edge && edgeTest < z[edge]) {
        --edge;
    }
}

void NextDistinctGenPart(std::vector<int> &z, int &boundary,
                         int &edge, int &pivot, int &tarDiff,
                         int lastCol, int lastElem) {

    int vertex = edge + 1;
    tarDiff = 3;

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

        while (boundary > 1 && (z[boundary] - z[boundary - 1]) < 2) {
            --boundary;
        }

        pivot = (z[lastCol] < lastElem) ? lastCol : boundary - 1;
    }

    if (vertex < boundary) {
        if (z[vertex] - z[vertex - 1] == 1)
            ++vertex;

        while (vertex < pivot) {
            --z[vertex];
            ++z[pivot];

            if (z[vertex] - z[vertex - 1] == 1)
                ++vertex;

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

    edge = boundary - 1;
    tarDiff = 3;

    while (edge && (z[boundary] - z[edge]) < tarDiff) {
        --edge;
        ++tarDiff;
    }
}

using nextPartsPtr = void (*const)(std::vector<int> &rpsCnt,
                           std::vector<int> &z, int &e, int &b, int &p,
                           int &tarDiff, int lastCol, int lastElem);

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
};

void NextRepGen(std::vector<int> &rpsCnt, std::vector<int> &z, int &e,
                int &b, int &p, int &tarDiff, int lastCol, int lastElem) {
    NextRepGenPart(z, b, e, p, lastCol, lastElem);
};

void NextDistinctGen(std::vector<int> &rpsCnt, std::vector<int> &z,
                     int &e, int &b, int &p, int &tarDiff, int lastCol,
                     int lastElem) {
    NextDistinctGenPart(z, b, e, p, tarDiff, lastCol, lastElem);
};

nextPartsPtr GetNextPartsPtr(PartitionType ptype, bool IsGen) {

    if (IsGen) {
        if (ptype == PartitionType::Multiset) {
            return(nextPartsPtr(NextMultisetGen));
        } else if (std::find(RepPTypeArr.cbegin(), RepPTypeArr.cend(),
                             ptype) != RepPTypeArr.cend()) {
            return(nextPartsPtr(NextRepGen));
        } else {
            return(nextPartsPtr(NextDistinctGen));
        }
    } else {
        if (std::find(RepPTypeArr.cbegin(), RepPTypeArr.cend(),
                      ptype) != RepPTypeArr.cend()) {
            return(nextPartsPtr(NextRep));
        } else {
            return(nextPartsPtr(NextDistinct));
        }
    }
}
