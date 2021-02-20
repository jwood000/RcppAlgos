#include "StandardCount.h"
#include "NextPartitions.h"
#include "Cpp14MakeUnique.h"

template <typename typeVector>
inline void PopulateVecPerm(int m, const std::vector<typeVector> &v, std::vector<int> &z,
                        int &count, int maxRows, std::vector<typeVector> &partitionsVec) {
    
    do {
        for (int k = 0; k < m; ++k)
            partitionsVec.push_back(v[z[k]]);
        
        ++count;
    } while (std::next_permutation(z.begin(), z.end()) && count < maxRows);
}

inline bool keepGoing(const std::vector<int> &rpsCnt, int lastElem,
                      const std::vector<int> &z, int edge, int boundary) {
    
    if (edge >= 0) {
        const int myDiff = z[boundary] - z[edge];
        
        if (myDiff < 2) {
            return false;
        } else if (myDiff == 2) {
            return (rpsCnt[z[edge] + 1] > 1);
        } else {
            return (rpsCnt[z[edge] + 1] && rpsCnt[z[boundary] - 1]);
        }
    } else {
        return false;
    }
}

template <typename typeVector>
void PartitionsMultiSet(int m, const std::vector<typeVector> &v, std::vector<int> &z,
                        int lastElem, int lastCol, int maxRows, bool IsComb,
                        const std::vector<int> &Reps, std::vector<typeVector> &partitionsVec) {
    
    std::vector<int> rpsCnt(Reps.cbegin(), Reps.cend());
    
    for (const auto ind: z)
        --rpsCnt[ind];
    
    // boundary is the greatest index that can be decremented
    // pivot is the greatest index that can be incremented
    // edge is the greatest index smaller than boundary that can be incremented
    int p, e, b = lastCol, count = 0;
    PrepareMultiSetPart(rpsCnt, z, b, p, e, lastCol, lastElem);
    
    if (IsComb) {
        while (keepGoing(rpsCnt, lastElem, z, e, b)) {
            for (int k = 0; k < m; ++k)
                partitionsVec.push_back(v[z[k]]);
            
            ++count;
            
            if (count >= maxRows)
                break;
            
            NextMultiSetGenPart(rpsCnt, z, e, b,  p, lastCol, lastElem);
        }
        
        if (count < maxRows)
            for (int k = 0; k < m; ++k)
                partitionsVec.push_back(v[z[k]]);
    } else {
        while (keepGoing(rpsCnt, lastElem, z, e, b)) {
            do {
                for (int k = 0; k < m; ++k)
                    partitionsVec.push_back(v[z[k]]);
                
                ++count;
            } while (std::next_permutation(z.begin(), z.end()) && count < maxRows);
            
            if (count >= maxRows)
                break;
            
            NextMultiSetGenPart(rpsCnt, z, e, b,  p, lastCol, lastElem);
        }
        
        if (count < maxRows) {
            do {
                for (int k = 0; k < m; ++k)
                    partitionsVec.push_back(v[z[k]]);
                
                ++count;
            } while (std::next_permutation(z.begin(), z.end()) && count < maxRows);
        }
    }
}

template <typename typeRcpp, typename typeVector>
void PartitionsRep(int m, const std::vector<typeVector> &v, std::vector<int> &z,
                   int lastElem, int lastCol, int maxRows, bool IsComb,
                   std::vector<typeVector> &partitionsVec, typeRcpp &matRcpp) {
    
    // smallest index such that z[boundary] == currMax
    int boundary = lastCol;
    
    while (boundary > 1 && z[boundary - 1] == z[lastCol])
        --boundary;
    
    // pivot is the greatest index that can be incremented.
    // We know that if z[boundary] < lastElem ==>> pivot = lastCol
    int pivot = (z[boundary] < lastElem) ? lastCol : boundary - 1;
    
    // edge is the greatest index such that z[boundary] - z[edge] >= 2
    // This is the index that will be be used as a starting point
    // to determine the next combination that meets the criteria
    int edge = boundary - 1;
    int edgeTest = z[boundary] - 2;
    
    while (edge && edgeTest < z[edge])
        --edge;
    
    if (IsComb) {
        for (int count = 0; count < maxRows; ++count) {
            for (int k = 0; k < m; ++k)
                matRcpp(count, k) = v[z[k]];
            
            NextRepGenPart(z, boundary, edge, pivot, lastCol, lastElem);
        }
    } else {
        int count = 0;
        
        while ((edge >= 0) && (z[boundary] - z[edge] >= 2)) {
            PopulateVecPerm(m, v, z, count, maxRows, partitionsVec);
            
            if (count >= maxRows)
                break;
            
            NextRepGenPart(z, boundary, edge, pivot, lastCol, lastElem);
        }
        
        if (count < maxRows)
            PopulateVecPerm(m, v, z, count, maxRows, partitionsVec);
    }
}

template <typename typeRcpp, typename typeVector>
void PartitionsDistinct(int m, const std::vector<typeVector> &v, std::vector<int> &z,
                        int lastElem, int lastCol, int maxRows, bool IsComb,
                        std::vector<typeVector> &partitionsVec, typeRcpp &matRcpp) {
    
    // Largest index such that z[boundary] - z[boundary - 1] > 1
    // This is the index that will be decremented at the same
    // time edge is incremented. See below.
    int boundary = lastCol;
    
    while (boundary > 1 && (z[boundary] - z[boundary - 1]) < 2)
        --boundary;
    
    // pivot is the greatest index that can be incremented...
    // Either z[pivot + 1] - z[pivot] > 1 or if z[lastCol] < lastElem
    // pivot = lastCol since incrementing z[lastCol] is possible
    int pivot = (z[lastCol] < lastElem) ? lastCol : boundary - 1;
    
    // edge is the greatest index such that when incremented
    // the result will be at least one less than its neighbor
    // even if its neighbor is decremented
    int edge = boundary - 1;
    int tarDiff = 3;
    
    while (edge && (z[boundary] - z[edge]) < tarDiff) {
        --edge;
        ++tarDiff;
    }
    
    if (IsComb) {
        for (int count = 0; count < maxRows; ++count) {
            for (int k = 0; k < m; ++k)
                matRcpp(count, k) = v[z[k]];
            
            NextDistinctGenPart(z, boundary, edge, pivot, tarDiff, lastCol, lastElem);
        } 
    } else  {
        int count = 0;
        
        while (edge >= 0 && (z[boundary] - z[edge]) >= tarDiff) {
            PopulateVecPerm(m, v, z, count, maxRows, partitionsVec);
            
            if (count >= maxRows)
                break;
            
            NextDistinctGenPart(z, boundary, edge, pivot, tarDiff, lastCol, lastElem);
        }
        
        if (count < maxRows)
            PopulateVecPerm(m, v, z, count, maxRows, partitionsVec);
    }
}

template void PartitionsDistinct(int, const std::vector<int>&, std::vector<int>&,
                                 int, int, int, bool, std::vector<int>&, Rcpp::IntegerMatrix&);

template void PartitionsDistinct(int, const std::vector<double>&, std::vector<int>&,
                                 int, int, int, bool, std::vector<double>&, Rcpp::NumericMatrix&);

template void PartitionsRep(int, const std::vector<int>&, std::vector<int>&,
                                 int, int, int, bool, std::vector<int>&, Rcpp::IntegerMatrix&);

template void PartitionsRep(int, const std::vector<double>&, std::vector<int>&,
                                 int, int, int, bool, std::vector<double>&, Rcpp::NumericMatrix&);

template void PartitionsMultiSet(int, const std::vector<int>&, std::vector<int>&,
                                 int, int, int, bool, const std::vector<int>&, std::vector<int>&);

template void PartitionsMultiSet(int, const std::vector<double>&, std::vector<int>&,
                                 int, int, int, bool, const std::vector<int>&, std::vector<double>&);
