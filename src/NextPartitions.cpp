#include <vector>

void NextDistinct(std::vector<int> &z, int &boundary,
                  int &edge, int &tarDiff, int lastCol) {
    
    if (z[boundary] - z[edge] != tarDiff)
        boundary = edge + 1;
    
    ++z[edge];
    --z[boundary];
    int myEdgePlus = z[edge] + boundary - edge;
    
    for (; boundary < lastCol; ++boundary, ++myEdgePlus) {
        z[lastCol] += (z[boundary] - myEdgePlus);
        z[boundary] = myEdgePlus;
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

void NextPartition(std::vector<int> &z, int &boundary, int &edge, int lastCol) {
    
    if (z[boundary] - z[edge] != 2)
        boundary = edge + 1;
    
    ++z[edge];
    --z[boundary];
    const int myEdge = z[edge];
    
    for (; boundary < lastCol; ++boundary) {
        z[lastCol] += (z[boundary] - myEdge);
        z[boundary] = myEdge;
    }
    
    const int currMax = z[boundary];
    
    while (boundary > 1 && z[boundary - 1] == currMax)
        --boundary;
    
    edge = boundary - 1;
    const int edgeTest = z[boundary] - 2;
    
    while (edge && edgeTest < z[edge])
        --edge;
}

// BndDecrementPossible, VtxDecrementPossible, EdgeIncrementPossible, GetPivotExtr
// and PivotDecrementPossible are all helper functions for PartitionsMultiSet
// ************************** Start Helper Functions ***************************
inline bool BndDecrementPossible(const std::vector<int> &rpsCnt,
                                 const std::vector<int> &z, int boundary) {
    
    if (boundary > 1) {
        if ((z[boundary] - z[boundary - 1]) < 2) {
            return ((z[boundary] != z[boundary - 1]) ? !rpsCnt[z[boundary] - 1] : true);
        } else {
            return false;
        }
    } else {
        return false;
    }
}

inline bool VtxDecrementPossible(const std::vector<int> &rpsCnt, int lastCol,
                                 const std::vector<int> &z, int vertex, int edge) {
    
    if (vertex < lastCol) {
        const int myDiff = z[vertex] - z[edge];
        
        if (myDiff < 2) {
            return true;
        } else  {
            return ((myDiff != 2 || rpsCnt[z[edge] + 1] <= 1) &&
                    (myDiff <= 2 || !rpsCnt[z[edge] + 1] || !rpsCnt[z[vertex] - 1]));
        }
    } else {
        return false;
    }
}

inline bool EdgeIncrementPossible(const std::vector<int> &rpsCnt,
                                  const std::vector<int> &z, int edge, int boundary) {
    
    if (edge) {
        const int myDiff = z[boundary] - z[edge];
        
        if (myDiff < 2) {
            return true;
        } else {
            return (myDiff != 2 || rpsCnt[z[edge] + 1] < 2) && (myDiff <= 2 || !rpsCnt[z[edge] + 1]);
        }
    } else {
        return false;
    }
}

inline int GetPivotExtr(const std::vector<int> &rpsCnt,
                        const std::vector<int> &z, int lastCol, int lastElem) {
    
    int res = lastCol - 1;
    
    while (res > 0 && z[res] == lastElem)
        --res;
    
    while (res > 0 && !rpsCnt[z[res] + 1])
        --res;
    
    return res;
}

inline bool PivotDecrementPossible(const std::vector<int> &rpsCnt, int lastElem,
                                   const std::vector<int> &z, int pivot, int vertex) {
    
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

void PrepareMultiSetPart(const std::vector<int> &rpsCnt, const std::vector<int> &z,
                         int &b, int &p, int &e, int lastCol, int lastElem) {
    
    while (BndDecrementPossible(rpsCnt, z, b))
        --b;
    
    p = (z[lastCol] < lastElem) ? lastCol : GetPivotExtr(rpsCnt, z, lastCol, lastElem);
    e = b - 1;
    
    while (EdgeIncrementPossible(rpsCnt, z, e, b))
        --e;
}

// **************************** End  Helper Functions ***************************

void NextMultiSetGenPart(std::vector<int> &rpsCnt, std::vector<int> &z,
                         int &e, int &b, int &p, int lastCol, int lastElem) {
    
    // vertex is the smallest index greater than edge that can be decremented
    int v = e + 1;
    
    while (VtxDecrementPossible(rpsCnt, lastCol, z, v, e))
        ++v;
    
    ++rpsCnt[z[e]];
    ++z[e];
    --rpsCnt[z[e]];
    
    ++rpsCnt[z[v]];
    --z[v];
    --rpsCnt[z[v]];
    
    if (v == b) {
        if (b < lastCol)
            ++b;
        
        while (BndDecrementPossible(rpsCnt, z, b))
            --b;
        
        p = (z[lastCol] < lastElem) ? lastCol : GetPivotExtr(rpsCnt, z, lastCol, lastElem);
    }
    
    while ((v < lastCol) && (z[v] == z[v - 1] || 
           z[v] == z[e] || (z[v] - z[v - 1] == 1 && !rpsCnt[z[v - 1]]))) {
        ++v;
    }
    
    while (v < p && rpsCnt[z[v] - 1] && rpsCnt[z[p] + 1]) {
        ++rpsCnt[z[v]];
        --z[v];
        --rpsCnt[z[v]];
        
        ++rpsCnt[z[p]];
        ++z[p];
        --rpsCnt[z[p]];
        
        while (z[v] == z[v - 1] || (z[v] - z[v - 1] == 1 && !rpsCnt[z[v - 1]]))
            ++v;
        
        while (PivotDecrementPossible(rpsCnt, lastElem, z, p, v))
            --p;
    }
    
    b = p;
    
    while ((b < lastCol) && ((z[b] == z[b + 1]) ||
           (z[b + 1] > z[b] && (rpsCnt[z[b + 1] - 1] || rpsCnt[z[b + 1]])))) {
        ++b;
    }
    
    while (BndDecrementPossible(rpsCnt, z, b))
        --b;
    
    e = b - 1;
    
    while (EdgeIncrementPossible(rpsCnt, z, e, b))
        --e;
}

void NextRepGenPart(std::vector<int> &z, int &boundary, int &edge, 
                    int &pivot, int lastCol, int lastElem) {
    
    int vertex = (z[boundary] - z[edge] == 2) ? boundary : edge + 1;
    
    ++z[edge];
    --z[vertex];
    
    if (vertex == boundary) {
        if (boundary < lastCol)
            ++boundary;
        
        const int currMax = z[boundary];
        
        while (boundary > 1 && z[boundary - 1] == currMax)
            --boundary;
        
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
    
    if (boundary < lastCol && z[boundary] < z[boundary + 1])
        ++boundary;
    
    const int currMax = z[boundary];
    
    while (boundary > 1 && z[boundary - 1] == currMax)
        --boundary;
    
    edge = boundary - 1;
    int edgeTest = z[boundary] - 2;
    
    while (edge && edgeTest < z[edge])
        --edge;
}

void NextDistinctGenPart(std::vector<int> &z, int &boundary, int &edge, 
                         int &pivot, int &tarDiff, int lastCol, int lastElem) {
    
    int vertex = edge + 1;
    tarDiff = 3;
    
    while (vertex < lastCol && (z[vertex] - z[edge]) < tarDiff) {
        ++vertex;
        ++tarDiff;
    }
    
    ++z[edge];
    --z[vertex];
    
    if (vertex == boundary) {
        if (boundary < lastCol)
            ++boundary;
        
        while (boundary > 1 && (z[boundary] - z[boundary - 1]) < 2)
            --boundary;
        
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
            
            if (z[pivot] == lastElem || (pivot < lastCol && z[pivot + 1] - z[pivot] == 1))
                --pivot;
        }
        
        boundary = pivot;
        
        if (boundary < lastCol && z[boundary + 1] - z[boundary] > 1)
            ++boundary;
    }
    
    edge = boundary - 1;
    tarDiff = 3;
    
    while (edge && (z[boundary] - z[edge]) < tarDiff) {
        --edge;
        ++tarDiff;
    }
}
