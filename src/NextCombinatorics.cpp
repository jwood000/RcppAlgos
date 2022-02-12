#include "Combinations/NextComboSection.h"
#include <algorithm>

using nextIterPtr = bool (*const)(const std::vector<int> &freqs,
                          std::vector<int> &z, int n1, int m1);

bool nextCombDistinct(const std::vector<int> &freqs,
                      std::vector<int> &z, int n1, int m1) {

    if (z.front() == (n1 - m1)) {
        return false;
    }

    if (z[m1] != n1) {
        ++z[m1];
    } else {
        nextCombSec(z, m1, n1 - m1);
    }

    return true;
}

bool nextCombRep(const std::vector<int> &freqs,
                 std::vector<int> &z, int n1, int m1) {

    if (z.front() == n1) {
        return false;
    }

    if (z[m1] != n1) {
        ++z[m1];
    } else {
        nextCombSecRep(z, m1, n1);
    }

    return true;
}

bool nextCombMulti(const std::vector<int> &freqs,
                   std::vector<int> &z, int n1, int m1) {

    bool res = false;

    for (int i = 0, k = freqs.size() - (m1 + 1); i <= m1; ++i, ++k) {
        if (z[i] != freqs[k]) {
            res = true;
            break;
        }
    }

    if (!res) {
        return false;
    }

    if (z[m1] != n1) {
        ++z[m1];
    } else {
        std::vector<int> zIndex(n1 + 1);

        for (int i = 0; i <= n1; ++i) {
            zIndex[i] = std::find(freqs.cbegin(),
                                  freqs.cend(), i) - freqs.cbegin();
        }

        nextCombSecMulti(freqs, zIndex, z, m1, freqs.size() - m1 - 1);
    }

    return true;
}

bool nextPermFull(const std::vector<int> &freqs,
                  std::vector<int> &z, int n1, int m1) {

    return std::next_permutation(z.begin(), z.end());
}

bool nextPermPartial(const std::vector<int> &freqs,
                     std::vector<int> &z, int n1, int m1) {

    bool res = false;

    if (freqs.size()) {
        for (int i = 0, k = freqs.size() - 1; i <= m1; ++i, --k) {
            if (z[i] != freqs[k]) {
                res = true;
                break;
            }
        }
    } else {
        for (int i = 0, j = n1; i <= m1; ++i, --j) {
            if (z[i] != j) {
                res = true;
                break;
            }
        }
    }

    if (!res) {
        return false;
    }

    int p1 = m1 + 1;

    while (p1 <= n1 && z[m1] >= z[p1]) {
        ++p1;
    }

    if (p1 <= n1) {
        std::swap(z[p1], z[m1]);
    } else {
        std::reverse(z.begin() + m1 + 1, z.end());
        p1 = m1;

        while (z[p1 + 1] <= z[p1]) {
            --p1;
        }

        int p2 = n1;

        while (z[p2] <= z[p1]) {
            --p2;
        }

        std::swap(z[p1], z[p2]);
        std::reverse(z.begin() + p1 + 1, z.end());
    }

    return true;
}

bool nextPermRep(const std::vector<int> &freqs,
                 std::vector<int> &z, int n1, int m1) {

    bool res = true;

    for (int i = 0; i <= m1; ++i) {
        if (z[i] != n1) {
            res = false;
            break;
        }
    }

    if (res) {
        return false;
    }

    for (int i = m1; i >= 0; --i) {
        if (z[i] != n1) {
            ++z[i];
            break;
        } else {
            z[i] = 0;
        }
    }

    return true;
}

nextIterPtr GetNextIterPtr(bool IsComb, bool IsMult,
                           bool IsRep, bool IsFull) {

    if (IsComb) {
        if (IsMult) {
            return(nextIterPtr(nextCombMulti));
        } else if (IsRep) {
            return(nextIterPtr(nextCombRep));
        } else {
            return(nextIterPtr(nextCombDistinct));
        }
    } else {
        if (IsRep) {
            return(nextIterPtr(nextPermRep));
        } else if (IsFull) {
            return(nextIterPtr(nextPermFull));
        } else {
            return(nextIterPtr(nextPermPartial));
        }
    }
}
