#include <vector>
#include <algorithm>

void NextCompositionZero(std::vector<int> &z, int lastCol) {

    if (z[lastCol] != 0) {
        --z[lastCol];
        ++z[lastCol - 1];
    } else {
        int j = lastCol - 1;

        while (z[j] == 0)
            --j;

        ++z[j - 1];
        std::reverse(z.begin() + j, z.end());
        --z[lastCol];
    }
}

void NextCompositionOne(std::vector<int> &z, int lastCol) {

    if (z[lastCol] != 1) {
        --z[lastCol];
        ++z[lastCol - 1];
    } else {
        int j = lastCol - 1;

        while (z[j] == 1)
            --j;

        ++z[j - 1];
        std::reverse(z.begin() + j, z.end());
        --z[lastCol];
    }
}
