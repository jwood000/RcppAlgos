#include <vector>
#include <algorithm>

template <int num>
void NextComposition(std::vector<int> &z, int lastCol) {

    if (z[lastCol] != num) {
        --z[lastCol];
        ++z[lastCol - 1];
    } else {
        int j = lastCol - 1;

        while (z[j] == num) {
            --j;
        }

        ++z[j - 1];
        std::reverse(z.begin() + j, z.end());
        --z[lastCol];
    }
}

template void NextComposition<0>(std::vector<int>&, int);
template void NextComposition<1>(std::vector<int>&, int);
