#include <algorithm>
#include <vector>

template <int one_or_zero>
void NextCompositionRep(std::vector<int> &z, int lastCol) {

    if (z[lastCol] != one_or_zero) {
        --z[lastCol];
        ++z[lastCol - 1];
    } else {
        int j = lastCol - 1;

        while (j > 0 && z[j] == one_or_zero) {
            --j;
        }

        if (j > 0) {
            ++z[j - 1];
            std::reverse(z.begin() + j, z.end());
            --z[lastCol];
        }
    }
}

void NextCompositionDistinct(
    std::vector<int> &z, int &idx_1, int &idx_2, int &final_idx, int lastCol
) {

    if (z[lastCol] != final_idx) {

        bool keepGoing = true;
        int last_idx = z[lastCol];
        int pent_idx = z[lastCol - 1];

        const int z_size = z.size();
        const int tar = last_idx + pent_idx;

        while (keepGoing) {
            const bool check_less = idx_1 < idx_2;

            while (z[idx_1] < pent_idx && idx_1 < z_size) {
                ++idx_1;
            }

            // Need distinct
            if (idx_1 < z_size && tar != (2 * z[idx_1])) {

                int gap = tar - z[idx_1];

                while (gap > z[idx_2] && idx_2 > lastCol) {
                    --idx_2;
                }

                keepGoing = idx_2 > lastCol && gap != z[idx_2];
            }
        }

        std::swap(z[lastCol - 1], z[idx_1]);
        std::swap(z[lastCol], z[idx_2]);
    } else {
        int j = lastCol - 1;

        while (j > 0 && z[j] == 0) {
            --j;
        }

        if (j > 0) {
            ++z[j - 1];
            std::reverse(z.begin() + j, z.end());
            --z[lastCol];
        }
    }
}

template void NextCompositionRep<0>(std::vector<int>&, int);
template void NextCompositionRep<1>(std::vector<int>&, int);
