#include <algorithm>
#include <vector>

template <typename T>
void PopulateVec(const std::vector<T> &v,
                 std::vector<T> &cnstrntVec,
                 std::vector<int> &z, int &count,
                 int m, int nRows, bool IsComb) {

    if (IsComb) {
        for (int k = 0; k < m; ++k) {
            cnstrntVec.push_back(v[z[k]]);
        }

        ++count;
    } else {
        do {
            for (int k = 0; k < m; ++k) {
                cnstrntVec.push_back(v[z[k]]);
            }

            ++count;
        } while (std::next_permutation(z.begin(), z.end()) && count < nRows);
    }
}

template void PopulateVec(const std::vector<int>&, std::vector<int>&,
                          std::vector<int>&, int&, int, int, bool);

template void PopulateVec(const std::vector<double>&, std::vector<double>&,
                          std::vector<int>&, int&, int, int, bool);
