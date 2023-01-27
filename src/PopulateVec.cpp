#include <algorithm>
#include <vector>

template <typename T>
void PopulateVec(const std::vector<T> &v,
                 std::vector<T> &cnstrntVec,
                 std::vector<int> &z, std::size_t &count,
                 std::size_t m, std::size_t nRows, bool IsComb) {

    if (IsComb) {
        for (std::size_t k = 0; k < m; ++k) {
            cnstrntVec.push_back(v[z[k]]);
        }

        ++count;
    } else {
        do {
            for (std::size_t k = 0; k < m; ++k) {
                cnstrntVec.push_back(v[z[k]]);
            }

            ++count;
        } while (std::next_permutation(z.begin(), z.end()) && count < nRows);
    }
}

template void PopulateVec(const std::vector<int>&, std::vector<int>&,
                          std::vector<int>&, std::size_t&,
                          std::size_t, std::size_t, bool);

template void PopulateVec(const std::vector<double>&, std::vector<double>&,
                          std::vector<int>&, std::size_t&,
                          std::size_t, std::size_t, bool);
