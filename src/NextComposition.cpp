#include <algorithm>
#include <numeric>
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

int NextDistinctSection(
    const std::vector<int> &v, std::vector<int> &idx, int target, int m
) {

    int lenV = v.size();
    const int testMax = std::accumulate(v.cend() - m, v.cend(), 0);

    // The length is too small
    if (testMax < target) {
        return -2;
    }

    const int testMin = std::accumulate(v.cbegin(), v.cbegin() + m, 0);

    // The length is too long
    if (testMin > target) {
        return -1;
    }

    std::iota(idx.begin(), idx.end(), 0);

    auto GetSum = [](const std::vector<int> &_v,
                     const std::vector<int> &_idx) {
        return std::accumulate(
            _idx.begin(), _idx.end(), 0, [&](int sum, int i) {
                return sum + _v[i];
            }
        );
    };

    int tempSum = GetSum(v, idx);
    bool keepGoing = tempSum != target;
    int partial = target - (tempSum - v[idx.back()]);
    int j = m - 1;

    while (keepGoing) {
        auto lower = std::lower_bound(v.begin() + j, v.end(), partial);

        if (lower != v.end() && *lower == partial) {
            idx[m - 1] = std::distance(v.begin(), lower);
            return 1;
        }

        int k = m - 2;
        int g = lenV - 2;

        while (k > 0 && idx[k] == g) {
            --k;
            --g;
        }

        if (k == 0 && idx.front() == (lenV - m)) {
            return 0;
        }

        bool impossible = true;

        while (impossible && k >= 0) {
            ++idx[k];

            for (int i = k + 1; i < m; ++i) {
                idx[i] = idx[i - 1] + 1;
            }

            tempSum = GetSum(v, idx);
            int maxSum = tempSum;

            for (int i = k; i < m; ++i) {
                maxSum -= v[idx[i]];
                maxSum += v[lenV - (m - i)];
            }

            impossible = (tempSum > target) ||
                (maxSum < target);
            --k;
        }

        if (impossible) {
            return 0;
        }

        keepGoing = tempSum != target;
        partial = target - (tempSum - v[idx.back()]);
        j = idx.back();
    }

    return 1;
}

bool NextRoutine(
    std::vector<int> &z, std::vector<int> &compliment,
    int &i1, int &i2, int lastCol, int lastIdx
) {

    auto&& ref_one = z[lastCol - 1];
    auto&& ref_two = z[lastCol];

    if (i1 < 0) {
        std::swap(ref_one, ref_two);
        return true;
    }

    const int target = ref_one + ref_two;
    bool less_flag = ref_one < ref_two;
    int check = compliment[i1] + compliment[i2];

    while (check != target) {
        if (check < target) {
            ++i1;
        } else {
            --i2;
        }

        check = compliment[i1] + compliment[i2];
    }

    if (less_flag) {
        if (i1 == i2) {
            std::swap(ref_one, ref_two);
            ++i1;
            --i2;
            return true;
        } else if (ref_two - ref_one < compliment[i1] - compliment[i2]) {
            std::swap(ref_one, ref_two);
            return true;
        }
    }

    std::swap(ref_one, compliment[i1]);
    std::swap(ref_two, compliment[i2]);
    std::sort(compliment.begin(), compliment.end());

    if (i1 < lastIdx) ++i1; else return false;
    if (i2 > 0) --i2; else return false;

    if (i1 == i2 && i1 < lastIdx && i2 > 0) {
        ++i1;
        --i2;
    }

    return true;
}

void NextCompositionDistinct(
    std::vector<int> &z, std::vector<int> &compliment,
    int i1, int i2, int myMax, int lastCol, int lastIdx
) {

    if (z[lastCol - 1] < myMax) {
        NextRoutine(z, compliment, i1, i2, lastCol, lastIdx);
    } else if (lastCol > 1) {
        int m = 2;
        int j = lastCol - m;
        const int target = std::accumulate(z.cbegin(), z.cend(), 0);

        compliment.insert(compliment.end(), z.end() - m, z.end());
        std::sort(compliment.begin(), compliment.end());
        bool keepGoing = true;

        while (j >= 0 && keepGoing) {
            int res = 0;
            std::vector<int> idx(m);

            while (res == 0) {
                auto lower = std::lower_bound(
                    compliment.begin(), compliment.end(), z[j]
                );

                if (lower != compliment.end()) {
                    std::swap(*lower, z[j]);
                    int partial = target -
                        std::accumulate(z.cbegin(), z.cend() - m, 0);
                    res = NextDistinctSection(compliment, idx, partial, m);
                } else {
                    res = 100;
                }
            }

            if (res == 1) {
                keepGoing = false;

                for (int i = m - 1, k = lastCol; i >= 0; --i, --k) {
                    z[k] = compliment[idx[i]];
                    compliment.erase(compliment.begin() + idx[i]);
                }

                myMax = z[lastCol];

                idx.resize(2);
                int lastTwo = z[lastCol - 1] + z[lastCol];
                int res2 = NextDistinctSection(compliment, idx, lastTwo, 2);

                if (res2 == 1) {
                    i1 = idx.front();
                    i2 = idx.back();
                } else {
                    i1 = -1;
                }
            } else if (res == 100) {
                keepGoing = true;
                compliment.push_back(z[j]);
                --j;
                ++m;
            }
        }
    }
}

template void NextCompositionRep<0>(std::vector<int>&, int);
template void NextCompositionRep<1>(std::vector<int>&, int);
