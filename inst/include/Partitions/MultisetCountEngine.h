#pragma once

#include <cstddef>
#include <vector>

template <bool IsComp>
class MultisetCountEngine {
public:
    double Count(
        int n, int m,
        const std::vector<int>& allowed,
        const std::vector<int>& Reps,
        int strtLen = 0
    );

private:
    std::vector<double> dp;
    std::vector<double> dp_old;

    // Flattened table: dp[j][s] is stored at dp[j * n1 + s].
    static constexpr std::size_t idx(
        std::size_t j, std::size_t s, std::size_t n1
    ) noexcept {
        return j * n1 + s;
    }
};
