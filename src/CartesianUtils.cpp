#include <cstdlib>
#include <numeric>
#include <vector>
#include <gmpxx.h>

double CartesianCount(const std::vector<int> &lenGrps) {
    return std::accumulate(lenGrps.begin(), lenGrps.end(),
                           1.0, std::multiplies<double>());
}

void CartesianCountGmp(mpz_class &result, const std::vector<int> &lenGrps) {

    result = 1;

    for (auto len: lenGrps) {
        result *= len;
    }
}

std::vector<int> nthProduct(double dblIdx, const std::vector<int> &lenGrp) {

    double index1 = dblIdx;
    const int m = lenGrp.size();

    std::vector<int> res(m);
    double temp = CartesianCount(lenGrp);

    for (int k = 0; k < m; ++k) {
        temp /= lenGrp[k];
        int j = static_cast<int>(index1 / temp);
        res[k] = j;
        index1 -= (temp * j);
    }

    for (auto &v_i: res) {
        v_i *= m;
    }

    return res;
}

std::vector<int> nthProductGmp(const mpz_class &mpzIdx,
                               const std::vector<int> &lenGrp) {

    mpz_class index1(mpzIdx);
    const int m = lenGrp.size();

    std::vector<int> res(m);
    mpz_class temp;
    mpz_class temp2;
    CartesianCountGmp(temp, lenGrp);

    for (int k = 0; k < m; ++k) {
        mpz_divexact_ui(temp.get_mpz_t(), temp.get_mpz_t(), lenGrp[k]);
        temp2 = index1 / temp;
        int j = temp2.get_si();
        res[k] = j;
        index1 -= (temp * j);
    }

    for (auto &v_i: res) {
        v_i *= m;
    }

    return res;
}

bool nextProduct(const std::vector<int> &lenGrps,
                 std::vector<int> &z, int m) {

    if (z.back() < lenGrps.back()) {
        z.back() += m;
        return true;
    } else {
        z.back() = 0;

        for (int i = m - 2; i >= 0; --i) {
            if (z[i] < lenGrps[i]) {
                z[i] += m;
                return true;
            } else {
                z[i] = 0;
            }
        }
    }

    return false;
}
