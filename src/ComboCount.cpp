#include "cpp11/R.hpp"
#include <algorithm> // std::min
#include <gmpxx.h>
#include <numeric>   // std::accumulate
#include <vector>
#include <deque>
#include <cmath>     // std::round

// Returns number of k-combinations from n elements.
// Mathematically speaking, we have: n! / (m! * (n - m)!)
//
// Using gmp function to get rid of floating point errors.
// E.g. choose(92, 13) = 2232785266876081
// shoule be gmp::chooseZ(92, 13) = 2232785266876080
double nChooseK(int n, int m) {

    if (m == n || m == 0) {
        return 1.0;
    }

    mpz_class temp;
    mpz_bin_uiui(temp.get_mpz_t(), n, m);

    double res = temp.get_d();
    return res;
}

double NumCombsWithRep(int n, int m) {
    return nChooseK(n + m - 1, m);
}

// The resulting vector, "triangleVec" resembles triangle numbers. In
// fact, this vector is obtained in a very similar method as generating
// triangle numbers, albeit in a repeating fashion. Two things to keep
// in mind is that we can't guarantee the following:
//      1) the repetition of each element is greater than or equal to n
//      2) that the repetition of each element isn't the same
double MultisetCombRowNumFast(int n, int m, const std::vector<int> &Reps) {

    if (m < 1 || n <= 1) {
        return 1.0;
    }

    if (m == n && std::accumulate(Reps.begin(), Reps.end(), 0) == n) {
        return 1.0;
    }

    const int m1 = m + 1;
    std::vector<double> triangleVec(m1);
    std::vector<double> temp(m1);

    int myMax = std::min(m1, Reps[0] + 1);

    for (int i = 0; i < myMax; ++i) {
        triangleVec[i] = temp[i] = 1.0;
    }

    --myMax;
    int ind = 1;

    for (; myMax < m; ++ind) {
        int myMin = std::min(Reps[ind], m);

        for (int i = 1; i <= myMin; ++i) {
            triangleVec[i] += triangleVec[i - 1];
        }

        myMin = std::min(Reps[ind] + myMax, m);
        int j = 0;

        for (int i = (Reps[ind] + 1); i <= myMin; ++i, ++j) {
            triangleVec[i] += triangleVec[i - 1];
            triangleVec[i] -= temp[j];
            temp[j] = triangleVec[j];
        }

        for (; j <= myMin; ++j) {
            temp[j] = triangleVec[j];
        }

        myMax = myMin;
    }

    const int n1 = n - 1;

    for (; ind < n1; ++ind) {
        double t = triangleVec[m];
        const int s = std::min(Reps[ind], m);

        for (int i = 1; i <= s; ++i) {
            triangleVec[m] += triangleVec[m - i];
        }

        double mySum = triangleVec[m];

        for (int i = m - 1; i >= s; --i) {
            mySum -= t;
            t = triangleVec[i];
            mySum += triangleVec[i - s];
            triangleVec[i] = mySum;
        }

        for (int i = s - 1; i > 0; --i) {
            mySum -= t;
            t = triangleVec[i];
            triangleVec[i] = mySum;
        }
    }

    if (ind < n) {
        const int myMin2 = std::min(Reps[n1], m);

        for (int i = 1; i <= myMin2; ++i) {
            triangleVec[m] += triangleVec[m - i];
        }
    }

    return triangleVec[m];
}

// This function will be used in the main function to determine whether
// gmp analogs are needed as the fast algorithm above could potentially
// produce negative results because of issues with double precision
double MultisetCombRowNum(int n, int m, const std::vector<int> &Reps) {

    if (m < 1 || n <= 1) {
        return 1;
    }

    const int m1 = m + 1;
    const int limit = (m1 > Reps[0] + 1) ? Reps[0] + 1 : m1;

    std::vector<double> temp(m1);
    std::vector<double> triangleVec(m1);

    for (int i = 0; i < limit; ++i) {
        triangleVec[i] = 1;
    }

    temp = triangleVec;

    for (int k = 1; k < n; ++k) {
        for (int i = m; i > 0; --i) {
            int myMax = i - Reps[k];

            if (myMax < 0) {
                myMax = 0;
            }

            double tempSum = 0;

            for (int j = myMax; j <= i; ++j) {
                tempSum += triangleVec[j];
            }

            temp[i] = tempSum;
        }

        triangleVec = temp;
    }

    return triangleVec[m];
}

void ManageCountsVector(std::vector<int> &Counts, int &n1) {
    if (!Counts.empty()) {
        --Counts.front();

        if (Counts.front() == 0 && Counts.size() > 1) {
            --n1;
            Counts.erase(Counts.begin());
        }
    }
}

void ManageCountsDeque(std::deque<int> &Counts, int &n1) {
    if (!Counts.empty()) {
        --Counts.front();

        if (Counts.size() > 1 && Counts.front() == 0) {
            --n1;
            Counts.pop_front();
        }
    }
}
