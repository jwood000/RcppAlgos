#include <vector>
#include <algorithm> // std::min
#include <numeric>   // std::accumulate
#include <cmath>     // std::round

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

// Returns number of k-combinations from n elements.
// Mathematically speaking, we have: n!/(k!*(n-k)!)
double nChooseK(int n, int k) {

    if (k == n || k == 0) {
        return 1.0;
    }

    double nCk = 1;

    for (double i = (n - k + 1), d = 1; d <= k; ++i, ++d) {
        nCk *= i;
        nCk /= d;
    }

    return std::round(nCk);
}

double NumCombsWithRep(int n, int r) {
    return nChooseK(n + r - 1, r);
}

// The resulting vector, "triangleVec" resembles triangle numbers. In
// fact, this vector is obtained in a very similar method as generating
// triangle numbers, albeit in a repeating fashion. Two things to keep
// in mind is that we can't guarantee the following:
//      1) the repetition of each element is greater than or equal to n
//      2) that the repetition of each element isn't the same
double MultisetCombRowNumFast(int n, int r, const std::vector<int> &Reps) {

    if (r < 1 || n <= 1) {
        return 1.0;
    }

    if (r == n) {
        if (std::accumulate(Reps.cbegin(), Reps.cend(), 0) == n) {
            return 1.0;
        }
    }

    const int r1 = r + 1;
    std::vector<double> triangleVec(r1);
    std::vector<double> temp(r1);

    int myMax = r1;

    if (myMax > Reps[0] + 1) {
        myMax = Reps[0] + 1;
    }

    for (int i = 0; i < myMax; ++i) {
        triangleVec[i] = temp[i] = 1;
    }

    --myMax;
    int ind = 1;

    for (; myMax < r; ++ind) {
        int myMin = std::min(Reps[ind], r);

        for (int i = 1; i <= myMin; ++i) {
            triangleVec[i] += triangleVec[i - 1];
        }

        myMin = std::min(Reps[ind] + myMax, r);
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
        double t = triangleVec[r];
        const int s = std::min(Reps[ind], r);

        for (int i = 1; i <= s; ++i) {
            triangleVec[r] += triangleVec[r - i];
        }

        double mySum = triangleVec[r];

        for (int i = r - 1; i >= s; --i) {
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
        const int myMin2 = std::min(Reps[n1], r);

        for (int i = 1; i <= myMin2; ++i) {
            triangleVec[r] += triangleVec[r - i];
        }
    }

    return triangleVec[r];
}

// This function will be used in the main function to determine whether
// gmp analogs are needed as the fast algorithm above could potentially
// produce negative results because of issues with double precision
double MultisetCombRowNum(int n, int r, const std::vector<int> &Reps) {

    if (r < 1 || n <= 1) {
        return 1;
    }

    const int r1 = r + 1;
    const int limit = (r1 > Reps[0] + 1) ? Reps[0] + 1 : r1;

    std::vector<double> temp(r1);
    std::vector<double> triangleVec(r1);

    for (int i = 0; i < limit; ++i) {
        triangleVec[i] = 1;
    }

    temp = triangleVec;

    for (int k = 1; k < n; ++k) {
        for (int i = r; i > 0; --i) {
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

    return triangleVec[r];
}
