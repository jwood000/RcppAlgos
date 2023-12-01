// This algorithm is nearly identical to the one found in the
// standard algorithm library. Also, this and the following
// algo are used in functions found in Permutations.h as well
// as PermuteResults.h
#include <cstdlib>

void nextFullPerm(int *const myArray, int maxInd) {

    int p1 = maxInd - 1;
    int p2 = maxInd;

    while (myArray[p1 + 1] <= myArray[p1])
        --p1;

    while (myArray[p2] <= myArray[p1])
        --p2;

    int temp = myArray[p1];
    myArray[p1] = myArray[p2];
    myArray[p2] = temp;

    for (int k = p1 + 1, q = maxInd; k < q; ++k, --q) {
        int temp = myArray[k];
        myArray[k] = myArray[q];
        myArray[q] = temp;
    }
}

// This algorithm is the same as above except that since we are
// not using the entire vector, we have to first check that the
// mth element is the largest. If it is, we have to reverse all
// of the elements to the right of the rth position before
// finding the next permutation. This is so because if we did
// not, all of the next perms of the entire vector would produce
// many duplicate m-length perms. If it isn't the largest, we
// find the element to the right and swap them. We can then
// proceed to the next perm. We can do this because the standard
// algo would end up performing two unnecessary reversings.
void nextPartialPerm(int *const myArray, int m1, int maxInd) {

    int p1 = m1 + 1;

    while (p1 <= maxInd && myArray[m1] >= myArray[p1])
        ++p1;

    if (p1 <= maxInd) {
        int temp = myArray[p1];
        myArray[p1] = myArray[m1];
        myArray[m1] = temp;
    } else {
        for (int k = m1 + 1, q = maxInd; k < q; ++k, --q) {
            int temp = myArray[k];
            myArray[k] = myArray[q];
            myArray[q] = temp;
        }

        p1 = m1;
        while (myArray[p1 + 1] <= myArray[p1])
            --p1;

        int p2 = maxInd;

        while (myArray[p2] <= myArray[p1])
            --p2;

        int temp = myArray[p1];
        myArray[p1] = myArray[p2];
        myArray[p2] = temp;

        for (int k = p1 + 1, q = maxInd; k < q; ++k, --q) {
            int temp = myArray[k];
            myArray[k] = myArray[q];
            myArray[q] = temp;
        }
    }
}
