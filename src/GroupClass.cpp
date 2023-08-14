#include "ComboGroup/GroupClass.h"

void Group::situate(std::vector<int> &z, int idx1, int offset) const {
    int idx3 = idx1 + 1;
    std::sort(z.begin() + idx3, z.end());

    while (z[idx3] < z[idx1]) {
        ++idx3;
    }

    std::swap(z[idx3], z[idx1]);
    std::rotate(z.begin() + idx1 + 1,
                z.begin() + idx3 + 1,
                z.begin() + idx3 + offset - idx1);
}

bool Group::is_max(const std::vector<int> &z, int i) const {
    std::vector<int> v(z.begin() + lbound[i], z.end());
    std::sort(v.begin(), v.end(), std::greater<int>());

    int n_grps = 1;
    const int limit = static_cast<int>(grp.size() - 1);

    for (int j = i; j < limit && grp[j] == grp[j + 1]; ++j) {
        ++n_grps;
    }

    return z[lbound[i]] == v[(n_grps * grp[i]) - 1];
}

Group::Group(const std::vector<int> &_grp,
             const std::vector<int> &_ubound,
             const std::vector<int> &_lbound,
             const std::vector<bool> &_same) :
    ubound(_ubound), lbound(_lbound), same(_same), grp(_grp) {}

void Group::balance(std::vector<int> &z, int idx1,
                    int curr_bnd, int i) const {

    situate(z, idx1, curr_bnd + grp[i]);

    if (same[i] && z[lbound[i]] > z[lbound[i + 1L]]) {
        int j = i;

        while (grp[j] == grp[j + 1]) {
            ++j;
        }

        int len_rng    = lbound[j + 1] - lbound[i + 1];
        int idx4       = ubound[i] + 1L;
        const int edge = z[lbound[i]];

        while (z[idx4] < edge) {
            ++idx4;
        }

        std::rotate(z.begin() + ubound[i] + 1L,
                    z.begin() + idx4,
                    z.begin() + idx4 + len_rng);
    }
}

bool Group::require_external(const std::vector<int> &z, int i) const {
    if (!same[i])     return false;
    if (is_max(z, i)) return false;
    return grp[i] != grp.back();
}

bool Group::flip_external(std::vector<int> &z, int &idx1, int i) const {

    int j = i;

    while (grp[j] == grp[j + 1L]) {
        ++j;
    }

    idx1         = ubound[j - 1L];
    int idx2     = ubound[j + 1L];
    int grp_size = grp[i] * 2;
    int curr_bnd = lbound[i];

    for (int j = 0; idx1 > curr_bnd && z[idx2] < z[idx1]; ++j) {
        --idx1;

        if (j == grp[i]) {
            j = 0;
            grp_size += grp[i];
        }
    }

    if (z[idx1] < z[idx2]) {
        situate(z, idx1, curr_bnd + grp_size);
        return true;
    }

    return false;
}

void Group::step(int &idx1, int &idx2, int &curr_bnd, int i) const {
    idx2     -= grp[i + 1];
    curr_bnd -= grp[i - 1];
    idx1     -= same[i] ? 2 : 1;
}

int Group::get_low(int curr_bnd, int i) const {
    return same[i] ? curr_bnd + 1 : curr_bnd;
}

int Group::get_size() const {
    return grp.size();
}
