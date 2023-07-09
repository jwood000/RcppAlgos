#pragma once

#include <functional>
#include <algorithm>
#include <vector>
#include <cstddef>

class Group {
private:
    const std::vector<int> ubound;
    const std::vector<int> lbound;
    const std::vector<bool> same;

    void situate(std::vector<int> &z, int idx1, int offset) const;
    bool is_max(const std::vector<int> &z, int i) const;

public:

    const std::vector<int> grp;

    Group(const std::vector<int> &_grp,
          const std::vector<int> &_ubound,
          const std::vector<int> &_lbound,
          const std::vector<bool> &_same);

    void balance(std::vector<int> &z, int idx1,
                 int curr_bnd, int i) const;
    bool require_external(const std::vector<int> &z, int i) const;
    bool flip_external(std::vector<int> &z, int &idx1, int i) const;
    void step(int &idx1, int &idx2, int &curr_bnd, int i) const;
    int get_low(int curr_bnd, int i) const;
    int get_size() const;
};
