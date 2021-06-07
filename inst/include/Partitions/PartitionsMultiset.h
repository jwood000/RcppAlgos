#ifndef PARTITIONS_MULTISET_H
#define PARTITIONS_MULTISET_H

#include <vector>

bool keepGoing(const std::vector<int> &rpsCnt, int lastElem,
               const std::vector<int> &z, int edge, int boundary);

template <typename T>
void PartsGenMultiset(std::vector<T> &partitionsVec, const std::vector<T> &v,
                      const std::vector<int> &Reps, std::vector<int> &z,
                      int width, int lastElem, int lastCol, int maxRows);

template <typename T>
void PartsGenPermMultiset(std::vector<T> &partitionsVec,
                          const std::vector<T> &v,
                          const std::vector<int> &Reps, std::vector<int> &z,
                          int width, int lastElem, int lastCol, int maxRows);

#endif
