#include <unordered_set>
#include <gmpxx.h>
#include <cstdint>
#include <limits>
#include <vector>
#include <string>

void AddComb(std::unordered_set<std::uint64_t> &uintHash,
             std::vector<std::uint64_t> &uintKeyKeeper,
             const std::vector<int> &rowIdx, std::vector<int> &cartCombs,
             std::uint64_t &maxKey, std::uint64_t masterKey,
             int prime, int idx) {

    const std::uint64_t key = masterKey * static_cast<std::uint64_t>(prime);
    if (key > maxKey) maxKey = key;

    if (uintHash.find(key) == uintHash.end()) {
        uintHash.insert(key);
        cartCombs.insert(cartCombs.end(),
                         rowIdx.cbegin(),
                         rowIdx.cend());
        cartCombs.push_back(idx);
        uintKeyKeeper.push_back(key);
    }
}

void AddComb(std::unordered_set<std::string> &strHash,
             std::vector<std::string> &strKeyKeeper,
             const std::vector<int> &rowIdx,
             std::vector<int> &cartCombs,
             mpz_class key, mpz_class masterKey,
             int prime, int idx) {

    key = masterKey * prime;
    const std::string strKey = key.get_str();

    if (strHash.find(strKey) == strHash.end()) {
        strHash.insert(strKey);
        cartCombs.insert(cartCombs.end(),
                         rowIdx.cbegin(),
                         rowIdx.cend());
        cartCombs.push_back(idx);
        strKeyKeeper.push_back(strKey);
    }
}

void comboGrid(std::vector<int> &cartCombs,
               std::vector<int> &lastCol,
               std::vector<int> &lenGrps,
               const std::vector<std::vector<int>> &myVecs,
               const std::vector<int> &primes, bool IsRep) {

    std::unordered_set<std::uint64_t> uintHash;
    std::vector<std::uint64_t> uintKeyKeeper;
    std::uint64_t maxKey = 0;

    for (const auto idx: myVecs[0]) {
        const std::uint64_t key = primes[idx];
        if (key > maxKey) maxKey = key;

        if (uintHash.find(key) == uintHash.end()) {
            uintHash.insert(key);
            cartCombs.push_back(idx);
            uintKeyKeeper.push_back(key);
        }
    }

    std::size_t i = 1;
    bool NeedsMpz = false;

    for (std::size_t size = myVecs.size() - 1; i < size; ++i) {
        if ((std::numeric_limits<std::uint64_t>::max() / maxKey) <
            static_cast<std::uint64_t>(primes[myVecs[i].back()])) {
            NeedsMpz = true;
            break;
        }

        const std::vector<int> tempCombs = cartCombs;
        const std::vector<std::uint64_t> tempKeyKeeper = uintKeyKeeper;

        uintHash.clear();
        cartCombs.clear();
        uintKeyKeeper.clear();

        uintHash.reserve(tempKeyKeeper.size() * myVecs[i].size());
        uintKeyKeeper.reserve(tempKeyKeeper.size() * myVecs[i].size());
        cartCombs.reserve(tempKeyKeeper.size() * myVecs[i].size() * (i + 1));

        if (IsRep) {
            for (std::size_t j = 0, myStrt = 0, myEnd = i;
                 j < tempKeyKeeper.size(); ++j, myStrt = myEnd, myEnd += i) {

                const std::uint64_t masterKey = tempKeyKeeper[j];
                const std::vector<int> rowIdx(tempCombs.begin() + myStrt,
                                              tempCombs.begin() + myEnd);

                for (auto idx: myVecs[i]) {
                    AddComb(uintHash, uintKeyKeeper, rowIdx, cartCombs,
                            maxKey, masterKey, primes[idx], idx);
                }
            }
        } else {
            for (std::size_t j = 0, myStrt = 0, myEnd = i;
                 j < tempKeyKeeper.size(); ++j, myStrt = myEnd, myEnd += i) {

                const std::uint64_t masterKey = tempKeyKeeper[j];
                const std::vector<int> rowIdx(tempCombs.begin() + myStrt,
                                              tempCombs.begin() + myEnd);

                for (auto idx: myVecs[i]) {
                    if (masterKey % static_cast<std::uint64_t>(primes[idx]) != 0) {
                        AddComb(uintHash, uintKeyKeeper, rowIdx, cartCombs,
                                maxKey, masterKey, primes[idx], idx);
                    }
                }
            }
        }
    }

    if ((std::numeric_limits<std::uint64_t>::max() / maxKey) <
        static_cast<std::uint64_t>(primes[myVecs.back().back()])) {
        NeedsMpz = true;
    }
    if (!NeedsMpz && myVecs.size() > 1) {
        uintHash.clear();
        lenGrps.assign(uintKeyKeeper.size(), 0);
        lastCol.reserve(uintKeyKeeper.size() * myVecs.back().size());
        uintHash.reserve(uintKeyKeeper.size() * myVecs.back().size());

        if (IsRep) {
            for (std::size_t j = 0; j < uintKeyKeeper.size(); ++j) {
                const std::uint64_t masterKey = uintKeyKeeper[j];

                for (auto idx: myVecs.back()) {
                    if (uintHash.find(masterKey *
                        static_cast<std::uint64_t>(
                            primes[idx]
                        )) == uintHash.end()) {
                        uintHash.insert(masterKey *
                            static_cast<std::uint64_t>(primes[idx]));
                        lastCol.push_back(idx);
                        ++lenGrps[j];
                    }
                }
            }
        } else {
            for (std::size_t j = 0; j < uintKeyKeeper.size(); ++j) {
                const std::uint64_t masterKey = uintKeyKeeper[j];

                for (auto idx: myVecs.back()) {
                    if (masterKey % static_cast<std::uint64_t>(primes[idx]) != 0) {
                        if (uintHash.find(masterKey * static_cast<std::uint64_t>(primes[idx])) == uintHash.end()) {
                            uintHash.insert(masterKey * static_cast<std::uint64_t>(primes[idx]));
                            lastCol.push_back(idx);
                            ++lenGrps[j];
                        }
                    }
                }
            }
        }
    } else if (myVecs.size() == 1) {
        lastCol = cartCombs;
        cartCombs.clear();
        lenGrps.assign(lastCol.size(), 1);
    } else {
        std::vector<std::string> strKeyKeeper(uintKeyKeeper.size());

        for (std::size_t j = 0; j < uintKeyKeeper.size(); ++j) {
            strKeyKeeper[j] = std::to_string(uintKeyKeeper[j]);
        }

        mpz_class key;
        mpz_class masterKey;
        std::unordered_set<std::string> strHash;

        for (std::size_t size = myVecs.size() - 1; i < size; ++i) {
            const std::vector<int> tempCombs = cartCombs;
            const std::vector<std::string> tempKeyKeeper = strKeyKeeper;

            cartCombs.clear();
            strKeyKeeper.clear();
            strHash.clear();

            cartCombs.reserve(tempKeyKeeper.size() * myVecs[i].size() * (i + 1));
            strKeyKeeper.reserve(tempKeyKeeper.size() * myVecs[i].size());
            strHash.reserve(tempKeyKeeper.size() * myVecs[i].size());

            if (IsRep) {
                for (std::size_t j = 0, myStrt = 0, myEnd = i;
                     j < tempKeyKeeper.size(); ++j, myStrt = myEnd, myEnd += i) {

                    mpz_set_str(masterKey.get_mpz_t(),
                                tempKeyKeeper[j].c_str(), 10);
                    const std::vector<int> rowIdx(tempCombs.begin() + myStrt,
                                                  tempCombs.begin() + myEnd);

                    for (auto idx: myVecs[i]) {
                        AddComb(strHash, strKeyKeeper, rowIdx, cartCombs,
                                key, masterKey, primes[idx], idx);
                    }
                }
            } else {
                for (std::size_t j = 0, myStrt = 0, myEnd = i;
                     j < tempKeyKeeper.size(); ++j, myStrt = myEnd, myEnd += i) {

                    mpz_set_str(masterKey.get_mpz_t(),
                                tempKeyKeeper[j].c_str(), 10);
                    const std::vector<int> rowIdx(tempCombs.begin() + myStrt,
                                                  tempCombs.begin() + myEnd);

                    for (auto idx: myVecs[i]) {
                        if (mpz_divisible_ui_p(masterKey.get_mpz_t(),
                                               primes[idx]) == 0) {
                            AddComb(strHash, strKeyKeeper, rowIdx, cartCombs,
                                    key, masterKey, primes[idx], idx);
                        }
                    }
                }
            }
        }

        strHash.clear();
        lenGrps.assign(strKeyKeeper.size(), 0);
        lastCol.reserve(strKeyKeeper.size() * myVecs.back().size());
        strHash.reserve(strKeyKeeper.size() * myVecs.back().size());

        if (IsRep) {
            for (std::size_t j = 0; j < strKeyKeeper.size(); ++j) {
                mpz_set_str(masterKey.get_mpz_t(),
                            strKeyKeeper[j].c_str(), 10);

                for (auto idx: myVecs.back()) {
                    key = masterKey * primes[idx];
                    const std::string strKey = key.get_str();

                    if (strHash.find(strKey) == strHash.end()) {
                        strHash.insert(strKey);
                        lastCol.push_back(idx);
                        ++lenGrps[j];
                    }
                }
            }
        } else {
            for (std::size_t j = 0; j < strKeyKeeper.size(); ++j) {
                mpz_set_str(masterKey.get_mpz_t(),
                            strKeyKeeper[j].c_str(), 10);

                for (auto idx: myVecs.back()) {
                    if (mpz_divisible_ui_p(masterKey.get_mpz_t(),
                                           primes[idx]) == 0) {
                        key = masterKey * primes[idx];
                        const std::string strKey = key.get_str();

                        if (strHash.find(strKey) == strHash.end()) {
                            strHash.insert(strKey);
                            lastCol.push_back(idx);
                            ++lenGrps[j];
                        }
                    }
                }
            }
        }
    }
}
