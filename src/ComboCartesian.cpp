#include "Cpp14MakeUnique.h"
#include <unordered_set>
#include <vector>
#include <string>
#include <gmp.h>

void comboGrid(std::vector<int> &cartCombs, bool IsRep,
               const std::vector<std::vector<int>> &myVecs, 
               const std::vector<int> &primes) {

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

    for (; i < myVecs.size(); ++i) {
        if ((std::numeric_limits<std::uint64_t>::max() / maxKey) <
            primes[myVecs[i].back()]) {
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
                    const std::uint64_t key = masterKey * primes[idx];
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
            }
        } else {
            for (std::size_t j = 0, myStrt = 0, myEnd = i;
                 j < tempKeyKeeper.size(); ++j, myStrt = myEnd, myEnd += i) {

                const std::uint64_t masterKey = tempKeyKeeper[j];
                const std::vector<int> rowIdx(tempCombs.begin() + myStrt,
                                              tempCombs.begin() + myEnd);

                for (auto idx: myVecs[i]) {
                    if (masterKey % primes[idx] != 0) {
                        const std::uint64_t key = masterKey * primes[idx];
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
                }
            }
        }
    }

    if (NeedsMpz) {
        std::vector<std::string> strKeyKeeper(uintKeyKeeper.size());

        for (std::size_t j = 0; j < uintKeyKeeper.size(); ++j) {
            strKeyKeeper[j] = std::to_string(uintKeyKeeper[j]);
        }

        mpz_t key;
        mpz_t masterKey;

        mpz_init(key);
        mpz_init(masterKey);
        std::unordered_set<std::string> strHash;

        for (; i < myVecs.size(); ++i) {
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

                    mpz_set_str(masterKey, tempKeyKeeper[j].c_str(), 10);
                    const std::vector<int> rowIdx(tempCombs.begin() + myStrt,
                                                  tempCombs.begin() + myEnd);

                    for (auto idx: myVecs[i]) {
                        mpz_mul_si(key, masterKey, primes[idx]);
                        auto buffer = FromCpp14::make_unique<char[]>(
                            mpz_sizeinbase(key, 10) + 2
                        );
                        mpz_get_str(buffer.get(), 10, key);
                        const std::string strKey = buffer.get();

                        if (strHash.find(strKey) == strHash.end()) {
                            strHash.insert(strKey);
                            cartCombs.insert(cartCombs.end(),
                                             rowIdx.cbegin(),
                                             rowIdx.cend());
                            cartCombs.push_back(idx);
                            strKeyKeeper.push_back(strKey);
                        }
                    }
                }
            } else {
                for (std::size_t j = 0, myStrt = 0, myEnd = i;
                     j < tempKeyKeeper.size(); ++j, myStrt = myEnd, myEnd += i) {

                    mpz_set_str(masterKey, tempKeyKeeper[j].c_str(), 10);
                    const std::vector<int> rowIdx(tempCombs.begin() + myStrt,
                                                  tempCombs.begin() + myEnd);

                    for (auto idx: myVecs[i]) {
                        if (mpz_divisible_ui_p(masterKey, primes[idx]) == 0) {
                            mpz_mul_si(key, masterKey, primes[idx]);
                            auto buffer = FromCpp14::make_unique<char[]>(
                                mpz_sizeinbase(key, 10) + 2
                            );
                            mpz_get_str(buffer.get(), 10, key);
                            const std::string strKey = buffer.get();

                            if (strHash.find(strKey) == strHash.end()) {
                                strHash.insert(strKey);
                                cartCombs.insert(cartCombs.end(),
                                                 rowIdx.cbegin(),
                                                 rowIdx.cend());
                                cartCombs.push_back(idx);
                                strKeyKeeper.push_back(strKey);
                            }
                        }
                    }
                }
            }
        }

        mpz_clear(masterKey);
        mpz_clear(key);
    }
}
