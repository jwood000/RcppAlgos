#include <Rcpp.h>
#include <gmp.h>
#include <unordered_set>
#include "Cpp14MakeUnique.h"

constexpr int base10 = 10;

void comboGrid(std::vector<int> &cartCombs, bool IsRep,
               const std::vector<std::vector<int>> &myVec, 
               const std::vector<int> &primes) {
    
    std::unordered_set<uint64_t> uintHash;
    std::vector<uint64_t> uintKeyKeeper;
    uint64_t maxKey = 0;
    
    for (const auto ind: myVec[0]) {
        uint64_t key = primes[ind];
        
        if (key > maxKey) {
            maxKey = key;
        }
        
        const auto it = uintHash.find(key);
        
        if (it == uintHash.end()) {
            uintHash.insert(key);
            cartCombs.push_back(ind);
            uintKeyKeeper.push_back(key);
        }
    }
    
    std::size_t i = 1;
    bool NeedsMpz = false;

    for (; i < myVec.size(); ++i) {
        if ((std::numeric_limits<uint64_t>::max() / maxKey) < primes[myVec[i].back()]) {
            NeedsMpz = true;
            break;
        }

        const std::vector<int> tempCombs = cartCombs;
        const std::vector<uint64_t> tempKeyKeeper = uintKeyKeeper;

        cartCombs.clear();
        uintKeyKeeper.clear();
        uintHash.clear();
        
        cartCombs.reserve(tempKeyKeeper.size() * myVec[i].size() * (i + 1));
        uintKeyKeeper.reserve(tempKeyKeeper.size() * myVec[i].size());
        uintHash.reserve(tempKeyKeeper.size() * myVec[i].size());
        
        if (IsRep) {
            for (std::size_t j = 0, myStrt = 0, myEnd = i;
                 j < tempKeyKeeper.size(); ++j, myStrt = myEnd, myEnd += i) {
    
                const uint64_t masterKey = tempKeyKeeper[j];
                const std::vector<int> masterVec(tempCombs.begin() + myStrt, tempCombs.begin() + myEnd);
    
                for (const auto ind: myVec[i]) {
                    uint64_t key = masterKey * primes[ind];
                    const auto it = uintHash.find(key);
                    
                    if (key > maxKey) {
                        maxKey = key;
                    }
    
                    if (it == uintHash.end()) {
                        std::vector<int> vecVal = masterVec;
                        vecVal.push_back(ind);
                        uintHash.insert(key);
                        cartCombs.insert(cartCombs.end(), vecVal.cbegin(), vecVal.cend());
                        uintKeyKeeper.push_back(key);
                    }
                }
            }
        } else {
            for (std::size_t j = 0, myStrt = 0, myEnd = i;
                 j < tempKeyKeeper.size(); ++j, myStrt = myEnd, myEnd += i) {
                
                const uint64_t masterKey = tempKeyKeeper[j];
                const std::vector<int> masterVec(tempCombs.begin() + myStrt, tempCombs.begin() + myEnd);
                
                for (const auto ind: myVec[i]) {
                    if (masterKey % primes[ind] != 0) {
                        uint64_t key = masterKey * primes[ind];
                        const auto it = uintHash.find(key);
                        
                        if (key > maxKey) {
                            maxKey = key;
                        }
                        
                        if (it == uintHash.end()) {
                            std::vector<int> vecVal = masterVec;
                            vecVal.push_back(ind);
                            uintHash.insert(key);
                            cartCombs.insert(cartCombs.end(), vecVal.cbegin(), vecVal.cend());
                            uintKeyKeeper.push_back(key);
                        }
                    }
                }
            }
        }
    }
    
    if (NeedsMpz) {
        std::vector<std::string> strKeyKeeper(uintKeyKeeper.size());
        
        for (std::size_t j = 0; j < uintKeyKeeper.size(); ++j)
            strKeyKeeper[j] = std::to_string(uintKeyKeeper[j]);
        
        std::unordered_set<std::string> strHash;
        
        mpz_t masterKey, key;
        mpz_init(masterKey); mpz_init(key);
        
        for (; i < myVec.size(); ++i) {
            const std::vector<int> tempCombs = cartCombs;
            const std::vector<std::string> tempKeyKeeper = strKeyKeeper;
            
            cartCombs.clear();
            strKeyKeeper.clear();
            strHash.clear();
            
            cartCombs.reserve(tempKeyKeeper.size() * myVec[i].size() * (i + 1));
            strKeyKeeper.reserve(tempKeyKeeper.size() * myVec[i].size());
            strHash.reserve(tempKeyKeeper.size() * myVec[i].size());
            
            if (IsRep) {
                for (std::size_t j = 0, myStrt = 0, myEnd = i;
                     j < tempKeyKeeper.size(); ++j, myStrt = myEnd, myEnd += i) {
                    
                    mpz_set_str(masterKey, tempKeyKeeper[j].c_str(), base10);
                    const std::vector<int> masterVec(tempCombs.begin() + myStrt, tempCombs.begin() + myEnd);
                    
                    for (const auto ind: myVec[i]) {
                        mpz_mul_si(key, masterKey, primes[ind]);
                        auto buffer = FromCpp14::make_unique<char[]>(mpz_sizeinbase(key, base10) + 2);
                        mpz_get_str(buffer.get(), base10, key);
                        const std::string strKey = buffer.get();
                        const auto it = strHash.find(strKey);
    
                        if (it == strHash.end()) {
                            std::vector<int> vecVal = masterVec;
                            vecVal.push_back(ind);
                            strHash.insert(strKey);
                            cartCombs.insert(cartCombs.end(), vecVal.cbegin(), vecVal.cend());
                            strKeyKeeper.push_back(strKey);
                        }
                    }
                }
            } else {
                for (std::size_t j = 0, myStrt = 0, myEnd = i;
                     j < tempKeyKeeper.size(); ++j, myStrt = myEnd, myEnd += i) {
                    
                    mpz_set_str(masterKey, tempKeyKeeper[j].c_str(), base10);
                    const std::vector<int> masterVec(tempCombs.begin() + myStrt, tempCombs.begin() + myEnd);
                    
                    for (const auto ind: myVec[i]) {
                        if (mpz_divisible_ui_p(masterKey, primes[ind]) == 0) {
                            mpz_mul_si(key, masterKey, primes[ind]);
                            auto buffer = FromCpp14::make_unique<char[]>(mpz_sizeinbase(key, base10) + 2);
                            mpz_get_str(buffer.get(), base10, key);
                            const std::string strKey = buffer.get();
                            const auto it = strHash.find(strKey);
                            
                            if (it == strHash.end()) {
                                std::vector<int> vecVal = masterVec;
                                vecVal.push_back(ind);
                                strHash.insert(strKey);
                                cartCombs.insert(cartCombs.end(), vecVal.cbegin(), vecVal.cend());
                                strKeyKeeper.push_back(strKey);
                            }
                        }
                    }
                }
            }
        }
        
        mpz_clear(masterKey); mpz_clear(key);
    }
}
