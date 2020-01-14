#include "PollardRhoUtils.h"
#include "CleanConvert.h"

template <typename typeReturn>
void FactorTrialDivision(std::int64_t &t, std::vector<typeReturn> &factors) {
    
    int p = 3;
    
    while ((t & 1) == 0) {
        factors.push_back(2);
        t >>= 1;
    }
    
    for (std::size_t i = 1; i < primesDiffPR.size();) {
        if ((t % p) != 0) {
            p += primesDiffPR[i++];
            if (t < (p * p))
                break;
        } else {
            t /= p;
            factors.push_back(static_cast<typeReturn>(p));
        }
    }
    
    if ((t % p) == 0) {
        t /= p;
        factors.push_back(static_cast<typeReturn>(p));
    }
}

std::int64_t PositiveMod(std::int64_t i, std::int64_t n) {
    return ((i % n) + n) % n;
}

std::int64_t myGCD(std::int64_t u, std::int64_t v) {
    
    u = PositiveMod(u, v);
    
    while (v != 0) {
        std::int64_t r = u % v;
        u = v;
        v = r;
    }
    
    return u;
}

std::int64_t ProdBigMod(std::int64_t x1_i64, std::int64_t x2_i64, std::int64_t p_i64) {
    
    double prodX = static_cast<double>(x1_i64) * static_cast<double>(x2_i64);
    std::int64_t result = 0;
    
    if (prodX < static_cast<double>(p_i64)) {
        result = prodX;
    } else if (p_i64 < Sqrt63Max || prodX < my63Max) {
        result = (x1_i64 * x2_i64) % p_i64;
    } else {
        std::int64_t nChunks = 1;
        std::int64_t chunk;
        double part1 = my63Max;
        
        while (part1 >= my63Max) {
            std::int64_t cSize = static_cast<std::int64_t>(my63Max / x1_i64);
            chunk = (x1_i64 * cSize) % p_i64;
            nChunks = x2_i64 / cSize;
            std::int64_t part2 = ((x2_i64 - cSize * nChunks) * x1_i64) % p_i64;
            part1 = static_cast<double>(nChunks) * static_cast<double>(chunk);
            x1_i64 = chunk;
            x2_i64 = nChunks;
            result = (result + part2) % p_i64;
        }
        
        std::int64_t part1_i64 = (nChunks * chunk) % p_i64;
        result = (part1_i64 + result) % p_i64;
    }
    
    return result;
}

std::int64_t ExpBySquaring(std::int64_t x, std::int64_t n, std::int64_t p) {
    
    std::int64_t result;
    
    if (n == 1) {
        result = PositiveMod(x, p);
    } else if (n % 2 == 0) {
        result = ExpBySquaring(ProdBigMod(x, x, p), n / 2, p);
    } else {
        result = ProdBigMod(x, ExpBySquaring(ProdBigMod(x, x, p), (n - 1) / 2, p), p);
    }
    
    return result;
}

int MillerRabin(std::int64_t n, std::int64_t nm1, std::int64_t x,
                std::int64_t& y, std::int64_t q, std::uint64_t k) {
    
    y = ExpBySquaring(x, q, n);
    if (y == 1 || y == nm1)
        return 1;
    
    for (std::size_t i = 1; i < k; ++i) {
        y = ExpBySquaring(y, 2, n);
        if (y == nm1)
            return 1;
        if (y == 1)
            return 0;
    }
    
    return 0;
}

int IsPrime(std::int64_t n) {
    
    int k, primeTestReturn;
    std::int64_t nm1, q, tmp, a;
    
    std::vector<std::int64_t> factors;
    
    if (n < 2)
        return 0;
    
    /* We have already casted out small primes. */
    if (n < FirstOmittedPrime * FirstOmittedPrime)
        return 1;
    
    /* Precomputation for Miller-Rabin.  */
    q = nm1 = n - 1;
    
    /* Find q and k, where q is odd and n = 1 + 2**k * q.  */
    k = 0;
    while ((q & 1) == 0) {
        q /= 2;
        ++k;
    }
    
    a = 2;
    
    /* Perform a Miller-Rabin test, finds most composites quickly.  */
    if (!MillerRabin (n, nm1, a, tmp, q, (std::uint64_t) k)) {
        primeTestReturn = 0;
        goto ret2;
    }
    
    /* Factor n-1 for Lucas.  */
    tmp = nm1;
    getPrimeFactors(tmp, factors);
    
    /* Loop until Lucas proves our number prime, or Miller-Rabin proves our
    number composite.  */
    for (std::size_t r = 0; r < primesDiffPR.size(); ++r) {
        
        primeTestReturn = 1;
        
        for (std::size_t i = 0; i < factors.size() && primeTestReturn; ++i) {
            tmp = nm1 / factors[i];
            tmp = ExpBySquaring(a, tmp, n);
            primeTestReturn = (tmp != 1);
        }
        
        if (primeTestReturn)
            goto ret1;
        
        a += primesDiffPR[r];	/* Establish new base.  */
            
        if (!MillerRabin(n, nm1, a, tmp, q, (std::uint64_t) k)) {
            primeTestReturn = 0;
            goto ret1;
        }
    }
    
    Rcpp::stop("Lucas prime test failure. This should not happen");
    ret1:
        factors.resize(0);
        
    ret2:
        return primeTestReturn;
}

void PollardRhoMpzT(mpz_t n, std::size_t a, std::vector<double> &factors) {
    
    mpz_t x, z, y, P, t;
    std::size_t k, q;
    
    mpz_init(t);
    mpz_init_set_si(y, 2);
    mpz_init_set_si(x, 2);
    mpz_init_set_si(z, 2);
    mpz_init_set_ui(P, 1);
    k = q = 1;
    
    while(mpz_cmp_ui(n, 1) != 0) {
        for (;;) {
            do {
                mpz_mul(t, x, x);
                mpz_mod(x, t, n);
                mpz_add_ui(x, x, a);
                
                mpz_sub(t, z, x);
                mpz_mul(t, P, t);
                mpz_mod(P, t, n);
                
                if (k % 32 == 1) {
                    mpz_gcd(t, P, n);
                    if (mpz_cmp_ui(t, 1) != 0)
                        goto factor_found;
                    mpz_set(y, x);
                }
                
            } while (--k != 0);
            
            mpz_set(z, x);
            k = q;
            q <<= 1;
            
            for (std::size_t i = 0; i < k; ++i) {
                mpz_mul(t, x, x);
                mpz_mod(x, t, n);
                mpz_add_ui(x, x, a);
            }
            
            mpz_set(y, x);
        }
        
        factor_found:
            do {
                mpz_mul(t, y, y);
                mpz_mod(y, t, n);
                mpz_add_ui(y, y, a);
                
                mpz_sub(t, z, y);
                mpz_gcd(t, t, n);
                
            } while (mpz_cmp_ui(t, 1) == 0);
        
        mpz_divexact(n, n, t);	/* divide by t, before t is overwritten */

        if (mpz_probab_prime_p(t, MR_REPS) == 0) {
            PollardRhoMpzT(t, a + 1, factors);
        } else {
            double dblT = mpz_get_d(t);
            factors.push_back(dblT);
            
            while (mpz_divisible_p(n, t)) {
                mpz_divexact(n, n, t);
                factors.push_back(dblT);
            }
        }
        
        if (mpz_probab_prime_p(n, MR_REPS) != 0) {
            factors.push_back(mpz_get_d(n));
            break;
        }
        
        mpz_mod(x, x, n);
        mpz_mod(z, z, n);
        mpz_mod(y, y, n);
    }
    
    mpz_clear(P);
    mpz_clear(t);
    mpz_clear(z);
    mpz_clear(x);
    mpz_clear(y);
}

void PollardRho(std::int64_t n, std::int64_t a, std::vector<int>& factors) {
    
    std::int64_t x, z, y, P, t;
    std::size_t k = 1u;
    std::size_t q = 1u;
    
    y = x = z = 2;
    P = 1;
    
    while (n != 1) {
        for (;;) {
            do {
                x *= x;
                x %= n;
                x += a;
                
                t = z - x;
                t = PositiveMod(t, n);
                P *= t;
                P %= n;
                
                if (k % 32 == 1) {
                    t = myGCD(P, n);
                    if (t != 1)
                        goto factor_found;
                    y = x;
                }
                
            } while (--k != 0);
            
            z = x;
            k = q;
            q <<= 1;
            
            for (std::size_t i = 0; i < k; ++i) {
                x *= x;
                x %= n;
                x += a;
            }
            
            y = x;
        }
        
        factor_found:
            do {
                y *= y;
                y %= n;
                y += a;
                t = myGCD(z - y, n);
                
            } while (t == 1);
        
        n /= t;	/* divide by t, before t is overwritten */
        
        if (!IsPrime(t)) {
            PollardRho(t, a + 1, factors);
        } else {
            factors.push_back(static_cast<int>(t));
            
            while ((n % t) == 0) {
                n /= t;
                factors.push_back(static_cast<int>(t));
            }
        }
        
        if (IsPrime(n)) {
            factors.push_back(static_cast<int>(n));
            break;
        }
        
        x %= n;
        z %= n;
        y %= n;
    }
}

template <typename typeReturn>
void getPrimeFactors(std::int64_t& t, std::vector<typeReturn>& factors) {
    FactorTrialDivision(t, factors);
    
    if (t > 1) {
        if (t < std::numeric_limits<int>::max()) {
            if (IsPrime(t)) {
                factors.push_back(t);
            } else {
                std::vector<int> intFactors;
                PollardRho(t, 1, intFactors);
                factors.insert(factors.end(), intFactors.cbegin(), intFactors.cend());
            }
        } else {
            mpz_t bigT;
            mpz_init(bigT);
            mpz_set_d(bigT, static_cast<double>(t));
            
            if (mpz_probab_prime_p(bigT, MR_REPS)) {
                factors.push_back(t);
            } else {
                std::vector<double> dblFactors;
                PollardRhoMpzT(bigT, 1, dblFactors);
                factors.insert(factors.end(),
                               std::make_move_iterator(dblFactors.cbegin()),
                               std::make_move_iterator(dblFactors.cend()));
            }
            
            mpz_clear(bigT);
        }
    }
    
    std::sort(factors.begin(), factors.end());
}

template void getPrimeFactors(std::int64_t&, std::vector<int>&);
template void getPrimeFactors(std::int64_t&, std::vector<double>&);

