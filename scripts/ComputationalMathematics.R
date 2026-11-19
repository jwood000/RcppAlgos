reprex::reprex({
    #'
    #' This document serves as an overview for solving problems common in [Computational Mathematics](https://en.wikipedia.org/wiki/Computational_mathematics). Of note, `primeSieve` and `primeCount` are based on the excellent work by [Kim Walisch](<https://github.com/kimwalisch>).
    #'
    #' ***
    #'
    #' ## `primeSieve`
    #'
    #' The primeSieve function is based on the [Segmented Sieve of Eratosthenes](https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes#Segmented_sieve). As stated in the linked article, the sieve itself is already very efficient. The problem from an efficiency standpoint, is due to the memory requirements. The segmented version overcomes this by only sieving small sections at a time, which greatly facilitates use of the cache.
    #'

    library(RcppAlgos)
    library(microbenchmark)
    options(width = 90)

    microbenchmark(primeSieve(1e6))

    ## Single threaded primes under a billion!!!
    system.time(primeSieve(10^9))

    ## Using 8 threads we can get under 0.5 seconds!!!
    system.time(primeSieve(10^9, nThreads = 8))

    ## Quickly generate large primes over small interval. N.B. The
    ## order for the bounds does not matter.
    options(scipen = 50)
    system.time(myPs <- primeSieve(10^13 + 10^3, 10^13))

    myPs

    ## Object created is small
    object.size(myPs)

    #'
    #' ### Larger primes
    #'
    #' Since version `2.3.0`, we are implementing the cache-friendly improvements for larger primes originally developed by [Tomás Oliveira](<https://sweet.ua.pt/tos/software/prime_sieve.html>).
    #'

    ## Version <= 2.2.0.. i.e. older versions
    system.time(old <- RcppAlgos220::primeSieve(1e15, 1e15 + 1e9))

    invisible(gc())
    ## v2.3.0+ is faster
    system.time(a <- primeSieve(1e15, 1e15 + 1e9))

    invisible(gc())
    ## And using nThreads is much much faster
    system.time(b <- primeSieve(1e15, 1e15 + 1e9, nThreads = 8))

    identical(old, a)

    identical(a, b)

    #'
    #' ## `primeCount`
    #'
    #' `primeCount` gives the number of primes less than or equal to `n`, without constructing a vector of all those primes. It supports Legendre, Meissel, and Lehmer methods, all of which return the same exact count throughout the supported range.
    #'
    #' The Legendre implementation and shared partial-sieve routines build on the excellent work by [Kim Walisch](<https://github.com/kimwalisch/primecount>), which requires OpenMP. One of the original goals was to make parallel prime counting available without requiring OpenMP, including with the default macOS toolchain. For this reason, `RcppAlgos` uses C++ threads for its parallel calculations.
    #'
    #' By default, `method = "auto"` selects a method based on `n` and the number of threads it can use. You can also select a method explicitly, which is useful when comparing them on your own inputs. See `?primeCount` for the selection rules, counting formulas, and memory requirements.
    #'
    #' Calculations are serial by default. Use `nThreads` for parallel computation when the problem is large enough; the actual number of threads is limited by the available hardware and useful work.
    #'

    ## Count the primes less than or equal to a trillion
    system.time(underOneTrillion <- primeCount(10^12))

    underOneTrillion

    ## Select a method explicitly
    primeCount(10^6, method = "legendre")
    primeCount(10^6, method = "meissel")
    primeCount(10^6, method = "lehmer")

    ## Count the primes less than or equal to ten billion
    ## and rely on auto choice for the method
    microbenchmark(primeCount(10^10))

    ## Use multiple threads for larger inputs
    ## again relying on auto choice for the method
    system.time(underOneHundredTrillion <- primeCount(1e14, nThreads = 8))

    underOneHundredTrillion

    system.time(primeCount(1e14, nThreads = 8, method = "legendre"))

    system.time(primeCount(1e14, nThreads = 8, method = "meissel"))

    system.time(primeCount(1e14, nThreads = 8, method = "lehmer"))

    ## Still not as fast as Kim Walisch's primecount library:
    cat(paste(
        system("primecount 1e14 --legendre --time", intern = TRUE),
        collapse = "\n"
    ))

    cat(paste(
        system("primecount 1e14 --meissel --time", intern = TRUE),
        collapse = "\n"
    ))

    ## For some reason, lehmer's method is not as efficient as the other
    ## methods in the primecount library
    cat(paste(
        system("primecount 1e14 --lehmer --time", intern = TRUE),
        collapse = "\n"
    ))

    #'
    #' ## `nthPrime`
    #'
    #' `nthPrime` finds the prime at a given position, starting with `nthPrime(1)` returning 2. For example, `nthPrime(25)` returns 97, whereas `primeCount(25)` returns 9.
    #'
    #' For larger positions, `nthPrime` combines prime counting with a segmented sieve, without storing all preceding primes. Like `primeCount`, it runs serially by default and accepts `nThreads` for parallel computation when useful. See `?nthPrime` for the supported range.
    #'

    nthPrime(25)

    ## Find the millionth prime
    nthPrime(10^6)

    ## Confirm a prime's position
    p <- nthPrime(1000)
    primeCount(p)
    primeCount(p - 1)

    #'
    #' While `nthPrime` is convenient for finding individual primes by position, generating consecutive primes this way repeats the counting and search work for every prime. A better approach is to find the two endpoints with `nthPrime`, then generate the primes between them with `primeSieve`:
    #'

    primes_by_position <- function(first, last) {
        sapply(seq.int(first, last), nthPrime)
    }

    primes_by_interval <- function(first, last) {
        lower <- nthPrime(first)
        upper <- nthPrime(last)
        primeSieve(lower, upper)
    }

    first <- 1e9
    last <- first + 1000

    system.time(by_position <- primes_by_position(first, last))
    system.time(by_interval <- primes_by_interval(first, last))

    identical(by_position, by_interval)

    #'
    #' ## Other Sieving Functions
    #'
    #' `RcppAlgos` comes equipped with several functions for quickly generating essential components for problems common in computational mathematics. All functions below can be executed in parallel by using the argument `nThreads`.
    #'
    #' The following sieving functions (`primeFactorizeSieve`, `divisorsSieve`, `numDivisorSieve`, & `eulerPhiSieve`) are very useful and flexible. Generate components up to a number or between two bounds.
    #'

    ## get the number of divisors for every number from 1 to n
    numDivisorSieve(20)

    ## If you want the complete factorization from 1 to n, use divisorsList
    system.time(allFacs <- divisorsSieve(10^5, namedList = TRUE))

    allFacs[c(4339, 15613, 22080)]


    ## Between two bounds
    primeFactorizeSieve(10^12, 10^12 + 5)


    ## Creating a named object
    eulerPhiSieve(20, namedVector = TRUE)


    system.time(a <- eulerPhiSieve(1e12, 1e12 + 1e7))

    ## Using nThreads for greater efficiency
    system.time(b <- eulerPhiSieve(1e12, 1e12 + 1e7, nThreads = 8))

    identical(a, b)

    #'
    #' ## Vectorized Functions
    #'
    #' There are three very fast vectorized functions for general factoring (e.g. all divisors of number), primality testing, as well as prime factoring (`divisorsRcpp`, `isPrimeRcpp`, `primeFactorize`).
    #'

    ## get result for individual numbers
    primeFactorize(123456789)


    ## or for an entire vector
    ## N.B. The R Version you are using... random sampling
    ## has changed throughout the years
    R.version[["version.string"]]
    set.seed(100)

    myVec <- sample(-100000000:100000000, 5)
    divisorsRcpp(myVec, namedList = TRUE)


    ## Creating a named object
    isPrimeRcpp(995:1000, namedVector = TRUE)

    system.time(a <- primeFactorize(1e12:(1e12 + 1e5)))

    ## Using nThreads for greater efficiency
    system.time(b <- primeFactorize(1e12:(1e12 + 1e5), nThreads = 8))

    identical(a, b)
}, advertise = FALSE, venue = "r", html_preview = FALSE, wd = ".")
