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
    #' Since version `2.3.0`, we are implementing the cache-friendly improvements for larger primes originally developed by [Tomás Oliveira](<http://sweet.ua.pt/tos/software/prime_sieve.html>).
    #'

    ## Version <= 2.2.0.. i.e. older versions
    system.time(old <- RcppAlgos220::primeSieve(1e15, 1e15 + 1e9))

    invisible(gc())
    ## v2.3.0+ is faster
    system.time(a <- primeSieve(1e15, 1e15 + 1e9))

    invisible(gc())
    ## And nThreads we much faster
    system.time(b <- primeSieve(1e15, 1e15 + 1e9, nThreads = 8))

    identical(a, b)

    #'
    #' ## `primeCount`
    #'
    #' The library by Kim Walisch relies on [OpenMP](<https://en.wikipedia.org/wiki/OpenMP>) for parallel computation with [Legendre's Formula](<http://mathworld.wolfram.com/LegendresFormula.html>). Currently, the default compiler on `macOS` is `clang`, which does not support `OpenMP`. James Balamuta (a.k.a. TheCoatlessProfessor... well at least [we think so](<https://thecoatlessprofessor.com/about/>)) has written a great article on this topic, which you can find here: <https://thecoatlessprofessor.com/programming/openmp-in-r-on-os-x/>. One of the goals of `RcppAlgos` is to be accessible by all users. With this in mind, we set out to count primes in parallel _without_ `OpenMP`.
    #'
    #' At first glance, this seems trivial as we have a function in `Primes.cpp` called `phiWorker` that counts the primes up to `x`. If you look in [phi.cpp](<https://github.com/kimwalisch/primecount/blob/master/src/phi.cpp>) in the `primecount` library by Kim Walisch, we see that `OpenMP` does its magic on a for loop that makes repeated calls to `phi` (which is what `phiWorker` is based on). All we need to do is break this loop into _n_ intervals where _n_ is the number of threads. Simple, right?
    #'
    #' We can certainly do this, but what you will find is that _n - 1_ threads will complete very quickly and the _n<sup>th</sup>_ thread will be left with a heavy computation. In order to alleviate this unbalanced load, we divide the loop mentioned above into smaller intervals. The idea is to completely calculate `phi` up to a limit _m_ using all _n_ threads and then gradually increase _m_. The advantage here is that we are benefiting greatly from the caching done by the work of the previous _n_ threads.
    #'
    #' With this is mind, here are some results:
    #'

    ## Enumerate the number of primes below trillion
    system.time(underOneTrillion <- primeCount(10^12))

    underOneTrillion

    ## Enumerate the number of primes below ten billion in 2 milliseconds
    microbenchmark(primeCount(10^10))

    system.time(underOneHundredTrillion <- primeCount(1e14, nThreads = 8))

    underOneHundredTrillion

    ## Still not as fast as Kim Walisch's primecount library:
    system("primecount 1e14 --legendre --time")

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

