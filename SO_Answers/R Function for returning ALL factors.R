reprex::reprex({
    #' A lot has changed in the R language since this question was originally asked. In version `0.6-3` of the `numbers` package, the function `divisors` was included that is very useful for getting all of the factors of a number. It will meet the needs of most users, however if you are looking for raw speed or you are working with larger numbers, you will need an alternative method. I have authored two packages (partially inspired by this question, I might add) that contain highly optimized functions aimed at problems just like this. The first one is `RcppAlgos` and the other is `RcppBigIntAlgos` (formerly called `bigIntegerAlgos`).
    #'
    #' `RcppAlgos`
    #' =========
    #' `RcppAlgos` contains two functions for obtaining divisors of numbers less than `2^53 - 1` : `divisorsRcpp` (a vectorized function for quickly obtaining the complete factorization of many numbers) & `divisorsSieve` (quickly generates the complete factorization over a range). First up, we factor many random numbers using `divisorsRcpp`:

    library(gmp)  ## for all_divisors by @GeorgeDontas
    library(RcppAlgos)
    library(numbers)
    options(scipen = 999)
    set.seed(42)
    testSamp <- sample(10^10, 10)

    all_divisors <- function(x) {
        sort_listz <- function(z) z <- z[order(as.numeric(z))]
        mult_listz <- function(x,y) do.call('c', lapply(y, function(i) i*x))

        if (abs(x)<=1) return(x)
        else {
            factorsz <- as.bigz(factorize(as.bigz(x)))
            factorsz <- sort_listz(factorsz)
            prime_factorsz <- unique(factorsz)
            prime_ekt <- vapply(prime_factorsz, function(i) sum(factorsz==i), integer(1), USE.NAMES=FALSE)
            spz <- vector()
            all <-1
            n <- length(prime_factorsz)
            for (i in 1:n) {
                pr <- prime_factorsz[i]
                pe <- prime_ekt[i]
                all <- all*(pe+1)

                prz <- as.bigz(pr)
                pse <- vector(mode="raw",length=pe+1)
                pse <- c( as.bigz(1), prz)

                if (pe>1) {
                    for (k in 2:pe) {
                        prz <- prz*pr
                        pse[k+1] <- prz
                    }
                }

                spz <- if (i>1) mult_listz (spz, pse) else pse;
            }
            spz <- sort_listz (spz)
            return (spz)
        }
    }

    ## vectorized so you can pass the entire vector as an argument
    testRcpp <- divisorsRcpp(testSamp)
    testDontas <- lapply(testSamp, all_divisors)

    identical(lapply(testDontas, as.numeric), testRcpp)

    #' And now, factor many numbers over a range using `divisorsSieve`:

    identical(lapply(testDontas, as.numeric), testRcpp)

    system.time(testSieve <- divisorsSieve(10^13, 10^13 + 10^5))

    system.time(testDontasSieve <- lapply((10^13):(10^13 + 10^5), all_divisors))

    identical(lapply(testDontasSieve, asNumeric), testSieve)

    #' Both `divisorsRcpp` and `divisorsSieve` are nice functions that are flexible and efficient, however they are limited to `2^53 - 1`.
    #'
    #' `RcppBigIntAlgos`
    #' =============
    #' The `RcppBigIntAlgos` package (formerly called `bigIntegerAlgos` prior to version 0.2.0) links directly to the [C library gmp](https://gmplib.org/) and features `divisorsBig` which is designed for very large numbers.

    library(RcppBigIntAlgos)
    ## testSamp is defined above... N.B. divisorsBig is not quite as
    ## efficient as divisorsRcpp. This is so because divisorsRcpp
    ## can take advantage of more efficient data types.
    testBig <- divisorsBig(testSamp)

    identical(testDontas, testBig)

    #' Functions for reproducing results:

    FUN <- function(x) {
        x <- as.integer(x)
        div <- seq_len(abs(x))
        factors <- div[x %% div == 0L]
        factors <- list(neg = -factors, pos = factors)
        return(factors)
    }

    get_all_factors <- function(n) {
        prime_factor_tables <- lapply(
            setNames(n, n),
            function(i) {
                if(i == 1) return(data.frame(x = 1L, freq = 1L))
                plyr::count(as.integer(gmp::factorize(i)))
            }
        )
        lapply(
            prime_factor_tables,
            function(pft) {
                powers <- plyr::alply(pft, 1, function(row) row$x ^ seq.int(0L, row$freq))
                power_grid <- do.call(expand.grid, powers)
                sort(unique(apply(power_grid, 1, prod)))
            }
        )
    }

    #' And here are the benchmark as defined in my original post (N.B. `MyFactors` is replaced by `divisorsRcpp` and `divisorsBig`).

    library(rbenchmark)
    set.seed(199)
    samp <- sample(10^9, 10^5)
    benchmark(RcppAlgos=divisorsRcpp(samp),
              RcppBigIntAlgos=divisorsBig(samp),
              DontasDivs=lapply(samp, all_divisors),
              replications=10,
              columns = c("test", "replications", "elapsed", "relative"),
              order = "relative")

    set.seed(97)
    samp <- sample(10^6, 10^4)
    benchmark(
        RcppAlgos       = divisorsRcpp(samp),
        RcppBigIntAlgos = divisorsBig(samp),
        numbers         = lapply(samp, divisors), ## From the numbers package
        DontasDivs      = lapply(samp, all_divisors),
        CottonDivs      = lapply(samp, get_all_factors),
        ChaseDivs       = lapply(samp, FUN),
        replications    = 5,
        columns         = c("test", "replications", "elapsed", "relative"),
        order           = "relative"
    )

    #' The next benchmarks demonstrate the true power of the underlying algorithm in the `divisorsBig` function. The number being factored is a power of `10`, so the prime factoring step can almost be completely ignored (e.g. `system.time(factorize(pow.bigz(10,30)))` registers `0` on my machine). Thus, the difference in timing is due solely to how quickly the prime factors can be combined to produce all factors.

    library(microbenchmark)
    powTen <- pow.bigz(10, 30)

    microbenchmark(
        algos  = divisorsBig(powTen),
        Dontas = all_divisors(powTen),
        unit   = "relative"
    )

    ## Negative numbers show an even greater increase in efficiency
    negPowTen <- powTen * -1

    microbenchmark(
        algos  = divisorsBig(negPowTen),
        Dontas = all_divisors(negPowTen),
        unit   = "relative"
    )

    #' Very Large Numbers
    #' ----------
    #' With `divisorsBig`, obtaining the complete factorization with very large inputs is no problem. The algorithm dynamically adjusts based off of the input and applies different algorithms in different situations. We can also take advantage of multithreading if [Lenstra's Elliptic Curve method](https://en.wikipedia.org/wiki/Lenstra_elliptic-curve_factorization) or the [Quadratic Sieve](https://en.wikipedia.org/wiki/Quadratic_sieve) is utilized.
    #'
    #' Here are some examples using `n5` and `n9` defined in [this answer](https://stackoverflow.com/a/36537823/4408538).

    n5 <- as.bigz("94968915845307373740134800567566911")
    system.time(print(divisorsBig(n5)))

    n9 <- prod(nextprime(urand.bigz(2, 82, 42)))
    system.time(print(divisorsBig(n9, nThreads = 4)))

    #' Here is an example provided by @Dontas with one large prime and one smaller prime:

    x <- pow.bigz(2, 256) + 1
    divisorsBig(x, showStats = TRUE, nThreads = 8)

    #' Compare this to finding the prime factorization using `gmp::factorize`:

    system.time(factorize(x))

    #' Lastly, here is an example with a large semiprime (N.B. since we know it's a semiprime, we skip the extended Pollard's rho algorithm as well as Lentra's elliptic curve method).

    ## https://members.loria.fr/PZimmermann/records/rsa.html
    rsa79 <- as.bigz("7293469445285646172092483905177589838606665884410340391954917800303813280275279")
    divisorsBig(
        rsa79, nThreads = 8, showStats = TRUE,
        skipPolRho = TRUE, skipECM = TRUE
    )
}, advertise = FALSE, venue = "r", html_preview = FALSE, wd = ".")
