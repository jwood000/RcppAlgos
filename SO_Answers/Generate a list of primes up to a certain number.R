reprex::reprex({

    #' ## Prime Numbers in R
    #'
    #' The OP asked to generate all prime numbers below one billion. All of the answers provided thus far are either not capable of doing this, will take a long a time to execute, or currently not available in R (see the [answer](https://stackoverflow.com/a/3790109/4408538) by @Charles). The package `RcppAlgos` (I am the author) is capable of generating the requested output in just under `1 second` using only one thread.  It is based off of the segmented sieve of Eratosthenes by [Kim Walisch](https://github.com/kimwalisch/primesieve).
    #'
    #' ### `RcppAlgos`

    ser <- system.time(RcppAlgos::primeSieve(1e9))  ## using 1 thread
    ser

    #' ### Using Multiple Threads
    #'
    #' And in recent versions (i.e. `>= 2.3.0`), we can utilize multiple threads for even faster generation. For example, now we can generate the primes up to **1 billion in under a quarter of a second!**

    par <- system.time(RcppAlgos::primeSieve(10^9, nThreads = 8))
    par

    #' ### Summary of Available Packages in R for Generating Primes
    #'
    #'   1. `schoolmath`
    #'   2. `primefactr`
    #'   3. `sfsmisc`
    #'   4. `primes`
    #'   5. `numbers`
    #'   6. `spuRs`
    #'   7. `randtoolbox`
    #'   8. `matlab`
    #'   9. `sieve` from @John (see below.. we removed `if(n > 1e8) stop("n too large")`)

    sieve <- function(n) {
        n <- as.integer(n)
        primes <- rep(TRUE, n)
        primes[1] <- FALSE
        last.prime <- 2L
        fsqr <- floor(sqrt(n))
        while (last.prime <= fsqr) {
            primes[seq.int(2L*last.prime, n, last.prime)] <- FALSE
            sel <- which(primes[(last.prime+1):(fsqr+1)])
            if(any(sel)){
                last.prime <- last.prime + min(sel)
            }else last.prime <- fsqr+1
        }
        which(primes)
    }

    #' Before we begin, we note that the problems pointed out by @Henrik in `schoolmath` still exists. Observe:

    ## 1 is NOT a prime number
    schoolmath::primes(start = 1, end = 20)

    ## This should return 1, however it is saying that 52
    ##  "prime" numbers less than 10^4 are divisible by 7!!
    sum(schoolmath::primes(start = 1, end = 10^4) %% 7L == 0)

    #' The point is, don't use `schoolmath` for generating primes at this point (I have filed an issue with the maintainer).
    #'
    #' Let's look at `randtoolbox` as it appears to be incredibly efficient. Observe:

    library(microbenchmark)
    options(digits = 4)
    options(width = 90)

    ## the argument for get.primes is for how many prime numbers you need
    ## whereas most packages get all primes less than a certain number
    microbenchmark(randtoolbox = randtoolbox::get.primes(78498),
                   RcppAlgos = RcppAlgos::primeSieve(10^6),
                   unit = "relative")

    #' A closer look reveals that it is essentially a lookup table (found in the file `randtoolbox.c` from the source code).

    #include "primes.h"
    #
    #     void reconstruct_primes()
    #     {
    #         int i;
    #         if (primeNumber[2] == 1)
    #             for (i = 2; i < 100000; i++)
    #                 primeNumber[i] = primeNumber[i-1] + 2*primeNumber[i];
    #     }

    #' Where `primes.h` is a header file that contains an array of _"halves of differences between prime numbers"_. Thus, you will be limited by the number of elements in that array for generating primes (i.e. the first one hundred thousand primes). If you are only working with smaller primes (less than `1,299,709` (i.e. the 100,000th prime)) and you are working on a project that requires the `nth` prime, `randtoolbox` is not a bad option.
    #'
    #' Below, we perform benchmarks on the rest of the packages.
    #'
    #' ### Primes up to One Million

    million <- microbenchmark(
        RcppAlgos = RcppAlgos::primeSieve(10^6),
        numbers = numbers::Primes(10^6),
        spuRs = spuRs::primesieve(c(), 2:10^6),
        primes = primes::generate_primes(1, 10^6),
        primefactr = primefactr::AllPrimesUpTo(10^6),
        sfsmisc = sfsmisc::primes(10^6),
        matlab = matlab::primes(10^6),
        JohnSieve = sieve(10^6),
        unit = "relative"
    )

    million

    #' ### Primes up to Ten Million

    ten_million <- microbenchmark(
        RcppAlgos = RcppAlgos::primeSieve(10^7),
        numbers = numbers::Primes(10^7),
        spuRs = spuRs::primesieve(c(), 2:10^7),
        primes = primes::generate_primes(1, 10^7),
        primefactr = primefactr::AllPrimesUpTo(10^7),
        sfsmisc = sfsmisc::primes(10^7),
        matlab = matlab::primes(10^7),
        JohnSieve = sieve(10^7),
        unit = "relative",
        times = 20
    )

    ten_million

    #' ### Primes up to One Hundred Million
    #'
    #' For the next two benchmarks, we only consider `RcppAlgos`, `sfsmisc`, `primes`, and the `sieve` function by @John.

    hundred_million <- microbenchmark(
        RcppAlgos = RcppAlgos::primeSieve(10^8),
        sfsmisc = sfsmisc::primes(10^8),
        primes = primes::generate_primes(1, 10^8),
        JohnSieve = sieve(10^8),
        unit = "relative",
        times = 20
    )

    hundred_million

    #' ### Primes up to One Billion
    #'
    #' N.B. We must remove the condition `if(n > 1e8) stop("n too large")` in the `sieve` function.

    ## See top section
    ## system.time(primeSieve(10^9))

    invisible(gc())
    pm_1e9 <- system.time(primes::generate_primes(1, 10^9))
    pm_1e9

    invisible(gc())
    sieve_1e9 <- system.time(sieve(10^9))
    sieve_1e9

    invisible(gc())
    sfs_1e9 <- system.time(sfsmisc::primes(10^9))
    sfs_1e9

    #' From these comparison, we see that `RcppAlgos` scales much better as _n_ gets larger.

    suppressWarnings(suppressMessages(library(dplyr)))
    suppressWarnings(suppressMessages(library(purrr)))

    billion <- tibble(
        expr = c("RcppAlgos", "primes", "sfsmisc", "JohnSieve"),
        time = c(ser["elapsed"], pm_1e9["elapsed"],
                  sfs_1e9["elapsed"], sieve_1e9["elapsed"])
    )

    my_scale <- \(x) x / min(x, na.rm = TRUE)

    time_table <- map(
        list(million, ten_million, hundred_million, billion),
        ~ .x %>% group_by(expr) %>% summarise(med = median(time))
    ) %>%
        reduce(left_join, by = join_by(expr)) %>%
        rename(`1e6` = 2, `1e7` = 3, `1e8` = 4, `1e9` = 5) %>%
        mutate(across(2:5, ~ my_scale(.x)))

    knitr::kable(
        time_table %>%
            mutate(across(2:5, ~ round(my_scale(.x), 2)))
    )

    #' The difference is even more dramatic when we utilize multiple threads:

    algos_time <- list(
        algos_1e6 = microbenchmark(
            ser = RcppAlgos::primeSieve(1e6),
            par = RcppAlgos::primeSieve(1e6, nThreads = 8), unit = "relative"
        ),
        algos_1e7 = microbenchmark(
            ser = RcppAlgos::primeSieve(1e7),
            par = RcppAlgos::primeSieve(1e7, nThreads = 8), unit = "relative"
        ),
        algos_1e8 = microbenchmark(
            ser = RcppAlgos::primeSieve(1e8),
            par = RcppAlgos::primeSieve(1e8, nThreads = 8),
            unit = "relative", times = 20
        ),
        algos_1e9 = microbenchmark(
            ser = RcppAlgos::primeSieve(1e9),
            par = RcppAlgos::primeSieve(1e9, nThreads = 8),
            unit = "relative", times = 10
        )
    )

    ser_vs_par <- map(
        algos_time, ~ .x %>% group_by(expr) %>% summarise(med = median(time))
    ) %>%
        reduce(left_join, by = join_by(expr)) %>%
        rename(`1e6` = 2, `1e7` = 3, `1e8` = 4, `1e9` = 5) %>%
        mutate(across(2:5, ~ my_scale(.x))) %>%
        filter(expr == "ser") %>%
        unlist()

    #' And multiplying the table above by the respective median times for the serial results:

    time_parallel <- bind_cols(
        expr = time_table[, 1],
        map(2:5, \(x) time_table[, x] * ser_vs_par[x]),
    ) %>%
        add_row(expr = "RcppAlgos-Par",
                `1e6` = 1, `1e7` = 1, `1e8` = 1, `1e9` = 1) %>%
        arrange(`1e6`)

    knitr::kable(
        time_parallel %>%
            mutate(across(2:5, ~ round(my_scale(.x), 2)))
    )

    #' ### Primes Over a Range

    microbenchmark(priRcppAlgos = RcppAlgos::primeSieve(10^9, 10^9 + 10^6),
                   priNumbers = numbers::Primes(10^9, 10^9 + 10^6),
                   priPrimes = primes::generate_primes(10^9, 10^9 + 10^6),
                   unit = "relative", times = 20)

    #' ### Primes up to 10 billion in Under 3 Seconds

    ##  primes less than 10 billion
    system.time(tenBillion <- RcppAlgos::primeSieve(10^10, nThreads = 8))

    length(tenBillion)

    ## Warning!!!... Large object created
    tenBillionSize <- object.size(tenBillion)
    print(tenBillionSize, units = "Gb")

    rm(tenBillionSize)
    invisible(gc())

    #' ### Primes Over a Range of Very Large Numbers:
    #'
    #' Prior to version `2.3.0`, we were simply using the same algorithm for numbers of every magnitude. This is okay for smaller numbers when most of the sieving primes have at least one multiple in each segment (Generally, the segment size is limited by the size of `L1 Cache ~32KiB`). However, when we are dealing with larger numbers, the sieving primes will contain many numbers that will have fewer than one multiple per segment. This situation creates a lot of overhead, as we are performing many worthless checks that pollutes the cache. Thus, we observe much slower generation of primes when the numbers are very large. If you want to test yourself, see [Installing older version of R package](https://stackoverflow.com/a/17082609/4408538)).
    #'
    #' In later versions (>= `2.3.0`), we are using the cache friendly improvement originally developed by [Tomás Oliveira](<http://sweet.ua.pt/tos/software/prime_sieve.html>). The improvements are drastic:

    ## Over 3x faster than older versions
    system.time(cacheFriendly <- RcppAlgos::primeSieve(1e15, 1e15 + 1e9))

    ##  Over 8x faster using multiple threads
    system.time(RcppAlgos::primeSieve(1e15, 1e15 + 1e9, nThreads = 8))

    #' ### Take Away
    #'
    #' 1. There are many great packages available for generating primes
    #' 2. If you are looking for speed in general, there is no match to `RcppAlgos::primeSieve`, especially for larger numbers.
    #' 3. If you need primes in a range, the packages `numbers`, `primes`, & `RcppAlgos` are the way to go.
    #' 4. The importance of good programming practices cannot be overemphasized (e.g. vectorization, using correct data types, etc.). This is most aptly demonstrated by the pure base R solution provided by @John. It is concise, clear, and very efficient.

}, advertise = FALSE, venue = "r", html_preview = FALSE, wd = ".")