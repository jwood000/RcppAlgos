# Run from the repository root:
# Rscript --vanilla .github/development/prime-count-dispatch-check.R [temporary_directory]
# Exercises the real R wrapper, native adapter, converters and cpp11 boundary.
# Counting kernels are replaced by recorded sentinels; this checks dispatch,
# argument forwarding and result conversion, not mathematics or worker safety.
main <- function() {
    repository <- normalizePath(".", winslash = "/")
    arguments <- commandArgs(TRUE)
    if (length(arguments) > 1L) stop("Expected at most one temporary output directory")
    output <- if (length(arguments)) arguments[[1L]] else tempfile("prime-count-dispatch-")
    dir.create(output, recursive = TRUE, showWarnings = FALSE)
    output <- normalizePath(output, winslash = "/")

    sources <- file.path(repository, "src", c(
        "PrimeCount.cpp", "StdConvert.cpp", "GmpConvert.cpp", "ConvertUtils.cpp"
    ))
    copies <- file.path(output, basename(sources))
    stopifnot(all(file.copy(sources, copies, overwrite = TRUE)))
    wrapper <- file.path(output, "NumberTheory.R")
    stopifnot(file.copy(file.path(repository, "R/NumberTheory.R"), wrapper,
                       overwrite = TRUE))
    write.csv(data.frame(file = c(copies, wrapper),
                         md5 = unname(tools::md5sum(c(copies, wrapper)))),
              file.path(output, "source-hashes.csv"), row.names = FALSE)

    previous <- Sys.getenv(c("CPATH", "PKG_LIBS"), unset = NA_character_)
    on.exit(for (name in names(previous)) {
        if (is.na(previous[[name]])) Sys.unsetenv(name)
        else do.call(Sys.setenv, setNames(list(previous[[name]]), name))
    }, add = TRUE)
    includes <- file.path(repository, "inst/include")
    if (!is.na(previous[["CPATH"]]) && nzchar(previous[["CPATH"]])) {
        includes <- c(includes, previous[["CPATH"]])
    }
    Sys.setenv(CPATH = paste(includes, collapse = .Platform$path.sep))
    libraries <- if (is.na(previous[["PKG_LIBS"]])) "" else previous[["PKG_LIBS"]]
    Sys.setenv(PKG_LIBS = paste(libraries, "-lgmp"))

    stubs <- file.path(output, "dispatch-stubs.cpp")
    writeLines(c(
        '#include "cpp11.hpp"',
        '#include "NumbersUtils/PrimeCount.h"',
        '#include <cstdint>',
        '#include <stdexcept>',
        'namespace {',
        '    int selected = 0, requested = 0, maximum = 0, calls = 0;',
        '    std::int64_t input = 0, forcedResult = -1;',
        '    bool fail = false;',
        '    std::int64_t record(int method, std::int64_t n, int threads, int cap) {',
        '        selected = method; input = n; requested = threads; maximum = cap; ++calls;',
        '        if (fail) throw std::runtime_error("injected dispatch failure");',
        '        return forcedResult >= 0 ? forcedResult : 100 + method;',
        '    }',
        '}',
        'namespace PrimeCounting {',
        '    std::int64_t PiPrime(std::int64_t n) { return record(1, n, 0, 0); }',
        '    std::int64_t PiLegendre(std::int64_t n, int t, int m) { return record(2, n, t, m); }',
        '    std::int64_t PiMeissel(std::int64_t n, int t, int m) { return record(3, n, t, m); }',
        '    std::int64_t PiLehmer(std::int64_t n, int t, int m) { return record(4, n, t, m); }',
        '}',
        '[[cpp11::register]]',
        'void reset_dispatch_state(double result, bool fault) {',
        '    selected = requested = maximum = calls = 0; input = 0;',
        '    forcedResult = static_cast<std::int64_t>(result); fail = fault;',
        '}',
        '[[cpp11::register]]',
        'cpp11::doubles dispatch_trace() {',
        '    cpp11::writable::doubles result(5);',
        '    result[0] = selected; result[1] = static_cast<double>(input);',
        '    result[2] = requested; result[3] = maximum; result[4] = calls;',
        '    return result;',
        '}'
    ), stubs)
    compiled <- new.env(parent = globalenv())
    build <- tempfile("build-", tmpdir = output)
    cpp11::cpp_source(c(copies, stubs), env = compiled, clean = FALSE,
                      quiet = FALSE, cxx_std = "CXX11", dir = build)

    # Keep the wrapper byte-for-byte intact. Only this isolated environment's
    # .Call binding changes; generated cpp11 functions retain their normal
    # environment and therefore still cross the actual R/C++ boundary.
    sandbox <- new.env(parent = globalenv())
    sys.source(wrapper, envir = sandbox)
    sandbox$pkgEnv <- list(maxThreads = 3L)
    assign("_RcppAlgos_PrimeCountCpp", "dispatch-fixture", envir = sandbox)
    sandbox$.Call <- function(symbol, ...) {
        stopifnot(identical(symbol, "dispatch-fixture"))
        compiled$PrimeCountCpp(...)
    }

    checks <- 0L
    reset <- function(result = -1, fault = FALSE) {
        compiled$reset_dispatch_state(result, fault)
    }
    check <- function(expression, method, n = 0, threads = 0, maximum = 0,
                      expected = as.integer(100 + method)) {
        reset()
        result <- eval(expression, envir = sandbox)
        trace <- if (method == 0) rep(0, 5) else c(method, n, threads, maximum, 1)
        stopifnot(identical(result, expected),
                  identical(as.numeric(compiled$dispatch_trace()), trace))
        checks <<- checks + 1L
    }
    check_error <- function(expression, pattern = NULL) {
        reset()
        message <- tryCatch({ force(expression); NA_character_ }, error = conditionMessage)
        stopifnot(!is.na(message), all(compiled$dispatch_trace() == 0))
        if (!is.null(pattern)) stopifnot(grepl(pattern, message, fixed = TRUE))
        checks <<- checks + 1L
    }

    # Expected automatic dispatch: retain tiny lookup/sieve, then select from
    # the input and usable thread budget. Sentinels distinguish actual driver
    # selection even though the mathematical algorithms return equal counts.
    tiny <- c(0L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 4L)
    for (n in seq_along(tiny)) {
        check(call("primeCount", n), 0, expected = tiny[[n]])
    }
    for (n in c(10, 99999)) check(call("primeCount", n), 1, n)
    inputs <- c(100000, 100001, 1e10 - 1, 1e10, 1e10 + 1,
                1e12 - 1, 1e12, 1e12 + 1, 2^53 - 1)
    serial <- c(2, 2, 2, 3, 3, 3, 3, 3, 3)
    parallel <- c(2, 2, 2, 3, 3, 3, 3, 3, 3)
    for (i in seq_along(inputs)) {
        n <- inputs[[i]]
        check(call("primeCount", n), serial[[i]], n, 1, 3)
        check(call("primeCount", n, 1L, "auto"), serial[[i]], n, 1, 3)
        check(call("primeCount", n, 2L, "auto"), parallel[[i]], n, 2, 3)
    }
    sandbox$pkgEnv$maxThreads <- 8L
    inputs <- c(1e12 - 1, 1e12, 1e12 + 1, 2^53 - 1)
    sixThreads <- c(3, 4, 4, 4)
    for (i in seq_along(inputs)) {
        n <- inputs[[i]]
        check(call("primeCount", n, 5L, "auto"), 3, n, 5, 8)
        check(call("primeCount", n, 6L, "auto"), sixThreads[[i]], n, 6, 8)
        check(call("primeCount", n, 8L, "auto"), sixThreads[[i]], n, 8, 8)
    }
    sandbox$pkgEnv$maxThreads <- 3L

    # Explicit methods retain their native driver for every input, including
    # small inputs and both sides of the automatic method boundaries.
    ids <- c(legendre = 2, meissel = 3, lehmer = 4)
    for (method in names(ids)) {
        for (n in c(1, 9, 10, 99999)) {
            check(call("primeCount", n, NULL, method), ids[[method]], n, 1, 3)
            check(call("primeCount", n, "ignored", method), ids[[method]], n, 1, 3)
        }
        for (n in c(100000, 1e10 - 1, 1e10, 1e12 - 1, 1e12, 2^53 - 1)) {
            check(call("primeCount", n, NULL, method), ids[[method]], n, 1, 3)
            check(call("primeCount", n, 1000L, method), ids[[method]], n, 1000, 3)
        }
    }
    check(quote(primeCount(9, "ignored")), 0, expected = 4L)
    check(quote(primeCount(99999, "ignored")), 1, 99999)
    check(quote(primeCount(100000, 2L)), 2, 100000, 2, 3)
    check(quote(primeCount(100000, method = "lehmer", nThreads = 4L)),
          4, 100000, 4, 3)
    sandbox$pkgEnv$maxThreads <- 7L
    check(quote(primeCount(100000, 100L, "meissel")), 3, 100000, 100, 7)

    # Selection uses the usable budget; drivers still receive the original
    # request and cap separately. NULL remains a request for one thread.
    sandbox$pkgEnv$maxThreads <- 1L
    check(quote(primeCount(1e12, 8L)), 3, 1e12, 8, 1)
    check(quote(primeCount(2^53 - 1, .Machine$integer.max)),
          3, 2^53 - 1, .Machine$integer.max, 1)
    check(quote(primeCount(1e12, NULL)), 3, 1e12, 1, 1)
    sandbox$pkgEnv$maxThreads <- 2L
    check(quote(primeCount(1e12, 1L)), 3, 1e12, 1, 2)
    check(quote(primeCount(1e12, 2L)), 3, 1e12, 2, 2)
    check(quote(primeCount(1e12, .Machine$integer.max)),
          3, 1e12, .Machine$integer.max, 2)
    check(quote(primeCount(1e12 - 1, .Machine$integer.max)),
          3, 1e12 - 1, .Machine$integer.max, 2)
    check(quote(primeCount(1e12, NULL)), 3, 1e12, 1, 2)
    sandbox$pkgEnv$maxThreads <- 5L
    check(quote(primeCount(1e12, .Machine$integer.max)),
          3, 1e12, .Machine$integer.max, 5)
    check(quote(primeCount(2^53 - 1, 6L)), 3, 2^53 - 1, 6, 5)
    sandbox$pkgEnv$maxThreads <- 6L
    check(quote(primeCount(1e12, 5L)), 3, 1e12, 5, 6)
    check(quote(primeCount(1e12, 6L)), 4, 1e12, 6, 6)
    check(quote(primeCount(1e12, .Machine$integer.max)),
          4, 1e12, .Machine$integer.max, 6)
    check(quote(primeCount(1e12 - 1, 6L)), 3, 1e12 - 1, 6, 6)
    check(quote(primeCount(1e12, NULL)), 3, 1e12, 1, 6)
    check(quote(primeCount(1e12, 6L, NULL)), 4, 1e12, 6, 6)
    check(quote(primeCount(1e12, 6L, "a")), 4, 1e12, 6, 6)
    check(quote(primeCount(1e12, 6L, "leg")), 2, 1e12, 6, 6)
    check(quote(primeCount(1e12, 6L, "m")), 3, 1e12, 6, 6)
    check(quote(primeCount(1e12, 6L, "leh")), 4, 1e12, 6, 6)
    sandbox$pkgEnv$maxThreads <- 3L
    check(quote(primeCount(1e12, 2L, NULL)), 3, 1e12, 2, 3)
    check(quote(primeCount(1e12, 2L, "a")), 3, 1e12, 2, 3)
    check(quote(primeCount(1e12, 2L, "leg")), 2, 1e12, 2, 3)
    check(quote(primeCount(1e12, 2L, "m")), 3, 1e12, 2, 3)
    check(quote(primeCount(1e12, 2L, "leh")), 4, 1e12, 2, 3)

    for (method in c("auto", names(ids))) {
        for (n in list(0, -1, 1.5, NA_real_, NaN, Inf, 2^53, "10", c(1, 2))) {
            check_error(sandbox$primeCount(n, method = method))
        }
        for (n in c(100000, 1e12)) {
            for (threads in list(0, -1, 1.5, NA_real_, Inf, "2", c(1L, 2L))) {
                check_error(sandbox$primeCount(n, threads, method))
            }
        }
        for (value in c(.Machine$integer.max, .Machine$integer.max + 1)) {
            reset(value)
            result <- sandbox$primeCount(1e10, method = method)
            expected <- if (value == .Machine$integer.max) as.integer(value) else value
            stopifnot(identical(result, expected))
            checks <- checks + 1L
        }
    }
    for (method in list("unknown", "le", NA_character_, 1, character(0))) {
        check_error(sandbox$primeCount(100000, method = method))
    }
    check_error(compiled$PrimeCountCpp(100000, NULL, 3L, "unknown"),
                "Unknown prime-counting method")
    check_error(compiled$PrimeCountCpp(1e12, "invalid", 3L, "unknown"),
                "Unknown prime-counting method")

    # This checks cpp11 translation on the calling thread only. There are no
    # worker threads in the sentinel fixture and no worker-safety claim here.
    reset(fault = TRUE)
    message <- tryCatch({ sandbox$primeCount(100000, method = "lehmer"); NA_character_ },
                        error = conditionMessage)
    stopifnot(identical(message, "injected dispatch failure"),
              identical(as.numeric(compiled$dispatch_trace()), c(4, 100000, 1, 3, 1)))
    checks <- checks + 1L
    check(quote(primeCount(100000, method = "lehmer")), 4, 100000, 1, 3)
    sandbox$pkgEnv$maxThreads <- 6L
    for (i in 1:2) {
        n <- c(1e10, 1e12)[[i]]
        threads <- c(1L, 6L)[[i]]
        method <- c(3, 4)[[i]]
        reset(fault = TRUE)
        message <- tryCatch({ sandbox$primeCount(n, threads); NA_character_ },
                            error = conditionMessage)
        stopifnot(identical(message, "injected dispatch failure"),
                  identical(as.numeric(compiled$dispatch_trace()),
                            c(method, n, threads, 6, 1)))
        checks <- checks + 1L
        check(call("primeCount", n, threads), method, n, threads, 6)
    }
    sandbox$pkgEnv$maxThreads <- 3L
    cat(sprintf("PASS dispatch/type/validation checks=%d build_directory=%s\n", checks, build))
}
main()
