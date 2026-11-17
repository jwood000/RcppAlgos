# Run from the repository root:
# Rscript --vanilla .github/development/prime-sieve-thread-check.R [temporary_directory]
# Only temporary files are generated; production and registration stay untouched.
main <- function() {
    repository <- normalizePath(".", winslash = "/")
    checker <- file.path(repository, ".github/development/prime-sieve-thread-check.cpp")
    preparer <- file.path(repository, ".github/development/prepare-prime-sieve-thread-check.py")
    arguments <- commandArgs(TRUE)
    if (length(arguments) > 1L) stop("Expected at most one temporary output directory")
    output <- if (length(arguments)) arguments[[1L]] else tempfile("rcppalgos-sieve-thread-r-")
    dir.create(output, recursive = TRUE, showWarnings = FALSE)
    output <- normalizePath(output, winslash = "/")
    status <- system2("python3", c(shQuote(preparer), shQuote(output)))
    if (status != 0L) stop("Preparing the instrumented source failed")

    previous_cpath <- Sys.getenv("CPATH", unset = NA_character_)
    on.exit(if (is.na(previous_cpath)) Sys.unsetenv("CPATH") else
        Sys.setenv(CPATH = previous_cpath), add = TRUE)
    includes <- c(file.path(repository, "inst/include"), output)
    if (!is.na(previous_cpath) && nzchar(previous_cpath)) includes <- c(includes, previous_cpath)
    Sys.setenv(CPATH = paste(includes, collapse = .Platform$path.sep))

    source <- c(
        '#include "cpp11.hpp"',
        '#define main standalone_sieve_thread_check_main',
        sprintf('#include "%s"', checker),
        '#undef main',
        '[[cpp11::register]]',
        'double invoke_sieve_thread_fault(int fixture, int mode) {',
        '    using namespace SieveThreadFaults;',
        '    const int maximum = fixture == 1 ? 450000 : 1000000;',
        '    const int threads = fixture == 2 ? 1 : (fixture == 1 ? 8 : 4);',
        '    const int workers = fixture == 2 ? 1 : (fixture == 1 ? 3 : 4);',
        '    Fault fault = Fault::None;',
        '    if (mode == 1 || mode == 3) fault = Fault::Standard;',
        '    if (mode == 4) fault = Fault::NonStandard;',
        '    if (mode == 5) fault = Fault::BadAlloc;',
        '    if (mode == 8) fault = Fault::Ordered;',
        '    int launch = -1;',
        '    if (mode == 2 || mode == 3) launch = 1;',
        '    if (mode == 6) launch = 0;',
        '    if (mode == 7) launch = workers - 1;',
        '    Configure(Scenario(fault, launch, mode == 4 ? workers - 1 : 0, false, mode == 3));',
        '    return static_cast<double>(SieveThreadChecks::Invoke(maximum, threads).size());',
        '}',
        '[[cpp11::register]]',
        'std::string expected_sieve_thread_fault(int mode) {',
        '    if (mode == 2 || mode == 3 || mode == 6 || mode == 7) return SieveThreadFaults::LaunchError().what();',
        '    if (mode == 4) return "C++ error (unknown cause)";',
        '    if (mode == 5) return std::bad_alloc().what();',
        '    return SieveThreadFaults::WorkerError(0).what();',
        '}',
        '[[cpp11::register]]',
        'cpp11::integers sieve_thread_fault_counters() {',
        '    const auto state = SieveThreadFaults::Inspect();',
        '    SieveThreadChecks::CheckJoined(state);',
        '    cpp11::writable::integers result(10);',
        '    result[0] = state.attempts;',
        '    result[1] = state.started;',
        '    result[2] = state.bodies;',
        '    result[3] = state.exited;',
        '    result[4] = state.joined;',
        '    result[5] = state.launchFailed;',
        '    result[6] = state.throwOrder.size();',
        '    result[7] = state.throwOrder.empty() ? -1 : state.throwOrder.front();',
        '    result[8] = state.throwOrder.empty() ? -1 : state.throwOrder.back();',
        '    result[9] = state.joinAttempts;',
        '    return result;',
        '}'
    )
    wrapper <- file.path(output, "prime-sieve-thread-r-bridge.cpp")
    writeLines(source, wrapper)
    functions <- new.env(parent = globalenv())
    build_directory <- tempfile("r-bridge-build-", tmpdir = output)
    cpp11::cpp_source(wrapper, env = functions, clean = FALSE, quiet = FALSE,
                      cxx_std = "CXX11", dir = build_directory)

    checks <- 0L
    for (fixture in 0:2) {
        modes <- if (fixture == 2L) c(1L, 4L, 5L) else 1:8
        workers <- c(4L, 3L, 0L)[fixture + 1L]
        count <- c(78498, 37706, 78498)[fixture + 1L]
        for (mode in modes) {
            expected_message <- functions$expected_sieve_thread_fault(mode)
            message <- tryCatch({ functions$invoke_sieve_thread_fault(fixture, mode); NA_character_ },
                                error = conditionMessage)
            stopifnot(identical(message, expected_message))
            state <- functions$sieve_thread_fault_counters()
            launch_failure <- mode %in% c(2L, 3L, 6L, 7L)
            launched <- if (mode %in% c(2L, 3L)) 1L else
                if (mode == 6L) 0L else if (mode == 7L) workers - 1L else workers
            throws <- if (mode %in% c(2L, 6L, 7L)) 0L else if (mode == 8L) 2L else 1L
            first <- if (throws == 0L) -1L else if (mode == 8L) 1L else
                if (mode == 4L && fixture != 2L) workers - 1L else 0L
            stopifnot(state[1L] == launched + as.integer(launch_failure),
                      all(state[2:5] == launched), state[6L] == as.integer(launch_failure),
                      state[7L] == throws, state[8L] == first,
                      state[10L] == launched)
            recovered <- functions$invoke_sieve_thread_fault(fixture, 0L)
            stopifnot(identical(recovered, count))
            recovered_state <- functions$sieve_thread_fault_counters()
            stopifnot(all(recovered_state[2:5] == workers))
            checks <- checks + 1L
            cat(sprintf("PASS fixture=%d mode=%d joined=%d message=%s recovery=%.0f\n",
                        fixture, mode, state[5L], message, recovered))
        }
    }
    cat(sprintf("PASS R_errors=%d recovered_sieves=%d build_directory=%s\n",
                checks, checks, build_directory))
}
main()
