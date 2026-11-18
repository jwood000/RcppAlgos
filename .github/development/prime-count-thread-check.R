# Run from the repository root:
# Rscript --vanilla .github/development/prime-count-thread-check.R [temporary_directory]
# Builds only in a temporary directory; production and registration files stay untouched.
main <- function() {
    repository <- normalizePath(".", winslash = "/")
    checker <- file.path(repository, ".github/development/prime-count-thread-check.cpp")
    preparer <- file.path(repository, ".github/development/prepare-prime-count-thread-check.py")
    stopifnot(file.exists(checker), file.exists(preparer))
    arguments <- commandArgs(TRUE)
    if (length(arguments) > 1L) stop("Expected at most one temporary output directory")
    output <- if (length(arguments)) arguments[[1L]] else tempfile("rcppalgos-thread-r-")
    dir.create(output, recursive = TRUE, showWarnings = FALSE)
    output <- normalizePath(output, winslash = "/")
    status <- system2("python3", c(shQuote(preparer), shQuote(output)))
    if (status != 0L) stop("Preparing the instrumented source failed")

    # cpp_source supplies cpp11's include directory; CPATH adds only local sources.
    previous_cpath <- Sys.getenv("CPATH", unset = NA_character_)
    on.exit(if (is.na(previous_cpath)) Sys.unsetenv("CPATH") else
        Sys.setenv(CPATH = previous_cpath), add = TRUE)
    includes <- c(file.path(repository, "inst/include"), output)
    if (!is.na(previous_cpath) && nzchar(previous_cpath)) includes <- c(includes, previous_cpath)
    Sys.setenv(CPATH = paste(includes, collapse = .Platform$path.sep))

    source <- c(
        '#include "cpp11.hpp"',
        '#define main standalone_thread_check_main',
        sprintf('#include "%s"', checker),
        '#undef main',
        '[[cpp11::register]]',
        'double invoke_thread_fault(int region, int mode) {',
        '    using namespace ThreadFaults;',
        '    const Region selected = static_cast<Region>(region);',
        '    WorkerFault fault = WorkerFault::None;',
        '    if (mode == 1 || mode == 3) fault = WorkerFault::Standard;',
        '    if (mode == 4) fault = WorkerFault::NonStandard;',
        '    if (mode == 5) fault = WorkerFault::BadAlloc;',
        '    if (mode == 0) Observe();',
        '    else Configure(Scenario(selected, fault, (mode == 2 || mode == 3) ? 1 : -1,',
        '                            false, mode == 3));',
        '    const std::int64_t n = selected == Region::Small ? 10000000000LL :',
        '        (selected == Region::Final ? 999999999999LL : 10000000001LL);',
        '    return static_cast<double>(PrimeCounting::PiLegendre(n, 4, 4));',
        '}',
        '[[cpp11::register]]',
        'std::string expected_thread_fault(int mode) {',
        '    if (mode == 2 || mode == 3) return ThreadFaults::LaunchError().what();',
        '    if (mode == 4) return "C++ error (unknown cause)";',
        '    if (mode == 5) return std::bad_alloc().what();',
        '    return ThreadFaults::WorkerError(0).what();',
        '}',
        '[[cpp11::register]]',
        'cpp11::integers thread_fault_counters() {',
        '    const auto snapshot = ThreadFaults::Inspect();',
        '    if (!snapshot.protocolError.empty()) throw std::runtime_error(snapshot.protocolError);',
        '    cpp11::writable::integers result(13);',
        '    result[0] = snapshot.target;',
        '    result[1] = static_cast<int>(snapshot.batches.size());',
        '    int started = 0, bodies = 0, exited = 0, joined = 0;',
        '    for (const auto& batch : snapshot.batches) {',
        '        started += batch.started; bodies += batch.bodies;',
        '        exited += batch.exited; joined += batch.joined;',
        '        if (batch.started != batch.bodies || batch.started != batch.exited ||',
        '            batch.started != batch.joined) throw std::runtime_error("unjoined batch");',
        '    }',
        '    result[2] = started; result[3] = bodies; result[4] = exited; result[5] = joined;',
        '    if (snapshot.target >= 0) {',
        '        const auto& batch = snapshot.batches[snapshot.target];',
        '        result[6] = batch.started; result[7] = batch.attempts;',
        '        result[8] = batch.exited; result[9] = batch.joined;',
        '    } else for (int i = 6; i <= 9; ++i) result[i] = 0;',
        '    result[10] = snapshot.launchFailed;',
        '    result[11] = static_cast<int>(snapshot.throwOrder.size());',
        '    result[12] = snapshot.throwOrder.empty() ? -1 : snapshot.throwOrder.front();',
        '    return result;',
        '}'
    )
    wrapper <- file.path(output, "prime-count-thread-r-bridge.cpp")
    writeLines(source, wrapper)
    functions <- new.env(parent = globalenv())
    # Always build afresh: the included checker and generated core are external files.
    build_directory <- tempfile("r-bridge-build-", tmpdir = output)
    cpp11::cpp_source(c(wrapper, file.path(repository, "src/Phi.cpp"),
                       file.path(repository, "src/Eratosthenes.cpp")),
        env = functions, clean = FALSE, quiet = FALSE, cxx_std = "CXX11",
        dir = build_directory)

    regions <- c("first", "middle", "final", "small")
    modes <- c("worker", "launch", "combined", "nonstandard", "allocation")
    counts <- c(455052511, 455052511, 37607912018, 455052511)
    sizes <- c(4L, 4L, 2L, 4L)
    checks <- 0L
    for (region in 0:3) for (mode in 1:5) {
        expected_message <- functions$expected_thread_fault(mode)
        message <- tryCatch({ functions$invoke_thread_fault(region, mode); NA_character_ },
                            error = conditionMessage)
        stopifnot(identical(message, expected_message))
        state <- functions$thread_fault_counters()
        launched <- if (mode %in% c(2L, 3L)) 1L else sizes[region + 1L]
        stopifnot(state[1L] >= 0L, state[2L] == state[1L] + 1L,
                  all(state[3:6] == state[3L]), state[7L] == launched,
                  state[8L] == launched + as.integer(mode %in% c(2L, 3L)),
                  state[9L] == launched, state[10L] == launched,
                  state[11L] == as.integer(mode %in% c(2L, 3L)),
                  state[12L] == as.integer(mode != 2L),
                  state[13L] == if (mode == 2L) -1L else 0L)
        recovered <- functions$invoke_thread_fault(region, 0L)
        stopifnot(identical(recovered, counts[region + 1L]))
        recovered_state <- functions$thread_fault_counters()
        stopifnot(recovered_state[1L] == -1L,
                  all(recovered_state[3:6] == recovered_state[3L]))
        checks <- checks + 1L
        cat(sprintf("PASS region=%s fault=%s joined=%d message=%s recovery=%.0f\n",
            regions[region + 1L], modes[mode], state[6L], message, recovered))
    }
    cat(sprintf("PASS R_errors=%d recovered_counts=%d regions=4 build_directory=%s\n",
                checks, checks, build_directory))
}
main()
