# Run from the repository root in a fresh R process:
# Rscript .github/development/prime-sieve-benchmarks/run.R LIBRARY OUTPUT [REPEATS]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Supply the installed development library and output directory")
library(RcppAlgos, lib.loc = normalizePath(args[1]))
out <- args[2]
repeats <- if (length(args) > 2L) as.integer(args[3]) else 7L
stopifnot(is.finite(repeats), repeats >= 3L)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(out, "runs.csv"))) stop("Choose a fresh output directory")
set.seed(20261101)
inspect <- RcppAlgos:::primeSieveInspect
policy <- get("pkgEnv", asNamespace("RcppAlgos"))
bigCut <- (policy$nCores * 6 * 1024^2)^2
limit <- 2^53 - 1
if (bigCut >= 8.5e15) stop("No useful serial Big domain on this host; design an explicitly fixed-policy experiment")
bases <- if (bigCut < 4e15) c(4e15, 8e15) else {
    floor(bigCut + c(0.2, 0.75) * (limit - bigCut))
}
threads <- c(1L, 2L, 4L)
threads <- threads[threads <= policy$maxThreads]
budgets <- c(default = NA_real_, wheel = 3754, `16KiB` = 16 * 1024,
             `64KiB` = 64 * 1024, `256KiB` = 256 * 1024,
             `1MiB` = 1024^2, `4MiB` = 4 * 1024^2)
budget_arg <- function(x) if (is.na(x)) NULL else x

metadata <- c(
    paste("Started:", format(Sys.time(), tz = "America/New_York", usetz = TRUE)),
    paste("Revision:", system2("git", c("rev-parse", "HEAD"), stdout = TRUE)),
    paste("Policy physical cores:", policy$nCores),
    paste("Policy max threads:", policy$maxThreads),
    paste("Serial Big cutoff:", format(bigCut, scientific = FALSE)),
    paste("Random seed:", 20261101), paste("Screening repeats:", repeats),
    "No affinity, frequency locking, cache flushing, or hardware counters.",
    "Each timed call freshly generates its base primes; repeated process is warm.",
    capture.output(sessionInfo()), "Git status:",
    system2("git", c("status", "--short"), stdout = TRUE)
)
if (Sys.info()[["sysname"]] == "Darwin") {
    keys <- c("machdep.cpu.brand_string", "hw.model", "hw.physicalcpu", "hw.logicalcpu",
              "hw.cachelinesize", "hw.l1dcachesize", "hw.l2cachesize",
              "hw.perflevel0.name", "hw.perflevel0.physicalcpu",
              "hw.perflevel0.l1dcachesize", "hw.perflevel0.l2cachesize",
              "hw.perflevel0.cpusperl2", "hw.perflevel1.name",
              "hw.perflevel1.physicalcpu", "hw.perflevel1.l1dcachesize",
              "hw.perflevel1.l2cachesize", "hw.perflevel1.cpusperl2")
    metadata <- c(metadata, "Hardware reports (missing keys are explicit):",
        unlist(lapply(keys, function(key) {
            value <- suppressWarnings(system2("sysctl", key, stdout = TRUE, stderr = TRUE))
            if (!is.null(attr(value, "status"))) paste(key, "unavailable") else value
        })))
}
writeLines(metadata, file.path(out, "metadata.txt"))
sources <- c("src/Eratosthenes.cpp", "src/PrimeSieve.cpp", "src/CacheInspector.cpp",
             "src/cpp11.cpp", "inst/include/NumbersUtils/Eratosthenes.h",
             "inst/include/NumbersUtils/CacheInspector.h", "R/CacheInspector.R",
             list.files(".github/development/prime-sieve-benchmarks", "\\.R$", full.names = TRUE),
             file.path(find.package("RcppAlgos"), "libs", "RcppAlgos.so"))
write.csv(data.frame(path = sources, md5 = unname(tools::md5sum(sources))),
          file.path(out, "source-hashes.csv"), row.names = FALSE)

source(".github/development/prime-sieve-benchmarks/validate.R")
validate_inspector(inspect, out, bases[1])

rows <- list()
counter <- 0L
save_rows <- function() write.csv(do.call(rbind, rows), file.path(out, "runs.csv"), row.names = FALSE)
record <- function(result, stage, case, config, lower, upper, requested, bytes, round, wall) {
    s <- result$summary; k <- result$kernels; p <- result$phases
    row <- data.frame(stage, case, config, lower, upper, width = upper - lower,
        requested, actual = s$effective_threads, budget_bytes = bytes, round,
        kernel = paste(unique(k$kernel), collapse = "+"),
        segment_span = if (nrow(k)) max(k$segment_span) else 0,
        segments = sum(k$segments), prime_count = length(result$primes),
        wall_ms = wall * 1000,
        native_ms = 1000 * sum(unlist(s[c("base_seconds", "policy_seconds", "sieve_seconds", "result_seconds")])),
        base_ms = 1000 * s$base_seconds, policy_ms = 1000 * s$policy_seconds,
        sieve_ms = 1000 * s$sieve_seconds, result_ms = 1000 * s$result_seconds,
        sieve_bytes = sum(k$sieve_bytes), offset_bytes = sum(k$offset_bytes),
        bucket_entries = sum(k$bucket_entries), bucket_capacity_bytes = sum(k$bucket_capacity_bytes),
        output_capacity_bytes = sum(k$output_capacity_bytes),
        shared_bytes = sum(unlist(result$shared)))
    for (phase in c("setup", "offsets", "buckets", "dense_mark", "bucket_mark", "extract", "reset")) {
        row[[paste0(phase, "_worker_ms")]] <- 1000 * sum(p$seconds[p$phase == phase])
    }
    counter <<- counter + 1L
    rows[[counter]] <<- row
}

run_case <- function(stage, case, lower, upper, configurations, rounds) {
    cat(stage, case, "configurations", nrow(configurations), "rounds", rounds, "\n")
    flush.console()
    expected <- primeSieve(lower, upper, nThreads = 1L)
    # Warm every configuration and verify complete results before measurement.
    for (i in sample(seq_len(nrow(configurations)))) {
        cfg <- configurations[i, ]
        result <- inspect(lower, upper, nThreads = cfg$threads,
                          segmentBytes = budget_arg(cfg$bytes))
        stopifnot(identical(result$primes, expected))
    }
    rm(result); gc()
    for (round in seq_len(rounds)) {
        for (i in sample(seq_len(nrow(configurations)))) {
            cfg <- configurations[i, ]
            start <- proc.time()[["elapsed"]]
            result <- inspect(lower, upper, nThreads = cfg$threads,
                              segmentBytes = budget_arg(cfg$bytes))
            wall <- proc.time()[["elapsed"]] - start
            record(result, stage, case, cfg$config, lower, upper,
                   cfg$threads, cfg$bytes, round, wall)
            # Kept outside both timers. Do not retain prime vectors in the CSV.
            stopifnot(identical(result$primes, expected))
            rm(result)
        }
    }
    save_rows()
}

serial <- data.frame(config = names(budgets), bytes = unname(budgets), threads = 1L)
for (base in bases) for (width in c(1e3, 1e6, 1e8)) {
    lower <- base + 17
    case <- paste(format(base, scientific = TRUE), format(width, scientific = TRUE), sep = "_")
    run_case("serial", case, lower, lower + width, serial, repeats)
}

# Select two fixed budgets by equal-weight geometric mean of per-case ratios.
# Screening data only nominates finalists; later stages collect fresh samples.
d <- do.call(rbind, rows)
med <- aggregate(native_ms ~ case + config, d[d$stage == "serial", ], median)
reference <- med[med$config == "default", c("case", "native_ms")]
names(reference)[2] <- "default_ms"
med <- merge(med, reference, by = "case")
med$ratio <- med$native_ms / med$default_ms
scores <- aggregate(log(ratio) ~ config, med[med$config != "default", ], mean)
best <- scores$config[order(scores[[2]])][1:2]
selected <- c("default", best)
writeLines(selected, file.path(out, "selected-configurations.txt"))
parallel <- expand.grid(config = selected, threads = threads, stringsAsFactors = FALSE)
parallel$bytes <- unname(budgets[parallel$config])
for (base in bases) for (width in c(1e6, 1e8)) {
    lower <- base + 17
    case <- paste(format(base, scientific = TRUE), format(width, scientific = TRUE), sep = "_")
    run_case("parallel", case, lower, lower + width, parallel, repeats)
}

# Fresh serial confirmation on shifted absolute intervals, common to all sizes.
finalists <- serial[serial$config %in% selected, ]
for (base in bases) for (width in c(1e3, 1e6, 1e8)) {
    lower <- base + 1234567
    case <- paste(format(base, scientific = TRUE), format(width, scientific = TRUE), sep = "_")
    run_case("confirmation", case, lower, lower + width, finalists, 15L)
}

# Exercise existing policy where requested threads may change kernels or fall back.
controls <- data.frame(config = paste0("default_t", threads), bytes = NA_real_, threads = threads)
for (width in c(1e3, 1e6, 1e8)) {
    run_case("policy", paste0("1e15_", format(width, scientific = TRUE)),
             1e15 + 17, 1e15 + 17 + width, controls, repeats)
}
cat("Finished", counter, "timed calls; all complete outputs matched.\n")
writeLines(c(metadata, paste("Completed:", format(Sys.time(), tz = "America/New_York", usetz = TRUE)),
             paste("Timed calls:", counter)), file.path(out, "metadata.txt"))
