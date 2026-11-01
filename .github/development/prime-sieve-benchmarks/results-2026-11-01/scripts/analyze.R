# Run after run.R finishes:
# Rscript .github/development/prime-sieve-benchmarks/analyze.R RESULTS_DIRECTORY
summarize_controls <- function(output_dir) {
    input <- file.path(output_dir, "controls.csv")
    if (!file.exists(input)) return(invisible(NULL))
    if (!file.exists(file.path(output_dir, "controls-metadata.txt"))) {
        message("Control run is not marked complete; skipping its analysis")
        return(invisible(NULL))
    }
    runs <- read.csv(input, stringsAsFactors = FALSE)
    groups <- split(runs, interaction(runs[c("case", "requested", "mode")], drop = TRUE))
    summary <- do.call(rbind, lapply(groups, function(d) {
        q <- unname(quantile(d$per_call_ms, c(0.25, 0.5, 0.75)))
        data.frame(case = d$case[1L], requested = d$requested[1L], mode = d$mode[1L],
                   rounds = nrow(d), calls = sum(d$calls),
                   per_call_ms_q1 = q[1L], per_call_ms_median = q[2L], per_call_ms_q3 = q[3L])
    }))
    summary$ratio_to_ordinary <- NA_real_
    for (i in seq_len(nrow(summary))) {
        reference <- summary$case == summary$case[i] &
            summary$requested == summary$requested[i] & summary$mode == "ordinary"
        stopifnot(sum(reference) == 1L)
        summary$ratio_to_ordinary[i] <- summary$per_call_ms_median[i] /
            summary$per_call_ms_median[reference]
    }
    phases_path <- file.path(output_dir, "control-phases.csv")
    if (file.exists(phases_path)) {
        phases <- read.csv(phases_path, stringsAsFactors = FALSE)
        behavior <- unique(phases[c("case", "requested", "actual", "kernel")])
        summary <- merge(summary, behavior, by = c("case", "requested"), all.x = TRUE)
    }
    write.csv(summary, file.path(output_dir, "controls-summary.csv"), row.names = FALSE)
    inspected <- summary[summary$mode == "inspected", ]
    scores <- do.call(rbind, lapply(split(inspected, inspected$requested), function(d) {
        data.frame(requested = d$requested[1L], cases = nrow(d),
                   inspected_ordinary_ratio_geomean = exp(mean(log(d$ratio_to_ordinary))),
                   best_ratio = min(d$ratio_to_ordinary), worst_ratio = max(d$ratio_to_ordinary))
    }))
    write.csv(scores, file.path(output_dir, "controls-scores.csv"), row.names = FALSE)
    cat("\nFull R-call control comparisons:\n")
    print(inspected, row.names = FALSE)
    invisible(list(summary = summary, scores = scores))
}

summarize_native_confirmation <- function(output_dir) {
    input <- file.path(output_dir, "native-confirmation.csv")
    if (!file.exists(input)) return(invisible(NULL))
    metadata <- file.path(output_dir, "native-confirmation-metadata.json")
    if (!file.exists(metadata) ||
        !any(grepl('"status": "completed"', readLines(metadata), fixed = TRUE))) {
        message("Native confirmation is not marked complete; skipping its analysis")
        return(invisible(NULL))
    }
    runs <- read.csv(input, stringsAsFactors = FALSE)
    stopifnot(all(runs$matches_reference))
    groups <- split(runs, interaction(runs[c("base", "width", "requested", "config")], drop = TRUE))
    summary <- do.call(rbind, lapply(groups, function(d) {
        row <- d[1L, c("base", "lower", "upper", "width", "requested", "config",
                        "budget_bytes", "segment_span", "prime_count")]
        q <- unname(quantile(d$native_ms, c(0.25, 0.5, 0.75)))
        row$samples <- nrow(d)
        row$actual_min <- min(d$actual)
        row$actual_max <- max(d$actual)
        row$native_ms_q1 <- q[1L]
        row$native_ms_median <- q[2L]
        row$native_ms_q3 <- q[3L]
        row
    }))
    summary$ratio_to_default <- NA_real_
    summary$ratio_to_serial_default <- NA_real_
    for (i in seq_len(nrow(summary))) {
        group <- summary$base == summary$base[i] & summary$width == summary$width[i]
        reference <- group & summary$requested == summary$requested[i] & summary$config == "default"
        serial <- group & summary$requested == 1L & summary$config == "default"
        stopifnot(sum(reference) == 1L, sum(serial) == 1L)
        summary$ratio_to_default[i] <- summary$native_ms_median[i] / summary$native_ms_median[reference]
        summary$ratio_to_serial_default[i] <- summary$native_ms_median[i] / summary$native_ms_median[serial]
    }
    summary$native_iqr_percent <- 100 *
        (summary$native_ms_q3 - summary$native_ms_q1) / summary$native_ms_median
    write.csv(summary, file.path(output_dir, "native-summary.csv"), row.names = FALSE)
    score_groups <- split(summary, interaction(summary[c("requested", "config")], drop = TRUE))
    scores <- do.call(rbind, lapply(score_groups, function(d) {
        data.frame(requested = d$requested[1L], config = d$config[1L], cases = nrow(d),
            default_ratio_geomean = exp(mean(log(d$ratio_to_default))),
            serial_default_ratio_geomean = exp(mean(log(d$ratio_to_serial_default))),
            best_ratio = min(d$ratio_to_default), worst_ratio = max(d$ratio_to_default))
    }))
    scores <- scores[order(scores$requested, scores$default_ratio_geomean), ]
    write.csv(scores, file.path(output_dir, "native-scores.csv"), row.names = FALSE)
    cat("\nMinimally instrumented native confirmation:\n")
    print(scores, row.names = FALSE)
    invisible(list(summary = summary, scores = scores))
}

draw_serial_sizing <- function(summary) {
    par(mfrow = c(2, 3), mar = c(4.7, 4.2, 3, 1), oma = c(1.5, 0.5, 3, 0.5))
    serial <- summary[summary$stage == "serial", ]
    for (case_name in unique(serial$case)) {
        d <- serial[serial$case == case_name, ]
        d <- d[order(d$sieve_bytes_median), ]
        x <- d$sieve_bytes_median / 1024
        near <- format(d$lower[1L], scientific = TRUE, digits = 1)
        width <- format(d$width[1L], scientific = TRUE, digits = 1)
        plot(x, d$native_ms_median, log = "x", type = "b", pch = 16,
             ylim = range(c(d$native_ms_q1, d$native_ms_q3)),
             xlab = "Reported sieve payload (KiB)", ylab = "Native elapsed (ms)",
             main = paste("Near", near, "/ width", width), col = "#29648A")
        arrows(x, d$native_ms_q1, x, d$native_ms_q3, angle = 90,
               code = 3, length = 0.03, col = "#29648A")
        baseline <- d$config == "default"
        points(x[baseline], d$native_ms_median[baseline], pch = 21,
               bg = "#B64926", col = "#B64926", cex = 1.4)
    }
    # One figure-wide key keeps all data and error bars unobstructed.
    par(fig = c(0, 1, 0, 1), mar = rep(0, 4), oma = rep(0, 4), new = TRUE)
    plot.new()
    legend("top", inset = 0.01, c("Median with interquartile range", "Current default"),
           col = c("#29648A", "#B64926"), pch = c(16, 16), bty = "n", horiz = TRUE,
           cex = 0.9)
}

analyze_inspector <- function(output_dir) {
    metadata <- readLines(file.path(output_dir, "metadata.txt"))
    if (!any(grepl("^Completed:", metadata))) {
        stop("Benchmark run is not marked complete; wait before analyzing")
    }
    runs <- read.csv(file.path(output_dir, "runs.csv"), stringsAsFactors = FALSE)
    keys <- c("stage", "case", "config", "requested")
    measures <- c("wall_ms", "native_ms", "base_ms", "policy_ms", "sieve_ms",
                  "result_ms", "setup_worker_ms", "offsets_worker_ms",
                  "buckets_worker_ms", "dense_mark_worker_ms", "bucket_mark_worker_ms",
                  "extract_worker_ms", "reset_worker_ms", "sieve_bytes", "offset_bytes",
                  "bucket_entries", "bucket_capacity_bytes", "output_capacity_bytes",
                  "shared_bytes", "segments")
    groups <- split(runs, interaction(runs[keys], drop = TRUE))
    summary <- do.call(rbind, lapply(groups, function(d) {
        row <- d[1L, c(keys, "lower", "upper", "width", "budget_bytes",
                        "segment_span", "prime_count", "kernel")]
        row$samples <- nrow(d)
        row$actual_min <- min(d$actual)
        row$actual_max <- max(d$actual)
        for (measure in measures) {
            quartiles <- unname(quantile(d[[measure]], c(0.25, 0.5, 0.75)))
            row[[paste0(measure, "_q1")]] <- quartiles[1L]
            row[[paste0(measure, "_median")]] <- quartiles[2L]
            row[[paste0(measure, "_q3")]] <- quartiles[3L]
        }
        row
    }))
    rownames(summary) <- NULL

    # Fixed-size comparisons use the default with the same requested workers.
    # The policy stage has only defaults, named default_t1/default_t2/etc.;
    # there its reference is the one-requested-worker default.
    summary$native_ratio_to_default <- NA_real_
    summary$wall_ratio_to_default <- NA_real_
    summary$native_ratio_to_serial_default <- NA_real_
    for (i in seq_len(nrow(summary))) {
        group <- summary$stage == summary$stage[i] & summary$case == summary$case[i]
        same_default <- group & summary$config == "default" &
            summary$requested == summary$requested[i]
        serial_default <- group & summary$config == "default" & summary$requested == 1L
        if (summary$stage[i] == "policy") {
            same_default <- serial_default <- group & summary$config == "default_t1"
        }
        if (sum(same_default) != 1L || sum(serial_default) != 1L) {
            stop("Missing or duplicated comparison baseline")
        }
        summary$native_ratio_to_default[i] <- summary$native_ms_median[i] /
            summary$native_ms_median[same_default]
        summary$wall_ratio_to_default[i] <- summary$wall_ms_median[i] /
            summary$wall_ms_median[same_default]
        summary$native_ratio_to_serial_default[i] <- summary$native_ms_median[i] /
            summary$native_ms_median[serial_default]
    }
    summary$native_iqr_percent <- 100 *
        (summary$native_ms_q3 - summary$native_ms_q1) / summary$native_ms_median
    summary <- summary[order(summary$stage, summary$case, summary$requested,
                             summary$native_ms_median), ]
    write.csv(summary, file.path(output_dir, "summary.csv"), row.names = FALSE)

    score_groups <- split(summary, interaction(summary[c("stage", "config", "requested")],
                                                 drop = TRUE))
    scores <- do.call(rbind, lapply(score_groups, function(d) {
        data.frame(stage = d$stage[1L], config = d$config[1L],
            requested = d$requested[1L], cases = nrow(d),
            native_ratio_geomean = exp(mean(log(d$native_ratio_to_default))),
            wall_ratio_geomean = exp(mean(log(d$wall_ratio_to_default))),
            serial_default_ratio_geomean = exp(mean(log(d$native_ratio_to_serial_default))),
            worst_native_ratio = max(d$native_ratio_to_default),
            best_native_ratio = min(d$native_ratio_to_default))
    }))
    scores <- scores[order(scores$stage, scores$requested, scores$native_ratio_geomean), ]
    write.csv(scores, file.path(output_dir, "balanced-scores.csv"), row.names = FALSE)

    best_groups <- split(summary, interaction(summary[c("stage", "case", "requested")],
                                                drop = TRUE))
    best <- do.call(rbind, lapply(best_groups, function(d) {
        d[which.min(d$native_ms_median), c(keys, "actual_min", "actual_max", "kernel",
            "budget_bytes", "native_ms_median", "native_ms_q1", "native_ms_q3",
            "wall_ms_median", "native_ratio_to_default", "native_ratio_to_serial_default")]
    }))
    write.csv(best, file.path(output_dir, "best-per-case.csv"), row.names = FALSE)

    workers <- unique(summary[c(keys, "actual_min", "actual_max", "kernel",
                                "segment_span", "segments_median", "sieve_bytes_median")])
    write.csv(workers, file.path(output_dir, "worker-behavior.csv"), row.names = FALSE)

    phase_names <- c("setup", "offsets", "buckets", "dense_mark", "bucket_mark", "extract", "reset")
    phases <- do.call(rbind, lapply(phase_names, function(phase) {
        d <- summary[c(keys, "actual_min", "actual_max", "kernel")]
        d$phase <- phase
        d$worker_ms_q1 <- summary[[paste0(phase, "_worker_ms_q1")]]
        d$worker_ms_median <- summary[[paste0(phase, "_worker_ms_median")]]
        d$worker_ms_q3 <- summary[[paste0(phase, "_worker_ms_q3")]]
        # Only serial phase times may be compared directly with total native
        # elapsed time. Concurrent worker times sum overlapping intervals.
        d$percent_native_serial_only <- ifelse(summary$actual_max == 1L,
            100 * d$worker_ms_median / summary$native_ms_median, NA_real_)
        d
    }))
    write.csv(phases, file.path(output_dir, "phase-summary.csv"), row.names = FALSE)

    pdf(file.path(output_dir, "serial-sizing.pdf"), width = 12, height = 7)
    draw_serial_sizing(summary)
    dev.off()
    # Native Quartz avoids an optional XQuartz dependency in macOS R builds.
    png_type <- if (Sys.info()[["sysname"]] == "Darwin" && capabilities("aqua")) {
        "quartz"
    } else if (capabilities("cairo")) "cairo" else getOption("bitmapType")
    png(file.path(output_dir, "serial-sizing.png"), width = 1600, height = 1000, res = 140,
        type = png_type)
    draw_serial_sizing(summary)
    dev.off()

    print(scores, row.names = FALSE)
    cat("\nRatios below 1 are faster; geometric means weight cases equally.\n",
        "Worker phase times are summed elapsed times, not CPU time or cache counters.\n",
        "No timing curve establishes cache residency or a cache-size boundary.\n", sep = "")
    controls <- summarize_controls(output_dir)
    native <- summarize_native_confirmation(output_dir)
    invisible(list(summary = summary, scores = scores, best = best, phases = phases,
                   controls = controls, native = native))
}

if (sys.nframe() == 0L) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) != 1L) stop("Supply the completed benchmark results directory")
    analyze_inspector(args[1L])
}
