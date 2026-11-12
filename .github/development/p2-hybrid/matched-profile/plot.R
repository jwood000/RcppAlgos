# Base graphics only; representative measured traces, not a schematic.
root <- commandArgs(trailingOnly = TRUE)[1]
options(warn = 2)
d <- read.csv(file.path(root, "plot-traces.csv"))
palette <- c(prefix = "#3869a1", tail = "#e39537", boundary = "#808080", stream = "#22806b")
png(file.path(root, "worker-timelines.png"), width = 1560, height = 1000,
    res = 150, type = "quartz")
par(mfrow = c(2, 2), mar = c(4, 5.5, 3.5, 1.5), oma = c(3, 0, 2.5, 0),
    family = "sans", las = 1, cex = 0.9)
for (cutoff in c("zero", "meissel")) {
    xmax <- max(d$total_ms[d$cutoff == cutoff]) * 1.05
    for (mode in c("baseline", "hybrid")) {
        s <- d[d$cutoff == cutoff & d$mode == mode, ]
        hybrid <- mode == "hybrid"
        plot(0, 0, type = "n", xlim = c(0, xmax), ylim = c(-0.55, 3.55),
             xaxt = "n", yaxt = "n", xlab = "Elapsed time (ms)", ylab = "", bty = "n",
             main = paste(if (cutoff == "zero") "Zero cutoff" else "Meissel cutoff",
                          if (hybrid) "- Phi then stream" else "- one Phi batch"))
        axis(1)
        axis(2, at = if (hybrid) c(3, 2, 1, 0) else c(3, 2), tick = FALSE,
             labels = if (hybrid) c("Phi 1", "Phi 2", "Stream 1", "Stream 2") else c("Worker 1", "Worker 2"))
        abline(v = axTicks(1), col = "#eeeeee", lwd = 1)
        for (j in seq_len(nrow(s))) {
            y <- 3 - 2 * s$group[j] - s$slot[j]
            rect(s$begin_ms[j], y - 0.25, s$end_ms[j], y + 0.25,
                 col = palette[s$kind[j]], border = NA)
        }
        if (!hybrid) {
            text(xmax * 0.02, 0.7, if (cutoff == "zero")
                 "The tail finishes while the first chunk is still running." else
                 "The balanced prefix finishes before the tail.", adj = 0, cex = 0.85)
        }
        mtext(sprintf("Selected warmed call: %.3f ms", unique(s$total_ms)), side = 3, line = 0.15, cex = 0.76)
    }
}
mtext("Matched P2 worker timelines: x = 8e10, two threads, table limit = 339410",
      side = 3, outer = TRUE, line = 0.6, font = 2)
mtext("Blue: Phi prefix   Orange: Phi tail   Gray: boundary chunk   Green: stream worker (including seed/setup)",
      side = 1, outer = TRUE, line = 0.4, cex = 0.78)
mtext("Separate instrumented traces; row scales differ. Bars show elapsed work intervals, not CPU utilization.",
      side = 1, outer = TRUE, line = 1.7, cex = 0.75)
invisible(dev.off())
