# Native Quartz bitmap output, no extra R packages.
options(warn = 2)
root <- commandArgs(trailingOnly = TRUE)[1]
q <- read.csv(file.path(root, "plot-query.csv"))
d <- read.csv(file.path(root, "timings-summary.csv"))
d <- d[d$mode == "p2" & d$cutoff == "zero" & d$x == 8e10 & d$state == "warm", ]
variants <- c("baseline", "head1", "head2", "head4", "head8", "head32", "hybrid16")
colors <- c("#303b47", "#007f66", "#368cb0", "#8b6caf", "#c89523", "#939ba2", "#cf5c52")
png(file.path(root, "initial-chunks.png"), width = 1600, height = 760, res = 150, type = "quartz")
par(mfrow = c(1, 2), mar = c(5, 4.5, 4, 1), oma = c(2, 0, 2, 0), las = 1, family = "sans")
barplot(q$phi_ms, names.arg = seq_len(nrow(q)), col = "#3869a1", border = NA,
        xlab = "Query position in the first chunk", ylab = "Phi time (ms)",
        main = "Several early queries carry substantial work", cex.names = 0.65,
        ylim = c(0, max(q$phi_ms) * 1.15))
mtext("One representative warmed original-scheduler trace", side = 3, line = 0.3, cex = 0.75)
plot(0, 0, type = "n", xlim = c(1, 8), ylim = c(0, max(d$primary_ms) * 1.30),
     xaxt = "n", xlab = "Requested threads", ylab = "Warmed P2 time (ms)",
     main = "Smaller initial chunks improve scaling", bty = "l")
axis(1, at = c(1, 2, 4, 8)); abline(h = axTicks(2), col = "#eeeeee")
for (i in seq_along(variants)) {
    s <- d[d$variant == variants[i], ]; s <- s[order(s$threads), ]
    lines(s$threads, s$primary_ms, type = "b", pch = 15 + ((i - 1) %% 5),
          col = colors[i], lty = if (i == 6) 2 else 1, lwd = if (i == 2) 2.5 else 1.5)
}
legend("topright", legend = c("Original 32", "Initial 1", "Initial 2", "Initial 4", "Initial 8", "Queue control 32", "Unchanged hybrid 16L"),
       col = colors, lty = c(1, 1, 1, 1, 1, 2, 1), lwd = 2, bty = "n", cex = 0.72)
mtext("Zero-cutoff P2 at x = 8e10, table limit = 339410", side = 3, outer = TRUE, line = 0.5, font = 2)
mtext("Left: diagnostic query windows. Right: separate uninstrumented timings, median of process medians; all later chunks remain 32 queries.",
      side = 1, outer = TRUE, line = 0.5, cex = 0.76)
invisible(dev.off())
