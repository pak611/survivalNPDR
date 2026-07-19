#!/usr/bin/env Rscript

library(survival)

# ── Paths (portable: resolve relative to this script's location) ──────────────
script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("--file=", "", script_args[grep("--file=", script_args)])
if (length(script_file) > 0) {
  root <- normalizePath(file.path(dirname(script_file), "../..")
) } else {
  root <- getwd()
}
input_path <- file.path(root, "data", "GSE9893", "GSE9893_clinicalData.txt")
output_dir <- file.path(root, "figures")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Read clinical data (tab-delimited, preserve column names)
clin <- read.delim(input_path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# Columns of interest
col_time <- "Follow-up period (months)"
col_status <- "State of health"

if (!all(c(col_time, col_status) %in% colnames(clin))) {
  stop("Required columns not found in clinical data.")
}

# Clean and parse time/status
raw_time <- trimws(clin[[col_time]])
raw_status <- trimws(clin[[col_status]])

time_months <- suppressWarnings(as.numeric(raw_time))
status_clean <- tolower(raw_status)

# Event: deceased; Censored: alive / alive with meta / others
status_event <- ifelse(grepl("^deceased", status_clean), 1,
                       ifelse(status_clean == "" | is.na(status_clean), NA, 0))

# Drop rows with missing time or status
valid_idx <- !is.na(time_months) & !is.na(status_event)
clin_valid <- clin[valid_idx, , drop = FALSE]
time_valid <- time_months[valid_idx]
status_valid <- status_event[valid_idx]

# Summary statistics
n_total <- length(time_valid)
n_events <- sum(status_valid == 1)
n_censored <- sum(status_valid == 0)

summary_stats <- data.frame(
  n_total = n_total,
  n_events = n_events,
  n_censored = n_censored,
  min_time = min(time_valid),
  q1_time = as.numeric(quantile(time_valid, 0.25)),
  median_time = median(time_valid),
  mean_time = mean(time_valid),
  q3_time = as.numeric(quantile(time_valid, 0.75)),
  max_time = max(time_valid),
  sd_time = sd(time_valid)
)

cat("GSE9893 Survival Summary\n")
cat(sprintf("  Total samples: %d\n", n_total))
cat(sprintf("  Events (deceased): %d\n", n_events))
cat(sprintf("  Censored (alive/other): %d\n", n_censored))
cat(sprintf("  Time (months): min=%.2f, Q1=%.2f, median=%.2f, mean=%.2f, Q3=%.2f, max=%.2f, sd=%.2f\n",
            summary_stats$min_time, summary_stats$q1_time, summary_stats$median_time,
            summary_stats$mean_time, summary_stats$q3_time, summary_stats$max_time,
            summary_stats$sd_time))

# Kaplan-Meier fit
km_fit <- survfit(Surv(time_valid, status_valid) ~ 1)

# Save survival curve and histogram
curve_path <- file.path(output_dir, "gse9893_survival_curve.png")
hist_path <- file.path(output_dir, "gse9893_survival_time_hist.png")
combo_path <- file.path(output_dir, "gse9893_survival_figure.png")
overlay_path <- file.path(output_dir, "gse9893_survival_overlay.png")

png(curve_path, width = 900, height = 700)
plot(km_fit, xlab = "Months", ylab = "Survival probability",
     main = "GSE9893 Kaplan-Meier Survival Curve",
     conf.int = TRUE, col = "#2C7FB8", lwd = 2)
grid()
dev.off()

# Overlay figure: KM curve with histogram (density) on secondary axis
hist_obj <- hist(time_valid, breaks = 30, plot = FALSE)
density_vals <- hist_obj$density
max_density <- max(density_vals)

png(overlay_path, width = 1000, height = 700)
par(mar = c(6, 5, 5, 5), cex = 1.3, cex.axis = 1.1, cex.lab = 1.2, cex.main = 1.3)
plot(km_fit, xlab = "Months", ylab = "Survival probability",
  main = "",
  conf.int = TRUE, col = "#2C7FB8", lwd = 2)
grid()
title(main = "GSE9893 Survival: KM + Time Distribution", line = 2)
legend("topright",
  legend = c("Kaplan–Meier survival", "Survival time density"),
  col = c("#2C7FB8", rgb(0.49, 0.80, 0.73, 0.8)),
  lty = c(1, NA), lwd = c(2, NA),
  pch = c(NA, 15), pt.cex = 2,
  bty = "n", cex = 1.0)
mtext("Bottom: Survival time histogram (density)", side = 1, line = 4.2, cex = 1.0)
par(new = TRUE)
plot(hist_obj$mids, density_vals, type = "h", lwd = 8,
  xlim = range(hist_obj$breaks), ylim = c(0, max_density * 1.05),
  axes = FALSE, xlab = "", ylab = "", col = rgb(0.49, 0.80, 0.73, 0.45))
axis(4)
mtext("Density", side = 4, line = 3, cex = 1.1)
dev.off()

png(hist_path, width = 900, height = 700)
hist(time_valid, breaks = 30, col = "#7FCDBB", border = "white",
     main = "GSE9893 Survival Time Distribution",
     xlab = "Months")
lines(density(time_valid), col = "#2C7FB8", lwd = 2)
grid()
dev.off()

# Combined, presentation-ready figure (KM + histogram)
png(combo_path, width = 1400, height = 650)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))
plot(km_fit, xlab = "Months", ylab = "Survival probability",
  main = "GSE9893 Kaplan-Meier",
  conf.int = TRUE, col = "#2C7FB8", lwd = 2)
grid()
hist(time_valid, breaks = 30, col = "#7FCDBB", border = "white",
  main = "Survival Time Distribution",
  xlab = "Months")
lines(density(time_valid), col = "#2C7FB8", lwd = 2)
grid()
dev.off()

# Write summary stats to CSV
write.csv(summary_stats, file.path(output_dir, "gse9893_survival_summary.csv"), row.names = FALSE)

cat(sprintf("\nSaved survival curve to: %s\n", curve_path))
cat(sprintf("Saved survival time histogram to: %s\n", hist_path))
cat(sprintf("Saved combined figure to: %s\n", combo_path))
cat(sprintf("Saved overlay figure to: %s\n", overlay_path))
cat(sprintf("Saved summary stats to: %s\n", file.path(output_dir, "gse9893_survival_summary.csv")))
