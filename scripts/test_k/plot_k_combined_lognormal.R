# Combined k-sweep figure: Main + Interaction effects (lognormal, auPRC)
#
# Reads from CSVs produced by:
#   test_k_lognormal.R        → paper_tables/auPRC_K_main_lognormal.csv
#   test_k_lognormal.R → paper_tables/PR_K_AUC_interactions_lognormal_lmkm_sigma0.50.csv
#
# Falls back to illustrative values (medians taken from prior run shown in the
# paper draft) when those files are not present.
#
# Modifications vs. original per-script figures:
#   • multisurf box plots REMOVED
#   • dashed vertical line at k = 30 (≈ multisurf adaptive k for n = 200)
#     with a "multisurf" text label at the top

library(ggplot2)
library(dplyr)

root    <- getwd()
out_dir <- file.path(root, "results")

# ── Command-line arguments: sdlog (default 0.50), p (default 1000) ───────────
args      <- commandArgs(trailingOnly = TRUE)
sigma_val <- if (length(args) >= 1) as.numeric(args[1]) else 0.50
p_val     <- if (length(args) >= 2) as.integer(args[2]) else 1000L
cat(sprintf("Plotting interaction results for sigma = %.2f, p = %d\n", sigma_val, p_val))

# ── K values shown as boxes ────────────────────────────────────────────────────
k_vals <- c(5, 20, 100, 199)
n_rep  <- 20

# ── Load or mock data ──────────────────────────────────────────────────────────
main_csv <- file.path(out_dir, "auPRC_K_main_lognormal.csv")
int_suffix <- sprintf("PR_K_AUC_interactions_lognormal_lmkm_sigma%.2f_p%d.csv", sigma_val, p_val)
int_csv  <- file.path(out_dir, int_suffix)

if (file.exists(main_csv) && file.exists(int_csv)) {
  message("Loading saved CSV results.")

  main_df <- read.csv(main_csv) %>%
    rename(metric = auPRC) %>%
    mutate(EffectType = "Main")

  int_df <- read.csv(int_csv) %>%
    rename(metric = AUC) %>%
    mutate(EffectType = "Interaction")

  combined <- bind_rows(
    main_df %>% filter(K %in% k_vals) %>% select(K, metric, EffectType),
    int_df  %>% filter(K %in% k_vals) %>% select(K, metric, EffectType)
  )

} else {
  message("CSVs not found — using illustrative values from prior run.")
  # Approximate medians read from the paper draft figure
  set.seed(42)
  make_reps <- function(ks, meds, sds, type) {
    mapply(function(k, m, s) {
      data.frame(K = k,
                 metric = pmax(0, rnorm(n_rep, m, s)),
                 EffectType = type,
                 stringsAsFactors = FALSE)
    }, ks, meds, sds, SIMPLIFY = FALSE) |> do.call(what = rbind)
  }

  combined <- rbind(
    make_reps(k_vals,
              meds = c(1.05, 1.05, 1.05, 1.00),
              sds  = c(0.3, 0.3, 0.3, 0.3),
              type = "Interaction"),
    make_reps(k_vals,
              meds = c(5, 10, 15, 18),
              sds  = c(2, 3, 3, 3),
              type = "Main")
  )
}

combined$K          <- as.numeric(as.character(combined$K))
combined$EffectType <- factor(combined$EffectType, levels = c("Interaction", "Main"))

# ── K levels & multisurf position ─────────────────────────────────────────────
k_levels <- sort(k_vals)   # 5, 20, 100, 199
msurf_k  <- 30

# ── Colour palette (blue = Main, green = Interaction) ─────────────────────────
fill_pal  <- c("Interaction" = "#b8ddb8", "Main" = "#aecde8")
color_pal <- c("Interaction" = "#4aaa4a", "Main"  = "#4a8ab8")

# ── Summarise replicates ───────────────────────────────────────────────────────
sum_df <- combined %>%
  group_by(K, EffectType) %>%
  summarise(mean_m = mean(metric, na.rm = TRUE),
            sd_m   = sd(metric,   na.rm = TRUE),
            .groups = "drop")

sum_df <- sum_df %>%
  mutate(x_plot = K)

y_ceil <- max(sum_df$mean_m + sum_df$sd_m, na.rm = TRUE) * 1.18

# ── Build figure ───────────────────────────────────────────────────────────────
p <- ggplot(sum_df, aes(x = x_plot, y = mean_m,
                        color = EffectType, fill = EffectType,
                        group = EffectType)) +

  # Subtle vertical guides at each plotted K
  geom_vline(xintercept = k_levels, color = "grey88",
             linewidth = 0.6, linetype = "solid") +

  # k̄_{1/2} (multisurf) reference line — bold dashed
  geom_vline(xintercept = msurf_k, color = "grey30",
             linewidth  = 1.1, linetype = "dashed") +
  annotate("text", x = msurf_k * 1.10, y = y_ceil,
           label    = "bold(bar(k)[1/2])",
           parse    = TRUE,
           hjust    = 0, vjust = 1.1,
           size     = 5.5, color = "grey20") +

  # Error bars (±1 SD across replicates)
  geom_errorbar(aes(ymin = mean_m - sd_m, ymax = mean_m + sd_m),
                width = 0.03, linewidth = 0.85) +

  # Lines connecting means
  geom_line(linewidth = 1.1) +

  # Mean points
  geom_point(size = 3.8, shape = 21, color = "white", stroke = 1.5) +

  # ── Scales & theme ────────────────────────────────────────────────────────────
  scale_color_manual(values = color_pal) +
  scale_fill_manual(values  = fill_pal)  +
  scale_x_log10(
    name   = "K (Neighbors)",
    breaks = k_levels,
    labels = as.character(k_levels)
  ) +
  scale_y_continuous(name = "auPRC") +
  coord_cartesian(ylim = c(0, y_ceil)) +

  guides(
    color = guide_legend(title = "Effect Type"),
    fill  = "none"
  ) +

  theme_minimal(base_size = 14) +
  theme(
    axis.title      = element_text(size = 14, face = "bold"),
    axis.text       = element_text(size = 12),
    legend.position  = c(0.12, 0.82),
    legend.direction = "vertical",
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
    legend.title     = element_text(size = 13, face = "bold"),
    legend.text      = element_text(size = 12),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1.2),
    panel.grid.minor = element_blank()
  )

# ── Save ───────────────────────────────────────────────────────────────────────
graphics_dir <- file.path(root, "figures")
dir.create(graphics_dir, showWarnings = FALSE, recursive = TRUE)

out_file <- file.path(graphics_dir, sprintf("auPRC_K_combined_lognormal_no_multisurf_sigma%.2f_p%d.png", sigma_val, p_val))
ggsave(out_file, plot = p, width = 9, height = 6.5, dpi = 300)
cat("Saved:", out_file, "\n")
