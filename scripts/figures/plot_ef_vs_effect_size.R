suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ggrepel)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
tables_dir   <- "results"
graphics_dir <- "figures"
dir.create(graphics_dir, showWarnings = FALSE, recursive = TRUE)

int_csv  <- file.path(tables_dir, "interaction_effect_size_sweep_rsurv_results_summary.csv")
main_csv <- file.path(tables_dir, "main_effect_size_sweep_rsurv_results_summary.csv")

int_df  <- read.csv(int_csv,  stringsAsFactors = FALSE)
main_df <- read.csv(main_csv, stringsAsFactors = FALSE)

# ── Fold-level data for SE = SD / sqrt(n_folds) and Wilcoxon tests ───────────
int_fold  <- read.csv(file.path(tables_dir, "interaction_effect_size_sweep_rsurv_fold_metrics.csv"),
                       stringsAsFactors = FALSE)
main_fold <- read.csv(file.path(tables_dir, "main_effect_size_sweep_rsurv_fold_metrics.csv"),
                       stringsAsFactors = FALSE)
int_n  <- aggregate(fold ~ method + beta_int,  data = int_fold,  FUN = length)
main_n <- aggregate(fold ~ method + beta_main, data = main_fold, FUN = length)
colnames(int_n)[3]  <- "n_folds"
colnames(main_n)[3] <- "n_folds"
int_df  <- merge(int_df,  int_n,  by = c("method", "beta_int"),  sort = FALSE)
main_df <- merge(main_df, main_n, by = c("method", "beta_main"), sort = FALSE)
int_df$se_ef_1pct   <- int_df$sd_ef_at_1pct   / sqrt(int_df$n_folds)
int_df$se_ef_5pct   <- int_df$sd_ef_at_5pct   / sqrt(int_df$n_folds)
int_df$se_ef_10pct  <- int_df$sd_ef_at_10pct  / sqrt(int_df$n_folds)
main_df$se_ef_1pct  <- main_df$sd_ef_at_1pct  / sqrt(main_df$n_folds)
main_df$se_ef_5pct  <- main_df$sd_ef_at_5pct  / sqrt(main_df$n_folds)
main_df$se_ef_10pct <- main_df$sd_ef_at_10pct / sqrt(main_df$n_folds)

# ── Method ordering and palette ───────────────────────────────────────────────
method_order <- c("Cox", "SVM", "Ranger", "LASSO-Cox",
                  "sNPDR-signed", "sNPDR-absolute", "LASSO-sNPDR")

palette <- c(
  "Cox"           = "#E41A1C",
  "SVM"           = "#FF7F00",
  "Ranger"        = "#4DAF4A",
  "LASSO-Cox"     = "#F781BF",
  "sNPDR-signed"  = "#1F78B4",
  "sNPDR-absolute"= "#984EA3",
  "LASSO-sNPDR"   = "#A65628"
)

shapes <- c(16, 17, 15, 4, 18, 3, 8)
names(shapes) <- method_order

# ── Base theme ────────────────────────────────────────────────────────────────
theme_paper <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_line(colour = "grey92"),
    strip.background   = element_rect(fill = "grey94", colour = "grey60"),
    strip.text         = element_text(face = "bold", size = 10),
    axis.title         = element_text(size = 11),
    axis.text          = element_text(size = 10),
    legend.key.size    = unit(0.45, "cm"),
    legend.text        = element_text(size = 9),
    legend.title       = element_text(size = 9, face = "bold"),
    plot.title         = element_text(face = "bold", size = 11),
    plot.subtitle      = element_text(size = 9, colour = "grey40")
  )

# ── Filter helpers ────────────────────────────────────────────────────────────
prep <- function(df) {
  df <- df[df$method %in% method_order, ]
  df$method <- factor(df$method, levels = method_order)
  df
}

int_methods  <- method_order
main_methods <- method_order

int_df  <- prep(int_df[int_df$method  %in% int_methods, ])
main_df <- prep(main_df[main_df$method %in% main_methods, ])

# ── Helper: pivot EF columns to long form ─────────────────────────────────────
make_ef_long <- function(df, beta_col) {
  df %>%
    select(method, beta = all_of(beta_col),
           `EF@1%`  = mean_ef_at_1pct,  sd_1  = sd_ef_at_1pct,  se_1  = se_ef_1pct,
           `EF@5%`  = mean_ef_at_5pct,  sd_5  = sd_ef_at_5pct,  se_5  = se_ef_5pct,
           `EF@10%` = mean_ef_at_10pct, sd_10 = sd_ef_at_10pct, se_10 = se_ef_10pct) %>%
    pivot_longer(
      cols      = c(`EF@1%`, `EF@5%`, `EF@10%`),
      names_to  = "k",
      values_to = "mean_ef"
    ) %>%
    mutate(
      se_ef = case_when(
        k == "EF@1%"  ~ se_1,
        k == "EF@5%"  ~ se_5,
        k == "EF@10%" ~ se_10
      ),
      k = factor(k, levels = c("EF@1%", "EF@5%", "EF@10%"))
    ) %>%
    select(method, beta, k, mean_ef, se_ef)
}

# ── Helper: build one faceted EF line plot ────────────────────────────────────
plot_ef <- function(long_df, x_label, title, subtitle, random_ef) {
  ggplot(long_df,
         aes(x = beta, y = mean_ef,
             colour = method, shape = method, group = method)) +
    geom_hline(yintercept = 1, linetype = "dashed",
               colour = "grey55", linewidth = 0.5) +
    geom_hline(aes(yintercept = random_ef),
               linetype = "dotted", colour = "grey70", linewidth = 0.5) +
    geom_ribbon(aes(ymin = pmax(mean_ef - 1.96 * se_ef, 0),
                    ymax = mean_ef + 1.96 * se_ef,
                    fill = method),
                alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.85) +
    geom_point(size = 2.8) +
    facet_wrap(~k, nrow = 1, scales = "free_y") +
    scale_colour_manual(values = palette, name = "Method") +
    scale_fill_manual(values   = palette, name = "Method") +
    scale_shape_manual(values  = shapes,  name = "Method") +
    scale_x_continuous(breaks = unique(long_df$beta)) +
    scale_y_continuous(limits = c(0, NA),
                       expand = expansion(mult = c(0.02, 0.12))) +
    labs(title = title, subtitle = subtitle,
         x = x_label, y = "Enrichment Factor") +
    theme_paper +
    theme(legend.position = "right")
}

# ── Interaction panel ─────────────────────────────────────────────────────────
r_int    <- int_df[1, ]
# EF random baseline: n_functional / (p * k_fraction) normalised to 1x by definition,
# but the scripts report EF = observed / expected where expected = k * n_func / p.
# Random = 1 (already shown as dashed line); no extra dotted line needed.
int_long <- make_ef_long(int_df, "beta_int")
sub_int  <- sprintf(
  "Lognormal AFT, n=%d, p=%d, %d interaction features, main effect size=%.1f, ~20%% censoring",
  r_int$n, r_int$p, r_int$n_int, r_int$beta_main
)

p_int <- plot_ef(
  int_long,
  x_label  = "Effect Size",
  title    = "Interaction effects — EF vs. interaction effect size",
  subtitle = sub_int,
  random_ef = 1
)

# ── Main effects panel ────────────────────────────────────────────────────────
r_main    <- main_df[1, ]
main_long <- make_ef_long(main_df, "beta_main")
sub_main  <- sprintf(
  "Lognormal AFT, n=%d, p=%d, %d causal features, interaction effect size=0, ~20%% censoring",
  r_main$n, r_main$p, r_main$n_main
)

p_main <- plot_ef(
  main_long,
  x_label  = "Effect Size",
  title    = "Main effects — EF vs. main effect size",
  subtitle = sub_main,
  random_ef = 1
)

# ── Combined figure ───────────────────────────────────────────────────────────
fig <- (p_int / p_main) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 13))

ggsave(
  file.path(graphics_dir, "ef_vs_effect_size.png"),
  plot = fig, width = 14, height = 9, dpi = 300
)
ggsave(
  file.path(graphics_dir, "ef_vs_effect_size.pdf"),
  plot = fig, width = 14, height = 9, device = cairo_pdf
)

cat("EF vs effect size figure saved to", graphics_dir, "\n")

# ── Wilcoxon tests: EF at each beta vs. beta=0 null ──────────────────────────
run_wilcox_ef <- function(fold_df, beta_col, ef_col, label, k_label) {
  methods <- sort(unique(fold_df$method))
  betas   <- sort(unique(fold_df[[beta_col]]))
  beta0   <- min(betas)
  cat(sprintf("\n=== %s effects: Wilcoxon %s vs. %s=%.2f null (one-sided, greater) ===\n",
              label, k_label, beta_col, beta0))
  cat(sprintf("%-20s", "Method"))
  for (b in betas[betas != beta0]) cat(sprintf("  %s=%.2f", beta_col, b))
  cat("\n")
  for (m in methods) {
    cat(sprintf("%-20s", m))
    null_vals <- fold_df[fold_df$method == m & abs(fold_df[[beta_col]] - beta0) < 1e-9, ef_col]
    for (b in betas[betas != beta0]) {
      alt_vals <- fold_df[fold_df$method == m & abs(fold_df[[beta_col]] - b) < 1e-9, ef_col]
      pval <- tryCatch(
        wilcox.test(alt_vals, null_vals, alternative = "greater")$p.value,
        error = function(e) NA_real_)
      sig <- if (is.na(pval)) "  " else if (pval < 0.001) "***" else if (pval < 0.01) "** " else if (pval < 0.05) "*  " else "   "
      cat(sprintf("  %6.4f%s", ifelse(is.na(pval), NaN, pval), sig))
    }
    cat("\n")
  }
}

kmap <- c(ef_at_1pct = "EF@1%", ef_at_5pct = "EF@5%", ef_at_10pct = "EF@10%")
for (kv in names(kmap)) {
  run_wilcox_ef(int_fold,  "beta_int",  kv, "Interaction", kmap[kv])
  run_wilcox_ef(main_fold, "beta_main", kv, "Main",        kmap[kv])
}
cat("\nSignificance: *** p<0.001  ** p<0.01  * p<0.05\n")

# ── Single-panel EF(0.01) lineplots (median ± IQR across simulations) ─────────
plot_ef1pct <- function(summ_df, beta_col, x_label) {
  ord <- c("Cox", "SVM", "Ranger", "LASSO-Cox",
           "sNPDR-signed", "sNPDR-absolute", "LASSO-sNPDR")
  pal <- c(
    "Cox"            = "#E41A1C",
    "SVM"            = "#FF7F00",
    "Ranger"         = "#984EA3",
    "LASSO-Cox"      = "#F781BF",
    "sNPDR-signed"   = "#1F78B4",
    "sNPDR-absolute" = "#33A02C",
    "LASSO-sNPDR"    = "#006D77"
  )
  shp <- c(Cox=16, SVM=17, Ranger=15, "LASSO-Cox"=25,
            "sNPDR-signed"=18, "sNPDR-absolute"=3, "LASSO-sNPDR"=6)

  df <- summ_df
  df$method <- as.character(df$method)
  df <- df[df$method %in% ord, ]
  df$method <- factor(df$method, levels = ord)
  betas     <- sort(unique(df[[beta_col]]))
  n_sims    <- round(mean(df$n_sims, na.rm = TRUE))
  last_df   <- df[df[[beta_col]] == max(df[[beta_col]]), ]

  ggplot(df, aes_string(x = beta_col, y = "median_ef_at_1pct",
                        colour = "method", fill = "method",
                        group = "method", shape = "method")) +
    geom_ribbon(aes(ymin = pmax(q1_ef_at_1pct, 0),
                    ymax = q3_ef_at_1pct),
                alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(
      data          = last_df,
      aes(label     = method),
      direction     = "y",
      hjust         = 0,
      nudge_x       = diff(range(betas)) * 0.04,
      segment.size  = 0.3,
      segment.color = "grey50",
      size          = 5.0,
      fontface      = "bold",
      force         = 3,
      force_pull    = 0.5,
      max.overlaps  = Inf,
      colour        = pal[as.character(last_df$method)],
      show.legend   = FALSE
    ) +
    scale_colour_manual(values = pal, name = "Method", breaks = ord) +
    scale_fill_manual(values   = pal, name = "Method", breaks = ord) +
    scale_shape_manual(values  = shp, name = "Method", breaks = ord) +
    scale_x_continuous(breaks = betas,
                       expand = expansion(mult = c(0.02, 0.30))) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.02, 0.12))) +
    labs(x = x_label,
         y = sprintf("EF(0.01) (median \u00b1 IQR, %d simulations)", n_sims)) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92"),
      axis.title       = element_text(size = 14, face = "bold"),
      axis.text        = element_text(size = 12),
      legend.text      = element_text(size = 12),
      legend.title     = element_text(size = 13, face = "bold"),
      legend.position  = "none",
      plot.title       = element_text(size = 15, face = "bold"),
      plot.subtitle    = element_text(size = 10, colour = "grey40"),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 1)
    )
}

r_main <- main_df[1, ]
p_main_1pct <- plot_ef1pct(
  summ_df  = main_df,
  beta_col = "beta_main",
  x_label  = "Main Effect Size"
)
ggsave(
  file.path(graphics_dir, "main_ef1pct_sweep_rsurv_lineplot.png"),
  plot = p_main_1pct, width = 12, height = 8.5, dpi = 300
)
cat("Main EF(0.01) lineplot saved.\n")

r_int <- int_df[1, ]
p_int_1pct <- plot_ef1pct(
  summ_df  = int_df,
  beta_col = "beta_int",
  x_label  = "Interaction Effect Size"
)
ggsave(
  file.path(graphics_dir, "interaction_ef1pct_sweep_rsurv_lineplot.png"),
  plot = p_int_1pct, width = 12, height = 8.5, dpi = 300
)
cat("Interaction EF(0.01) lineplot saved.\n")
