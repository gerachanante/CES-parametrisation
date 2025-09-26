options(scipen = 999) #avoids scientific notation unless necessary
# install.packages(c("micEconCES","dplyr","readr","purrr","ggplot2","parallel",
#"ggpmisc","pheatmap","GGally","ggcorrplot","ggridges","ggpubr"))
#install.packages("ggpubr")
#install.packages("fmsb")
library(dplyr)
library(tidyr)
library(readr)
library(viridis)
library(ggplot2)
library(ggpmisc)
library(ggrepel)
library(pheatmap)
library(patchwork)
library(tibble)
library(GGally)
library(ggcorrplot)
library(ggridges)
library(readxl)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggpubr)
library(rstatix)
library(fmsb)
library(grid)
library(gridExtra)


### Settings
setwd("C:/Users/escami_g/OneDrive - Paul Scherrer Institut/05.Models/MERGE updates/CES-parametrisation/finest grid")
infile <- "MERGE macro.csv"

# Load original data and scale
df <- read_csv(infile, show_col_types = TRUE)

dfS <- df %>%
  group_by(r) %>%
  mutate(
    Ybase = Y[t == 2022][1],
    Kbase = K[t == 2022][1],
    Lbase = L[t == 2022][1],
    Ebase = E[t == 2022][1]
  ) %>%
  mutate(
    Ys = Y/Ybase,
    Ks = K/Kbase,
    Ls = L/Lbase,
    Es = E/Ebase
  ) %>%
  ungroup()

# Load previous results
results_grid <- readRDS("results_run1.rds")
best_methods <- read_csv("CES_best_methods.csv")
aic_weights <- read_csv("CES_AICc_weights.csv")
iam_table <- read_csv("IAM_params.csv")
# Adding validity to the runs, economically feasible and converged from the solver
add_validity <- function(df) {
  df %>%
    mutate(
      across(any_of(c("delta_KL","delta_VAE","gamma","nu","lambda")),
             ~ suppressWarnings(as.numeric(.))),
      valid =
        (conv %in% TRUE) &
        is.finite(delta_KL) & delta_KL >= 0 & delta_KL <= 1 &
        is.finite(delta_VAE) & delta_VAE >= 0 & delta_VAE <= 1 &
        is.finite(gamma) & gamma > 0.2 & gamma < 5 &
        is.finite(nu) & nu > 0.2 & nu < 5 &
        is.finite(lambda)
    )
}
results_grid <- results_grid %>% select(-any_of("valid")) %>% add_validity()
results_grid_valid <- results_grid %>% filter(valid)
results_grid_invalid <- results_grid %>% filter(!valid)


# PReset theme for graphs
theme_nat <- function(base_size = 11){
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "#e9e9e9"),
      axis.title  = element_text(colour = "#111111"),
      axis.text   = element_text(colour = "#4d4d4d"),
      plot.title  = element_text(face = "bold", size = base_size+2, colour = "#111111",
                                 margin = margin(b = 6)),
      plot.subtitle = element_text(size = base_size, colour = "#4d4d4d", margin = margin(b=8)),
      plot.caption = element_text(size = base_size-1, colour = "#4d4d4d"),
      legend.position = "right",
      legend.title = element_text(colour = "#4d4d4d"),
      legend.text  = element_text(colour = "#4d4d4d"),
      strip.text   = element_text(face = "bold", colour = "#111111")
    )
}

NAVY <- "#143d66"
BLUE <- "#2E75B6"
BBLUE <- "#DAE3F3"
LGREY <- "#c7c7c7"
MGREY <- "#A6A6A6"
DGREY <- "#4d4d4d"
YELLOW <- "#FFC000"
RED <- "#C00000"
ORANGE <- "#FFDFC5"

# Dynamic captions so symbols match each plot’s content
cap_none      <- ""
cap_elast     <- "Symbols: σK-L, σVA–E = substitution elasticities."
cap_rho       <- "Symbols: ρKL, ρVA–E = CES exponents on the K–L and VA–E nests."
cap_fit       <- expression(paste("Symbols: ", R^2, " = coefficient of determination; ε = residual."))
cap_params    <- "Symbols: γ, λ, ν = scale/growth/curvature; δK-L, δVA–E = share parameters."


pct <- function(x) scales::percent(x, accuracy = 1)
clip01 <- function(x) pmax(pmin(x, 1), 0)
near_   <- function(x, a, tol = 1e-12) abs(x - a) <= tol

# Numbers
lab_si <- scales::label_number_si(accuracy = 0.1)
fmt_n <- function(x) ifelse(x < 1000,
                            scales::number(x, accuracy = 1, big.mark = ","),
                            lab_si(x))

# MAth labels
lab_sigma_KL  <- bquote(sigma[K - L])
lab_sigma_VAE <- bquote(sigma[VA - E])
lab_rho_KL    <- bquote(rho[KL])
lab_rho_VAE   <- bquote(rho[VAE])

# Consistent method ordering (fastest & most stable first) = by highest convergence share,
# break ties by lowest median runtime (mins)
method_perf <- results_grid_valid %>%
  distinct(r, method, conv, runtime_total) %>%
  group_by(method) %>%
  summarise(
    share_conv = mean(conv, na.rm = TRUE),
    med_rt_min = median(runtime_total, na.rm = TRUE) / 60,
    .groups = "drop"
  ) %>% arrange(desc(share_conv), med_rt_min)
method_order <- method_perf$method

# Global medians for best-method elasticities
med_KL  <- median(best_methods$sigma_KL,  na.rm = TRUE)
med_VAE <- median(best_methods$sigma_VAE, na.rm = TRUE)

# Axis formatting
sigma_axis_breaks <- function(r_lim, s_grid = c(10, 5, 2, 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01)) {
  r_from_sigma <- (1 / s_grid) - 1
  keep <- r_from_sigma >= min(r_lim, na.rm = TRUE) & r_from_sigma <= max(r_lim, na.rm = TRUE)
  list(breaks = r_from_sigma[keep], labels = s_grid[keep])
}

# Load mapping
region_map <- read_excel("MERGE regions proposal.xlsx", sheet = 1) %>%
  rename(r = MERGE) %>% mutate(r = as.character(r))
world <- ne_countries(scale = "medium", returnclass = "sf")
best_map <- best_methods %>% inner_join(region_map, by = "r")
world_best <- world %>% left_join(best_map, by = c("iso_a3" = "ISO3"))


# Small helpers used by plots
plot_bar <- function(df, x, y, title, subtitle="", ylab="", fillcol=LGREY, label_fmt=scales::percent){
  ggplot(df, aes(x=reorder(.data[[x]], .data[[y]]), y=.data[[y]])) +
    geom_col(width=.7, fill=fillcol, colour="white") +
    geom_text(aes(label=if(!is.null(label_fmt)) label_fmt(.data[[y]]) else round(.data[[y]],2)),
              hjust=-0.1, size=3, colour=DGREY) +
    coord_flip(ylim=c(0, max(df[[y]], na.rm=TRUE)*1.1)) +
    labs(title=title, subtitle=subtitle, x=NULL, y=ylab) +
    theme_nat()
}
plot_elasticity_scatter <- function(df, med_KL, med_VAE){
  ggplot(df, aes(sigma_KL, sigma_VAE)) +
    geom_vline(xintercept=med_KL, linetype="dashed", colour=NAVY, linewidth=.4) +
    geom_hline(yintercept=med_VAE, linetype="dashed", colour=NAVY, linewidth=.4) +
    geom_point(aes(shape=method, fill=method), size=2.8, colour=DGREY, alpha=.9) +
    scale_fill_manual(values=scales::hue_pal()(length(unique(df$method)))) +
    scale_shape_manual(values=21:25) +
    ggrepel::geom_text_repel(aes(label=r), colour=DGREY, size=3, max.overlaps=25, box.padding=.3, segment.alpha=.3) +
    labs(title="Best-method elasticities per region (quadrant view)",
         subtitle="Dashed lines = global medians",
         x=expression(sigma[K-L]), y=expression(sigma[VA-E])) +
    theme_nat()
}



### 1. Diagnostics of all runs
rg <- results_grid %>% filter(!is.na(method) & nzchar(method))

# Coverage per (region, method)
cov_region <- rg %>%
  group_by(r, method) %>%
  summarise(
    share_converged = mean(conv,  na.rm = TRUE),
    share_valid     = mean(valid, na.rm = TRUE),
    .groups = "drop"
  )

# Respect your method ordering if it exists; else order by median valid share
method_order_use <- if (exists("method_order")) method_order else
  cov_region %>%
    group_by(method) %>%
    summarise(med_valid = median(share_valid, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(med_valid)) %>% pull(method)

# Build median + IQR across regions for both statuses
cov_stats <- cov_region %>%
  tidyr::pivot_longer(c(share_converged, share_valid),
                      names_to = "which", values_to = "share") %>%
  mutate(
    which  = factor(which,
                    levels = c("share_converged","share_valid"),
                    labels = c("Converged","Valid")),
    method = factor(method, levels = method_order_use)
  ) %>%
  group_by(method, which) %>%
  summarise(med = median(share, na.rm = TRUE),
            q1  = quantile(share, 0.25, na.rm = TRUE),
            q3  = quantile(share, 0.75, na.rm = TRUE),
            .groups = "drop")


# Plot (points = medians; bars = IQR), dodged so the two series don’t sit on top of each other
fig1_1 <- ggplot(cov_stats, aes(x = method, y = med, colour = which, shape = which)) +
  geom_errorbar(aes(ymin = q1, ymax = q3),
                width = .5, position = position_dodge(width = .55), colour = DGREY) +
  geom_point(size = 2.8, position = position_dodge(width = .55)) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_colour_manual(values = c("Converged" = DGREY, "Valid" = NAVY)) +
  scale_shape_manual(values  = c("Converged" = 16,    "Valid" = 18)) +
  labs(
    title    = expression(rho*"-grid coverage by method (across regions)"),
    subtitle = "Points = median share; bars = IQR across regions",
    x = NULL, y = "Share of grid", colour = NULL, shape = NULL, caption = cap_rho
  ) +
  theme_nat() +
  theme(legend.position = "top", legend.direction = "horizontal")
print(fig1_1)

# 1.2 Region × Method convergence (share of methods converged per region)
if (!exists("method_order")) method_order <- sort(unique(results_grid$method))

rm_share <- results_grid %>%
  group_by(r, method) %>%
  summarise(
    p_grid_solved = mean(is.finite(rss)),            # share of rho-grid with an RSS
    .groups = "drop"
  ) %>%
  # make sure every region × method shows (missing -> 0)
  complete(r, method, fill = list(p_grid_solved = 0)) %>%
  mutate(method = factor(method, levels = method_order))

# Order regions by average coverage (just for a tidy display)
ord_r <- rm_share %>%
  group_by(r) %>% summarise(mu = mean(p_grid_solved), .groups = "drop") %>%
  arrange(desc(mu), r) %>% pull(r)
rm_share <- rm_share %>% mutate(r = factor(r, levels = ord_r))

# ---- Fig 1.2 (revised): Convergence by region × method (fraction of ρ-grid solved)
fig1_2 <- ggplot(rm_share, aes(x = method, y = r, fill = p_grid_solved)) +
  geom_tile(width = 0.98, height = 0.98, colour = "white", linewidth = 0.25) +
  # label only the extremes to reduce clutter
  geom_text(
    data = subset(rm_share, p_grid_solved <= .05 | p_grid_solved >= .95),
    aes(label = scales::percent(p_grid_solved, accuracy = 1)),
    colour = "white", size = 2.6, fontface = "bold"
  ) +
  scale_fill_gradientn(
    colours = c(RED, "#F4B9B9", LGREY, "#9BB3C9", NAVY),
    values  = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)),
    limits  = c(0, 1),
    labels  = scales::percent
  ) +
  labs(
    title    = "Convergence by region × method",
    subtitle = "Fill = share of ρ-grid points that produced an RSS (i.e., actually solved)",
    x = "Method", y = "Region", fill = "Solved share", caption = cap_rho
  ) +
  theme_nat() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(fig1_2)


# 1.3 Best-fit rho medians & IQR by method (valid runs) + secondary axis for sigma
if (!exists("method_order")) method_order <- sort(unique(results_grid$method))

rho_long <- results_grid %>%
  filter(valid) %>%                                   # valid runs only
  select(method, rho_KL, rho_VAE) %>%
  tidyr::pivot_longer(c(rho_KL, rho_VAE), names_to = "rho_type", values_to = "rho") %>%
  mutate(
    rho_type = factor(rho_type,                       # avoid recode() masking
                      levels = c("rho_KL","rho_VAE"),
                      labels = c("ρKL","ρVA–E")),
    method   = factor(method, levels = method_order)
  ) %>%
  filter(is.finite(rho))

rho_stats <- rho_long %>%
  group_by(method, rho_type) %>%
  summarise(
    med = median(rho, na.rm = TRUE),
    q1  = quantile(rho, 0.25, na.rm = TRUE),
    q3  = quantile(rho, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# y-limits with a little padding so the sec.axis looks nice
y_min <- min(rho_stats$q1, rho_stats$med, na.rm = TRUE)
y_max <- max(rho_stats$q3, rho_stats$med, na.rm = TRUE)
pad   <- (y_max - y_min) * 0.08
y_lim <- c(y_min - pad, y_max + pad)

fig1_3 <- ggplot(rho_stats, aes(x = method, y = med, colour = rho_type)) +
  geom_errorbar(aes(ymin = q1, ymax = q3),
                width = 0.25, linewidth = 0.6, colour = DGREY,
                position = position_dodge(width = 0.55)) +
  geom_point(position = position_dodge(width = 0.55), size = 2.9) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = LGREY) +  # ρ=0 ⇒ σ=1
  coord_flip() +
  scale_colour_manual(values = c("ρKL" = DGREY, "ρVA–E" = NAVY), name = NULL) +
  scale_y_continuous(
    limits = y_lim,
    name   = expression(rho),
    sec.axis = sec_axis(~ 1/(1 + .), name = expression(sigma))
  ) +
  labs(
    title    = "Best-fit ρ medians with IQR by method",
    subtitle = "Points = medians; bars = IQR (valid runs only). Right axis shows σ = 1/(1+ρ).",
    x = NULL,
    caption  = cap_rho
  ) +
  theme_nat()

print(fig1_3)

    .groups = "drop"

# nice limits so the σ axis maps well (add a little padding)
y_min <- min(rho_stats$q1, rho_stats$med, na.rm = TRUE)
y_max <- max(rho_stats$q3, rho_stats$med, na.rm = TRUE)
pad   <- diff(range(c(y_min, y_max))) * 0.08
y_lim <- c(y_min - pad, y_max + pad)

fig1_3 <- ggplot(rho_stats, aes(x = method, y = med, colour = rho_type)) +
  geom_errorbar(aes(ymin = q1, ymax = q3),
                width = .25, linewidth = .6, colour = DGREY,
                position = position_dodge(width = .55)) +
  geom_point(position = position_dodge(width = .55), size = 2.9) +
  # reference lines: ρ = 0 (⇒ σ = 1)
  geom_hline(yintercept = 0, linetype = "dotted", colour = LGREY) +
  coord_flip() +
  scale_colour_manual(values = c("ρKL" = DGREY, "ρVA–E" = NAVY), name = NULL) +
  scale_y_continuous(
    limits = y_lim,
    name   = expression(rho),
    sec.axis = sec_axis(~ 1 / (1 + .),
                        name = expression(sigma))
  ) +
  labs(
    title    = "Best-fit ρ medians with IQR by method",
    subtitle = "Points = medians; bars = IQR (valid runs only). Right axis shows σ = 1/(1+ρ).",
    x = NULL,
    caption  = cap_rho
  ) +
  theme_nat() +
  theme(axis.text.x = element_text())  # after flip, this is the numeric axis

print(fig1_3)




# 1.4 Runtime by method (minutes)
rt_stats <- results_grid %>%
  group_by(method) %>%
  summarise(med_min = median(runtime_total, na.rm = TRUE)/60, .groups = "drop") %>%
  arrange(med_min)
fig1_4 <- results_grid %>%
  mutate(method = factor(method, levels = rt_stats$method)) %>%
  ggplot(aes(x = method, y = runtime_total/60)) +
  geom_boxplot(width = .72, fill= LGREY, colour = DGREY, outlier.size = 1.4, outlier.alpha = .6) +
  stat_summary(fun = median, geom = "point", size = 3, shape = 21, fill = NAVY, colour = NAVY) +
  coord_flip() +
  labs(title = "Runtime by method", subtitle = "Minutes; navy = method median",
       x = NULL, y = "Runtime (minutes)", caption = cap_none) + theme_nat()
print(fig1_4)

# 1.5 Runtime vs fit quality
rt_quality <- results_grid %>%
  filter(conv, is.finite(R2), R2 > 0) %>%
  mutate(method = factor(method, levels = method_order))

fig1_5 <- ggplot(rt_quality, aes(x = runtime_total/60, y = R2, colour = method)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_y_continuous(limits = c(0,1)) +
  labs(title = "Runtime vs fit quality",
       subtitle = "Each point = converged run; x = runtime (log minutes), y = R²",
       x = "Runtime (minutes, log scale)", y = expression(R^2)) +
  theme_nat() + theme(legend.position = "top")
print(fig1_5)


# 1.5 RSS vs iterations (global quadratic)
results_grid <- results_grid %>%
  mutate(
    method = trimws(as.character(method))
  )

fit_df <- results_grid %>%
  filter(conv, is.finite(iter), iter > 0, is.finite(rss), rss > 0, !is.na(method), nzchar(method))

method_levels <- unique(c(if (exists("method_order")) method_order else character(0),
                          sort(unique(fit_df$method))))
fit_df <- fit_df %>% mutate(method_f = factor(method, levels = method_levels))

pal_methods <- scales::hue_pal()(length(method_levels))
names(pal_methods) <- method_levels

x_cap <- stats::quantile(fit_df$iter, 0.99, na.rm = TRUE)
y_cap <- stats::quantile(fit_df$rss,  0.99, na.rm = TRUE)

fig1_5 <- ggplot(fit_df, aes(x = iter, y = rss, colour = method_f)) +
  geom_point(alpha = 0.7, size = 1.9) +
  scale_colour_manual(values = pal_methods, drop = FALSE, name = "") +
  coord_cartesian(xlim = c(0, x_cap), ylim = c(0, y_cap)) +
  labs(
    title = "More iterations in the solvers don't necessarily achieve better fits",
    subtitle = "Converged runs RSS per number of iterations by solver",
    x = "Iterations", y = "RSS", caption = cap_fit
  ) +
  theme_nat() +
  theme(legend.position = "top")

print(fig1_5)

  
# 1.6 ΔAICc surfaces over the rho-grid for the chosen (r, method)
grid_bestmethod <- results_grid %>%
  inner_join(best_methods %>% select(r, method), by = c("r","method"))

aicc_surface <- grid_bestmethod %>%
  group_by(r, method) %>%
  mutate(dAICc_grid = AICc_plusRho - min(AICc_plusRho, na.rm = TRUE)) %>%
  ungroup()

edge_min <- aicc_surface %>%
  group_by(r, method) %>%
  slice_min(order_by = dAICc_grid, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(on_edge = (on_edge_KL | on_edge_VAE))

edge_counts <- edge_min %>%
  group_by(method) %>%
  summarise(share_edge = mean(on_edge, na.rm = TRUE), .groups = "drop")

fig1_6 <- ggplot(edge_counts, aes(x = reorder(method, share_edge), y = share_edge)) +
  geom_col(fill = LGREY, colour = "white", width = 0.65) +
  geom_text(aes(label = paste0(round(share_edge*100,1),"%")),
            hjust = -0.05, colour = DGREY, size = 3.1) +
  coord_flip(ylim = c(0,1)) +
  labs(title = "Edge solutions by method",
       subtitle = "Share of ΔAICc minima at grid boundaries",
       x = NULL, y = "Share of edge minima") +
  theme_nat()
print(fig1_6)




# 1.9 Failure reasons by method (invalid runs only)
# 1) Group non-convergence messages
msg_group_map <- function(x){
  x <- tolower(paste0(x))
  dplyr::case_when(
    grepl("false convergence", x)                 ~ "False convergence",
    grepl("max ?(iter|imum)", x)                  ~ "Max iterations",
    grepl("reduction", x)                         ~ "Reduction criterion",
    grepl("bounds|box|ftol|xtol|step|line", x)    ~ "Bounds/tolerance",
    TRUE                                          ~ "Unspecified"
  )
}

# 2) Diagnose which econ rule failed (first hit wins, purely for labeling)
econ_fail_reason <- function(df){
  with(df, dplyr::case_when(
    !is.finite(delta_KL)  | delta_KL  < 0 | delta_KL  > 1 ~ "δK–L out of [0,1]",
    !is.finite(delta_VAE) | delta_VAE < 0 | delta_VAE > 1 ~ "δVA–E out of [0,1]",
    !is.finite(gamma)     | gamma <= 0.2 | gamma >= 5     ~ "γ out of (0.2,5)",
    !is.finite(nu)        | nu    <= 0.2 | nu    >= 5     ~ "ν out of (0.2,5)",
    !is.finite(lambda)                                   ~ "λ not finite",
    TRUE                                                 ~ "other econ bound"
  ))
}

# 3) Final per-run status
runs_status <- results_grid %>%
  mutate(
    status = dplyr::case_when(
      !conv                    ~ msg_group_map(msg),            # non-convergence first
      conv & !valid            ~ econ_fail_reason(cur_data()),  # converged but econ-invalid
      conv &  valid            ~ "Valid",                        # converged & passes econ screen
      TRUE                     ~ "Unspecified"
    ),
    status = factor(
      status,
      levels = c(
        # show NON-CONV reasons first in legend/stack
        "False convergence","Max iterations","Reduction criterion",
        "Bounds/tolerance","Unspecified",
        # then econ-invalid reasons for converged runs
        "δK–L out of [0,1]","δVA–E out of [0,1]","γ out of (0.2,5)",
        "ν out of (0.2,5)","λ not finite","other econ bound",
        # finally the valid ones
        "Valid"
      )
    )
  )

# 4) Summarise to counts (or switch to shares if you prefer)
error_summary2 <- runs_status %>%
  count(method, status, name = "n") %>%
  tidyr::complete(method, status, fill = list(n = 0)) %>%
  group_by(method) %>%
  mutate(total = sum(n)) %>%
  ungroup()

method_order <- error_summary2 %>%
  group_by(method) %>%
  summarise(nonvalid = sum(n[status != "Valid"]), .groups = "drop") %>%
  arrange(desc(nonvalid)) %>%
  pull(method)

error_summary2 <- error_summary2 %>%
  mutate(method = factor(method, levels = method_order))

# pretty count labels (e.g. 1.2K)
lab_short <- scales::label_number(
  accuracy  = 0.1,
  big.mark  = ",",
  trim      = TRUE,
  scale_cut = scales::cut_short_scale()  # "", "k", "M", …
)

# 5) Single stacked bar: full scope of all runs, blue = valid, others split by reason
fig1_9 <- ggplot(error_summary2, aes(x = method, y = n, fill = status)) +
  geom_col(width = 0.9, colour = "white") +
  geom_text(
    aes(label = ifelse(n > 0, lab_short(n), "")),
    position = position_stack(vjust = 0.5),
    size = 2.8, colour = "black", fontface = "bold"
  ) +
  coord_flip() +
  scale_y_continuous(labels = lab_short, expand = expansion(mult = c(0, .05))) +
  scale_fill_manual(
    breaks = levels(error_summary2$status),
    values = c(
      # non-conv
      "False convergence"   = "#D5D8FF",
      "Max iterations"      = "#CDFFF5",
      "Reduction criterion" = "#FFDFC5",
      "Bounds/tolerance"    = "#B7E5FF",
      "Unspecified"         = "#D9D9D9",
      # econ-invalid (converged)
      "δK–L out of [0,1]"   = "#B3B3B3",
      "δVA–E out of [0,1]"  = "#A0A0A0",
      "γ out of (0.2,5)"    = "#8D8D8D",
      "ν out of (0.2,5)"    = "#7A7A7A",
      "λ not finite"        = "#676767",
      "other econ bound"    = "#5A5A5A",
      # valid
      "Valid"               = "#00B0F0"
    ),
    drop = FALSE
  ) +
  labs(
    title    = "Newton, L-BFGS-B and BFGS produce valid results",
    subtitle = "Invalid results are due to lack of solver convergence or economic parameters out of bounds",
    x = NULL, y = "Number of runs", fill = NULL
  ) +
  theme_nat() +
  theme(
    legend.position  = "top",
    legend.direction = "horizontal",
    legend.box       = "horizontal",
    legend.margin    = margin(t = 2, b = 2),
    legend.spacing.x = grid::unit(6, "pt")
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, keywidth = grid::unit(.9, "cm")))

print(fig1_9)




# 1.10 Radar chart
# Compute metrics per method
method_summary <- results_grid %>%
  group_by(method) %>%
  summarise(
    conv_share   = mean(conv, na.rm = TRUE),
    med_R2       = median(R2, na.rm = TRUE),
    med_runtime  = median(runtime_total, na.rm = TRUE),
    n_invalid    = sum(!conv, na.rm = TRUE),
    n_total      = n(),
    .groups = "drop"
  ) %>%
  mutate(
    fail_rate = n_invalid / n_total,
    inv_runtime = 1 / (1 + med_runtime/60),  # scale-inverted runtime
    inv_fail    = 1 - fail_rate
  ) %>%
  select(method, conv_share, med_R2, inv_runtime, inv_fail)

# Normalise to [0,1]
scaled <- as.data.frame(method_summary)
scaled[,-1] <- lapply(scaled[,-1], scales::rescale)

# Prepare radar data: max/min rows
radar_data <- rbind(rep(1, ncol(scaled)-1),
                    rep(0, ncol(scaled)-1),
                    scaled[,-1])
colnames(radar_data) <- colnames(scaled)[-1]
rownames(radar_data) <- c("max","min", scaled$method)

# Radar plot
fig1_10 <- fmsb::radarchart(radar_data,
                           axistype = 1,
                           pcol = rainbow(nrow(scaled)),
                           plwd = 2, plty = 1,
                           cglcol = "grey", cglty = 1,
                           axislabcol = "grey30", caxislabels = c("0","0.5","1"),
                           vlcex = 0.8,
                           title = "Method performance radar")
legend("topright", legend = scaled$method,
       col = rainbow(nrow(scaled)), lty = 1, cex = 0.7, bty = "n")
print(fig1_10)


### 2. Valid runs diagnostics
# 2.1 Method performance summary (n, median R², median runtime)
fig2_1 <- results_grid_valid %>%
  group_by(method) %>%
  summarise(n = n(),
            med_R2 = median(R2, na.rm = TRUE),
            med_runtime = median(runtime_total, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(med_R2)) %>%
  mutate(method = factor(method, levels = method)) %>%
  ggplot(aes(x = method, y = n)) +
  geom_col(width = .70, fill = LGREY, colour = "white") +
  geom_text(aes(label = paste0("R²~", round(med_R2, 2), " | t~", round(med_runtime/60, 1), " min")),
            vjust = -0.25, size = 3.1, colour = DGREY) +
  coord_flip() +
  labs(title = "Valid runs by method",
       subtitle = "Bar = count; label = median R² and runtime (min)",
       x = NULL, y = "Valid runs", caption = cap_fit) + theme_nat()
print(fig2_1)

# 2.1b Valid share per method ---
valid_share <- results_grid %>%
  group_by(method) %>%
  summarise(valid_frac = mean(conv, na.rm = TRUE), n = n(), .groups = "drop") %>%
  arrange(desc(valid_frac))

fig2_1b <- ggplot(valid_share, aes(x = reorder(method, valid_frac), y = valid_frac)) +
  geom_col(width = .7, fill = LGREY, colour = "white") +
  geom_text(aes(label = paste0(round(valid_frac*100,1), "% (n=", n, ")")),
            hjust = -0.05, colour = DGREY, size = 3.2) +
  coord_flip(ylim = c(0, 1.05)) +
  labs(title = "Valid share by method",
       subtitle = "Share of runs that converged successfully",
       x = NULL, y = "Share valid", caption = cap_none) +
  theme_nat()
print(fig2_1b)

# 2.2 Elasticity distributions overlaid (same x-axis), by method and overall
fig2_2 <- results_grid_valid %>%
  select(r, method, sigma_KL, sigma_VAE) %>%
  pivot_longer(c(sigma_KL, sigma_VAE), names_to = "elasticity", values_to = "sigma") %>%
  filter(is.finite(sigma)) %>%
  mutate(elasticity = recode(elasticity,
                             sigma_KL = "σK–L",
                             sigma_VAE = "σVA–E"),
         sigma_plot = pmin(sigma, 10)) %>%
  ggplot(aes(x = sigma_plot, y = method, fill = elasticity, colour = elasticity)) +
  ggridges::stat_density_ridges(alpha = 0.5, scale = 1,
                                rel_min_height = 0.01,
                                position = "identity",
                                kernel = "gaussian", adjust = 1.5, n = 512) +
  scale_fill_manual(values = c("σK–L" = LGREY, "σVA–E" = NAVY)) +
  scale_colour_manual(values = c("σK–L" = DGREY, "σVA–E" = DGREY)) +
  labs(title = "Substitution elasticities across methods (valid runs)",
       subtitle = "Overlayed densities by elasticity type",
       x = expression(sigma), y = "Method", caption = cap_elast,
       fill = NULL, colour = NULL) +
  theme_nat() +
  theme(legend.position = "top", legend.direction = "horizontal")

print(fig2_2)





# 2.3 Two bars per method (median σ + IQR as errorbars), sorted by σVA–E
elas_long <- results_grid_valid %>%
  select(method, sigma_KL, sigma_VAE) %>%
  pivot_longer(c(sigma_KL, sigma_VAE), names_to = "which", values_to = "sigma") %>%
  mutate(which = recode(which, "sigma_KL"="σK-L","sigma_VAE"="σVA–E"),
         method = factor(method, levels = method_order)) %>%
  filter(is.finite(sigma))
elas_stats <- elas_long %>%
  group_by(method, which) %>%
  summarise(med = median(sigma, na.rm = TRUE),
            q1  = quantile(sigma, .25, na.rm = TRUE),
            q3  = quantile(sigma, .75, na.rm = TRUE), .groups = "drop")
order_by_vae <- elas_stats %>% filter(which == "σVA–E") %>% arrange(med) %>% pull(method)
fig2_3 <- ggplot(elas_stats, aes(x = factor(method, levels = order_by_vae), y = med, fill = which)) +
  geom_col(position = position_dodge(width = .7), width = .6, colour = "white") +
  geom_errorbar(aes(ymin = q1, ymax = q3),
                position = position_dodge(width = .7), width = .2, colour = DGREY) +
  geom_text(aes(label = paste0("median ", round(med,3), "\n[", round(q1,2), ", ", round(q3,2), "]")),
            position = position_dodge(width = .7), hjust = -0.1, vjust = 0.5, size = 3.5, colour = DGREY) +
  coord_flip() +
  scale_fill_manual(values = c("σK-L" = DGREY, "σVA–E" = NAVY)) +
  labs(title = "Median elasticities by method with IQR",
       subtitle = "Bars = median, errorbars = IQR; text = values",
       x = NULL, y = expression(sigma), fill = NULL, caption = cap_elast) + theme_nat()
print(fig2_3)



# 2.4 Core parameter distributions (best methods)
best_params_long <- best_methods %>%
  select(r, gamma, lambda, nu, delta_KL, delta_VAE) %>%
  pivot_longer(-r, names_to = "param", values_to = "val") %>%
  mutate(
    param = recode(param,
                   gamma = "γ",
                   lambda = "λ",
                   nu = "ν",
                   delta_KL = "δK-L",
                   delta_VAE = "δVA–E"),
    # clamp shares to [0,1]
    val = ifelse(param %in% c("δK-L","δVA–E"), pmin(pmax(val, 0), 1), val)
  ) %>%
  filter(is.finite(val))

fig2_4 <- ggplot(best_params_long, aes(x = param, y = val)) +
  geom_violin(fill = LGREY, colour = DGREY, trim = TRUE, width = .8, alpha = .7) +
  geom_boxplot(width = .15, colour = NAVY, fill = "white", outlier.shape = NA) +
  stat_summary(fun = median, geom = "text",
               aes(label = round(..y.., 2)), vjust = -0.6, colour = NAVY, size = 2.8) +
  facet_wrap(~ param, scales = "free") +
  labs(
    title = "Core parameter distributions (best methods)",
    subtitle = "Violin + box; text = median",
    x = NULL, y = NULL, caption = cap_params
  ) +
  theme_nat()

print(fig2_4)



# 2.5 Parameter correlations (valid runs)
cor_df <- results_grid_valid %>%
  select(gamma, lambda, nu, delta_KL, delta_VAE, sigma_KL, sigma_VAE)

# drop columns with < 3 finite values or zero variance
keep <- vapply(cor_df, function(x) sum(is.finite(x)), integer(1)) >= 3 &
  vapply(cor_df, function(x) stats::sd(x, na.rm = TRUE) > 0, logical(1))
cor_df <- cor_df[, keep, drop = FALSE]

if (ncol(cor_df) >= 2) {
  cor_mat <- suppressWarnings(cor(cor_df, use = "pairwise.complete.obs"))
  cor_mat[!is.finite(cor_mat)] <- 0
  diag(cor_mat) <- 1
  
  fig2_5 <- ggcorrplot::ggcorrplot(
    cor_mat, hc.order = FALSE, type = "lower", lab = TRUE, outline.color = "white",
    colors = c("white", "#9BB3C9", NAVY)
  ) +
    labs(title = "Parameter correlations (valid runs)",
         subtitle = "Lower triangle; non-finite handled; clustering disabled",
         caption  = cap_params) +
    theme_nat()
  
  print(fig2_5)
}

# 2.6 Residual SD by region — SKIP if time-level file not available
if (file.exists("CES_results_time.csv")) {
  results_table_time <- read_csv("CES_results_time.csv", show_col_types = FALSE)
  res_by_reg <- results_table_time %>%
    group_by(r) %>% summarise(sd_res = sd(residual, na.rm = TRUE), n = sum(is.finite(residual)), .groups = "drop") %>%
    arrange(desc(sd_res))
  fig2_6 <- results_table_time %>%
    mutate(r = factor(r, levels = res_by_reg$r)) %>%
    ggplot(aes(x = r, y = residual)) +
    geom_boxplot(fill= LGREY, colour = DGREY, outlier.size = 1.2, outlier.alpha = .5) +
    geom_text(data = res_by_reg, aes(x = r, y = Inf, label = paste0("n=", n)),
              inherit.aes = FALSE, vjust = 1.3, colour = DGREY, size = 3) +
    coord_flip() +
    labs(title = "Residuals by region (valid fits)",
         subtitle = "Sorted by residual standard deviation; labels show number of time points",
         x = NULL, y = "Residual (log)", caption = cap_fit) + theme_nat()
  print(fig2_6)
} else {
  message("Skipping residual plots: CES_results_time.csv not found.")
}

# 2.7 Fit quality by method (valid runs)
rss_by_method <- results_grid_valid %>%
  group_by(method) %>%
  summarise(med_rss = median(rss, na.rm = TRUE),
            iqr_rss = IQR(rss, na.rm = TRUE), .groups = "drop")

fig2_7 <- ggplot(rss_by_method, aes(x = reorder(method, med_rss), y = med_rss)) +
  geom_col(fill = LGREY, colour = "white", width = 0.65) +
  geom_errorbar(aes(ymin = med_rss - iqr_rss/2, ymax = med_rss + iqr_rss/2),
                width = .2, colour = DGREY) +
  geom_text(aes(label = paste0("median ", round(med_rss,2),
                               " | IQR=", round(iqr_rss,2))),
            hjust = -0.05, vjust = -1, size = 3.5, colour = DGREY) +
  coord_flip() +
  labs(title = "Fit quality by method (RSS)",
       subtitle = "Median ± IQR of log(RSS) for valid runs",
       x = NULL, y = "Median log(RSS)") + theme_nat()
print(fig2_7)







### 3. BEST METHODS
# 3.1 Best-method count & share
fig3_1 <- best_methods %>% count(method) %>% mutate(share = n/sum(n)) %>% arrange(desc(share)) %>%
  ggplot(aes(x = reorder(method, share), y = share)) +
  geom_col(width = .70, fill= LGREY, colour = "white") +
  geom_text(aes(label = paste0(round(share*100,1), "% (n=", n, ")")),
            hjust = -0.05, size = 3.2, colour = DGREY) +
  coord_flip(ylim = c(0, 1.05)) +
  labs(title = "Best method by region", x = NULL, y = "Share of regions") + theme_nat()
print(fig3_1)

# 3.2 R² by region (best method), with global median
fig3_2a <- ggplot(best_methods, aes(x = reorder(r, R2), y = R2)) +
  geom_col(width = .70, fill= LGREY, colour = "white") +
  geom_hline(yintercept = median(best_methods$R2, na.rm = TRUE),
             colour = NAVY, linewidth = .4, linetype = "dashed") +
  coord_flip() +
  labs(title = expression("Fit quality by region ("*R^2*") — best methods"),
       subtitle = "Dashed = global median", x = NULL, y = expression(R^2)) +
  theme_nat()
print(fig3_2a)

fig3_2b <- ggplot(world_best) +
  geom_sf(aes(fill = R2), colour = "#f3f3f3", linewidth = .15) +
  scale_fill_gradient(limits = c(0, 1), low = "white", high = NAVY, na.value = "grey90") +
  labs(title = expression("Map of fit quality ("*R^2*") — best methods"),
       fill = expression(R^2)) +
  theme_nat()
print(fig3_2b)

# 3.3 Best-method elasticities per region — quadrant view
fig3_3 <- plot_elasticity_scatter(best_methods, med_KL, med_VAE)
print(fig3_3)

# 3.4 Combined-method regional view (cloud + IQR) with best marked
elas_reg <- results_grid_valid %>%
  select(r, method, sigma_KL, sigma_VAE) %>%
  pivot_longer(c(sigma_KL, sigma_VAE), names_to = "which", values_to = "sigma") %>%
  mutate(which = recode(which, sigma_KL = "σK-L", sigma_VAE = "σVA–E")) %>%
  filter(is.finite(sigma))
elas_iqr <- elas_reg %>%
  group_by(r, which) %>%
  summarise(q1 = quantile(sigma, 0.25, na.rm = TRUE),
            q3 = quantile(sigma, 0.75, na.rm = TRUE),
            med = median(sigma, na.rm = TRUE), .groups = "drop")
best_long <- best_methods %>%
  select(r, sigma_KL, sigma_VAE) %>%
  pivot_longer(c(sigma_KL, sigma_VAE), names_to = "which", values_to = "sigma") %>%
  mutate(which = recode(which, sigma_KL = "σK-L", sigma_VAE = "σVA–E"))
ord_reg <- elas_iqr %>% group_by(r) %>%
  summarise(mu = mean(med, na.rm = TRUE), .groups = "drop") %>% arrange(mu) %>% pull(r)
fig3_4 <- ggplot() +
  geom_linerange(
    data = elas_iqr,
    aes(y = factor(r, levels = ord_reg), xmin = q1, xmax = q3),
    linewidth = 2.2,              # <- was size
    colour = LGREY, alpha = .7
  ) +
  geom_point(data = elas_reg,
             aes(x = sigma, y = factor(r, levels = ord_reg)),
             size = 1.5, alpha = .5, colour = DGREY) +
  geom_point(data = best_long,
             aes(x = sigma, y = factor(r, levels = ord_reg)),
             shape = 8, size = 2.6, colour = NAVY) +
  facet_wrap(~ which, scales = "free_x") +
  labs(
    title = "Elasticities by region",
    subtitle = "Cloud = all methods; bar = IQR; star = best",
    x = expression(sigma), y = "Region"
  ) +
  theme_nat()
print(fig3_4)



# 3.5 Observed (y) vs Fitted (x), facets ordered by R²
if (file.exists("CES_results_time.csv")) {
  results_table_time <- read_csv("CES_results_time.csv", show_col_types = FALSE)
  obs_fit <- results_table_time %>%
    inner_join(best_methods %>% select(r, method, R2), by = c("r","method")) %>%
    mutate(Y_obs = fitted * exp(residual)) %>%
    filter(is.finite(Y_obs), is.finite(fitted))
  ord_facets <- best_methods %>% arrange(desc(R2)) %>% pull(r)
  fig3_5 <- obs_fit %>%
    mutate(r = factor(r, levels = ord_facets)) %>%
    ggplot(aes(x = fitted, y = Y_obs)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = DGREY) +
    geom_point(size = 1.6, alpha = .85, colour = DGREY) +
    facet_wrap(~ r, scales = "free") +
    labs(title = "Observed vs fitted (best methods)",
         subtitle = "Observed reconstructed as fitted × exp(residual); dashed = 1:1",
         x = "Fitted (scaled)", y = "Observed (scaled)", caption = cap_fit) + theme_nat()
  print(fig3_5)
} else {
  message("Skipping observed-vs-fitted facets: CES_results_time.csv not found.")
}

# 3.6 Best K-L elasticities map
fig3_6 <- ggplot(world_best) +
  geom_sf(aes(fill = pmin(sigma_KL,10)), colour = "white", linewidth = .1) +
  scale_fill_gradient(low="#FFC55D", high=NAVY, na.value="grey90") +
  labs(title = expression(paste(sigma[K-L], " by region")),
       fill = expression(sigma[K-L]),
       caption = "Note: elasticities capped at 10 for readability.") +
  theme_nat()
print(fig3_6)

# 3.7 Best VA-E elasticities map
fig3_7 <- ggplot(world_best) +
  geom_sf(aes(fill = pmin(sigma_VAE,10)), colour = "white", linewidth = .1) +
  scale_fill_gradient(low="#FFC55D", high=NAVY, na.value="grey90") +
  labs(title = expression(paste(sigma[VA-E], " by region")),
       fill = expression(sigma[VA-E]),
       caption = "Note: elasticities capped at 10 for readability.") +
  theme_nat()
print(fig3_7)

# 3.8 Selection certainty (Δ AICc weight best − runner-up) + distribution
per_method_best <- results_grid %>%
  filter(valid, !on_edge_KL, !on_edge_VAE) %>%                 # <- optional edges filter
  group_by(r, method) %>%
  slice_min(AICc_plusRho, with_ties = FALSE) %>%               # best (ρKL, ρVA–E) for this method
  ungroup()

aic_weights_methods <- per_method_best %>%
  group_by(r) %>%
  mutate(
    dAICc = AICc_plusRho - min(AICc_plusRho, na.rm = TRUE),
    wAICc = exp(-0.5 * dAICc),
    wAICc = wAICc / sum(wAICc, na.rm = TRUE)
  ) %>%
  ungroup()

aic_top2 <- aic_weights_methods %>%
  arrange(r, desc(wAICc)) %>%
  group_by(r) %>%
  summarise(
    method_best   = first(method),
    best          = first(wAICc),
    method_runner = if (n() >= 2) nth(method, 2) else NA_character_,
    runner        = if (n() >= 2) nth(wAICc, 2) else 0,
    delta         = best - runner,
    .groups = "drop"
  )
aic_top2_plot <- aic_top2 %>%
  mutate(
    pair = ifelse(is.na(method_runner),
                  paste(method_best, "(only)"),
                  paste(method_best, "vs", method_runner))
  )

# Order legend by frequency of pairs (most common first)
pair_order <- aic_top2_plot %>%
  count(pair, sort = TRUE) %>%
  pull(pair)

aic_top2_plot <- aic_top2_plot %>%
  mutate(
    pair   = factor(pair, levels = pair_order),
    r_plot = reorder(r, delta)
  )

# Palette for the pairs
pal_pairs <- setNames(
  scales::hue_pal()(length(pair_order)),
  pair_order
)

# Nicely formatted delta labels (outside the bars)
lab_delta <- scales::label_number(accuracy = 0.01)

fig3_8a <- ggplot(aic_top2_plot,
                  aes(x = r_plot, y = delta, fill = pair)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  geom_hline(yintercept = 0.2, colour = NAVY, linetype = "dotted", linewidth = 0.4) +
  geom_text(aes(label = lab_delta(delta)),
            hjust = -0.1, size = 3.0, colour = DGREY) +
  coord_flip(ylim = c(0, max(aic_top2_plot$delta, na.rm = TRUE) * 1.08)) +
  scale_fill_manual(values = pal_pairs, name = "Best vs runner-up") +
  labs(
    title    = "Selection certainty across regions",
    subtitle = "Δ weight = wAICc(best) − wAICc(runner-up), using each method’s best ρ-grid cell",
    x = NULL, y = "Δ AICc weight"
  ) +
  theme_nat() +
  theme(
    legend.position  = "top",
    legend.direction = "horizontal",
    legend.title     = element_text(),
    legend.box       = "horizontal",
    legend.margin    = margin(b = 4),
    legend.key.height= unit(10, "pt"),
    legend.key.width = unit(16, "pt")
  )

print(fig3_8a)

fig3_8b <- ggplot(aic_top2, aes(x = delta)) +
  geom_histogram(binwidth = .05, fill = LGREY, colour = "white") +
  labs(
    title = "Distribution of selection certainty (Δ weight)",
    x = "Δ AICc weight", y = "Count"
  ) +
  theme_nat()

print(fig3_8b)

# 3.9 CI width diagnostics for best methods (uncertainty summary)
ciw <- best_methods %>%
  transmute(r, method,
            ciw_gamma  = ci_hi_gamma  - ci_lo_gamma,
            ciw_lambda = ci_hi_lambda - ci_lo_lambda,
            ciw_dKL    = ci_hi_delta_KL  - ci_lo_delta_KL,
            ciw_dVAE   = ci_hi_delta_VAE - ci_lo_delta_VAE,
            ciw_nu     = ci_hi_nu - ci_lo_nu)
ciw_long <- ciw %>% pivot_longer(-c(r, method), names_to = "param", values_to = "ciw") %>%
  mutate(param = recode(param, ciw_gamma="γ", ciw_lambda="λ", ciw_dKL="δK-L", ciw_dVAE="δVA–E", ciw_nu="ν"))
fig3_9 <- ggplot(ciw_long, aes(x = param, y = ciw)) +
  geom_violin(fill = LGREY, colour = DGREY, trim = TRUE, width = .8, alpha = .7) +
  geom_boxplot(width = .15, colour = NAVY, fill = "white", outlier.shape = NA) +
  labs(title = "Confidence-interval widths (best methods)",
       x = NULL, y = "Width (95% approx = 1.96·SE*2)") + theme_nat()
print(fig3_9)

# 3.10 P-value distributions (best methods)
p_long2 <- best_methods %>%
  select(p_gamma, p_lambda, p_delta_KL, p_delta_VAE, p_nu) %>%
  pivot_longer(everything(), names_to = "param", values_to = "p") %>%
  mutate(param = recode(param, p_gamma="γ", p_lambda="λ", p_delta_KL="δK-L", p_delta_VAE="δVA–E", p_nu="ν")) %>%
  filter(is.finite(p), p >= 0, p <= 1)
fig3_10a <- ggplot(p_long2, aes(x = p)) +
  geom_histogram(binwidth = 0.05, fill = LGREY, colour = "white") +
  facet_wrap(~ param, ncol = 3) +
  labs(title = "P-value distributions (best methods)",
       x = "p-value", y = "Count") + theme_nat()
print(fig3_10a)

# 3.11 Parameter significance (best methods): p-value tiles
p_long <- best_methods %>%
  transmute(r,
            `γ` = p_gamma, `λ` = p_lambda, `δK-L` = p_delta_KL, `δVA–E` = p_delta_VAE, `ν` = p_nu) %>%
  pivot_longer(-r, names_to = "parameter", values_to = "p") %>%
  mutate(class = case_when(is.na(p) ~ "NA", p < 0.01 ~ "<0.01", p < 0.05 ~ "<0.05", TRUE ~ "ns"))
fig3_11 <- ggplot(p_long, aes(x = parameter, y = r, fill = class)) +
  geom_tile(width = .9, height = .9, colour = "white") +
  geom_text(aes(label = ifelse(is.na(p), "n/a", scales::pvalue(p, accuracy = .001))),
            size = 2.6, colour = DGREY) +
  scale_fill_manual(values = c("<0.01" = NAVY, "<0.05" = "#6B86A3", "ns"= LGREY, "NA" = "white")) +
  labs(title = "Parameter significance in best-method fits",
       subtitle = "Tiles show p-value classes with numeric values where available",
       x = "Parameter", y = "Region", fill = "p-value", caption = cap_params) + theme_nat()
print(fig3_11)

# 3.12 Stability of elasticities across grid (per region IQR) — SELF-CONTAINED
grid_bestmethod <- results_grid %>%
  mutate(valid = as.logical(valid)) %>%
  filter(valid) %>%
  inner_join(best_methods %>% select(r, method), by = c("r","method"))

elas_spread <- grid_bestmethod %>%
  group_by(r) %>%
  summarise(
    iqr_sigma_KL  = IQR(sigma_KL,  na.rm = TRUE),
    iqr_sigma_VAE = IQR(sigma_VAE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(-r, names_to = "which", values_to = "iqr") %>%
  mutate(which = recode(which, iqr_sigma_KL = "σK-L", iqr_sigma_VAE = "σVA–E"))

fig3_12 <- ggplot(elas_spread, aes(x = reorder(r, iqr), y = iqr, fill = which)) +
  geom_col(position = position_dodge(width = .7), width = .65, colour = "white") +
  coord_flip() +
  scale_fill_manual(values = c("σK-L" = DGREY, "σVA–E" = NAVY)) +
  labs(title = "Elasticity spread across the rho-grid (best method per region)",
       subtitle = "IQR of σ within the (ρ1,ρ) grid",
       x = NULL, y = "IQR(σ)") +
  theme_nat()
print(fig3_12)

# 3.13 Parameter estimation with uncertainty and significance
param_summary <- best_methods %>%
  select(r, method,
         gamma, se_gamma, p_gamma,
         lambda, se_lambda, p_lambda,
         nu, se_nu, p_nu,
         delta_KL, se_delta_KL, p_delta_KL,
         delta_VAE, se_delta_VAE, p_delta_VAE) %>%
  pivot_longer(-c(r, method),
               names_to = c("stat", "param"),
               names_pattern = "(se|p)?_?(gamma|lambda|nu|delta_KL|delta_VAE)",
               values_to = "val") %>%
  mutate(stat = ifelse(is.na(stat) | stat == "", "est", stat)) %>%
  pivot_wider(names_from = stat, values_from = val,
              values_fn = ~ mean(.x, na.rm = TRUE)) %>%
  filter(is.finite(est))

# Optional: pretty param symbols
param_summary <- param_summary %>%
  mutate(param = recode(param,
                        gamma     = "γ",
                        lambda    = "λ",
                        nu        = "ν",
                        delta_KL  = "δK–L",
                        delta_VAE = "δVA–E"))

# Plot
fig3_13 <- ggplot(param_summary,
                  aes(x = param, y = r, fill = est)) +
  geom_tile(colour = "white") +
  geom_errorbarh(aes(xmin = est - se, xmax = est + se),
                 colour = "black", height = 0.3, na.rm = TRUE) +
  geom_text(aes(label = paste0(round(est,2),
                               "\n(p=", scales::pvalue(p, accuracy=0.01),")")),
            size = 2.3, colour = DGREY) +
  scale_fill_gradient2(low = DGREY, mid = "white", high = NAVY, midpoint = 0) +
  labs(title = "Parameter estimates by region (best methods)",
       subtitle = "Tile = estimate (colour); errorbar = ±SE; text = value + p-value",
       x = "Parameter", y = "Region", fill = "Estimate") +
  theme_nat()
print(fig3_13)

# 3.14 Regional parameter distributions (valid runs, one graph per parameter)
# Parameters of interest
param_vars <- c("gamma","lambda","nu","delta_KL","delta_VAE","sigma_KL","sigma_VAE")

# Long format for valid runs
valid_long <- results_grid_valid %>%
  select(r, all_of(param_vars)) %>%
  pivot_longer(-r, names_to = "param", values_to = "val") %>%
  filter(is.finite(val))

# Median per region/param
median_vals <- valid_long %>%
  group_by(r, param) %>%
  summarise(median_val = median(val, na.rm = TRUE), .groups = "drop")

# Best method values
best_long <- best_methods %>%
  select(r, all_of(param_vars)) %>%
  pivot_longer(-r, names_to = "param", values_to = "best_val")

elas_long <- results_grid_valid %>%
  select(r, sigma_KL, sigma_VAE) %>%
  pivot_longer(c(sigma_KL, sigma_VAE), names_to = "param", values_to = "val") %>%
  mutate(param = recode(param, sigma_KL="σK–L", sigma_VAE="σVA–E")) %>%
  filter(is.finite(val))

best_elas <- best_methods %>%
  select(r, sigma_KL, sigma_VAE) %>%
  pivot_longer(c(sigma_KL, sigma_VAE), names_to = "param", values_to = "best_val") %>%
  mutate(param = recode(param, sigma_KL="σK–L", sigma_VAE="σVA–E"))

median_elas <- elas_long %>%
  group_by(r, param) %>%
  summarise(median_val = median(val, na.rm = TRUE), .groups = "drop")

# Pretty labels
param_labels <- c(gamma="γ", lambda="λ", nu="ν", delta_KL="δK–L", delta_VAE="δVA–E")
valid_long   <- valid_long   %>% mutate(param = recode(param, !!!param_labels))
median_vals  <- median_vals  %>% mutate(param = recode(param, !!!param_labels))
best_long    <- best_long    %>% mutate(param = recode(param, !!!param_labels))

# Function: violin plot for one parameter (enhanced visibility + legend + labels)
plot_param_group <- function(valid_df, median_df, best_df, params, title, colours) {
  v_sub <- valid_df %>% filter(param %in% params)
  m_sub <- median_df %>% filter(param %in% params)
  b_sub <- best_df   %>% filter(param %in% params)
  
  # Order regions by average best value
  ord <- b_sub %>% group_by(r) %>% summarise(mu = mean(best_val, na.rm=TRUE)) %>%
    arrange(mu) %>% pull(r)
  
  max_val <- max(c(v_sub$val, b_sub$best_val), na.rm = TRUE)
  min_val <- min(v_sub$val, na.rm = TRUE)
  x_offset <- (max_val - min_val) * 0.08
  
  ggplot(v_sub, aes(y = factor(r, levels = ord), x = val, fill = param, colour = param)) +
    geom_violin(alpha = 0.4, width = 0.8, scale = "width") +
    geom_point(data = m_sub,
               aes(x = median_val, y = factor(r, levels = ord), shape = "Median"),
               size = 2.5, inherit.aes = FALSE, colour = LGREY) +
    geom_point(data = b_sub,
               aes(x = best_val, y = factor(r, levels = ord), shape = "Best method"),
               size = 2.8, inherit.aes = FALSE, colour = NAVY) +
    geom_text(data = b_sub,
              aes(x = max_val + x_offset, y = factor(r, levels = ord),
                  label = round(best_val, 2), colour = param),
              inherit.aes = FALSE, hjust = 0, size = 3, fontface = "bold") +
    scale_fill_manual(values = colours) +
    scale_colour_manual(values = colours) +
    scale_shape_manual(values = c("Median" = 16, "Best method" = 18)) +
    labs(title = title,
         x = "Estimate", y = "Region", fill = "Parameter", shape = NULL, colour = "Parameter") +
    theme_nat() +
    theme(legend.position = "top", legend.direction = "horizontal",
          plot.margin = margin(5, 60, 5, 5)) +
    coord_cartesian(xlim = c(min_val, max_val + 3 * x_offset))
}


# Generate figures for each parameter
# γ + ν
fig3_14a <- plot_param_group(
  valid_df   = valid_long,
  median_df  = median_vals,
  best_df    = best_long,
  params     = c("γ","ν"),
  title      = "Regional distributions of γ and ν",
  colours    = c("γ" = DGREY, "ν" = NAVY)
)

# λ alone (keep violin style)
fig3_14b <- plot_param_group(
  valid_df   = valid_long,
  median_df  = median_vals,
  best_df    = best_long,
  params     = c("λ"),
  title      = "Regional distributions of λ",
  colours    = c("λ" = NAVY)
)

# δK–L + δVA–E
fig3_14c <- plot_param_group(
  valid_df   = valid_long,
  median_df  = median_vals,
  best_df    = best_long,
  params     = c("δK–L","δVA–E"),
  title      = "Regional distributions of δ parameters",
  colours    = c("δK–L" = DGREY, "δVA–E" = NAVY)
)

# σK–L + σVA–E
fig3_14d <- plot_param_group(
  valid_df   = elas_long,
  median_df  = median_elas,
  best_df    = best_elas,
  params     = c("σK–L","σVA–E"),
  title      = "Regional distributions of elasticities",
  colours    = c("σK–L" = DGREY, "σVA–E" = NAVY)
)

# Print them (can also arrange with patchwork if needed)
print(fig3_14a)
print(fig3_14b)
print(fig3_14c)
print(fig3_14d)


# 3.15 Distribution of rho among best methods
# Long best-ρ table
best_rho_long <- best_methods %>%
  dplyr::select(rho_KL, rho_VAE) %>%
  tidyr::pivot_longer(everything(), names_to = "rho_type", values_to = "rho") %>%
  dplyr::mutate(rho_type = factor(rho_type, levels = c("rho_KL","rho_VAE"),
                                  labels = c("ρKL","ρVA–E"))) %>%
  dplyr::filter(is.finite(rho))

# sample sizes & medians
rho_stats <- best_rho_long %>%
  dplyr::group_by(rho_type) %>%
  dplyr::summarise(n = sum(is.finite(rho)),
                   med = median(rho, na.rm = TRUE), .groups = "drop")

# density y at each median (for label placement)
dens_y_at <- function(xx, vec, adjust = 1.4, n = 2048){
  d <- stats::density(vec[is.finite(vec)], adjust = adjust, n = n)
  stats::approx(d$x, d$y, xout = xx, rule = 2)$y
}
rho_stats$y_med <- vapply(seq_len(nrow(rho_stats)), function(i){
  dens_y_at(rho_stats$med[i], best_rho_long$rho[best_rho_long$rho_type == rho_stats$rho_type[i]])
}, numeric(1))

# legend shows n & median (keeps canvas clean)
legend_labs <- setNames(
  paste0(levels(best_rho_long$rho_type)),
  levels(best_rho_long$rho_type))

n_regions <- dplyr::n_distinct(best_methods$r)
cap_rho_text <- paste0("Symbols: ρKL, ρVA–E = CES exponents on the K–L and VA–E nests.  n = ",
                       n_regions, " regions")

x_lim <- range(best_rho_long$rho, na.rm = TRUE)

fig3_15a <- ggplot(best_rho_long, aes(x = rho, fill = rho_type)) +
  # smooth filled densities (no outlines)
  stat_density(geom = "area", position = "identity",
               alpha = 0.35, adjust = 1.4, n = 2048, colour = NA) +
  # slim median lines
  geom_segment(data = rho_stats,
               aes(x = med, xend = med, y = 0, yend = y_med),
               inherit.aes = FALSE, linewidth = 0.8, linetype = "dashed", colour = DGREY) +
  # median labels next to the dashed lines (slight nudge up/right)
  geom_text(data = rho_stats,
            aes(x = med, y = y_med, label = paste0("median = ", round(med, 2))),
            inherit.aes = FALSE, nudge_y = max(rho_stats$y_med, na.rm = TRUE)*0.04,
            hjust = -0.05, size = 3.5, colour = DGREY) +
  scale_fill_manual(values = c("ρKL" = MGREY, "ρVA–E" = NAVY),
                    breaks = names(legend_labs),
                    labels = unname(legend_labs), name = NULL) +
  scale_x_continuous(limits = x_lim, name = expression(rho)) +
  labs(
    title    = "ρ exponents are denser around 0 for K-L and more evenly distributed for VA-E",
    subtitle = "Distributions of best-method ρ exponents",
    y        = "Kernel density (relative frequency)",
    caption  = cap_rho_text
  ) +
  theme_nat() +
  theme(
    legend.position  = "top",
    legend.direction = "horizontal",
    panel.grid.minor = element_blank()
  )

print(fig3_15a)


# -- B) 2D grid heatmap: where best-method cells fall in the full first-run grid
# Use the FIRST-RUN grid from results_grid to draw the full lattice (context),
# then fill with counts of regions whose best (ρKL, ρVA–E) fell on each cell.
# Base grid + counts (same as you already have)
# Build full first-run grid from the results you have
# ---------- B) 2D grid of best cells (counts per (rho_KL, rho_VAE)) ----------
# ----- Counts per (rho_KL, rho_VAE) cell selected as "best" -----
# --- Extremes across best methods (one value per axis) ---
ext_vals <- best_methods %>%
  dplyr::summarise(
    min_KL  = min(rho_KL,  na.rm = TRUE),
    max_KL  = max(rho_KL,  na.rm = TRUE),
    min_VAE = min(rho_VAE, na.rm = TRUE),
    max_VAE = max(rho_VAE, na.rm = TRUE)
  ) %>% as.list()

ext_cells <- grid_counts %>%
  dplyr::filter(
    n > 0 &
      (rho_KL %in% c(ext_vals$min_KL, ext_vals$max_KL) |
         rho_VAE %in% c(ext_vals$min_VAE, ext_vals$max_VAE))
  )

# σ tick hints (no sec.axis)
sigma_ticks <- c(10, 5, 2, 1, 0.5, 0.2, 0.1)
sx_df <- tibble::tibble(rho = 1/sigma_ticks - 1, lab = sigma_ticks) %>%
  dplyr::filter(rho >= min(k_vals), rho <= max(k_vals))
sy_df <- tibble::tibble(rho = 1/sigma_ticks - 1, lab = sigma_ticks) %>%
  dplyr::filter(rho >= min(v_vals), rho <= max(v_vals))

x_exp <- (max(k_vals) - min(k_vals)) * 0.06
y_exp <- (max(v_vals) - min(v_vals)) * 0.06

# helper: numeric labels in plotmath strings
lab_num <- function(x) sprintf("==%s", signif(x, 3))

fig3_15b <- ggplot() +
  geom_tile(data = grid_base,
            aes(x = rho_KL, y = rho_VAE),
            fill = "white", colour = "#EAEAEA",
            linewidth = 0.25, width = 0.98, height = 0.98) +
  geom_tile(data = subset(grid_counts, edge_KL | edge_VAE),
            aes(x = rho_KL, y = rho_VAE),
            fill = NA, colour = "#BDBDBD",
            linewidth = 0.5, width = 0.98, height = 0.98) +
  geom_tile(data = grid_counts,
            aes(x = rho_KL, y = rho_VAE, fill = n),
            colour = "white", linewidth = 0.35, width = 0.98, height = 0.98) +
  geom_text(data = subset(grid_counts, n > 0),
            aes(x = rho_KL, y = rho_VAE, label = n, colour = lab_col),
            size = 3.1, fontface = "bold", show.legend = FALSE) +
  scale_colour_identity() +
  scale_fill_gradientn(
    colours = c(BBLUE, BLUE, NAVY),
    limits  = c(0, max(grid_counts$n, na.rm = TRUE)),
    breaks  = scales::breaks_pretty(n = 5),
    name    = "Regions"
  ) +
  # highlight extreme cells actually chosen
  geom_tile(data = ext_cells,
            aes(x = rho_KL, y = rho_VAE),
            fill = NA, colour = YELLOW, linewidth = 1.0, width = 0.98, height = 0.98) +
  # global min/max guides
  geom_vline(xintercept = c(ext_vals$min_KL, ext_vals$max_KL),
             linetype = "dashed", colour = DGREY, linewidth = 0.6) +
  geom_hline(yintercept = c(ext_vals$min_VAE, ext_vals$max_VAE),
             linetype = "dashed", colour = DGREY, linewidth = 0.6) +
  # labels for guides — NOTE: use plotmath (rho), not Unicode ρ
  annotate("label", x = ext_vals$min_KL, y = max(v_vals) + y_exp*0.7,
           label = paste0("min~rho[K-L]", lab_num(ext_vals$min_KL)),
           parse = TRUE, size = 3, fill = "white", colour = DGREY) +
  annotate("label", x = ext_vals$max_KL, y = max(v_vals) + y_exp*0.7,
           label = paste0("max~rho[K-L]", lab_num(ext_vals$max_KL)),
           parse = TRUE, size = 3, fill = "white", colour = DGREY) +
  annotate("label", x = max(k_vals) + x_exp*0.5, y = ext_vals$min_VAE,
           label = paste0("min~rho[VA-E]", lab_num(ext_vals$min_VAE)),
           parse = TRUE, size = 3, fill = "white", colour = DGREY) +
  annotate("label", x = max(k_vals) + x_exp*0.5, y = ext_vals$max_VAE,
           label = paste0("max~rho[VA-E]", lab_num(ext_vals$max_VAE)),
           parse = TRUE, size = 3, fill = "white", colour = DGREY) +
  # axes (ρ only) + σ tick hints around panel
  scale_x_continuous(name = expression(rho[KL]),
                     expand = expansion(mult = c(0.00, 0.08))) +
  scale_y_continuous(name = expression(rho[VA-E]),
                     expand = expansion(mult = c(0.00, 0.10))) +
  geom_text(data = sx_df,
            aes(x = rho, y = max(v_vals) + y_exp, label = lab),
            inherit.aes = FALSE, vjust = 0, size = 3, colour = DGREY) +
  annotate("text", x = mean(range(k_vals)), y = max(v_vals) + y_exp*1.9,
           label = "sigma[K-L]", parse = TRUE, size = 3.2, colour = DGREY) +
  geom_text(data = sy_df,
            aes(x = max(k_vals) + x_exp, y = rho, label = lab),
            inherit.aes = FALSE, hjust = 0, size = 3, colour = DGREY) +
  annotate("text", x = max(k_vals) + x_exp*1.2, y = mean(range(v_vals)),
           label = "sigma[VA-E]", parse = TRUE, angle = -90, size = 3.2, colour = DGREY) +
  coord_cartesian(clip = "off") +
  labs(
    title    = "Best-method ρ cells on the first-run grid",
    subtitle = "Fill = # regions per (ρ[K–L], ρ[VA–E]); dashed lines = global min/max ρ; yellow outline = extreme cell chosen by any region.",
    caption  = "Symbols: ρ = CES exponent; σ tick hints shown outside the panel (σ = 1/(1+ρ))."
  ) +
  theme_nat() +
  theme(
    legend.position = "right",
    panel.grid      = element_blank(),
    plot.margin     = margin(10, 40, 10, 10)
  )

print(fig3_15b)




### 4. Insights from IAM export
# 4.1 TFP over time (faceted by region)
fig4_1a <- best_methods %>%
  mutate(region = r) %>%
  ggplot(aes(x = reorder(region, lambda), y = lambda)) +
  geom_col(width = .70, fill= LGREY, colour = "white") +
  coord_flip() +
  labs(title = "Estimated TFP growth rate (λ) by region — best methods",
       x = NULL, y = expression(lambda), caption = cap_params) + theme_nat()
print(fig4_1a)

fig4_1b <- iam_table %>%
  group_by(region) %>%
  mutate(TFP_norm = total_factor_productivity / first(total_factor_productivity)) %>%
  ungroup() %>%
  ggplot(aes(x = year, y = TFP_norm, group = region)) +
  geom_line(colour = DGREY) +
  facet_wrap(~ region, scales = "free_y") +
  labs(title = "IAM input: TFP trajectories normalised to base",
       x = "Year", y = "TFP / TFP[base]", caption = cap_params) + theme_nat()
print(fig4_1b)

