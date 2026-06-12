options(scipen = 999) #avoids scientific notation unless necessary

###### PACKAGES ######

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
library(forcats)



###### SETTINGS ######
setwd("C:/Users/escami_g/OneDrive - Paul Scherrer Institut/05.Models/MERGE updates/CES-parametrisation/stage2")
infile <- "MERGE macro.csv"

# Preferred: bundled objects if available. Otherwise: legacy files.
STAGE_TAG <- "stage2"  
ces_file  <- "results_run1.rds"


###### LOAD DATA ######
df <- read_csv(infile, show_col_types = TRUE)

dfS <- df %>%
  group_by(r) %>%
  mutate(
    Ybase = Y[t == 2023][1],
    Kbase = K[t == 2023][1],
    Lbase = L[t == 2023][1],
    Ebase = E[t == 2023][1]
  ) %>%
  mutate(
    Ys = Y/Ybase,
    Ks = K/Kbase,
    Ls = L/Lbase,
    Es = E/Ebase
  ) %>%
  ungroup()

# Load previous results: bundled objects or legacy files
if (file.exists(ces_file)) {
  message("Loading objects from ", ces_file)
  ces_objects <- readRDS(ces_file)
  list2env(ces_objects, .GlobalEnv)
} else {
  message("Bundled file not found, using individual files (legacy mode).")
  
  # Minimal set needed for this analysis script
  results_grid <- readRDS("results_run1.rds")
  best_methods <- read_csv("CES_best_methods.csv", show_col_types = FALSE)
  iam_table    <- read_csv("IAM_params.csv", show_col_types = FALSE)
}

###### HELPER FUNCTIONS ######
# Preset theme for graphs
theme_nat <- function(base_size = 12){
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
cap_params    <- "Symbols: γ, λ, ν = scale/growth/curvature; δK–VA, δVA–Y = share parameters."


pct <- function(x) scales::percent(x, accuracy = 1)
clip01 <- function(x) pmax(pmin(x, 1), 0)
near_   <- function(x, a, tol = 1e-12) abs(x - a) <= tol

# Numbers
# Compact labels for large numbers (k, M, B, …)
lab_si <- scales::label_number(
  accuracy  = 0.1,
  big.mark  = ",",
  scale_cut = scales::cut_short_scale()   # "", k, M, B, …
)

fmt_n <- function(x) {
  ifelse(
    x < 1000,
    scales::number(x, accuracy = 1, big.mark = ","),
    lab_si(x)
  )
}


# MAth labels
lab_sigma_KL  <- bquote(sigma[K - L])
lab_sigma_VAE <- bquote(sigma[VA - E])
lab_rho_KL    <- bquote(rho[KL])
lab_rho_VAE   <- bquote(rho[VAE])

# Adding validity to the runs, economically feasible and converged from the solver
add_validity <- function(df) {
  df %>%
    mutate(
      across(any_of(c("delta_KVA","delta_VAY","gamma","nu","lambda")),
             ~ suppressWarnings(as.numeric(.))),
      valid =
        (conv == TRUE) &
        is.finite(rss) & rss > 0 &
        is.finite(R2)  & R2 > 0 &
        between(delta_KVA, 0, 1) &
        between(delta_VAY, 0, 1) &
        between(gamma,   0.5, 3) &
        between(nu,      0.7, 1.3) &
        between(lambda, -0.05, 0.05)
      
    )
}

# Only recompute validity if it isn't already present (backwards compatible)
if (!"valid" %in% names(results_grid)) {
  results_grid <- add_validity(results_grid)
}

# Standardise method column for downstream joins/factors
results_grid <- results_grid %>%
  mutate(method = trimws(as.character(method)))

results_grid_valid   <- results_grid %>% filter(valid)
results_grid_invalid <- results_grid %>% filter(!valid)

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

# Function that adds validity to the runs, economically feasible and converged from the solver
add_validity <- function(df) {
  df %>%
    mutate(
      across(
        any_of(c("delta_KVA","delta_VAY","gamma","nu","lambda")),
        ~ suppressWarnings(as.numeric(.)) # Force to be numeric
      ),
      
      # Solver level convergence
      converged = isTRUE(conv) | (is.logical(conv) & conv),
      
      # Validity based on statistical values of the grid point model
      valid_stat = is.finite(rss) & rss > 0 &
        is.finite(R2) & R2 > 0 & R2 <= 1,
      
      # Economically plausible validity
      valid_econ = 
        is.finite(delta_KVA) & between(delta_KVA, 0, 1) &
        is.finite(delta_VAY) & between(delta_VAY, 0, 1) &
        is.finite(gamma) & between(gamma, 0.5, 3) &
        is.finite(nu) & between(nu, 0.7, 1.3) &
        is.finite(lambda) & between(lambda, -0.05, 0.05),
      
      # Validity based on the convergence, statistical value validity and economically plausible parameter validity
      valid = converged & valid_stat & valid_econ,
      
      # Strict range economically plausible validity
      valid_strict =
        valid &
        between(delta_KVA, 0.2, 0.8) &
        between(delta_VAY, 0.2, 0.8) &
        !on_edge_KL &
        !on_edge_VAE,
      
      # high-level “which filter failed”
      solver_issue = !converged,
      stat_issue = converged & !valid_stat,
      econ_issue = converged & valid_stat & !valid_econ,
      
      # Structured solver_reason from 'msg'
      solver_reason = case_when(
        converged %in% TRUE ~ "Valid",
        is.na(msg) | msg == "" ~ "Unspecified",
        grepl("false|relative|singular", msg, ignore.case = TRUE) %in% TRUE ~ "False convergence",
        grepl("max", msg, ignore.case = TRUE) %in% TRUE ~ "Max iterations",
        grepl("tol|bounds", msg, ignore.case = TRUE) %in% TRUE ~ "Bounds/tolerance",
        grepl("reduction", msg, ignore.case = TRUE) %in% TRUE ~ "Reduction criterion",
        TRUE ~ "Unspecified"
      ),
      
      # Translating to readable generic validity reasons
      valid_reason = pmap_chr( # readable column with grid-values of why something is invalid
        list(converged, rss, R2, delta_KVA, delta_VAY, gamma, nu, lambda), # list of checked parameters and statistics
        function(conv, rss, R2, delta_KVA, delta_VAY, gamma, nu, lambda) { 
          reasons <- c()
          
          if (!isTRUE(conv)) reasons <- c(reasons, "Solver did not converge")
          if (!is.finite(rss) || rss <= 0) reasons <- c(reasons, "RSS invalid")
          if (!is.finite(R2) || R2 <= 0 || R2 > 1) reasons <- c(reasons, "R2 invalid")
          if (!is.finite(delta_KVA) || delta_KVA < 0 || delta_KVA > 1) reasons <- c(reasons, "dK-VA out of [0,1]")
          if (!is.finite(delta_VAY) || delta_VAY < 0 || delta_VAY > 1) reasons <- c(reasons, "dVA-Y out of [0,1]")
          if (!is.finite(gamma) || gamma < 0.5 || gamma > 3) reasons <- c(reasons, "gamma out of [0.5,3]")
          if (!is.finite(nu) || nu < 0.7 || nu > 1.3) reasons <- c(reasons, "v out of [0.7,1.3]")
          if (!is.finite(lambda) || lambda < -0.05 || lambda > 0.05) reasons <- c(reasons, "lambda out of [-0.05,0.05]")
          
          if (length(reasons) == 0) "OK" else paste(reasons, collapse = "; ") 
        }
      ),
      
      # Status combines solver and validity info into one
      status = case_when(
        valid ~ "Valid",
        # solver failed, why?
        !converged & solver_reason != "Valid" ~ solver_reason,
        # solver ok, statistical failure, why?
        grepl("RSS invalid", valid_reason) ~ "RSS invalid",
        grepl("R2 invalid", valid_reason) ~ "R2 invalid",
        # solver & statistics ok, economic failure, why?
        grepl("dK-VA", valid_reason) ~ "dK-VA out of [0,1]",
        grepl("dVA-Y", valid_reason) ~ "dVA-Y out of [0,1]",
        grepl("gamma out", valid_reason) ~ "gamma out of [0.5,3]",
        grepl("v out", valid_reason) ~ "v out of [0.7,1.3]",
        grepl("lambda out", valid_reason) ~ "lambda out of [-0.05,0.05]",
        TRUE ~ "Unspecified"
      ),
      
      status = factor(
        status,
        levels = c(
          # solver-level reasons
          "False convergence",
          "Max iterations",
          "Bounds/tolerance",
          "Reduction criterion",
          # statistical failures
          "RSS invalid",
          "R2 invalid",
          # economic failures
          "dK-VA out of [0,1]",
          "dVA-Y out of [0,1]", 
          "gamma out of [0.5,3]",
          "v out of [0.7,1.3]",
          "lambda out of [-0.05,0.05]",
          # residual
          "Unspecified",
          "Valid"
        )
      )
      
    )
}


###### STAGE 1: Diagnostics of all runs ######
rg <- results_grid %>% filter(!is.na(method) & nzchar(method))

# Coverage per (region, method)
cov_region <- rg %>%
  group_by(r, method) %>%
  summarise(
    share_converged = mean(conv,  na.rm = TRUE),
    share_valid     = mean(valid, na.rm = TRUE),
    .groups = "drop"
  )

# Stage-1 method ordering: based on all runs
method_order <- cov_region %>%
  group_by(method) %>%
  summarise(med_valid = median(share_valid, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med_valid)) %>%
  pull(method)

method_order_use <- method_order

# Fig S1-1: coverage by method (median + IQR across regions)
rg <- results_grid %>% filter(!is.na(method) & nzchar(method))

cov_stats <- cov_region %>%
  pivot_longer(c(share_converged, share_valid),
               names_to = "which", values_to = "share") %>%
  mutate(
    which  = factor(which,
                    levels = c("share_converged", "share_valid"),
                    labels = c("Converged", "Valid")),
    method = factor(method, levels = method_order_use)
  ) %>%
  group_by(method, which) %>%
  summarise(
    med = median(share, na.rm = TRUE),
    q1  = quantile(share, 0.25, na.rm = TRUE),
    q3  = quantile(share, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

fig_s1_1 <- ggplot(cov_stats, aes(x = method, y = med, colour = which, shape = which)) +
  geom_errorbar(
    aes(ymin = q1, ymax = q3),
    width = .5, position = position_dodge(width = .55), colour = DGREY
  ) +
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
print(fig_s1_1)


# Fig S1-2: region × method convergence heatmap (share solved)
rm_share <- results_grid %>%
  group_by(r, method) %>%
  summarise(p_grid_solved = mean(is.finite(rss)), .groups = "drop") %>%
  complete(r, method, fill = list(p_grid_solved = 0)) %>%
  mutate(method = factor(method, levels = method_order_use))

ord_r <- rm_share %>%
  group_by(r) %>%
  summarise(mu = mean(p_grid_solved), .groups = "drop") %>%
  arrange(desc(mu), r) %>%
  pull(r)

rm_share <- rm_share %>% mutate(r = factor(r, levels = ord_r))

fig_s1_2 <- ggplot(rm_share, aes(x = method, y = r, fill = p_grid_solved)) +
  geom_tile(width = 0.98, height = 0.98, colour = "white", linewidth = 0.25) +
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
    subtitle = "Fill = share of ρ-grid points that produced an RSS",
    x = "Method", y = "Region", fill = "Solved share", caption = cap_rho
  ) +
  theme_nat() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(fig_s1_2)


# Fig S1-3: rho medians + IQR by method (valid runs) with sigma axis
rho_long_s1 <- results_grid %>%
  filter(valid) %>%
  select(method, rho_KL, rho_VAE) %>%
  pivot_longer(c(rho_KL, rho_VAE), names_to = "rho_type", values_to = "rho") %>%
  mutate(
    rho_type = factor(rho_type, levels = c("rho_KL", "rho_VAE"), labels = c("ρKL", "ρVA–E")),
    method   = factor(method, levels = method_order_use)
  ) %>%
  filter(is.finite(rho))

rho_stats_s1 <- rho_long_s1 %>%
  group_by(method, rho_type) %>%
  summarise(
    med = median(rho, na.rm = TRUE),
    q1  = quantile(rho, 0.25, na.rm = TRUE),
    q3  = quantile(rho, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

y_min <- min(rho_stats_s1$q1, rho_stats_s1$med, na.rm = TRUE)
y_max <- max(rho_stats_s1$q3, rho_stats_s1$med, na.rm = TRUE)
pad   <- (y_max - y_min) * 0.08
y_lim <- c(y_min - pad, y_max + pad)

fig_s1_3 <- ggplot(rho_stats_s1, aes(x = method, y = med, colour = rho_type)) +
  geom_errorbar(
    aes(ymin = q1, ymax = q3),
    width = .25, linewidth = .6, colour = DGREY,
    position = position_dodge(width = .55)
  ) +
  geom_point(position = position_dodge(width = .55), size = 2.9) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = LGREY) +
  coord_flip() +
  scale_colour_manual(values = c("ρKL" = DGREY, "ρVA–E" = NAVY), name = NULL) +
  scale_y_continuous(
    limits = y_lim,
    name   = expression(rho),
    sec.axis = sec_axis(~ 1 / (1 + .), name = expression(sigma))
  ) +
  labs(
    title    = "Best-fit ρ medians with IQR by method",
    subtitle = "Points = medians; bars = IQR (valid runs only). Right axis shows σ = 1/(1+ρ).",
    x = NULL, caption = cap_rho
  ) +
  theme_nat()
print(fig_s1_3)

# Fig S1-4: runtime by method
rt_stats <- results_grid %>%
  group_by(method) %>%
  summarise(med_min = median(runtime_total, na.rm = TRUE) / 60, .groups = "drop") %>%
  arrange(med_min)

fig_s1_4 <- results_grid %>%
  mutate(method = factor(method, levels = rt_stats$method)) %>%
  ggplot(aes(x = method, y = runtime_total / 60)) +
  geom_boxplot(width = .72, fill = LGREY, colour = DGREY,
               outlier.size = 1.4, outlier.alpha = .6) +
  stat_summary(fun = median, geom = "point", size = 3, shape = 21,
               fill = NAVY, colour = NAVY) +
  coord_flip() +
  labs(
    title = "Runtime by method",
    subtitle = "Minutes; navy = method median",
    x = NULL, y = "Runtime (minutes)", caption = cap_none
  ) +
  theme_nat()
print(fig_s1_4)

# Fig S1-5: runtime vs fit quality (converged runs)
rt_quality <- results_grid %>%
  filter(conv, is.finite(R2), R2 > 0) %>%
  mutate(method = factor(method, levels = method_order_use))

fig_s1_5 <- ggplot(rt_quality, aes(x = runtime_total / 60, y = R2, colour = method)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title    = "Runtime vs fit quality",
    subtitle = "Each point = converged run; x = runtime (minutes), y = R²",
    x = "Runtime (minutes)", y = expression(R^2)
  ) +
  theme_nat() +
  theme(legend.position = "top")
print(fig_s1_5)


# Fig S1-6: RSS vs iterations (capped)
fit_df <- results_grid %>%
  filter(
    conv,
    is.finite(iter), iter > 0,
    is.finite(rss),  rss > 0,
    !is.na(method), nzchar(method)
  )

method_levels <- unique(c(method_order_use, sort(unique(fit_df$method))))
fit_df <- fit_df %>% mutate(method_f = factor(method, levels = method_levels))

pal_methods <- scales::hue_pal()(length(method_levels))
names(pal_methods) <- method_levels

x_cap <- stats::quantile(fit_df$iter, 0.99, na.rm = TRUE)
y_cap <- stats::quantile(fit_df$rss,  0.99, na.rm = TRUE)

fig_s1_6 <- ggplot(fit_df, aes(x = iter, y = rss, colour = method_f)) +
  geom_point(alpha = 0.7, size = 1.9) +
  scale_colour_manual(values = pal_methods, drop = FALSE, name = "") +
  coord_cartesian(xlim = c(0, x_cap), ylim = c(0, y_cap)) +
  labs(
    title    = "Iterations do not necessarily yield better fits",
    subtitle = "Converged runs: RSS vs solver iterations (axes capped at 99th percentile)",
    x = "Iterations", y = "RSS", caption = cap_fit
  ) +
  theme_nat() +
  theme(legend.position = "top")
print(fig_s1_6)


# Fig S1-7: edge minima share by method (best-method runs only)
grid_bestmethod <- results_grid %>%
  inner_join(best_methods %>% select(r, method), by = c("r", "method"))

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

fig_s1_7 <- ggplot(edge_counts, aes(x = reorder(method, share_edge), y = share_edge)) +
  geom_col(fill = LGREY, colour = "white", width = 0.65) +
  geom_text(
    aes(label = paste0(round(share_edge * 100, 1), "%")),
    hjust = -0.05, colour = DGREY, size = 3.1
  ) +
  coord_flip(ylim = c(0, 1)) +
  labs(
    title    = "Edge solutions by method",
    subtitle = "Share of ΔAICc minima at grid boundaries",
    x = NULL, y = "Share of edge minima"
  ) +
  theme_nat()
print(fig_s1_7)


# Fig S1-8: run status by method (non-convergence + econ bound failures + valid)
runs_status <- results_grid %>%
  mutate(
    status = solver_reason  # use the precomputed solver categories
  ) %>%
  mutate(
    status = factor(
      status,
      levels = c(
        "False convergence",
        "Max iterations",
        "Bounds/tolerance",
        "Reduction criterion",
        "Unspecified",
        "Valid"
      )
    )
  )


error_summary <- runs_status %>%
  count(method, status, name = "n") %>%
  complete(method, status, fill = list(n = 0)) %>%
  group_by(method) %>%
  mutate(total = sum(n)) %>%
  ungroup()

method_order_err <- error_summary %>%
  group_by(method) %>%
  summarise(nonvalid = sum(n[status != "Valid"]), .groups = "drop") %>%
  arrange(desc(nonvalid)) %>%
  pull(method)

error_summary <- error_summary %>%
  mutate(method = factor(method, levels = method_order_err))

lab_short <- scales::label_number(
  accuracy  = 0.1,
  big.mark  = ",",
  trim      = TRUE,
  scale_cut = scales::cut_short_scale()
)

fig_s1_8 <- ggplot(error_summary, aes(x = method, y = n, fill = status)) +
  geom_col(width = 0.9, colour = "white") +
  geom_text(
    aes(label = ifelse(n > 0, lab_short(n), "")),
    position = position_stack(vjust = 0.5),
    size = 2.8, colour = "black", fontface = "bold"
  ) +
  coord_flip() +
  scale_y_continuous(labels = lab_short, expand = expansion(mult = c(0, .05))) +
  scale_fill_manual(
    breaks = levels(error_summary$status),
    values = c(
      "False convergence"   = "#D5D8FF",
      "Max iterations"      = "#CDFFF5",
      "Reduction criterion" = "#FFDFC5",
      "Bounds/tolerance"    = "#B7E5FF",
      "Unspecified"         = "#D9D9D9",
      "δK–VA out of [0,1]"  = "#B3B3B3",
      "δVA–Y out of [0,1]"  = "#A0A0A0",
      "γ out of (0.2,5)"    = "#8D8D8D",
      "ν out of (0.2,5)"    = "#7A7A7A",
      "λ not finite"        = "#676767",
      "other econ bound"    = "#5A5A5A",
      "Valid"               = "#00B0F0"
    ),
    drop = FALSE
  ) +
  labs(
    title    = "Run status by method",
    subtitle = "Valid runs vs solver non-convergence and economic-bounds failures",
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
print(fig_s1_8)



###### STAGE 2: VALID RUNS (CONVERGED + ECON-VALID) ######

# Fig S2-1: method summary (n, median R², median runtime)
fig_s2_1 <- results_grid_valid %>%
  group_by(method) %>%
  summarise(
    n           = n(),
    med_R2      = median(R2, na.rm = TRUE),
    med_runtime = median(runtime_total, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(med_R2)) %>%
  mutate(method = factor(method, levels = method)) %>%
  ggplot(aes(x = method, y = n)) +
  geom_col(width = .70, fill = LGREY, colour = "white") +
  geom_text(
    aes(label = paste0("R²~", round(med_R2, 2), " | t~", round(med_runtime / 60, 1), " min")),
    vjust = -0.25, size = 3.1, colour = DGREY
  ) +
  coord_flip() +
  labs(
    title    = "Valid runs by method",
    subtitle = "Bar = count; label = median R² and runtime (min)",
    x = NULL, y = "Valid runs", caption = cap_fit
  ) +
  theme_nat()
print(fig_s2_1)


# Fig S2-2: valid share per method (conv + econ-valid)
valid_share <- results_grid %>%
  group_by(method) %>%
  summarise(
    valid_frac = mean(valid, na.rm = TRUE),
    n          = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(valid_frac))

fig_s2_2 <- ggplot(valid_share, aes(x = reorder(method, valid_frac), y = valid_frac)) +
  geom_col(width = .7, fill = LGREY, colour = "white") +
  geom_text(
    aes(label = paste0(round(valid_frac * 100, 1), "% (n=", n, ")")),
    hjust = -0.05, colour = DGREY, size = 3.2
  ) +
  coord_flip(ylim = c(0, 1.05)) +
  labs(
    title    = "Economically valid share by method",
    subtitle = "Share of runs that converged and passed economic bounds",
    x = NULL, y = "Share valid", caption = cap_none
  ) +
  theme_nat()
print(fig_s2_2)


# Fig S2-3: rho densities by method (valid runs)
rho_medians_s2 <- results_grid_valid %>%
  select(method, rho_KL, rho_VAE) %>%
  pivot_longer(c(rho_KL, rho_VAE), names_to = "rho_type", values_to = "rho") %>%
  filter(is.finite(rho)) %>%
  mutate(rho_type = recode(rho_type, rho_KL = "ρK–L", rho_VAE = "ρVA–E")) %>%
  group_by(method, rho_type) %>%
  summarise(rho_med = median(rho, na.rm = TRUE), .groups = "drop")

fig_s2_3 <- results_grid_valid %>%
  select(method, rho_KL, rho_VAE) %>%
  pivot_longer(c(rho_KL, rho_VAE), names_to = "rho_type", values_to = "rho") %>%
  filter(is.finite(rho)) %>%
  mutate(rho_type = recode(rho_type, rho_KL = "ρK–L", rho_VAE = "ρVA–E")) %>%
  ggplot(aes(x = rho, y = method, fill = rho_type, colour = rho_type)) +
  geom_point(
    data = rho_medians_s2,
    aes(x = rho_med, y = method, colour = rho_type),
    shape = "|", size = 8, stroke = 2,
    position = position_nudge(y = 0.1)
  ) +
  stat_density_ridges(
    alpha = 0.5, scale = 1, rel_min_height = 0.01,
    position = "identity", kernel = "gaussian", adjust = 1.5, n = 512
  ) +
  scale_fill_manual(values = c("ρK–L" = LGREY, "ρVA–E" = NAVY)) +
  scale_colour_manual(values = c("ρK–L" = DGREY, "ρVA–E" = DGREY)) +
  labs(
    title    = "ρ-grid posterior density by method (valid runs)",
    subtitle = "Ridge densities in ρ-space, after economic/fit filters",
    x = expression(rho), y = "Method", fill = NULL, colour = NULL
  ) +
  theme_nat() +
  theme(legend.position = "top", legend.direction = "horizontal")
print(fig_s2_3)


# Fig S2-4: rho medians + IQR by method
rho_stats_s2 <- rho_long_s2 %>%
  group_by(method, rho_type) %>%
  summarise(
    med = median(rho, na.rm = TRUE),
    q1  = quantile(rho, 0.25, na.rm = TRUE),
    q3  = quantile(rho, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# n(valid) per method (used for ordering + annotation)
n_method <- results_grid_valid %>%
  count(method, name = "n_valid")

rho_stats_plot <- rho_stats_s2 %>%
  left_join(n_method, by = "method") %>%
  mutate(method = fct_reorder(method, n_valid, .desc = TRUE))

y_min <- min(rho_stats_plot$q1, rho_stats_plot$med, na.rm = TRUE)
y_max <- max(rho_stats_plot$q3, rho_stats_plot$med, na.rm = TRUE)
pad   <- (y_max - y_min) * 0.08
y_lim <- c(y_min - pad, y_max + pad)

fig_s2_4 <- ggplot(rho_stats_plot, aes(x = method, y = med, colour = rho_type)) +
  geom_errorbar(
    aes(ymin = q1, ymax = q3),
    width = 0.25, linewidth = 1.2, colour = DGREY,
    position = position_dodge(width = 0.55)
  ) +
  geom_point(
    position = position_dodge(width = 0.55),
    shape = "|", size = 7, stroke = 1.4
  ) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = LGREY) +
  geom_text(
    data = n_method %>% mutate(method = fct_reorder(method, n_valid, .desc = TRUE)),
    aes(x = method, y = y_lim[2], label = paste0("n = ", fmt_n(n_valid))),
    inherit.aes = FALSE,
    hjust = -0.05, size = 3.2, colour = DGREY
  ) +
  coord_flip() +
  scale_colour_manual(values = c("ρKL" = DGREY, "ρVA–E" = NAVY), name = NULL) +
  scale_y_continuous(
    limits = y_lim,
    name = expression(rho),
    sec.axis = sec_axis(~ 1 / (1 + .), name = expression(sigma))
  ) +
  labs(
    title    = "Best-fit ρ medians with interquartile ranges by method",
    subtitle = "Bars = 25–75% range; ticks = medians; annotation = number of valid runs",
    x = NULL, caption = cap_rho
  ) +
  theme_nat()
print(fig_s2_4)



# Fig S2-5: elasticity densities by method (valid runs)
fig_s2_5 <- results_grid_valid %>%
  select(method, sigma_KL, sigma_VAE) %>%
  pivot_longer(c(sigma_KL, sigma_VAE), names_to = "elasticity", values_to = "sigma") %>%
  filter(is.finite(sigma)) %>%
  mutate(
    elasticity = recode(elasticity, sigma_KL = "σK–L", sigma_VAE = "σVA–E"),
    sigma_plot = pmin(sigma, 10)
  ) %>%
  ggplot(aes(x = sigma_plot, y = method, fill = elasticity, colour = elasticity)) +
  stat_density_ridges(
    alpha = 0.5, scale = 1, rel_min_height = 0.01,
    position = "identity", kernel = "gaussian", adjust = 1.5, n = 512
  ) +
  scale_fill_manual(values = c("σK–L" = LGREY, "σVA–E" = NAVY)) +
  scale_colour_manual(values = c("σK–L" = DGREY, "σVA–E" = DGREY)) +
  labs(
    title    = "Substitution elasticities across methods (valid runs)",
    subtitle = "Densities by elasticity type; elasticities capped at 10 for readability",
    x = expression(sigma), y = "Method", caption = cap_elast,
    fill = NULL, colour = NULL
  ) +
  theme_nat() +
  theme(legend.position = "top", legend.direction = "horizontal")
print(fig_s2_5)


###### STAGE 3: BEST METHODS AND UNCERTAINTIES ######
# Load mapping for map plots
region_map <- read_excel("MERGE regions proposal.xlsx", sheet = 1) %>%
  rename(r = MERGE) %>%
  mutate(r = as.character(r))

world <- ne_countries(scale = "medium", returnclass = "sf")
best_map <- best_methods %>% inner_join(region_map, by = "r")
world_best <- world %>% left_join(best_map, by = c("iso_a3" = "ISO3"))


# Fig S3-1: best-method share
fig_s3_1 <- best_methods %>%
  count(method) %>%
  mutate(share = n / sum(n)) %>%
  arrange(desc(share)) %>%
  ggplot(aes(x = reorder(method, share), y = share)) +
  geom_col(width = .70, fill = LGREY, colour = "white") +
  geom_text(
    aes(label = paste0(round(share * 100, 1), "% (n=", n, ")")),
    hjust = -0.05, size = 3.2, colour = DGREY
  ) +
  coord_flip(ylim = c(0, 1.05)) +
  labs(title = "Best method by region", x = NULL, y = "Share of regions") +
  theme_nat()
print(fig_s3_1)


# Fig S3-2a: R² by region 
fig_s3_2a <- ggplot(best_methods, aes(x = reorder(r, R2), y = R2)) +
  geom_col(width = .70, fill = LGREY, colour = "white") +
  geom_hline(
    yintercept = median(best_methods$R2, na.rm = TRUE),
    colour = NAVY, linewidth = .4, linetype = "dashed"
  ) +
  coord_flip() +
  labs(
    title    = expression("Fit quality by region ("*R^2*") — best methods"),
    subtitle = "Dashed = global median",
    x = NULL, y = expression(R^2)
  ) +
  theme_nat()
print(fig_s3_2a)

# Fig S3-2b: R² map
fig_s3_2b <- ggplot(world_best) +
  geom_sf(aes(fill = R2), colour = "#f3f3f3", linewidth = .15) +
  scale_fill_gradient(limits = c(0, 1), low = "white", high = NAVY, na.value = "grey90") +
  labs(title = expression("Map of fit quality ("*R^2*") — best methods"), fill = expression(R^2)) +
  theme_nat()
print(fig_s3_2b)


# Fig S3-3: best-method elasticities per region (quadrant)
fig_s3_3 <- plot_elasticity_scatter(best_methods, med_KL, med_VAE)
print(fig_s3_3)


# Fig S3-4: Elasticity dispersion across methods per region (cloud + IQR)
elas_reg <- results_grid_valid %>%
  select(r, method, sigma_KL, sigma_VAE) %>%
  pivot_longer(c(sigma_KL, sigma_VAE),
               names_to = "which", values_to = "sigma") %>%
  mutate(which = recode(
    which,
    sigma_KL  = "σK–L",
    sigma_VAE = "σVA–E"
  )) %>%
  filter(is.finite(sigma))

elas_iqr <- elas_reg %>%
  group_by(r, which) %>%
  summarise(
    q1  = quantile(sigma, 0.25, na.rm = TRUE),
    q3  = quantile(sigma, 0.75, na.rm = TRUE),
    med = median(sigma, na.rm = TRUE),
    .groups = "drop"
  )

best_long <- best_methods %>%
  select(r, sigma_KL, sigma_VAE) %>%
  pivot_longer(c(sigma_KL, sigma_VAE),
               names_to = "which", values_to = "sigma") %>%
  mutate(which = recode(
    which,
    sigma_KL  = "σK–L",
    sigma_VAE = "σVA–E"
  ))

ord_reg <- elas_iqr %>%
  group_by(r) %>%
  summarise(mu = mean(med, na.rm = TRUE), .groups = "drop") %>%
  arrange(mu) %>% pull(r)

fig_s3_4 <- ggplot() +
  geom_linerange(
    data = elas_iqr,
    aes(y = factor(r, levels = ord_reg), xmin = q1, xmax = q3),
    linewidth = 2.2, colour = LGREY, alpha = 0.8
  ) +
  geom_point(
    data = elas_reg,
    aes(x = sigma, y = factor(r, levels = ord_reg)),
    size = 1.5, alpha = 0.5, colour = DGREY
  ) +
  geom_point(
    data = best_long,
    aes(x = sigma, y = factor(r, levels = ord_reg)),
    shape = 8, size = 2.6, colour = NAVY
  ) +
  facet_wrap(~ which, scales = "free_x") +
  labs(
    title    = "Elasticities by region",
    subtitle = "Cloud = all valid methods; bar = IQR; star = best method",
    x = expression(sigma), y = "Region"
  ) +
  theme_nat()

print(fig_s3_4)


# Fig S3-5: Observed vs fitted output (best methods)
if (file.exists("CES_results_time.csv")) {
  
  results_time <- read_csv("CES_results_time.csv", show_col_types = FALSE)
  
  obs_fit <- results_time %>%
    inner_join(
      best_methods %>% select(r, method, R2),
      by = c("r", "method")
    ) %>%
    mutate(
      Y_obs = fitted * exp(residual)
    ) %>%
    filter(is.finite(Y_obs), is.finite(fitted))
  
  region_order_fit <- best_methods %>%
    arrange(desc(R2)) %>%
    pull(r)
  
  fig_s3_5 <- obs_fit %>%
    mutate(r = factor(r, levels = region_order_fit)) %>%
    ggplot(aes(x = fitted, y = Y_obs)) +
    geom_abline(
      slope     = 1,
      intercept = 0,
      linetype  = "dashed",
      colour    = DGREY
    ) +
    geom_point(
      size  = 1.6,
      alpha = 0.85,
      colour = DGREY
    ) +
    facet_wrap(~ r, scales = "free") +
    labs(
      title    = "Observed vs fitted output (best methods)",
      subtitle = "Observed reconstructed as fitted × exp(residual)",
      x = "Fitted output (scaled)",
      y = "Observed output (scaled)",
      caption = cap_fit
    ) +
    theme_nat()
  
  print(fig_s3_5)
  
} else {
  message("Skipping Fig S3-5: CES_results_time.csv not found.")
}


# Fig S3-6: Spatial distribution of best K–L elasticities
fig_s3_6 <- ggplot(world_best) +
  geom_sf(
    aes(fill = pmin(sigma_KL, 10)),
    colour   = "white",
    linewidth = 0.1
  ) +
  scale_fill_gradient(
    low      = "#FFC55D",
    high     = NAVY,
    na.value = "grey90"
  ) +
  labs(
    title   = expression(paste(sigma[K-L], " by region (best methods)")),
    fill    = expression(sigma[K-L]),
    caption = "Elasticities capped at 10 for readability."
  ) +
  theme_nat()

print(fig_s3_6)


# Fig S3-7: Spatial distribution of best VA–E elasticities
fig_s3_7 <- ggplot(world_best) +
  geom_sf(
    aes(fill = pmin(sigma_VAE, 10)),
    colour   = "white",
    linewidth = 0.1
  ) +
  scale_fill_gradient(
    low      = "#FFC55D",
    high     = NAVY,
    na.value = "grey90"
  ) +
  labs(
    title   = expression(paste(sigma[VA-E], " by region (best methods)")),
    fill    = expression(sigma[VA-E]),
    caption = "Elasticities capped at 10 for readability."
  ) +
  theme_nat()

print(fig_s3_7)


# Fig S3-8 Selection certainty: Δ AICc weight (best − runner-up)

per_method_best <- results_grid %>%
  filter(valid, !on_edge_KL, !on_edge_VAE) %>%
  group_by(r, method) %>%
  slice_min(AICc_plusRho, with_ties = FALSE) %>%
  ungroup()

aic_weights <- per_method_best %>%
  group_by(r) %>%
  mutate(
    dAICc = AICc_plusRho - min(AICc_plusRho),
    wAICc = exp(-0.5 * dAICc),
    wAICc = wAICc / sum(wAICc)
  ) %>%
  ungroup()

aic_top2 <- aic_weights %>%
  arrange(r, desc(wAICc)) %>%
  group_by(r) %>%
  summarise(
    best   = first(wAICc),
    runner = if (n() >= 2) nth(wAICc, 2) else 0,
    delta  = best - runner,
    .groups = "drop"
  )

fig_s3_8 <- ggplot(aic_top2, aes(x = reorder(r, delta), y = delta)) +
  geom_col(width = 0.7, fill = LGREY, colour = "white") +
  geom_hline(yintercept = 0.2, linetype = "dotted",
             colour = NAVY, linewidth = 0.5) +
  coord_flip() +
  labs(
    title    = "Selection certainty across regions",
    subtitle = "Δ weight = wAICc(best) − wAICc(runner-up)",
    x = NULL, y = "Δ AICc weight"
  ) +
  theme_nat()

print(fig_s3_8)


# Fig S3-9 CI width diagnostics (best methods)
ciw <- best_methods %>%
  transmute(
    r,
    γ  = ci_hi_gamma  - ci_lo_gamma,
    λ  = ci_hi_lambda - ci_lo_lambda,
    ν  = ci_hi_nu     - ci_lo_nu,
    `δK–VA` = ci_hi_delta_KVA - ci_lo_delta_KVA,
    `δVA–Y` = ci_hi_delta_VAY - ci_lo_delta_VAY
  ) %>%
  pivot_longer(-r, names_to = "param", values_to = "width")

fig_s3_9 <- ggplot(ciw, aes(x = param, y = width)) +
  geom_violin(fill = LGREY, colour = DGREY, alpha = 0.7) +
  geom_boxplot(width = 0.15, fill = "white",
               colour = NAVY, outlier.shape = NA) +
  labs(
    title = "Confidence-interval widths (best methods)",
    x = NULL, y = "Width"
  ) +
  theme_nat()

print(fig_s3_9)



# Fig S3-10 P-value distributions (best methods)
p_long2 <- best_methods %>%
  select(p_gamma, p_lambda, p_delta_KVA, p_delta_VAY, p_nu) %>%
  pivot_longer(everything(), names_to = "param", values_to = "p") %>%
  mutate(param = recode(
    param,
    p_gamma     = "γ",
    p_lambda    = "λ",
    p_delta_KVA = "δK–VA",
    p_delta_VAY = "δVA–Y",
    p_nu        = "ν"
  )) %>%
  filter(is.finite(p), p >= 0, p <= 1)

fig_s3_10 <- ggplot(p_long2, aes(x = p)) +
  geom_histogram(binwidth = 0.05, fill = LGREY, colour = "white") +
  facet_wrap(~ param, ncol = 3) +
  labs(title = "P-value distributions (best methods)", x = "p-value", y = "Count") +
  theme_nat()
print(fig_s3_10)


# S3-11: Parameter significance tiles
p_long <- best_methods %>%
  transmute(
    r,
    `γ`     = p_gamma,
    `λ`     = p_lambda,
    `δK–VA` = p_delta_KVA,
    `δVA–Y` = p_delta_VAY,
    `ν`     = p_nu
  ) %>%
  pivot_longer(-r, names_to = "parameter", values_to = "p") %>%
  mutate(class = case_when(
    is.na(p) ~ "NA",
    p < 0.01 ~ "<0.01",
    p < 0.05 ~ "<0.05",
    TRUE     ~ "ns"
  ))

fig_s3_11 <- ggplot(p_long, aes(x = parameter, y = r, fill = class)) +
  geom_tile(width = 0.9, height = 0.9, colour = "white") +
  geom_text(
    aes(label = ifelse(is.na(p), "n/a", scales::pvalue(p, accuracy = 0.001))),
    size = 2.6, colour = DGREY
  ) +
  scale_fill_manual(values = c("<0.01" = NAVY, "<0.05" = "#6B86A3", "ns" = LGREY, "NA" = "white")) +
  labs(
    title    = "Parameter significance in best-method fits",
    subtitle = "Tiles show p-value classes; numbers show p-values where available",
    x = "Parameter", y = "Region", fill = "p-value", caption = cap_params
  ) +
  theme_nat()
print(fig_s3_11)


# S3-12: Stability of elasticities across grid (IQR per region)
grid_bestmethod <- results_grid %>%
  mutate(valid = as.logical(valid)) %>%
  filter(valid) %>%
  inner_join(best_methods %>% select(r, method), by = c("r", "method"))

elas_spread <- grid_bestmethod %>%
  group_by(r) %>%
  summarise(
    iqr_sigma_KL  = IQR(sigma_KL,  na.rm = TRUE),
    iqr_sigma_VAE = IQR(sigma_VAE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(-r, names_to = "which", values_to = "iqr") %>%
  mutate(which = recode(
    which,
    iqr_sigma_KL  = "σK–L",
    iqr_sigma_VAE = "σVA–E"
  ))

fig_s3_12 <- ggplot(elas_spread, aes(x = reorder(r, iqr), y = iqr, fill = which)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65, colour = "white") +
  coord_flip() +
  scale_fill_manual(values = c("σK–L" = DGREY, "σVA–E" = NAVY)) +
  labs(
    title    = "Elasticity spread across the rho-grid (best method per region)",
    subtitle = "IQR of σ within the (ρKL, ρVA–E) grid for the selected method",
    x = NULL, y = "IQR(σ)", fill = NULL
  ) +
  theme_nat()
print(fig_s3_12)


# S3-13: Parameter estimates by region (tiles + SE + p)
param_summary <- best_methods %>%
  select(
    r, method,
    gamma,      se_gamma,      p_gamma,
    lambda,     se_lambda,     p_lambda,
    nu,         se_nu,         p_nu,
    delta_KVA,  se_delta_KVA,  p_delta_KVA,
    delta_VAY,  se_delta_VAY,  p_delta_VAY
  ) %>%
  pivot_longer(
    -c(r, method),
    names_to = c("stat", "param"),
    names_pattern = "(se|p)?_?(gamma|lambda|nu|delta_KVA|delta_VAY)",
    values_to = "val"
  ) %>%
  mutate(stat = ifelse(is.na(stat) | stat == "", "est", stat)) %>%
  pivot_wider(names_from = stat, values_from = val, values_fn = ~ mean(.x, na.rm = TRUE)) %>%
  filter(is.finite(est)) %>%
  mutate(param = recode(
    param,
    gamma     = "γ",
    lambda    = "λ",
    nu        = "ν",
    delta_KVA = "δK–VA",
    delta_VAY = "δVA–Y"
  ))

fig_s3_13 <- ggplot(param_summary, aes(x = param, y = r, fill = est)) +
  geom_tile(colour = "white") +
  geom_errorbarh(
    aes(xmin = est - se, xmax = est + se),
    colour = "black", height = 0.3, na.rm = TRUE
  ) +
  geom_text(
    aes(label = paste0(round(est, 2), "\n(p=", scales::pvalue(p, accuracy = 0.01), ")")),
    size = 2.3, colour = DGREY
  ) +
  scale_fill_gradient2(low = DGREY, mid = "white", high = NAVY, midpoint = 0) +
  labs(
    title    = "Parameter estimates by region (best methods)",
    subtitle = "Tile = estimate; errorbar = ±SE; text = estimate + p-value",
    x = "Parameter", y = "Region", fill = "Estimate"
  ) +
  theme_nat()
print(fig_s3_13)


# S3-14: Regional distributions (valid runs) + medians + best
param_vars <- c("gamma","lambda","nu","delta_KVA","delta_VAY")

valid_long <- results_grid_valid %>%
  select(r, all_of(param_vars)) %>%
  pivot_longer(-r, names_to = "param", values_to = "val") %>%
  filter(is.finite(val))

median_vals <- valid_long %>%
  group_by(r, param) %>%
  summarise(median_val = median(val, na.rm = TRUE), .groups = "drop")

best_long_params <- best_methods %>%
  select(r, all_of(param_vars)) %>%
  pivot_longer(-r, names_to = "param", values_to = "best_val")

param_labels <- c(
  gamma     = "γ",
  lambda    = "λ",
  nu        = "ν",
  delta_KVA = "δK–VA",
  delta_VAY = "δVA–Y"
)

valid_long       <- valid_long       %>% mutate(param = recode(param, !!!param_labels))
median_vals      <- median_vals      %>% mutate(param = recode(param, !!!param_labels))
best_long_params <- best_long_params %>% mutate(param = recode(param, !!!param_labels))

elas_long_all <- results_grid_valid %>%
  select(r, sigma_KL, sigma_VAE) %>%
  pivot_longer(c(sigma_KL, sigma_VAE), names_to = "param", values_to = "val") %>%
  mutate(param = recode(
    param,
    sigma_KL  = "σK–L",
    sigma_VAE = "σVA–E"
  )) %>%
  filter(is.finite(val))

best_elas <- best_methods %>%
  select(r, sigma_KL, sigma_VAE) %>%
  pivot_longer(c(sigma_KL, sigma_VAE), names_to = "param", values_to = "best_val") %>%
  mutate(param = recode(
    param,
    sigma_KL  = "σK–L",
    sigma_VAE = "σVA–E"
  ))

median_elas <- elas_long_all %>%
  group_by(r, param) %>%
  summarise(median_val = median(val, na.rm = TRUE), .groups = "drop")

plot_param_group <- function(valid_df, median_df, best_df, params, title, colours) {
  v_sub <- valid_df  %>% filter(param %in% params)
  m_sub <- median_df %>% filter(param %in% params)
  b_sub <- best_df   %>% filter(param %in% params)
  
  ord <- b_sub %>%
    group_by(r) %>%
    summarise(mu = mean(best_val, na.rm = TRUE), .groups = "drop") %>%
    arrange(mu) %>% pull(r)
  
  max_val  <- max(c(v_sub$val, b_sub$best_val), na.rm = TRUE)
  min_val  <- min(v_sub$val, na.rm = TRUE)
  x_offset <- (max_val - min_val) * 0.08
  
  ggplot(v_sub, aes(y = factor(r, levels = ord), x = val, fill = param, colour = param)) +
    geom_violin(alpha = 0.4, width = 0.8, scale = "width") +
    geom_point(
      data = m_sub,
      aes(x = median_val, y = factor(r, levels = ord), shape = "Median"),
      size = 2.5, inherit.aes = FALSE, colour = LGREY
    ) +
    geom_point(
      data = b_sub,
      aes(x = best_val, y = factor(r, levels = ord), shape = "Best method"),
      size = 2.8, inherit.aes = FALSE, colour = NAVY
    ) +
    geom_text(
      data = b_sub,
      aes(
        x = max_val + x_offset,
        y = factor(r, levels = ord),
        label = round(best_val, 2),
        colour = param
      ),
      inherit.aes = FALSE, hjust = 0, size = 3, fontface = "bold"
    ) +
    scale_fill_manual(values = colours) +
    scale_colour_manual(values = colours) +
    scale_shape_manual(values = c("Median" = 16, "Best method" = 18)) +
    labs(title = title, x = "Estimate", y = "Region",
         fill = "Parameter", shape = NULL, colour = "Parameter") +
    theme_nat() +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      plot.margin = margin(5, 60, 5, 5)
    ) +
    coord_cartesian(xlim = c(min_val, max_val + 3 * x_offset))
}

fig_s3_14a <- plot_param_group(
  valid_df  = valid_long,
  median_df = median_vals,
  best_df   = best_long_params,
  params    = c("γ", "ν"),
  title     = "Regional distributions of γ and ν",
  colours   = c("γ" = DGREY, "ν" = NAVY)
)

fig_s3_14b <- plot_param_group(
  valid_df  = valid_long,
  median_df = median_vals,
  best_df   = best_long_params,
  params    = c("λ"),
  title     = "Regional distributions of λ",
  colours   = c("λ" = NAVY)
)

fig_s3_14c <- plot_param_group(
  valid_df  = valid_long,
  median_df = median_vals,
  best_df   = best_long_params,
  params    = c("δK–VA", "δVA–Y"),
  title     = "Regional distributions of δ parameters (K–VA and VA–Y)",
  colours   = c("δK–VA" = DGREY, "δVA–Y" = NAVY)
)

fig_s3_14d <- plot_param_group(
  valid_df  = elas_long_all,
  median_df = median_elas,
  best_df   = best_elas,
  params    = c("σK–L", "σVA–E"),
  title     = "Regional distributions of elasticities",
  colours   = c("σK–L" = DGREY, "σVA–E" = NAVY)
)

print(fig_s3_14a)
print(fig_s3_14b)
print(fig_s3_14c)
print(fig_s3_14d)


# S3-15a Distribution of rho among best methods
best_rho_long <- best_methods %>%
  select(rho_KL, rho_VAE) %>%
  pivot_longer(everything(), names_to = "rho_type", values_to = "rho") %>%
  mutate(
    rho_type = factor(
      rho_type,
      levels = c("rho_KL", "rho_VAE"),
      labels = c("ρKL", "ρVA–E")
    )
  ) %>%
  filter(is.finite(rho))

rho_stats <- best_rho_long %>%
  group_by(rho_type) %>%
  summarise(
    n   = sum(is.finite(rho)),
    med = median(rho, na.rm = TRUE),
    .groups = "drop"
  )

dens_y_at <- function(xx, vec, adjust = 1.4, n = 2048) {
  d <- stats::density(vec[is.finite(vec)], adjust = adjust, n = n)
  stats::approx(d$x, d$y, xout = xx, rule = 2)$y
}

rho_stats$y_med <- vapply(
  seq_len(nrow(rho_stats)),
  function(i) dens_y_at(rho_stats$med[i],
                        best_rho_long$rho[best_rho_long$rho_type == rho_stats$rho_type[i]]),
  numeric(1)
)

n_regions <- dplyr::n_distinct(best_methods$r)
cap_rho_text <- paste0(
  "Symbols: ρKL, ρVA–E = CES exponents on the K–L and VA–E nests. n = ",
  n_regions, " regions"
)

x_lim <- range(best_rho_long$rho, na.rm = TRUE)

fig_s3_15a <- ggplot(best_rho_long, aes(x = rho, fill = rho_type)) +
  stat_density(
    geom = "area", position = "identity",
    alpha = 0.35, adjust = 1.4, n = 2048, colour = NA
  ) +
  geom_segment(
    data = rho_stats,
    aes(x = med, xend = med, y = 0, yend = y_med),
    inherit.aes = FALSE,
    linewidth = 0.8, linetype = "dashed", colour = DGREY
  ) +
  geom_text(
    data = rho_stats,
    aes(x = med, y = y_med, label = paste0("median = ", round(med, 2))),
    inherit.aes = FALSE,
    nudge_y = max(rho_stats$y_med, na.rm = TRUE) * 0.04,
    hjust = -0.05, size = 3.5, colour = DGREY
  ) +
  scale_fill_manual(values = c("ρKL" = MGREY, "ρVA–E" = NAVY), name = NULL) +
  scale_x_continuous(limits = x_lim, name = expression(rho)) +
  labs(
    title    = "Distributions of best-method ρ exponents",
    subtitle = "Kernel densities across regions",
    y        = "Kernel density (relative frequency)",
    caption  = cap_rho_text
  ) +
  theme_nat() +
  theme(legend.position = "top", legend.direction = "horizontal")
print(fig_s3_15a)


# S3-15b 2D grid heatmap: where best-method cells fall in the full first-run grid
if (all(c("grid_counts", "grid_base", "k_vals", "v_vals") %in% ls())) {
  
  ext_vals <- best_methods %>%
    summarise(
      min_KL  = min(rho_KL,  na.rm = TRUE),
      max_KL  = max(rho_KL,  na.rm = TRUE),
      min_VAE = min(rho_VAE, na.rm = TRUE),
      max_VAE = max(rho_VAE, na.rm = TRUE)
    ) %>% as.list()
  
  ext_cells <- grid_counts %>%
    filter(
      n > 0 &
        (rho_KL %in% c(ext_vals$min_KL, ext_vals$max_KL) |
           rho_VAE %in% c(ext_vals$min_VAE, ext_vals$max_VAE))
    )
  
  sigma_ticks <- c(10, 5, 2, 1, 0.5, 0.2, 0.1)
  sx_df <- tibble(rho = 1 / sigma_ticks - 1, lab = sigma_ticks) %>%
    filter(rho >= min(k_vals), rho <= max(k_vals))
  sy_df <- tibble(rho = 1 / sigma_ticks - 1, lab = sigma_ticks) %>%
    filter(rho >= min(v_vals), rho <= max(v_vals))
  
  x_exp <- (max(k_vals) - min(k_vals)) * 0.06
  y_exp <- (max(v_vals) - min(v_vals)) * 0.06
  
  lab_num <- function(x) sprintf("==%s", signif(x, 3))
  
  fig_s3_15b <- ggplot() +
    geom_tile(
      data = grid_base,
      aes(x = rho_KL, y = rho_VAE),
      fill = "white", colour = "#EAEAEA",
      linewidth = 0.25, width = 0.98, height = 0.98
    ) +
    geom_tile(
      data = subset(grid_counts, edge_KL | edge_VAE),
      aes(x = rho_KL, y = rho_VAE),
      fill = NA, colour = "#BDBDBD",
      linewidth = 0.5, width = 0.98, height = 0.98
    ) +
    geom_tile(
      data = grid_counts,
      aes(x = rho_KL, y = rho_VAE, fill = n),
      colour = "white", linewidth = 0.35, width = 0.98, height = 0.98
    ) +
    geom_text(
      data = subset(grid_counts, n > 0),
      aes(x = rho_KL, y = rho_VAE, label = n, colour = lab_col),
      size = 3.1, fontface = "bold", show.legend = FALSE
    ) +
    scale_colour_identity() +
    scale_fill_gradientn(
      colours = c(BBLUE, BLUE, NAVY),
      limits  = c(0, max(grid_counts$n, na.rm = TRUE)),
      breaks  = scales::breaks_pretty(n = 5),
      name    = "Regions"
    ) +
    geom_tile(
      data = ext_cells,
      aes(x = rho_KL, y = rho_VAE),
      fill = NA, colour = YELLOW,
      linewidth = 1.0, width = 0.98, height = 0.98
    ) +
    geom_vline(
      xintercept = c(ext_vals$min_KL, ext_vals$max_KL),
      linetype = "dashed", colour = DGREY, linewidth = 0.6
    ) +
    geom_hline(
      yintercept = c(ext_vals$min_VAE, ext_vals$max_VAE),
      linetype = "dashed", colour = DGREY, linewidth = 0.6
    ) +
    annotate(
      "label",
      x = ext_vals$min_KL, y = max(v_vals) + y_exp * 0.7,
      label = paste0("min~rho[K-L]", lab_num(ext_vals$min_KL)),
      parse = TRUE, size = 3, fill = "white", colour = DGREY
    ) +
    annotate(
      "label",
      x = ext_vals$max_KL, y = max(v_vals) + y_exp * 0.7,
      label = paste0("max~rho[K-L]", lab_num(ext_vals$max_KL)),
      parse = TRUE, size = 3, fill = "white", colour = DGREY
    ) +
    annotate(
      "label",
      x = max(k_vals) + x_exp * 0.5, y = ext_vals$min_VAE,
      label = paste0("min~rho[VA-E]", lab_num(ext_vals$min_VAE)),
      parse = TRUE, size = 3, fill = "white", colour = DGREY
    ) +
    annotate(
      "label",
      x = max(k_vals) + x_exp * 0.5, y = ext_vals$max_VAE,
      label = paste0("max~rho[VA-E]", lab_num(ext_vals$max_VAE)),
      parse = TRUE, size = 3, fill = "white", colour = DGREY
    ) +
    scale_x_continuous(
      name = expression(rho[KL]),
      expand = expansion(mult = c(0.00, 0.08))
    ) +
    scale_y_continuous(
      name = expression(rho[VA-E]),
      expand = expansion(mult = c(0.00, 0.10))
    ) +
    geom_text(
      data = sx_df,
      aes(x = rho, y = max(v_vals) + y_exp, label = lab),
      inherit.aes = FALSE, vjust = 0, size = 3, colour = DGREY
    ) +
    annotate(
      "text",
      x = mean(range(k_vals)), y = max(v_vals) + y_exp * 1.9,
      label = "sigma[K-L]", parse = TRUE,
      size = 3.2, colour = DGREY
    ) +
    geom_text(
      data = sy_df,
      aes(x = max(k_vals) + x_exp, y = rho, label = lab),
      inherit.aes = FALSE, hjust = 0, size = 3, colour = DGREY
    ) +
    annotate(
      "text",
      x = max(k_vals) + x_exp * 1.2, y = mean(range(v_vals)),
      label = "sigma[VA-E]", parse = TRUE,
      angle = -90, size = 3.2, colour = DGREY
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title    = "Best-method ρ cells on the first-run grid",
      subtitle = "Fill = # regions per (ρ[K–L], ρ[VA–E]); dashed = global min/max ρ; yellow outline = extreme cell used by any region.",
      caption  = "σ tick hints shown outside the panel (σ = 1/(1+ρ))."
    ) +
    theme_nat() +
    theme(
      legend.position = "right",
      panel.grid = element_blank(),
      plot.margin = margin(10, 40, 10, 10)
    )
  
  print(fig_s3_15b)
  
} else {
  message("Skipping Fig S3-15b: grid_counts/grid_base/k_vals/v_vals not found.")
}


###### STAGE 4: IAM-LEVEL CONSISTENCY CHECKS ######

# S4-1a: Estimated TFP growth (λ) by region 
fig_s4_1a <- best_methods %>%
  mutate(region = r) %>%
  ggplot(aes(x = reorder(region, lambda), y = lambda)) +
  geom_col(width = 0.70, fill = LGREY, colour = "white") +
  coord_flip() +
  labs(
    title   = "Estimated TFP growth rate (λ) by region — best methods",
    x       = NULL,
    y       = expression(lambda),
    caption = cap_params
  ) +
  theme_nat()
print(fig_s4_1a)


# S4-1b: IAM TFP trajectories (normalised)
fig_s4_1b <- iam_table %>%
  group_by(region) %>%
  mutate(TFP_norm = total_factor_productivity /
           first(total_factor_productivity)) %>%
  ungroup() %>%
  ggplot(aes(x = year, y = TFP_norm, group = region)) +
  geom_line(colour = DGREY, linewidth = 0.4) +
  facet_wrap(~ region, scales = "free_y") +
  labs(
    title   = "IAM input: TFP trajectories normalised to base year",
    x       = "Year",
    y       = "TFP / TFP(base)",
    caption = cap_params
  ) +
  theme_nat()
print(fig_s4_1b)


# S4-2: IAM TFP growth vs estimated λ
iam_lambda <- iam_table %>%
  group_by(region) %>%
  summarise(
    iam_lambda = coef(lm(log(total_factor_productivity) ~ year))[2],
    .groups = "drop"
  )

lambda_compare <- best_methods %>%
  select(r, lambda) %>%
  rename(region = r) %>%
  inner_join(iam_lambda, by = "region")

fig_s4_2 <- ggplot(lambda_compare,
                   aes(x = iam_lambda, y = lambda)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = DGREY) +
  geom_point(size = 2.4, colour = NAVY) +
  ggrepel::geom_text_repel(
    aes(label = region),
    size = 3, colour = DGREY,
    max.overlaps = 20
  ) +
  labs(
    title    = "Estimated vs IAM-implied TFP growth",
    subtitle = "IAM slope from log(TFP) ~ year; dashed = 1:1",
    x        = expression(lambda[IAM]),
    y        = expression(lambda[estimated])
  ) +
  theme_nat()
print(fig_s4_2)


# S4-3: Residual structure vs IAM growth
if (file.exists("CES_results_time.csv")) {
  
  results_time <- read_csv("CES_results_time.csv", show_col_types = FALSE)
  
  resid_growth <- results_time %>%
    inner_join(best_methods %>% select(r, method), by = c("r", "method")) %>%
    group_by(r) %>%
    summarise(
      sd_resid = sd(residual, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    inner_join(lambda_compare, by = "r")
  
  fig_s4_3 <- ggplot(resid_growth,
                     aes(x = lambda, y = sd_resid)) +
    geom_point(size = 2.6, colour = NAVY) +
    ggrepel::geom_text_repel(
      aes(label = r),
      size = 3, colour = DGREY,
      max.overlaps = 20
    ) +
    labs(
      title    = "Residual dispersion vs estimated TFP growth",
      subtitle = "Checks whether high λ regions systematically misfit",
      x        = expression(lambda),
      y        = "Residual SD"
    ) +
    theme_nat()
  print(fig_s4_3)
  
} else {
  message("Skipping Fig S4-3: CES_results_time.csv not found.")
}


# S4-4: IAM-consistency of factor shares
share_check <- best_methods %>%
  transmute(
    r,
    delta_KVA,
    delta_VAY,
    sum_shares = delta_KVA + delta_VAY
  )

fig_s4_4 <- ggplot(share_check,
                   aes(x = reorder(r, sum_shares), y = sum_shares)) +
  geom_col(width = 0.70, fill = LGREY, colour = "white") +
  geom_hline(yintercept = 1,
             linetype = "dashed", colour = NAVY) +
  coord_flip() +
  labs(
    title    = "Consistency of CES share parameters",
    subtitle = "δK–VA + δVA–Y should be ≤ 1 for IAM coherence",
    x        = NULL,
    y        = expression(delta[K-VA] + delta[VA-Y])
  ) +
  theme_nat()
print(fig_s4_4)


# S4-5: Elasticity vs IAM sectoral rigidity proxy
if ("sector_rigidity" %in% names(iam_table)) {
  
  elas_rigid <- best_methods %>%
    select(r, sigma_KL, sigma_VAE) %>%
    left_join(
      iam_table %>%
        group_by(region) %>%
        summarise(rigidity = mean(sector_rigidity, na.rm = TRUE),
                  .groups = "drop"),
      by = c("r" = "region")
    )
  
  fig_s4_5 <- ggplot(elas_rigid,
                     aes(x = rigidity, y = sigma_KL)) +
    geom_point(size = 2.4, colour = NAVY) +
    geom_smooth(method = "lm", se = FALSE,
                linetype = "dashed", colour = DGREY) +
    labs(
      title    = "Capital–labour elasticity vs IAM rigidity proxy",
      subtitle = "Illustrative consistency check",
      x        = "IAM rigidity proxy",
      y        = expression(sigma[K-L])
    ) +
    theme_nat()
  print(fig_s4_5)
  
} else {
  message("Skipping Fig S4-5: IAM rigidity proxy not available.")
}

