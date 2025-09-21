options(scipen = 999) #avoids scientific notation unless necessary
# install.packages(c("micEconCES","dplyr","readr","purrr","ggplot2","parallel",
#"ggpmisc","pheatmap","GGally","ggcorrplot","ggridges","ggpubr"))
#install.packages("ggpubr")
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



### Settings
setwd("C:/Users/escami_g/OneDrive - Paul Scherrer Institut/05.Models/MERGE updates/CES-parametrisation/fine grid")
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
results            <- readRDS("results_run1.rds")
results_table      <- read_csv("CES_region_method.csv")
results_table_valid<- read_csv("CES_region_method_valid.csv")
results_table_time <- read_csv("CES_region_method_year.csv")
results_grid       <- read_csv("CES_gridsearch.csv")
convergence_summary<- read_csv("CES_convergence_summary.csv")
best_methods       <- read_csv("CES_best_methods.csv")
aic_weights        <- read_csv("CES_AICc_weights.csv")
grid_conv_share    <- read_csv("CES_grid_convergence_share.csv")
iam_table          <- read_csv("IAM_params.csv")


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
BBLUE <- "#DAE3F3"
LGREY <- "#c7c7c7"
DGREY <- "#4d4d4d"
YELLOW <- "#FFC000"

# Dynamic captions so symbols match each plot’s content
cap_none      <- ""
cap_elast     <- "Symbols: σ[K–L], σ[VA–E] = substitution elasticities."
cap_rho       <- "Symbols: ρ[KL], ρ[VA–E] = CES exponents on the K–L and VA–E nests."
cap_fit       <- expression(paste("Symbols: ", R^2, " = coefficient of determination; ε = residual."))
cap_params    <- "Symbols: γ, λ, ν = scale/growth/curvature; δ[K–L], δ[VA–E] = share parameters."


pct <- function(x) scales::percent(x, accuracy = 1)
clip01 <- function(x) pmax(pmin(x, 1), 0)
near_   <- function(x, a, tol = 1e-12) abs(x - a) <= tol


# MAth labels
lab_sigma_KL  <- bquote(sigma[K - L])
lab_sigma_VAE <- bquote(sigma[VA - E])
lab_rho_KL    <- bquote(rho[KL])
lab_rho_VAE   <- bquote(rho[VAE])

# Consistent method ordering (fastest & most stable first) = by highest convergence share,
# break ties by lowest median runtime (mins)
method_perf <- results_table %>%
  distinct(r, method, conv) %>%
  group_by(method) %>%
  summarise(
    share_conv = mean(conv, na.rm = TRUE),
    med_rt_min = median(results_table$runtime[results_table$method == first(method)]/60, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(share_conv), med_rt_min)
method_order <- method_perf$method
shape_vals   <- c(21,22,23,24,25,3,4,1)[seq_along(method_order)]

# Global medians for best-method elasticities
med_KL  <- median(best_methods$sigma_KL,  na.rm = TRUE)
med_VAE <- median(best_methods$sigma_VAE, na.rm = TRUE)

# Load mapping
region_map <- read_excel("MERGE regions proposal.xlsx", sheet = 1)

region_map <- region_map %>%
  rename(r = MERGE) %>%
  mutate(r = as.character(r))

world <- ne_countries(scale = "medium", returnclass = "sf")

best_map <- best_methods %>%
  inner_join(region_map, by = "r")

world_best <- world %>%
  left_join(best_map, by = c("iso_a3" = "ISO3"))

### 1. Diagnostics of all runs
# 1.1 Convergence by method (% and counts)
grid_axis <- results_grid %>%
  mutate(conv = ifelse(is.null(convergence), !is.na(rss), as.logical(convergence))) %>%
  select(r, method, rho1, rho, conv)

kl_share <- grid_axis %>%
  group_by(method, rho1) %>% summarise(ok = any(conv, na.rm = TRUE), .groups = "drop") %>%
  group_by(method) %>% summarise(share = mean(ok), axis = "ρ[KL]", .groups = "drop")

vae_share <- grid_axis %>%
  group_by(method, rho) %>% summarise(ok = any(conv, na.rm = TRUE), .groups = "drop") %>%
  group_by(method) %>% summarise(share = mean(ok), axis = "ρ[VA–E]", .groups = "drop")

axis_conv <- bind_rows(kl_share, vae_share) %>%
  mutate(method = factor(method, levels = method_order),
         axis   = factor(axis, levels = c("ρ[KL]","ρ[VA–E]")))

fig1_1 <- ggplot(axis_conv, aes(x = method, y = share, fill = axis)) +
  geom_col(position = position_dodge(width = .7), width = .65, colour = "white") +
  geom_text(aes(label = pct(share)),
            position = position_dodge(width = .7), vjust = -0.15, size = 3, colour = DGREY) +
  coord_flip() +
  scale_fill_manual(values = c("ρ[KL]" = DGREY, "ρ[VA–E]" = NAVY)) +
  labs(
    title    = "Axis-wise convergence on ρ (grid coverage by method)",
    subtitle = "Two bars per method: share of ρ values with at least one converged partner on the other axis",
    x = NULL, y = "Share converged", fill = NULL,
    caption  = cap_rho
  ) +
  theme_nat()
print(fig1_1)

# 1.2 Region × Method convergence heatmap
region_conv <- results_table %>%
  distinct(r, method, conv) %>%
  group_by(r) %>%
  summarise(share_methods = mean(conv, na.rm = TRUE), n_methods = n(), .groups = "drop") %>%
  arrange(desc(share_methods))

fig1_2 <- ggplot(region_conv, aes(x = reorder(r, share_methods), y = share_methods)) +
  geom_col(width = .70, fill= LGREY, colour = "white") +
  geom_text(aes(label = pct(share_methods)), hjust = -0.05, colour = DGREY, size = 3) +
  coord_flip(ylim = c(0, 1.05)) +
  labs(
    title   = "Convergence by region (share of methods that converged)",
    subtitle= "Each bar = region; value = fraction of optimisation methods that reached convergence",
    x = NULL, y = "Share of methods converged", caption = cap_none
  ) +
  theme_nat()
print(fig1_2)

# 1.3 Best ρ medians & IQR by method (motivation for grid refinement)
rho_long <- results_table_valid %>%                     # use valid fits only
  select(method, rho_KL, rho_VAE) %>%
  pivot_longer(cols = c(rho_KL, rho_VAE),               # <- pivot only the rho columns
               names_to = "rho_type", values_to = "rho") %>%
  mutate(
    rho_type = recode(rho_type, rho_KL = "ρ[KL]", rho_VAE = "ρ[VA–E]"),
    method   = factor(method, levels = method_order)
  ) %>%
  filter(is.finite(rho))

rho_stats <- rho_long %>%
  group_by(method, rho_type) %>%
  summarise(
    med = median(rho, na.rm = TRUE),
    q1  = quantile(rho, .25, na.rm = TRUE),
    q3  = quantile(rho, .75, na.rm = TRUE),
    .groups = "drop"
  )

fig1_3 <- ggplot(rho_stats, aes(x = method, y = med, colour = rho_type)) +
  geom_errorbar(aes(ymin = q1, ymax = q3),
                width = .25, linewidth = .5, colour = DGREY,
                position = position_dodge(width = .5)) +
  geom_point(position = position_dodge(width = .5), size = 2.6) +
  coord_flip() +
  scale_colour_manual(values = c("ρ[KL]" = DGREY, "ρ[VA–E]" = NAVY)) +
  labs(
    title    = "Best-fit ρ medians with IQR by method",
    subtitle = "Points = medians; bars = interquartile range (valid runs only)",
    x = NULL, y = "ρ (exponent)", colour = NULL,
    caption  = "Symbols: ρ[KL], ρ[VA–E] = CES exponents on the K–L and VA–E nests."
  ) +
  theme_nat()
print(fig1_3)

# 1.4 Runtime by method (minutes), sorted fast→slow, with bigger outliers
rt_stats <- results_table %>%
  group_by(method) %>%
  summarise(med_min = median(runtime/60, na.rm = TRUE), .groups = "drop") %>%
  arrange(med_min)

fig1_4 <- results_table %>%
  mutate(method = factor(method, levels = rt_stats$method)) %>%
  ggplot(aes(x = method, y = runtime/60)) +
  geom_boxplot(width = .72, fill= LGREY, colour = DGREY, outlier.size = 1.4, outlier.alpha = .6) +
  stat_summary(fun = median, geom = "point", size = 2.2, shape = 21, fill = NAVY, colour = NAVY) +
  coord_flip() +
  labs(
    title    = "Runtime by method",
    subtitle = "Minutes; navy = method median",
    x = NULL, y = "Runtime (minutes)", caption = cap_none
  ) +
  theme_nat()
print(fig1_4)

# 1.5 RSS vs iterations
fit_df <- results_table %>%
  filter(is.finite(iter), iter > 0, is.finite(rss), rss > 0) %>%
  mutate(x = log(iter + 1), y = log(rss))

m_poly  <- lm(y ~ poly(x, 2, raw = TRUE), data = fit_df)
b       <- coef(m_poly)                   # (Intercept), x, x^2
R2_poly <- summary(m_poly)$r.squared

# Plotmath string (single character vector) for annotate()
eq_str <- sprintf(
  "log(RSS)==%.3f+%.3f*log(Iter+1)+%.3f*log(Iter+1)^2~~(R^2==%.3f)",
  b[1], b[2], b[3], R2_poly
)

pred_line <- tibble(x = seq(min(fit_df$x), max(fit_df$x), length.out = 200)) |>
  mutate(y = predict(m_poly, newdata = tibble(x = x)))

pal_methods <- colorRampPalette(c(DGREY, NAVY))(length(unique(results_table$method)))

fig1_5 <- ggplot(fit_df, aes(x = x, y = y, colour = factor(method, levels = method_order))) +
  geom_point(alpha = .7, size = 1.9) +
  geom_line(data = pred_line, aes(x = x, y = y),
            inherit.aes = FALSE, colour = NAVY, linetype = "dashed") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.02, vjust = 1.2,
           label = eq_str, parse = TRUE, colour = DGREY, size = 3.2) +
  scale_colour_manual(values = pal_methods, name = "Method") +
  labs(
    title    = "Fit quality vs iterations (global quadratic trend)",
    subtitle = "Points = runs coloured by method; dashed = pooled quadratic in log space",
    x = "log(Iterations+1)", y = "log(RSS)", caption = cap_fit
  ) +
  theme_nat()
print(fig1_5)

### 2. Valid runs diagnostics
# 2.1 Method performance summary (n, median R², median runtime)
valid_summary <- results_table_valid %>%
  distinct(r, method, .keep_all = TRUE) %>%
  group_by(method) %>%
  summarise(
    share_regions = n()/n_distinct(results_table$r),
    med_R2        = median(R2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(share_regions), desc(med_R2))

fig2_1 <- results_table_valid %>%
  group_by(method) %>%
  summarise(n = n(),
            med_R2 = median(R2, na.rm = TRUE),
            med_runtime = median(runtime, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(med_R2)) %>%
  mutate(method = factor(method, levels = method)) %>%
  ggplot(aes(x = method, y = n)) +
  geom_col(width = .70, fill = LGREY, colour = "white") +
  geom_text(aes(label = paste0("R²~", round(med_R2, 2),
                               " | t~", round(med_runtime/60, 1), " min")),
            vjust = -0.25, size = 3.1, colour = DGREY) +
  coord_flip() +
  labs(
    title = "Valid runs by method",
    subtitle = "Bar = count; label = median R² and runtime (min)",
    x = NULL, y = "Valid runs",
    caption = cap_fit
  ) +
  theme_nat()
print(fig2_1)

# 2.2 Elasticity distributions overlaid (same x-axis), by method and overall
fig2_2 <- results_table_valid %>%
  select(r, method, sigma_KL, sigma_VAE) %>%
  pivot_longer(c(sigma_KL, sigma_VAE), names_to = "elasticity", values_to = "sigma") %>%
  mutate(
    elasticity = recode(elasticity, sigma_KL = "σ[K-L]", sigma_VAE = "σ[VA-E]"),
    sigma_plot = pmin(sigma, 10)
  ) %>%
  ggplot(aes(x = sigma_plot, y = method)) +
  ggridges::geom_density_ridges(scale = 1.05, rel_min_height = 0.01,
                                fill = LGREY, colour = DGREY) +
  facet_wrap(~elasticity, scales = "free_x") +
  labs(
    title = "Substitution elasticities across methods (valid runs)",
    subtitle = "Distributions clipped at σ = 10",
    x = expression(sigma), y = "Method", caption = cap_elast
  ) +
  theme_nat()
print(fig2_2)

# 2.3 Two bars per method (median σ + IQR as errorbars), sorted by σ[VA–E]
elas_long <- results_table_valid %>%
  select(method, sigma_KL, sigma_VAE) %>%
  pivot_longer(c(sigma_KL, sigma_VAE),
               names_to = "which", values_to = "sigma") %>%
  mutate(
    which  = recode(which, sigma_KL = "σ[K–L]", sigma_VAE = "σ[VA–E]"),
    method = factor(method, levels = method_order)
  ) %>%
  filter(is.finite(sigma))

elas_stats <- elas_long %>%
  group_by(method, which) %>%
  summarise(med = median(sigma, na.rm = TRUE),
            q1  = quantile(sigma, .25, na.rm = TRUE),
            q3  = quantile(sigma, .75, na.rm = TRUE), .groups = "drop")

order_by_vae <- elas_stats %>% filter(which == "σ[VA–E]") %>% arrange(med) %>% pull(method)

fig2_3 <- ggplot(elas_stats, aes(x = factor(method, levels = order_by_vae), y = med, fill = which)) +
  geom_col(position = position_dodge(width = .7), width = .6, colour = "white") +
  geom_errorbar(aes(ymin = q1, ymax = q3),
                position = position_dodge(width = .7), width = .2, colour = DGREY) +
  coord_flip() +
  scale_fill_manual(values = c("σ[K–L]" = DGREY, "σ[VA–E]" = NAVY)) +
  labs(
    title    = "Median elasticities by method with IQR",
    subtitle = "Two bars per method (σ[K–L] and σ[VA–E])",
    x = NULL, y = expression(sigma), fill = NULL,
    caption  = cap_elast
  ) +
  theme_nat()
print(fig2_3)

# 2.4 Core parameter distributions (γ, λ, δ[K–L], δ[VA–E], ν)
best_params_long <- best_methods %>%
  select(r, gamma, lambda, nu, delta_KL, delta_VAE) %>%
  pivot_longer(-r, names_to = "param", values_to = "val") %>%
  mutate(param = recode(param,
                        gamma = "γ", lambda = "λ", nu = "ν",
                        delta_KL = "δ[K–L]", delta_VAE = "δ[VA–E]")) %>%
  mutate(val = case_when(
    param %in% c("δ[K–L]","δ[VA–E]") ~ pmin(pmax(val, 0), 1),   # clamp shares to [0,1]
    TRUE ~ val
  )) %>%
  filter(is.finite(val))

fig2_4 <- ggplot(best_params_long, aes(x = param, y = val)) +
  geom_violin(fill= LGREY, colour = DGREY, width = .8, alpha = .7, trim = TRUE) +
  geom_jitter(width = .12, height = 0, size = 1.1, alpha = .6, colour = DGREY) +
  stat_summary(fun = median, geom = "point", size = 2.2, shape = 21, fill = NAVY, colour = NAVY) +
  labs(
    title    = "Core parameter distributions (best methods)",
    subtitle = "Violin = distribution; dots = regions; navy = median",
    x = "Parameter", y = "Value", caption = cap_params
  ) +
  theme_nat()
print(fig2_4)

# 2.5 Parameter correlations
cor_df <- results_table_valid %>%
  select(gamma, lambda, nu, delta_KL, delta_VAE, sigma_KL, sigma_VAE) %>%
  select(where(~ sum(is.finite(.x)) > 1))
fig2_5 <- if (ncol(cor_df) >= 2) {
  cor_mat <- cor(cor_df, use = "pairwise.complete.obs")
  ggcorrplot::ggcorrplot(
    cor_mat, hc.order = TRUE, type = "lower", lab = TRUE, outline.color = "white",
    colors = c("white", "#9BB3C9", NAVY)
  ) +
    labs(title = "Parameter correlations (valid runs)",
         subtitle = "Lower triangle; hierarchical ordering",
         caption  = cap_params) +
    theme_nat()
} else NULL
print(fig2_5)

# 2.6 Residual SD distribution
res_by_reg <- results_table_time %>%
  group_by(r) %>%
  summarise(sd_res = sd(residual, na.rm = TRUE), n = sum(is.finite(residual)), .groups = "drop") %>%
  arrange(desc(sd_res))

fig2_6 <- results_table_time %>%
  mutate(r = factor(r, levels = res_by_reg$r)) %>%
  ggplot(aes(x = r, y = residual)) +
  geom_boxplot(fill= LGREY, colour = DGREY, outlier.size = 1.2, outlier.alpha = .5) +
  geom_text(data = res_by_reg, aes(x = r, y = Inf, label = paste0("n=", n)),
            inherit.aes = FALSE, vjust = 1.3, colour = DGREY, size = 3) +
  coord_flip() +
  labs(
    title    = "Residuals by region (valid fits)",
    subtitle = "Sorted by residual standard deviation; labels show number of time points",
    x = NULL, y = "Residual (log)", caption = cap_fit
  ) +
  theme_nat()
print(fig2_6)

# 2.7 Parameter significance (best methods): p-value tiles by region and parameter
p_long <- best_methods %>%
  transmute(r,
            `γ`  = p_gamma,
            `λ`  = p_lambda,
            `δ[K–L]`  = p_delta_KL,
            `δ[VA–E]` = p_delta_VAE,
            `ν`  = p_nu) %>%
  pivot_longer(-r, names_to = "parameter", values_to = "p") %>%
  mutate(class = case_when(
    is.na(p)            ~ "NA",
    p < 0.01            ~ "<0.01",
    p < 0.05            ~ "<0.05",
    TRUE                ~ "ns"
  ))

fig2_7 <- ggplot(p_long, aes(x = parameter, y = reorder(r, as.numeric(factor(class))), fill = class)) +
  geom_tile(width = .9, height = .9, colour = "white") +
  scale_fill_manual(values = c("<0.01" = NAVY, "<0.05" = "#6B86A3", "ns"= LGREY, "NA" = "white")) +
  labs(
    title    = "Parameter significance in best-method fits",
    subtitle = "Tiles by region × parameter; colour encodes p-value class",
    x = "Parameter", y = "Region", fill = "p-value",
    caption = cap_params
  ) +
  theme_nat()
print(fig2_7)


### 3. BEST METHODS
# 3.1 Best-method count and share
fig3_1 <- best_methods %>%
  count(method) %>%
  mutate(share = n/sum(n)) %>%
  arrange(desc(share)) %>%
  ggplot(aes(x = reorder(method, share), y = share)) +
  geom_col(width = .70, fill= LGREY, colour = "white") +
  geom_text(aes(label = paste0(pct(share), "  (n=", n, ")")),
            hjust = -0.05, size = 3.2, colour = DGREY) +
  coord_flip(ylim = c(0, 1.05)) +
  labs(title = "Best method by region (share)", x = NULL, y = "Share of regions") +
  theme_nat()
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
shape_vals <- rep(c(21,22,23,24,25), length.out = length(method_order))
fig3_3 <- best_methods %>%
  mutate(method = factor(method, levels = method_order)) %>%
  ggplot(aes(sigma_KL, sigma_VAE)) +
  geom_vline(xintercept = med_KL,  linetype = "dashed", colour = NAVY, linewidth = .4) +
  geom_hline(yintercept = med_VAE, linetype = "dashed", colour = NAVY, linewidth = .4) +
  geom_point(aes(shape = method, fill = method), size = 2.8, colour = DGREY, alpha = .95) +
  scale_shape_manual(values = shape_vals) +
  scale_fill_manual(values = scales::hue_pal()(length(method_order))) +
  ggrepel::geom_text_repel(aes(label = r), colour = DGREY, size = 3,
                           max.overlaps = 30, min.segment.length = 0.05,
                           box.padding = 0.3, segment.alpha = 0.3, show.legend = FALSE) +
  labs(title   = "Best-method elasticities per region (quadrant view)",
       subtitle= "Filled shapes by method; dashed lines = global medians",
       x = lab_sigma_KL, y = lab_sigma_VAE, shape = "Method", fill = "Method",
       caption = cap_elast) +
  theme_nat()
print(fig3_3)

# 3.4 Combined-method regional view (cloud + IQR) with best marked
elas_reg <- results_table_valid %>%
  select(r, method, sigma_KL, sigma_VAE) %>%
  pivot_longer(c(sigma_KL, sigma_VAE), names_to = "which", values_to = "sigma") %>%
  mutate(which = recode(which, sigma_KL = "σ[K–L]", sigma_VAE = "σ[VA–E]")) %>%
  filter(is.finite(sigma))

elas_iqr <- elas_reg %>%
  group_by(r, which) %>%
  summarise(q1 = quantile(sigma, 0.25, na.rm = TRUE),
            q3 = quantile(sigma, 0.75, na.rm = TRUE),
            med = median(sigma, na.rm = TRUE), .groups = "drop")

best_long <- best_methods %>%
  select(r, sigma_KL, sigma_VAE) %>%
  pivot_longer(c(sigma_KL, sigma_VAE), names_to = "which", values_to = "sigma") %>%
  mutate(which = recode(which, sigma_KL = "σ[K–L]", sigma_VAE = "σ[VA–E]"))

ord_reg <- elas_iqr %>% group_by(r) %>%
  summarise(mu = mean(med, na.rm = TRUE), .groups = "drop") %>%
  arrange(mu) %>% pull(r)

fig3_4 <- ggplot() +
  geom_linerange(data = elas_iqr %>% mutate(r = factor(r, levels = ord_reg)),
                 aes(y = r, xmin = q1, xmax = q3), size = 2.2, colour= LGREY, alpha = .7) +
  geom_point(data = elas_reg %>% mutate(r = factor(r, levels = ord_reg)),
             aes(x = sigma, y = r), size = 1.8, alpha = .55, colour = DGREY) +
  geom_point(data = best_long %>% mutate(r = factor(r, levels = ord_reg)),
             aes(x = sigma, y = r), shape = 8, size = 2.6, colour = NAVY) +
  facet_wrap(~ which, scales = "free_x") +
  labs(
    title    = "Elasticities by region: all methods (cloud) + IQR + best",
    subtitle = "Grey points = all valid methods; grey bar = IQR; navy star = best method",
    x = expression(sigma), y = "Region", caption = cap_elast
  ) +
  theme_nat()
print(fig3_4)

# 3.5 Observed (y) vs Fitted (x), facets ordered by R²
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
  labs(
    title    = "Observed vs fitted (best methods)",
    subtitle = "Observed reconstructed from CSV as fitted × exp(residual); dashed = 1:1",
    x = "Fitted (scaled)", y = "Observed (scaled)",
    caption = cap_fit
  ) +
  theme_nat()
print(fig3_5)

# 3.6 Best K-L elasticities map
fig3_6 <- ggplot(world_best) +
  geom_sf(aes(fill = pmin(sigma_KL, 10)), colour = LGREY, linewidth = .01) +
  scale_fill_gradient(low = BBLUE, high = NAVY, na.value = "grey90") +
  labs(title = expression(paste(sigma[K-L], " by region (best method)")),
       fill = bquote(.(lab_sigma_KL))) +
  theme_nat()
print(fig3_6)

# 3.7 Best VA-E elasticities map
fig3_7 <- ggplot(world_best) +
  geom_sf(aes(fill = pmin(sigma_VAE, 10)), colour = LGREY, linewidth = .01) +
  scale_fill_gradient(low = BBLUE, high = NAVY, na.value = "grey90") +
  labs(title = expression(paste(sigma[VA-E], " by region (best method)")),
       fill = bquote(.(lab_sigma_VAE))) +
  theme_nat()
print(fig3_7)

# 3.8 Selection certainty (Δ AICc weight best − runner-up) + distribution
aic_top2 <- aic_weights %>%
  arrange(r, desc(wAICc)) %>%
  group_by(r) %>% slice_head(n = 2) %>%
  summarise(best = first(wAICc), runner = dplyr::last(wAICc),
            delta = best - runner, .groups = "drop")

fig3_8a <- ggplot(aic_top2, aes(x = reorder(r, delta), y = delta)) +
  geom_col(width = .70, fill= LGREY, colour = "white") +
  geom_hline(yintercept = .2, colour = NAVY, linetype = "dotted", linewidth = .4) +
  coord_flip() +
  labs(
    title = "Selection certainty across regions",
    subtitle = "Δ weight (best – runner-up); dotted line at 0.2 ≈ clear preference",
    x = NULL, y = "Δ AICc weight"
  ) +
  theme_nat()
print(fig3_8a)

fig3_8b <- ggplot(aic_top2, aes(x = delta)) +
  geom_histogram(binwidth = .05, fill= LGREY, colour = "white") +
  labs(title = "Distribution of selection certainty (Δ weight)", x = "Δ AICc weight", y = "Count") +
  theme_nat()
print(fig3_8b)


### 4. Insights from IAM export
# 4.1 TFP over time (faceted by region)
fig4_1a <- best_methods %>%
  mutate(region = r) %>%
  ggplot(aes(x = reorder(region, lambda), y = lambda)) +
  geom_col(width = .70, fill= LGREY, colour = "white") +
  coord_flip() +
  labs(title = "Estimated TFP growth rate (λ) by region — best methods",
       x = NULL, y = expression(lambda), caption = cap_params) +
  theme_nat()
print(fig4_1a)

fig4_1b <- iam_table %>%
  group_by(region) %>%
  mutate(TFP_norm = total_factor_productivity / first(total_factor_productivity)) %>%
  ungroup() %>%
  ggplot(aes(x = year, y = TFP_norm, group = region)) +
  geom_line(colour = DGREY) +
  facet_wrap(~ region, scales = "free_y") +
  labs(title = "IAM input: TFP trajectories normalised to base",
       x = "Year", y = "TFP / TFP[base]", caption = cap_params) +
  theme_nat()
print(fig4_1b)

