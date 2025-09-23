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
setwd("C:/Users/escami_g/OneDrive - Paul Scherrer Institut/05.Models/MERGE updates/CES-parametrisation/test")
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
results_grid <- read_csv("CES_results_grid.csv")
best_methods       <- read_csv("CES_best_methods.csv")
aic_weights        <- read_csv("CES_AICc_weights.csv")
iam_table          <- read_csv("IAM_params.csv")
results_grid_valid <- read_csv("CES_results_grid_valid.csv") 
results_grid_invalid <- read_csv("CES_results_grid_invalid.csv") 


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
cap_elast     <- "Symbols: σK-L, σVA–E = substitution elasticities."
cap_rho       <- "Symbols: ρKL, ρVA–E = CES exponents on the K–L and VA–E nests."
cap_fit       <- expression(paste("Symbols: ", R^2, " = coefficient of determination; ε = residual."))
cap_params    <- "Symbols: γ, λ, ν = scale/growth/curvature; δK-L, δVA–E = share parameters."


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
# 1.1 Axis-wise convergence on rho values
grid_axis <- results_grid %>% select(r, method, rho_KL, rho_VAE, conv)

kl_share <- grid_axis %>%
  group_by(method, rho_KL) %>% summarise(ok = any(conv, na.rm = TRUE), .groups = "drop") %>%
  group_by(method) %>% summarise(share = mean(ok), axis = "ρKL", .groups = "drop")

vae_share <- grid_axis %>%
  group_by(method, rho_VAE) %>% summarise(ok = any(conv, na.rm = TRUE), .groups = "drop") %>%
  group_by(method) %>% summarise(share = mean(ok), axis = "ρVA–E", .groups = "drop")

axis_conv <- bind_rows(kl_share, vae_share) %>%
  mutate(method = factor(method, levels = method_order),
         axis   = factor(axis, levels = c("ρKL","ρVA–E")))

fig1_1 <- ggplot(axis_conv, aes(x = method, y = share, fill = axis)) +
  geom_col(position = position_dodge(width = .7), width = .65, colour = "white") +
  geom_text(aes(label = paste0(round(share*100,1), "%")),
            position = position_dodge(width = .7), vjust = -0.15, size = 3, colour = DGREY) +
  coord_flip() +
  scale_fill_manual(values = c("ρKL" = DGREY, "ρVA–E" = NAVY)) +
  labs(
    title    = "Axis-wise convergence on ρ (grid coverage by method)",
    subtitle = "Share of ρ values with at least one converged partner",
    x = NULL, y = "Share converged", fill = NULL, caption = cap_rho
  ) + theme_nat()
print(fig1_1)

# 1.2 Region × Method convergence (share of methods converged per region)
region_conv <- results_grid_valid %>%
  distinct(r, method, conv) %>%
  group_by(r) %>%
  summarise(share_methods = mean(conv, na.rm = TRUE), n_methods = n(), .groups = "drop") %>%
  arrange(desc(share_methods))
fig1_2 <- plot_bar(region_conv, "r", "share_methods",
                   "Convergence by region",
                   "Fraction of methods converged",
                   "Share", fillcol = NAVY,
                   label_fmt = function(x) paste0(round(x*100,1),"%"))
print(fig1_2)

# 1.3 Best-fit rho medians & IQR by method (valid runs)
rho_long <- results_grid %>%
  select(method, rho_KL, rho_VAE) %>%
  pivot_longer(c(rho_KL, rho_VAE), names_to = "rho_type", values_to = "rho") %>%
  mutate(rho_type = recode(rho_type, rho_KL="ρKL", rho_VAE="ρVA–E"),
         method = factor(method, levels = method_order)) %>%
  filter(is.finite(rho))
rho_stats <- rho_long %>%
  group_by(method, rho_type) %>%
  summarise(med = median(rho, na.rm = TRUE),
            q1  = quantile(rho, .25, na.rm = TRUE),
            q3  = quantile(rho, .75, na.rm = TRUE), .groups = "drop")
fig1_3 <- ggplot(rho_stats, aes(x = method, y = med, colour = rho_type)) +
  geom_errorbar(aes(ymin = q1, ymax = q3), width = .25, linewidth = .5, colour = DGREY,
                position = position_dodge(width = .5)) +
  geom_point(position = position_dodge(width = .5), size = 2.6) +
  coord_flip() +
  scale_colour_manual(values = c("ρKL" = DGREY, "ρVA–E" = NAVY)) +
  labs(title = "Best-fit ρ medians with IQR by method",
       subtitle = "Points = medians; bars = interquartile range (valid runs only)",
       x = NULL, y = "ρ (exponent)", colour = NULL,
       caption  = "Symbols: ρKL, ρVA–E = CES exponents on the K–L and VA–E nests.") +
  theme_nat()
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
  stat_summary(fun = median, geom = "point", size = 2.2, shape = 21, fill = NAVY, colour = NAVY) +
  coord_flip() +
  labs(title = "Runtime by method", subtitle = "Minutes; navy = method median",
       x = NULL, y = "Runtime (minutes)", caption = cap_none) + theme_nat()
print(fig1_4)

# 1.5 RSS vs iterations (global quadratic)
fit_df <- results_grid %>%
  filter(is.finite(iter), iter > 0, is.finite(rss), rss > 0) %>%
  mutate(x = log(iter + 1), y = log(rss))
if (nrow(fit_df) >= 5) {
  m_poly  <- lm(y ~ poly(x, 2, raw = TRUE), data = fit_df)
  b       <- coef(m_poly); R2_poly <- summary(m_poly)$r.squared
  eq_str <- sprintf("log(RSS)==%.3f+%.3f*log(Iter+1)+%.3f*log(Iter+1)^2~~(R^2==%.3f)", b[1], b[2], b[3], R2_poly)
  pred_line <- tibble(x = seq(min(fit_df$x), max(fit_df$x), length.out = 200)) |>
    mutate(y = predict(m_poly, newdata = tibble(x = x)))
  pal_methods <- scales::hue_pal()(length(unique(results_grid$method)))
  fig1_5 <- ggplot(fit_df, aes(x = x, y = y, colour = factor(method, levels = method_order))) +
    geom_point(alpha = .7, size = 1.9) +
    geom_line(data = pred_line, aes(x = x, y = y),
              inherit.aes = FALSE, colour = NAVY, linetype = "dashed") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.02, vjust = 1.2,
             label = eq_str, parse = TRUE, colour = DGREY, size = 3.2) +
    scale_colour_manual(values = pal_methods, name = "Method") +
    labs(title = "Fit quality vs iterations (global quadratic trend)",
         subtitle = "Points = runs coloured by method; dashed = pooled quadratic in log space",
         x = "log(Iterations+1)", y = "log(RSS)", caption = cap_fit) + theme_nat()
  print(fig1_5)
}

# 1.6 ΔAICc surfaces over the rho-grid for the chosen (r, method)
grid_bestmethod <- results_grid %>%
  inner_join(best_methods %>% select(r, method), by = c("r","method"))
aicc_surface <- grid_bestmethod %>%
  group_by(r, method) %>%
  mutate(dAICc_grid = AICc_plusRho - min(AICc_plusRho, na.rm = TRUE)) %>%
  ungroup()
# Facet surfaces (cap ΔAICc for readability)
fig1_6 <- ggplot(aicc_surface, aes(x = rho_KL, y = rho_VAE, fill = pmin(dAICc_grid, 50))) +
  geom_tile() +
  scale_fill_viridis(option = "C", name = expression(Delta*AIC[c]), na.value = "grey90") +
  facet_wrap(~ r, scales = "free", ncol = 4) +
  labs(title = expression("ΔAICc surfaces by region for selected best method"),
       x = expression(rho[KL]), y = expression(rho[VA-E])) + theme_nat()
print(fig1_6)

# Are best grid minima on the rho edges? (edge risk)
edge_min <- aicc_surface %>%
  group_by(r, method) %>%
  slice_min(order_by = dAICc_grid, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(on_edge = (on_edge_KL | on_edge_VAE))
fig1_7 <- edge_min %>% count(on_edge) %>%
  mutate(share = n/sum(n)) %>%
  ggplot(aes(x = c("Interior","Edge")[on_edge+1], y = share)) +
  geom_col(width = .6, fill = LGREY, colour = "white") +
  geom_text(aes(label = paste0(round(share*100,1),"% (n=", n, ")")), vjust = -0.3, colour = DGREY) +
  coord_cartesian(ylim = c(0, 1.05)) +
  labs(title = "Where do ΔAICc minima sit on the grid?",
       subtitle = "Share of minima lying on grid edges vs interior (best method per region)",
       x = NULL, y = "Share") + theme_nat()
print(fig1_7)

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

# 2.2 Elasticity distributions overlaid (same x-axis), by method and overall
fig2_2 <- results_grid_valid %>%
  select(r, method, sigma_KL, sigma_VAE) %>%
  pivot_longer(c(sigma_KL, sigma_VAE), names_to = "elasticity", values_to = "sigma") %>%
  mutate(elasticity = factor(elasticity, levels = c("sigma_KL","sigma_VAE"), labels = c("σK-L","σVA–E")),
         sigma_plot = pmin(sigma, 10)) %>%
  ggplot(aes(x = sigma_plot, y = method, fill = elasticity, colour = elasticity)) +
  ggridges::stat_density_ridges(alpha = 0.5, scale = 1, rel_min_height = 0.01,
                                position = "identity", kernel = "gaussian", adjust = 1.5, n = 512) +
  scale_fill_manual(values = c("σK-L" = LGREY, "σVA–E" = NAVY)) +
  scale_colour_manual(values = c("σK-L" = DGREY, "σVA–E" = DGREY)) +
  labs(title = "Substitution elasticities across methods (valid runs)",
       subtitle = "Overlayed densities",
       x = expression(sigma), y = "Method", caption = cap_elast, fill = NULL, colour = NULL) +
  theme_nat() + theme(legend.position = "top", legend.direction = "horizontal")
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
  coord_flip() +
  scale_fill_manual(values = c("σK-L" = DGREY, "σVA–E" = NAVY)) +
  labs(title = "Median elasticities by method with IQR",
       subtitle = "Two bars per method (σK-L and σVA–E)",
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
  geom_linerange(data=elas_iqr, aes(y=factor(r, levels = ord_reg), xmin=q1, xmax=q3), size=2.2, colour=LGREY, alpha=.7) +
  geom_point(data=elas_reg, aes(x=sigma, y=factor(r, levels = ord_reg)), size=1.5, alpha=.5, colour=DGREY) +
  geom_point(data=best_long, aes(x=sigma, y=factor(r, levels = ord_reg)), shape=8, size=2.6, colour=NAVY) +
  facet_wrap(~which, scales="free_x") +
  labs(title="Elasticities by region", subtitle="Cloud = all methods; bar = IQR; star = best",
       x=expression(sigma), y="Region") + theme_nat()
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
  scale_fill_gradient(low=BBLUE, high=NAVY, na.value="grey90") +
  labs(title = expression(paste(sigma[K-L], " by region")),
       fill = expression(sigma[K-L]),
       caption = "Note: elasticities capped at 10 for readability.") +
  theme_nat()
print(fig3_6)

# 3.7 Best VA-E elasticities map
fig3_7 <- ggplot(world_best) +
  geom_sf(aes(fill = pmin(sigma_VAE,10)), colour = "white", linewidth = .1) +
  scale_fill_gradient(low=BBLUE, high=NAVY, na.value="grey90") +
  labs(title = expression(paste(sigma[VA-E], " by region")),
       fill = expression(sigma[VA-E]),
       caption = "Note: elasticities capped at 10 for readability.") +
  theme_nat()
print(fig3_7)

# 3.8 Selection certainty (Δ AICc weight best − runner-up) + distribution
aic_top2 <- aic_weights %>%
  arrange(r, desc(wAICc)) %>% group_by(r) %>% slice_head(n = 2) %>%
  summarise(best = first(wAICc), runner = dplyr::last(wAICc), delta = best - runner, .groups = "drop")
fig3_8a <- ggplot(aic_top2, aes(x = reorder(r, delta), y = delta)) +
  geom_col(width = .70, fill= LGREY, colour = "white") +
  geom_hline(yintercept = .2, colour = NAVY, linetype = "dotted", linewidth = .4) +
  coord_flip() +
  labs(title = "Selection certainty across regions",
       subtitle = "Δ weight (best – runner-up); dotted line at 0.2 ≈ clear preference",
       x = NULL, y = "Δ AICc weight") + theme_nat()
fig3_8b <- ggplot(aic_top2, aes(x = delta)) +
  geom_histogram(binwidth = .05, fill= LGREY, colour = "white") +
  labs(title = "Distribution of selection certainty (Δ weight)",
       x = "Δ AICc weight", y = "Count") + theme_nat()
print(fig3_8a); print(fig3_8b)

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

