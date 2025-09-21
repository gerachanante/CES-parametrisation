options(scipen = 999) #avoids scientific notation unless necessary
# install.packages(c("micEconCES","dplyr","readr","purrr","ggplot2","parallel","ggpmisc","pheatmap","GGally","ggcorrplot","ggridges"))
#install.packages("ggridges")
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


### Settings
setwd("C:/Users/escami_g/OneDrive - Paul Scherrer Institut/05.Models/MERGE updates/CES-parametrisation/-0.99 to -0.95 by 0.04 and -0.95 to 5 by 0.05 and 6 to 80 by 1")



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

### Diagnostics of all runs
# Convergence by method
convergence_summary %>%
  ggplot(aes(x = method, y = count, fill = status)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = count), position = position_dodge(width=0.9),
            vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("Failed"="red","Converged"="darkgreen")) +
  theme_minimal(base_size=12) +
  labs(title = "Convergence by Method",
       y = "Number of Region-Method Runs", x = "Method")

# Proportion converged
results_table %>%
  distinct(r, method, conv) %>%
  group_by(method) %>%
  summarise(prop_converged = mean(conv, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = method, y = prop_converged, fill = method)) +
  geom_col(show.legend = FALSE) +
  coord_flip() + theme_minimal(base_size=12) +
  labs(title = "Proportion of Regions Converged by Method", x = "", y = "Proportion")

# Region-method convergence heatmap
conv_df <- results_table %>% distinct(r, method, conv)
heatmap_matrix <- conv_df %>%
  mutate(conv = as.integer(conv)) %>%
  pivot_wider(names_from = method, values_from = conv,
              values_fill = list(conv = 0)) %>%
  column_to_rownames("r") %>% as.matrix()
pheatmap(heatmap_matrix, color = c("red","green4"),
         main = "Convergence Heatmap\n(1 = converged, 0 = failed)")

# Grid convergence share
if (nrow(grid_conv_share) > 0) {
  grid_conv_share %>%
    ggplot(aes(method, share_converged, fill = method)) +
    geom_boxplot(show.legend = FALSE, outlier.size = 0.7) +
    theme_minimal(base_size = 12) +
    labs(title = "Share of Grid Converged (by Method)",
         x = "", y = "Share of ρ×ρ1 grid")
}

# Runtime distribution
ggplot(results_table, aes(x = method, y = runtime, fill = method)) +
  geom_boxplot(show.legend = FALSE, outlier.size = 0.5) +
  scale_y_log10() +
  theme_minimal(base_size = 12) +
  labs(title = "Runtime Distribution by Method",
       x = "Method", y = "Runtime (log-seconds)")

# RSS vs Iterations
ggplot(results_table, aes(x = iter, y = rss, color = method)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_smooth(method="lm", se=FALSE, formula=y~poly(log(x),2),
              linetype="dashed") +
  geom_text_repel(aes(label = r), size = 2, max.overlaps = 10) +
  scale_y_log10() + scale_x_log10() +
  theme_classic(base_size = 12) +
  labs(title   = "Iterations vs Fit Quality",
       x = "Iterations (log)", y = "RSS (log)")

### Valid runs diagnostics
# Correlation matrix
cor_df  <- results_table_valid %>%
  select(gamma, lambda, delta_KL, delta_VAE, nu, sigma_KL, sigma_VAE)
cor_df  <- cor_df[, colSums(!is.na(cor_df)) > 0, drop = FALSE]
cor_mat <- cor(cor_df, use = "pairwise.complete.obs")
if (ncol(cor_mat) >= 2) {
  pheatmap(cor_mat, main = "Parameter Correlations")
  ggcorrplot(cor_mat, hc.order = TRUE, type = "lower", lab = TRUE)
  ggpairs(cor_df)
}

# Histograms of γ and λ
ggplot(results_table_valid, aes(x = gamma)) +
  geom_histogram(binwidth = 0.4, fill = "purple", color = "white") +
  theme_classic(base_size=12) +
  labs(title="Distribution of γ (TFP intercept)", x=expression(gamma), y="Count")

ggplot(results_table_valid, aes(x = lambda)) +
  geom_histogram(binwidth = 0.002, fill = "green4", color = "white") +
  theme_classic(base_size=12) +
  labs(title="Distribution of λ (TFP growth rate)", x=expression(lambda), y="Count")

# Boxplot of parameters
results_table_valid %>%
  pivot_longer(cols = c(gamma, lambda, delta_KL, delta_VAE, nu),
               names_to = "param", values_to = "estimate") %>%
  filter(!is.na(estimate)) %>%
  ggplot(aes(x = param, y = estimate, fill = param)) +
  geom_boxplot(show.legend = FALSE, outlier.size = 1) +
  theme_minimal(base_size=12) +
  labs(title = "Distributions of Parameter Estimates", x = "Parameter", y = "Value")

# Elasticity distributions (σKL and σVAE, clipped for readability)
ggplot(filter(results_table_valid, is.finite(sigma_KL), sigma_KL < 10),
       aes(x = sigma_KL)) +
  geom_histogram(binwidth = 0.2, fill = "steelblue", color = "white") +
  theme_minimal(base_size=12) +
  labs(title="Distribution of σ[K-L] (clipped at 10)",
       x=expression(sigma[K-L]), y="Count")

ggplot(filter(results_table_valid, is.finite(sigma_VAE), sigma_VAE < 10),
       aes(x = sigma_VAE)) +
  geom_histogram(binwidth = 0.2, fill = "darkorange", color = "white") +
  theme_minimal(base_size=12) +
  labs(title="Distribution of σ[VA-E] (clipped at 10)",
       x=expression(sigma[VA-E]), y="Count")

# Elasticity density ridges by method
results_table_valid %>%
  pivot_longer(c(sigma_KL, sigma_VAE), names_to="elasticity", values_to="sigma") %>%
  ggplot(aes(x = sigma, y = method, fill = elasticity)) +
  geom_density_ridges(alpha=0.5) +
  facet_wrap(~elasticity, scales="free_x") +
  theme_minimal(base_size=12) +
  labs(title = "Distributions of Elasticities by Method",
       x = expression(sigma), y = "Method")

# Scatter of σKL vs σVAE
ggplot(filter(results_table_valid, is.finite(sigma_KL), is.finite(sigma_VAE)),
       aes(x = sigma_KL, y = sigma_VAE, color = method)) +
  geom_point(alpha=0.6) +
  theme_minimal(base_size=12) +
  labs(title = "Valid Runs: σ[K-L] vs σ[VA-E]",
       x = expression(sigma[K-L]), y = expression(sigma[VA-E]))

# Faceted histograms for σ[K-L]
ggplot(filter(results_table_valid, sigma_KL < 10),
       aes(x = sigma_KL, fill = method)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  facet_wrap(~r, scales = "free_y") +
  theme_minimal(base_size=12) +
  labs(title = expression("Distribution of σ[K-L] across Regions"),
       x = expression(sigma[K-L]), y = "Count")

# Ridge plot of γ across regions
results_table_valid %>%
  ggplot(aes(x = gamma, y = r, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01) +
  scale_fill_viridis_c() +
  theme_minimal(base_size=12) +
  labs(title = "Distribution of γ by Region",
       x = expression(gamma), y = "Region")

# TFP trajectories
results_table_time %>%
  inner_join(best_methods %>% select(r, method), by = c("r","method")) %>%
  ggplot(aes(x = t, y = TFP, color = method)) +
  geom_line() +
  facet_wrap(~r, scales = "free_y") +
  theme_minimal(base_size=12) +
  labs(title = "TFP Trajectories (Best Methods per Region)",
       x = "Year", y = "TFP")

# Residual distribution across all valid fits
results_table_time %>%
  ggplot(aes(x = residual, fill = method)) +
  geom_density(alpha = 0.4) +
  theme_minimal(base_size=12) +
  labs(title = "Distribution of Residuals (log-scale fit)",
       x = "Residual (log)", y = "Density")

### Best methods
# Count of best methods
best_methods %>%
  count(method) %>%
  ggplot(aes(x = method, y = n, fill = method)) +
  geom_col(show.legend = FALSE) + coord_flip() +
  theme_minimal(base_size=12) +
  labs(title = "Number of Regions Best Fit by Each Method",
       x = "Method", y = "Count of Regions")

# Elasticities by region
best_methods %>%
  ggplot(aes(x = reorder(r, sigma_KL), y = sigma_KL, fill = conv)) +
  geom_col() +
  scale_fill_manual(values=c("TRUE"="lightblue4","FALSE"="red")) +
  geom_hline(yintercept=c(0.5,1), linetype="dashed", color=c("blue","red")) +
  coord_flip() + theme_minimal(base_size=12) +
  labs(title="σ[K-L] by Region (Best Method)",
       x="Region", y=expression(sigma[K-L]))

best_methods %>%
  ggplot(aes(x = reorder(r, sigma_VAE), y = sigma_VAE, fill = conv)) +
  geom_col() +
  scale_fill_manual(values=c("TRUE"="lightblue4","FALSE"="red")) +
  geom_hline(yintercept=c(0.5,1), linetype="dashed", color=c("blue","red")) +
  coord_flip() + theme_minimal(base_size=12) +
  labs(title="σ[VA-E] by Region (Best Method)",
       x="Region", y=expression(sigma[VA-E]))

# σ[K-L] map
ggplot(world_best) +
  geom_sf(aes(fill = pmin(sigma_KL, 10))) +   # clip at 10
  scale_fill_viridis_c(option = "plasma", na.value = "grey90") +
  theme_minimal() +
  labs(title = expression("σ[K-L] by Region (Best Method)"),
       fill = expression(sigma[K-L]))

# σ[VA-E] map
ggplot(world_best) +
  geom_sf(aes(fill = pmin(sigma_VAE, 10))) +
  scale_fill_viridis_c(option = "viridis", na.value = "grey90") +
  theme_minimal() +
  labs(title = expression("σ[VA-E] by Region (Best Method)"),
       fill = expression(sigma[VA-E]))

# Quadrant analysis (one point per region)
median_KL  <- median(best_methods$sigma_KL, na.rm=TRUE)
median_VAE <- median(best_methods$sigma_VAE, na.rm=TRUE)
ggplot(filter(best_methods, is.finite(sigma_KL), is.finite(sigma_VAE)),
       aes(x = sigma_KL, y = sigma_VAE, color = r)) +
  geom_hline(yintercept = median_VAE, linetype="dashed", alpha=.5) +
  geom_vline(xintercept = median_KL,  linetype="dashed", alpha=.5) +
  geom_point(size=3, alpha=0.8) +
  theme_minimal(base_size=12) +
  labs(title="Elasticity Quadrants (Best Method)",
       x=expression(sigma[K-L]), y=expression(sigma[VA-E]))

### Insights from IAM export
# IAM parameter table preview
head(iam_table)

# RSS ridge plots
results_grid %>%
  filter(r %in% unique(best_methods$r)) %>%
  ggplot(aes(x = rho1, y = rho, z = rss)) +
  geom_contour_filled() +
  facet_grid(r ~ method) +
  scale_fill_viridis_d() +
  theme_minimal(base_size=12) +
  labs(title = "RSS Contour Surfaces by Region & Method",
       x = expression(rho[1]), y = expression(rho), fill = "RSS")

iam_table %>%
  ggplot(aes(x = year, y = total_factor_productivity, color = region)) +
  geom_line() +
  theme_minimal(base_size=12) +
  labs(title = expression("IAM Input: TFP over Time"),
       x = "Year", y = "TFP")


# AICc weights distribution (among valid runs)
aic_weights %>%
  filter(is.finite(wAICc)) %>%
  ggplot(aes(method, wAICc, fill = method)) +
  geom_violin(alpha=.4, show.legend=FALSE) +
  geom_boxplot(width=.2, outlier.size=.5, show.legend=FALSE) +
  theme_minimal(base_size=12) +
  labs(title="AICc Weights by Method (Valid Runs Only)",
       x="Method", y="AICc weight")