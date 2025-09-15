options(scipen = 999) #avoids scientific notation unless necessary
# install.packages(c("micEconCES","dplyr","readr","purrr","ggplot2","parallel","ggpmisc","pheatmap"))

library(micEconCES)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(furrr)
library(parallel)
library(viridis)
library(ggplot2)
library(ggpmisc)
library(ggrepel)
library(pheatmap)
library(patchwork)
library(tibble)

### Settings
setwd("C:/Users/escami_g/OneDrive - Paul Scherrer Institut/05.Models/MERGE updates/CES-parametrisation")
infile <- "MERGE macro.csv"

# Allows to use multiple CPU cores for faster solve time
plan(multisession, workers = parallel::detectCores() - 2)

### Load data
df <- read_csv(infile, show_col_types = TRUE)

# Normalisation of data by the average of historical values
dfS <- df %>%
  group_by(r) %>%
  mutate(
    Ys = Y / mean(Y, na.rm = TRUE),
    Ks = K / mean(K, na.rm = TRUE),
    Ls = L / mean(L, na.rm = TRUE),
    Es = E / mean(E, na.rm = TRUE)
  ) %>%
  ungroup()

# Search grid for rho values for the K-L and VA-E nests
rhoGrid_KL  <- seq(-1, 10, by = 1)
rhoGrid_VAE <- seq(-1, 10, by = 1)

# Helper function to extract rho
get_rho_val <- function(fit, which = c("rho1","rho")) {
  which <- match.arg(which)
  if (is.null(fit)) return(NA_real_)
  cf <- coef(fit)
  if (which %in% names(cf)) return(as.numeric(cf[which]))
  if (!is.null(fit$allRhoSum)) {
    best <- fit$allRhoSum %>%
      filter(rss == min(rss, na.rm=TRUE)) %>% slice(1)
    return(ifelse(which=="rho1", best$rho1, best$rho))
  }
  NA_real_
}

# --- estimation per region ---
estimate_region <- function(d, region_name) {
  message("\nEstimating region: ", region_name)
  d_num <- d %>% transmute(t, Ys, Ks, Ls, Es)
  
  methods <- c("LM", "NM", "Nelder-Mead", "BFGS", "PORT", "Newton", "DE", "CG", "L-BFGS-B", "SANN")
  
  fit_all   <- setNames(vector("list", length(methods)), methods)
  conv_all  <- setNames(rep(FALSE, length(methods)), methods)
  msg_all   <- setNames(rep(NA_character_, length(methods)), methods)
  times_all <- setNames(rep(NA_real_, length(methods)), methods)
  
  for (m in methods) {
    t0 <- Sys.time()
    
    # Setup method-specific control
    control_arg <- switch(m,
                          "DE"     = DEoptim.control(itermax = 1e4),
                          "PORT"   = list(eval.max = 1e4, iter.max = 1e4, reltol = 1e-8),
                          "L-BFGS-B" = list(maxit = 1e5),
                          list(maxiter = 1e6)
    )
    
    # Setup method-specific bounds
    if (m == "L-BFGS-B") {
      lower_arg <- c(gamma=0, delta=0, rho1=-1, rho=-1, nu=0)
      upper_arg <- c(gamma=Inf, delta=1, rho1=Inf, rho=Inf, nu=Inf)
    } else {
      lower_arg <- NULL
      upper_arg <- NULL
    }
    
    # Run cesEst appropriately for Newton vs others
    fit_try <- try(
      suppressWarnings({
        if (m == "Newton") {
          cesEst(
            yName = "Ys", xNames = c("Ks","Ls","Es"),
            tName = "t", data = d_num,
            vrs = TRUE, multErr = TRUE,
            method = m,
            rho1 = rhoGrid_KL, rho = rhoGrid_VAE,
            lower = lower_arg, upper = upper_arg
            # no control here
          )
        } else {
          cesEst(
            yName = "Ys", xNames = c("Ks","Ls","Es"),
            tName = "t", data = d_num,
            vrs = TRUE, multErr = TRUE,
            method = m,
            rho1 = rhoGrid_KL, rho = rhoGrid_VAE,
            lower = lower_arg, upper = upper_arg,
            control = control_arg
          )
        }
      }),
      silent = TRUE
    )
    
    runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    times_all[m] <- runtime
    
    # Determine success & messages safely
    success   <- inherits(fit_try, "cesEst")
    conv_flag <- if (success && !is.null(fit_try$convergence)) {
      as.logical(fit_try$convergence[1])
    } else {
      FALSE
    }
    
    msg_flag <- if (inherits(fit_try, "try-error")) {
      as.character(fit_try)[1]
    } else if (!inherits(fit_try, "cesEst")) {
      "Not a cesEst object"
    } else if (!is.null(fit_try$message)) {
      as.character(fit_try$message)[1]
    } else {
      ""  # empty message
    }
    
    fit_all[[m]] <- if (success) fit_try else NULL
    conv_all[m]  <- conv_flag
    msg_all[m]   <- msg_flag
    
    # Logging
    message(
      sprintf("  %s %-10s in %.1fs%s",
              if (conv_flag) "✓" else "✗",
              m,
              runtime,
              if (!conv_flag && nzchar(msg_flag)) paste0(" msg: ", msg_flag) else ""
      )
    )
  }
  
  list(
    fits  = fit_all,
    conv  = conv_all,
    msg   = msg_all,
    times = times_all,
    data  = d_num
  )
}




# Split into regions
splits <- dfS %>% group_split(r, .keep=TRUE)
region_names <- dfS %>% distinct(r) %>% pull(r)

# Run the estimation in parallel
fits_all <- future_map2(
  splits,
  region_names,
  estimate_region,
  .progress    = TRUE,
  .options     = furrr_options(
    packages = c("micEconCES","dplyr","tidyr","purrr","readr"),
    globals  = c("rhoGrid_KL","rhoGrid_VAE","get_rho_val"),
    seed     = TRUE
  ),
  .env_globals = globalenv()
)

# --- Extract results ---
extract_region <- function(region_name, region_fits) {
  if (is.null(region_fits)) return(NULL)
  
  diag_tbl     <- tibble()
  timevary_tbl <- tibble()
  grid_tbl     <- tibble()
  
  for (m in names(region_fits$fits)) {
    fit_obj   <- region_fits$fits[[m]]
    conv_flag <- region_fits$conv[[m]] %||% FALSE
    msg_flag  <- region_fits$msg[[m]] %||% NA_character_
    runtime   <- region_fits$times[[m]] %||% NA_real_
    
    # --- Default values
    p_gamma  <- p_delta1 <- p_delta <- p_lambda <- p_nu <- c(NA, NA, NA, NA)
    rho1 <- rhoT <- sigma_KL <- sigma_VAE <- R2_val <- NA_real_
    
    if (!is.null(fit_obj) && inherits(fit_obj, "cesEst")) {
      s        <- summary(fit_obj)
      coef_mat <- coef(s)
      
      get_par <- function(par) {
        if (par %in% rownames(coef_mat)) coef_mat[par, ] else c(NA, NA, NA, NA)
      }
      
      p_gamma  <- get_par("gamma")
      p_delta1 <- get_par("delta_1")
      p_delta  <- get_par("delta")
      p_lambda <- get_par("lambda")
      p_nu     <- get_par("nu")
      
      # --- rho handling (choose best grid result if present)
      if (!is.null(fit_obj$allRhoSum)) {
        best <- fit_obj$allRhoSum %>% slice_min(rss, n = 1)
        rho1 <- best$rho1
        rhoT <- best$rho
      } else {
        rho1 <- get_rho_val(fit_obj, "rho1")
        rhoT <- get_rho_val(fit_obj, "rho")
      }
      
      sigma_KL  <- ifelse(is.finite(1/(1 + rho1)), 1/(1 + rho1), NA_real_)
      sigma_VAE <- ifelse(is.finite(1/(1 + rhoT)), 1/(1 + rhoT), NA_real_)
      
      obs    <- log(region_fits$data$Ys + 1e-12)
      fitv   <- log(fit_obj$fitted.values + 1e-12)
      R2_val <- tryCatch(
        1 - sum((obs - fitv)^2)/sum((obs - mean(obs))^2),
        error = function(e) NA_real_
      )
      
      # --- time-varying
      timevary_tbl <- bind_rows(timevary_tbl, tibble(
        r        = region_name,
        method   = m,
        t        = region_fits$data$t,
        fitted   = as.numeric(fit_obj$fitted.values),
        residual = as.numeric(fit_obj$residuals),
        TFP      = p_gamma["Estimate"] * exp(p_lambda["Estimate"] * region_fits$data$t)
      ))
      
      # --- keep full grid search separately
      if (!is.null(fit_obj$allRhoSum)) {
        grid_tbl <- bind_rows(grid_tbl, fit_obj$allRhoSum %>% mutate(r = region_name, method = m))
      }
    }
    
    # --- always one row per region-method
    diag_tbl <- bind_rows(diag_tbl, tibble(
      r             = region_name,
      method        = m,
      rss           = if (!is.null(fit_obj) && inherits(fit_obj, "cesEst")) fit_obj$rss else NA_real_,
      R2            = R2_val,
      iter          = if (!is.null(fit_obj) && inherits(fit_obj, "cesEst")) fit_obj$iter else NA_integer_,
      conv          = conv_flag,
      message       = msg_flag,
      gamma         = p_gamma["Estimate"],
      se_gamma      = p_gamma["Std. Error"],
      t_gamma       = p_gamma["t value"],
      p_gamma       = p_gamma["Pr(>|t|)"],
      lambda        = p_lambda["Estimate"],
      se_lambda     = p_lambda["Std. Error"],
      t_lambda      = p_lambda["t value"],
      p_lambda      = p_lambda["Pr(>|t|)"],
      d_KL          = p_delta1["Estimate"],
      se_delta1     = p_delta1["Std. Error"],
      t_delta1      = p_delta1["t value"],
      p_delta1      = p_delta1["Pr(>|t|)"],
      d_VAE         = p_delta["Estimate"],
      se_delta      = p_delta["Std. Error"],
      t_delta       = p_delta["t value"],
      p_delta       = p_delta["Pr(>|t|)"],
      nu            = p_nu["Estimate"],
      se_nu         = p_nu["Std. Error"],
      t_nu          = p_nu["t value"],
      p_nu          = p_nu["Pr(>|t|)"],
      rho_KL        = rho1,
      rho_VAE       = rhoT,
      sigma_KL      = sigma_KL,
      sigma_VAE     = sigma_VAE,
      runtime       = runtime
    ))
  }
  
  list(diag = diag_tbl, timevary = timevary_tbl, grid = grid_tbl)
}


# combine results across regions
results <- map2(region_names, fits_all, extract_region)

# --- Build tables ---
results_table <- bind_rows(map(results,"diag"))
results_table_time <- bind_rows(map(results,"timevary"))
results_grid <- bind_rows(compact(map(results,"grid")))

# Best method per region by lowest RSS
best_methods <- results_table %>%
  group_by(r) %>%
  slice_min(order_by = rss, n = 1, with_ties = FALSE) %>%
  ungroup()

# Optimisation convergence summary for all methods
all_methods <- c("Kmenta","LM","NM","Nelder-Mead","PORT","BFGS","Newton","DE","CG","L-BFGS-B","SANN")

convergence_summary <- results_table %>%
  select(r, method, conv) %>%
  distinct() %>%                           # ensure one entry per region-method
  group_by(method) %>%
  summarise(
    n_tried   = n(),
    n_success = sum(conv, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(n_fail = n_tried - n_success) %>%
  pivot_longer(c(n_success, n_fail), names_to = "status", values_to = "count") %>%
  mutate(status = ifelse(status == "n_success", "Converged", "Failed"))


# --- IAM table ---
if (!"method" %in% names(results_table_time) || nrow(results_table_time) == 0) {
  stop("results_table_time is empty or missing 'method' column — cannot build IAM table.")
} else {
iam_table <- results_table_time %>%
  filter(method %in% best_methods$method & r %in% best_methods$r) %>%
  left_join(
    best_methods %>% select(r, method, sigma_KL, sigma_VAE, rho_KL, rho_VAE, d_KL, d_VAE, gamma, lambda),
    by = c("r", "method")
  ) %>%
  transmute(
    year                        = t,
    region                      = r,
    total_factor_productivity   = TFP,
    share_capital_valueadded    = d_KL,
    share_valueadded_output     = d_VAE,
    substitution_capital_labour = sigma_KL,
    substitution_valueadded_energy = sigma_VAE,
    valueadded_CES_exponent     = rho_KL / (rho_KL + 1),
    output_CES_exponent         = rho_VAE / (rho_VAE + 1),
    capital_productivity_factor = 1,
    labour_productivity_factor  = 1
  )
}

# --- Save CSVs ---
write_csv(results_table, "CES_region_method.csv")
write_csv(results_table_time, "CES_region_method_year.csv")
write_csv(results_grid, "CES_gridsearch.csv")
write_csv(iam_table, "IAM_params.csv")


# --- Plots ---
# =====================================================================
# PLOTS FOR SCIENTIFIC OUTPUT
# =====================================================================

# Histogram of σK-L
ggplot(filter(results_table, is.finite(sigma_KL)),
       aes(x = sigma_KL, fill = r)) +
  geom_histogram(binwidth = 0.1, color = "white", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Distribution of σK-L",
       x = "σK-L (Capital–Labour Elasticity)", y = "Count") +
  theme(legend.position = "bottom")

# Histogram of σVA-E
ggplot(filter(results_table, is.finite(sigma_VAE)),
       aes(x = sigma_VAE, fill = r)) +
  geom_histogram(binwidth = 0.1, color = "white", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Distribution of σVA-E",
       x = "σVA-E (Value-Added–Energy Elasticity)", y = "Count") +
  theme(legend.position = "bottom")

# Violin + box + mean/median labels for σK-L
ggplot(results_table, aes(x = method, y = sigma_KL, fill = method)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
  stat_summary(fun = median, geom = "point", shape = 21, size = 3, color = "black", fill = "white") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, color = "red", fill = "white") +
  stat_summary(fun = median, geom = "text", aes(label = sprintf("%.2f", ..y..)), color = "black", vjust = -1, size = 3) +
  stat_summary(fun = mean, geom = "text", aes(label = sprintf("%.2f", ..y..)), color = "red", vjust = 1.5, size = 3) +
  theme_minimal() +
  labs(title = "σK-L by Method (median = black, mean = red)", x = "Method", y = expression(sigma[K-L])) +
  theme(legend.position = "none")

# Violin + box + mean/median labels for σVA-E
ggplot(results_table, aes(x = method, y = sigma_VAE, fill = method)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
  stat_summary(fun = median, geom = "point", shape = 21, size = 3, color = "black", fill = "white") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, color = "red", fill = "white") +
  stat_summary(fun = median, geom = "text", aes(label = sprintf("%.2f", ..y..)), color = "black", vjust = -1, size = 3) +
  stat_summary(fun = mean, geom = "text", aes(label = sprintf("%.2f", ..y..)), color = "red", vjust = 1.5, size = 3) +
  theme_minimal() +
  labs(title = "σVA-E by Method (median = black, mean = red)", x = "Method", y = expression(sigma[VAE])) +
  theme(legend.position = "none")

# Convergence rate by method
ggplot(convergence_summary, aes(x = method, y = count, fill = status)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Failed" = "red", "Converged" = "darkgreen")) +
  theme_minimal() +
  labs(title = "Convergence by Method", x = "Method", y = "Number of Regions")

# Scatter of σK-L vs σVA-E
library(ggpmisc)

median_KL  <- median(results_table$sigma_KL, na.rm = TRUE)
median_VAE <- median(results_table$sigma_VAE, na.rm = TRUE)

ggplot(
  filter(results_table, is.finite(sigma_KL) & is.finite(sigma_VAE)),
  aes(x = sigma_KL, y = sigma_VAE, color = r)
) +
  geom_quadrant_lines(
    xintercept = median_KL,
    yintercept = median_VAE,
    colour = "grey50"
  ) +
  stat_quadrant_counts(
    xintercept = median_KL,
    yintercept = median_VAE,
    aes(label = after_stat(paste0("n=", count)), colour = NULL),
    colour = "grey20"
  ) +
  geom_point(size = 2, alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Elasticity Quadrants: σK-L vs σVA-E",
    subtitle = "Regions by their substitutabilities",
    x = expression(sigma[K-L]),
    y = expression(sigma[VAE])
  )

# RSS vs Iterations tradeoff
ggplot(results_table, aes(x = iter, y = rss, color = method)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_text_repel(aes(label = r), size = 2, max.overlaps = 10) +
  scale_y_log10() +
  theme_classic() +
  labs(title = "Tradeoff: Iterations vs Fit Quality (RSS)",
       x = "Iterations", y = "Residual Sum of Squares (log scale)")

# Grid RSS landscapes
if (nrow(results_grid) > 0) {
  ggplot(results_grid, aes(x = rho1, y = rho, z = rss, fill = rss)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_viridis_c(trans = "log") +
    facet_wrap(~r, labeller = label_wrap_gen(width=20)) + # cleaner region names
    theme_minimal() +
    labs(title = "RSS landscape (ρK-L vs ρVA-E)",
         x = "ρK-L", y = "ρVA-E", fill = "RSS (log)")
}

# Distribution of coefficients γ and λ
ggplot(results_table, aes(x = gamma)) +
  geom_histogram(binwidth = 0.05, fill = "purple", color = "white") +
  theme_classic() +
  labs(title = "Distribution of γ (TFP intercept)", x = "γ", y = "Count")

ggplot(results_table, aes(x = lambda)) +
  geom_histogram(binwidth = 0.01, fill = "green4", color = "white") +
  theme_classic() +
  labs(title = "Distribution of λ (TFP growth rate)", x = "λ", y = "Count")

# Distribution of R² across region-method
ggplot(results_table, aes(x = reorder(r, R2), y = R2, color = method)) +
  geom_point(size = 2) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "R² by Region and Method",
    x = "Region (sorted by R²)",
    y = expression(R^2),
    color = "Method"
  )

ggplot(results_table, aes(x = iter, y = R2, color = method)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Iterations vs R² by Method",
    x = "Iterations",
    y = expression(R^2),
    color = "Method"
  )

# P-value significance for dK-L
ggplot(
  filter(results_table, !is.na(p_delta1) & conv == TRUE),
  aes(x = method, y = d_KL, fill = (p_delta1 < 0.05))
) +
  geom_boxplot() +
  scale_fill_manual(
    values = c("FALSE" = "grey80", "TRUE" = "steelblue"),
    labels = c("Not significant", "p < 0.05")
  ) +
  theme_minimal() +
  labs(
    title = "dKL Across Methods with Significance Highlighted",
    x = "Method", y = expression(delta[K-L]), fill = "Significance"
  )

# P-value significance for dVA-E
ggplot(
  filter(results_table, !is.na(p_delta) & conv == TRUE),
  aes(x = method, y = d_VAE, fill = (p_delta < 0.05))
) +
  geom_boxplot() +
  scale_fill_manual(
    values = c("FALSE" = "grey80", "TRUE" = "steelblue"),
    labels = c("Not significant", "p < 0.05")
  ) +
  theme_minimal() +
  labs(
    title = "dVA-E Across Methods with Significance Highlighted",
    x = "Method", y = expression(delta[VA-E]), fill = "Significance"
  )

# Convergence Heatmap (Region × Method)
conv_df <- results_table %>%
  distinct(r, method, conv)

df_wide <- conv_df %>%
  mutate(conv = as.integer(conv)) %>%
  pivot_wider(
    names_from = method,
    values_from = conv,
    values_fill = list(conv = 0)
  )


heatmap_matrix <- conv_df %>%
  mutate(conv = as.integer(conv)) %>%
  pivot_wider(names_from = method, values_from = conv, values_fill = list(conv = 0)) %>%
  column_to_rownames("r") %>%
  as.matrix()

p_conv_heatmap <- pheatmap(
  heatmap_matrix,
  color = c("red","green"),
  main = "Convergence Heatmap\n(1 = converged, 0 = failed)",
  silent = TRUE
)

p_heatmap_grob <- p_conv_heatmap$gtable

# Proportion Converged by Method
p_conv_prop <- conv_df %>%
  distinct(r, method, conv) %>%
  group_by(method) %>%
  summarise(prop_converged = mean(conv, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = method, y = prop_converged, fill = method)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Proportion of Regions that Converged by Method",
    x = "", y = "Proportion"
  )

# Elasticity Distributions (σ_K-L, σ_VAE)
# Reference lines at σ = 0.5 and σ = 1
p_sigma_KL <- results_table %>%
  filter(is.finite(sigma_KL)) %>%
  ggplot(aes(x = method, y = sigma_KL, fill = method)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, show.legend = FALSE, outlier.shape = NA) +
  geom_hline(yintercept = c(0.5, 1), linetype = "dashed", color = c("blue","red")) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = expression("σ"[K-L] * " Distribution by Method"),
    x = "",
    y = expression(sigma[K-L])
  )

p_sigma_VAE <- results_table %>%
  filter(is.finite(sigma_VAE)) %>%
  ggplot(aes(x = method, y = sigma_VAE, fill = method)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, show.legend = FALSE, outlier.shape = NA) +
  geom_hline(yintercept = c(0.5, 1), linetype = "dashed", color = c("blue","red")) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = expression("σ"[VAE] * " Distribution by Method"),
    x = "",
    y = expression(sigma[VAE])
  )

# Parameter Distributions
p_params <- results_table %>%
  pivot_longer(cols = c("gamma", "lambda", "d_KL", "d_VAE", "nu"),
               names_to = "param", values_to = "estimate") %>%
  filter(!is.na(estimate)) %>%
  ggplot(aes(x = param, y = estimate, fill = param)) +
  geom_boxplot(show.legend = FALSE, outlier.size = 1) +
  theme_minimal() +
  labs(
    title = "Distributions of Key Parameter Estimates Across Regions",
    x = "Parameter", y = "Estimated Value"
  )

# Significance Frequency (p-value for d_KL)
p_signif_d_KL <- results_table %>%
  mutate(significance = case_when(
    is.na(p_delta1)          ~ "NA",
    p_delta1 < 0.05          ~ "< 0.05",
    TRUE                     ~ "≥ 0.05"
  )) %>%
  group_by(method, significance) %>%
  summarise(count = n(), .groups = "drop") %>%
  ggplot(aes(x = method, y = count, fill = significance)) +
  geom_col(position = "stack") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Count of Significance Levels for d_KL by Method",
    x = "Method", y = "Number of Fits", fill = "Significance"
  )

# RSS vs Iterations scatter with labels
p_rss_iter <- results_table %>%
  mutate(label = ifelse(conv, r, NA_character_)) %>%
  ggplot(aes(x = iter, y = rss, color = method)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(aes(label = label),
                  size = 2, max.overlaps = 15,
                  na.rm = TRUE) +
  scale_y_log10() +
  theme_classic() +
  labs(
    title = "RSS vs Iterations (Log RSS) — labels for converged regions",
    x = "Iterations", y = "RSS (log scale)"
  )

# H. Grid Search RSS Surfaces (sample regions/methods with grid data)
grid_samp <- results_grid %>%
  group_by(r, method) %>%
  filter(n() > 50) %>%
  slice(1:4)  # up to four region-method pairs

p_grid_rss <- grid_samp %>%
  ggplot(aes(x = rho1, y = rho, z = rss)) +
  geom_contour_filled() +
  facet_grid(r ~ method, scales = "free") +
  theme_minimal() +
  labs(
    title = "RSS Contour Surfaces (rho1 vs rho) for Sample Region-Method",
    x = expression(rho[1]), y = expression(rho)
  )



# K. Elasticity by Region (using best_methods)
p_sigma_best_KL <- best_methods %>%
  ggplot(aes(x = reorder(r, sigma_KL), y = sigma_KL)) +
  geom_col(fill = "darkblue") +
  geom_hline(yintercept = c(0.5, 1), linetype = "dashed", color = c("blue","red")) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "σ_K-L (Capital-Labour Elasticity) by Region (Best Fit Method)",
    x = "Region", y = expression(sigma[K-L])
  )

p_best_method_count <- best_methods %>%
  count(method) %>%
  ggplot(aes(x = method, y = n, fill = method)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Number of Regions best fit by Method",
    x = "Method", y = "Count of Regions"
  )

# Combine & Display Plots (in panels)
message("Plotting combined panels...")

# Top diagnostics
(p_conv_prop + p_iter + p_runtime) / (p_sigma_KL + p_sigma_VAE + p_params) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Diagnostics Panel")

# RSS vs Iterations
print(p_rss_iter)

# Grid RSS surfaces
print(p_grid_rss)

# Summary elasticity and best-method
(p_sigma_best_KL + p_best_method_count) +
  plot_annotation(title = "Elasticity Summary & Best-method per Region")

# Significance frequencies
print(p_signif_d_KL)

