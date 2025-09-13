# install.packages(c("micEconCES","dplyr","readr","purrr","ggplot2","parallel"))
library(micEconCES)
library(dplyr)
library(readr)
library(purrr)
library(ggplot2)
library(furrr)
library(parallel)

# --- Settings ---
setwd("C:/Users/escami_g/OneDrive - Paul Scherrer Institut/05.Models/MERGE updates/CES-parametrisation")
infile    <- "MERGE macro.csv"

# Allows R to work with multiple sessions using in parallel the CPU cores (faster runtime but less details while running)
plan(multisession, workers = parallel::detectCores() - 2)

# --- Load ---
df <- read_csv(infile, show_col_types = TRUE)

# --- Scale for stability ---
dfS <- df %>%
  group_by(r) %>%
  mutate(
    Ys = Y / mean(Y, na.rm = TRUE),
    Ks = K / mean(K, na.rm = TRUE),
    Ls = L / mean(L, na.rm = TRUE),
    Es = E / mean(E, na.rm = TRUE)
  ) %>%
  ungroup()

# --- Grid for rho values ---
rhoGrid <- seq(-1, 5, by = .05)

# --- helper: extract rho ---
get_rho_val <- function(fit, which = c("rho1","rho")) {
  which <- match.arg(which)
  if (is.null(fit)) return(NA_real_)
  cf <- coef(fit)
  if (!is.null(cf[which])) return(as.numeric(cf[which]))
  if (!is.null(fit$rssArray)) {
    idx <- which(fit$rssArray == min(fit$rssArray), arr.ind = TRUE)
    if (which=="rho1") return(fit$rho1Values[idx[1,1]])
    else               return(fit$rhoValues [idx[1,2]])
  }
  NA_real_
}

# --- estimation per region ---
estimate_region <- function(d, region_name) {
  message("Estimating: ", region_name)
  
  d_num <- d %>%
    transmute(
      t  = as.numeric(t),
      Ys = Y / mean(Y, na.rm=TRUE),
      Ks = K / mean(K, na.rm=TRUE),
      Ls = L / mean(L, na.rm=TRUE),
      Es = E / mean(E, na.rm=TRUE)
    )
  
  # methods to try
  methods <- c("Kmenta", "LM", "NM", "Nelder-Mead", "PORT", "BFGS", "Newton", "DE", "NM", "CG", "L-BFGS-B")
  fit_all <- list()
  times <- list()
  
  for (m in methods) {
    t0 <- Sys.time()
    fit_try <- try(
      cesEst(
        yName   = "Ys",
        xNames  = c("Ks","Ls","Es"),
        tName   = "t",
        data    = d_num,
        vrs     = TRUE,
        multErr = TRUE,
        method  = m,
        rho1    = rhoGrid,
        rho     = rhoGrid,
        control = list(maxiter = 1000),
        upper   = c(delta_1 = 1, delta = 1),
        lower   = c(delta_1 = 0, delta = 0)
      ),
      silent = TRUE
    )
    t1 <- Sys.time()
    times[[m]] <- as.numeric(difftime(t1, t0, units = "secs"))
    
    if (!(inherits(fit_try, "try-error") || is.null(fit_try))) {
      fit_all[[m]] <- fit_try
      used_method <- m
      message("✓ Success for ", region_name, " with method: ", m, " (", round(times[[m]],2),"s)")
    } else {
      message("✗ Failed with method: ", m)
    }
  }
  
  if (length(fit_all) == 0) {
    message("✗✗ All methods failed for ", region_name)
    return(NULL)
  } 

# After running all methods, choose the best fit
  rss_vals <- sapply(fit_all, function(f) f$rss)
  best_m <- names(which.min(rss_vals))
  best_fit <- fit_all[[best_m]]
  
  attr(best_fit, "method used") <- best_m
  attr(best_fit, "time used") <- times[[best_m]]
  attr(best_fit, "all methods") <- names(fit_all)
  attr(best_fit, "rss_all")  <- rss_vals
  attr(best_fit, "time_all") <- times
  
  return(best_fit)
}

# --- run for all regions ---
# split the data by region and assign region names
splits       <- dfS %>% group_split(r, .keep = TRUE)
region_names <- dfS %>% distinct(r) %>% pull(r)

# Run the CES estimation function for each region and estimation method
fits <- future_map2(splits, region_names, estimate_region, .progress = TRUE)


# --- Extract results ---
extract_fit <- function(fit, region) {
  if (is.null(fit)) {
    return(list(
      summary    = tibble(r = region, converged = FALSE),
      residuals  = NULL,
      fitted     = NULL,
      grid       = NULL
    ))
  }
  
  cf <- coef(fit)
  se <- try(sqrt(diag(fit$cov.unscaled)), silent = TRUE)
  se <- if(inherits(se, "try-error")) rep(NA, length(cf)) else se
  names(se) <- names(cf)
  
  rho1 <- get_rho_val(fit,"rho1")
  rhoT <- get_rho_val(fit,"rho")
  
  # compute R² (explained variance in log space)
  R2_val <- tryCatch(
    1 - sum(fit$residuals^2) / sum((log(fit$fitted.values + 1e-12) -
                                      mean(log(fit$fitted.values + 1e-12)))^2),
    error = function(e) NA
  )
  
  # compute t-stats and p-values
  t_stats <- cf / se
  p_vals <- 2 * (1 - pnorm(abs(t_stats)))
  
  summary_tbl <- tibble(
    r         = region,
    method    = attr(fit,"method used"),
    time      = attr(fit,"time used"),
    converged = fit$convergence,
    iter      = fit$iter,
    rss       = fit$rss,
    R2        = R2_val,
    
    gamma     = cf["gamma"],
    lambda    = cf["lambda"],
    d_KL      = pmin(pmax(cf["delta_1"], 0), 1),
    d_VAE     = pmin(pmax(cf["delta"], 0), 1),
    nu        = if ("nu" %in% names(cf)) cf["nu"] else NA,
    rho_KL    = rho1,
    rho_VAE   = rhoT,
    sigma_KL  = ifelse(is.finite(1/(1+rho1)), 1/(1+rho1), NA),
    sigma_VAE = ifelse(is.finite(1/(1+rhoT)), 1/(1+rhoT), NA),
    
    se_gamma  = se["gamma"],
    se_lambda = se["lambda"],
    se_delta1 = se["delta_1"],
    se_delta  = se["delta"],
    se_nu     = if ("nu" %in% names(se)) se["nu"] else NA,
    
    p_gamma   = p_vals["gamma"],
    p_lambda  = p_vals["lambda"],
    p_delta1  = p_vals["delta_1"],
    p_delta   = p_vals["delta"],
    p_nu      = if ("nu" %in% names(p_vals)) p_vals["nu"] else NA
  )

  diagnostics_tbl <- tibble(
    r      = region,
    method = names(attr(fit, "rss_all")),
    rss    = unlist(attr(fit, "rss_all")),
    time   = unlist(attr(fit, "time_all"))
  )
  
  residuals_tbl <- tibble(
    r = region,
    t = seq_along(fit$residuals),
    residuals = as.numeric(fit$residuals)
  )
  
  fitted_tbl <- tibble(
    r = region,
    t = seq_along(fit$fitted.values),
    fitted = as.numeric(fit$fitted.values)
  )
  
  tfp_tbl <- tibble(
    r = region,
    t = seq_along(fit$fitted.values),
    TFP = as.numeric(cf["gamma"]) * exp(as.numeric(cf["lambda"]) * seq_along(fit$fitted.values))
  )
  
  grid_tbl <- if (!is.null(fit$allRhoSum)) {
    fit$allRhoSum %>% mutate(r = region)
  } else {
    NULL
  }
  
  return(list(
    summary    = summary_tbl,
    diagnostics= diagnostics_tbl,
    residuals  = residuals_tbl,
    fitted     = fitted_tbl,
    grid       = grid_tbl
  ))
}


# Apply to all fits
results <- map2(fits, region_names, extract_fit)

# Combine summaries into one CSV
res_summary <- bind_rows(map(results, "summary"))
residuals_all <- bind_rows(map(results, "residuals"))
fitted_all    <- bind_rows(map(results, "fitted"))
grid_all      <- bind_rows(compact(map(results, "grid")))
IAM_params <- res_summary %>%
  select(r, d_KL, d_VAE, sigma_KL, sigma_VAE, nu, gamma, lambda, R2, method, time)

# Export CSV files
write_csv(res_summary, "CES_summary.csv")
write_csv(residuals_all, "CES_residuals.csv")
write_csv(fitted_all, "CES_fitted.csv")
write_csv(grid_all, "CES_gridsearch.csv")
write_csv(IAM_params, "IAM_parameters.csv")

# Build graphs from results
# 1. Histogram of elasticities of substitution
ggplot(res_summary, aes(x = sigma_VAE)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = "white") +
  theme_minimal() +
  labs(title = "Distribution of Energy–VA Substitution Elasticities (σ_KLE)",
       x = "σ_KLE", y = "Count")

ggplot(res_summary, aes(x = as.numeric(sigma_KL))) +
  geom_histogram(binwidth = 0.1, fill = "darkorange", color = "white") +
  theme_minimal() +
  labs(title = "Distribution of Capital–Labor Substitution Elasticities (σ_KL)",
       x = "σ_KL", y = "Count")

# 2. Regional substitutability profile
ggplot(res_summary, aes(x = sigma_KL, y = sigma_VAE, label = r)) +
  geom_point(color="firebrick") +
  geom_text(size=2.5, vjust=-0.5) +
  theme_minimal() +
  labs(title="Elasticities across regions",
       x="σ_KL (Capital–Labour)", y="σ_VAE (Energy–VA)")

# 3. Distribution of input shares
res_long <- res_summary %>%
  select(r, d_KL, d_VAE) %>%
  tidyr::pivot_longer(cols=c(d_KL,d_VAE), names_to="parameter", values_to="value")

ggplot(res_long, aes(x=reorder(r, value), y=value, fill=parameter)) +
  geom_col(position="dodge") +
  coord_flip() +
  theme_minimal() +
  labs(title="Share parameters by region",
       x="Region", y="Parameter value")

# 4. Runtime vs Fit quality (RSS)
# --- Combine diagnostics for plotting ---
res_diagnostics <- bind_rows(map(results, "diagnostics"))
write_csv(res_diagnostics, "CES_diagnostics.csv")

# --- Plot runtime vs fit quality (RSS) ---
ggplot(res_diagnostics, aes(x = time, y = rss, color = method, label = r)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text(size = 2.5, vjust = -0.8, check_overlap = TRUE) +
  scale_y_log10() +  # RSS can span orders of magnitude
  theme_minimal() +
  labs(title = "Tradeoff: Runtime vs Fit Quality",
       subtitle = "Across all methods and regions",
       x = "Runtime (seconds)",
       y = "Residual Sum of Squares (log scale)",
       color = "Method")

# 5. Observed vs Fitted values from the CES function
obs_vs_fit <- fitted_all %>%
  left_join(dfS %>% select(r, t, Ys), by = c("r","t"))

ggplot(obs_vs_fit, aes(x = Ys, y = fitted, color = r)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(title = "Observed vs Fitted Output",
       x = "Observed (scaled Ys)", y = "Fitted (scaled Ys)")

