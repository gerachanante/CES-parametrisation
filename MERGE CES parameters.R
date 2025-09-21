options(scipen = 999) # avoids scientific notation unless necessary
setTimeLimit(cpu = Inf, elapsed = Inf, transient = TRUE)

# ---- PACKAGES ----
# install.packages(c("micEconCES","dplyr","readr","purrr","ggplot2",
#                    "parallel","ggpmisc","pheatmap","ggrepel","patchwork","tibble","viridis"))

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
library(lmtest)
library(sandwich)

# ---- SETTINGS ----
setwd("C:/Users/escami_g/OneDrive - Paul Scherrer Institut/05.Models/MERGE updates/CES-parametrisation/fine grid")
infile <- "MERGE macro.csv"

RUN_ESTIMATION <- FALSE

rhoGrid_KL  <- c(seq(-0.9, 0.5, by = 0.02), seq(0.51, 3, by = 0.1), seq(3.5, 19.5, by = 0.5), seq(20, 40, by = 2))
rhoGrid_VAE <- c(seq(-0.8, 0.5, by = 0.02), seq(0.51, 3, by = 0.1), seq(3.5, 19.5, by = 0.5), seq(20, 40, by = 2))

# Helper functions
# AIC/BIC from RSS on the log scale (since multErr=TRUE -> we fit log(Y))
aic_bic_from_rss <- function(resid_log, k, add_k_rho = 0L) {
  n   <- length(resid_log)
  RSS <- sum(resid_log^2)
  k0  <- k + add_k_rho
  AIC <- n * log(RSS / n) + 2 * k0
  BIC <- n * log(RSS / n) + k0 * log(n)
  AICc <- if (n - k0 - 1 > 0) AIC + (2*k0*(k0+1))/(n - k0 - 1) else NA_real_
  list(AIC = AIC, BIC = BIC, AICc = AICc, RSS = RSS, n = n, k = k0)
}

# how many rhos to count as “parameters” (conservative)
rho_penalty <- as.integer(length(rhoGrid_KL)  > 1) +
  as.integer(length(rhoGrid_VAE) > 1)

`%||%` <- function(x, y) if (!is.null(x)) x else y


on_grid_edge <- function(val, grid, tol = 1e-12) {
  if (length(grid) == 0L) {
    return(rep(FALSE, length(val)))
  }
  g <- sort(unique(grid))
  gmin <- g[1]
  gmax <- g[length(g)]
  # compute whether each val is within tol of either edge
  near(val, gmin, tol = tol) | near(val, gmax, tol = tol)
}


# ---- ESTIMATION ----
if (isTRUE(RUN_ESTIMATION)) {
  # Run in multiple sessions for faster solve times
  plan(multisession, workers = parallel::detectCores() - 1)

  # ---- DATA ----
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
  
estimate_region <- function(d, region_name) {
  message("\nEstimating region: ", region_name)
  d_num <- d %>% transmute(t, Ys, Ks, Ls, Es)
  # Fast methods: NM, Nelder-Mead, Newton, L-BFGS-B
  # Slow methods: LM, PORT, BFGS
  # Terribly slow: SANN, DE, CG
  methods <- c("LM", "NM", "Nelder-Mead", "PORT", "Newton", "L-BFGS-B")
  #  methods <- c("LM", "NM", "Nelder-Mead", "PORT", "Newton", "L-BFGS-B")
  #  methods <- c("LM", "NM", "Nelder-Mead", "BFGS", "PORT", "Newton", "CG", "L-BFGS-B", "SANN", "DE")
  
  fit_all   <- setNames(vector("list", length(methods)), methods)
  conv_all  <- setNames(rep(FALSE, length(methods)), methods)
  msg_all   <- setNames(rep(NA_character_, length(methods)), methods)
  times_all <- setNames(rep(NA_real_, length(methods)), methods)
  
  for (m in methods) {
    t0 <- Sys.time()
    
    start_arg <- NULL; lower_arg <- NULL; upper_arg <- NULL; control_arg <- NULL
    
    if (m == "Newton") {
      start_arg <- c(
        gamma     = runif(1, 0.9, 1.1),
        lambda    = runif(1, -0.005, 0.005),
        delta_KL  = runif(1, 0.4, 0.6),
        delta_VAE = runif(1, 0.4, 0.6),
        nu        = runif(1, 0.9, 1.1)
      )
      # no lower/upper/control
    }
    if (m == "L-BFGS-B") {
      start_arg <- c(gamma=1, lambda=0.001, delta_KL=0.5, delta_VAE=0.5, nu=1)
      lower_arg <- c(gamma=1e-6, lambda=-0.1, delta_KL=1e-6, delta_VAE=1e-6, nu=1e-6)
      upper_arg <- c(gamma=10, lambda=0.1, delta_KL=1-1e-6, delta_VAE=1-1e-6, nu=10)
      control_arg <- list(maxit = 1000, factr = 1e9)
    }
    if (m == "PORT") control_arg <- list(eval.max=1e5, iter.max=2e4, reltol=1e-8)
    if (m == "BFGS") {
      start_arg <- c(gamma=1, lambda=0.001, delta_KL=0.5, delta_VAE=0.5, nu=1)
      lower_arg <- c(delta_KL=1e-6, delta_VAE=1e-6)
      upper_arg <- c(delta_KL=1-1e-6, delta_VAE=1-1e-6)
      control_arg <- list(maxit = 1000, reltol = 1e-8)
    }
    if (m == "CG") {
      start_arg <- c(gamma=1, lambda=0.001, delta_KL=0.5, delta_VAE=0.5, nu=1)
      lower_arg <- c(delta_KL=1e-6, delta_VAE=1e-6)
      upper_arg <- c(delta_KL=1-1e-6, delta_VAE=1-1e-6)
      control_arg <- list(maxit = 500, reltol = 1e-8)
    }    
    if (m %in% c("NM","Nelder-Mead")) {
      lower_arg <- c(delta_KL=1e-6, delta_VAE=1e-6)
      upper_arg <- c(delta_KL=1-1e-6, delta_VAE=1-1e-6)
      control_arg <- list(maxit=5000, reltol=1e-8)
    }
    if (m == "LM") {
      lower_arg <- c(delta_KL=1e-6, delta_VAE=1e-6)
      upper_arg <- c(delta_KL=1-1e-6, delta_VAE=1-1e-6)
      control_arg <- list(maxiter=10000, ftol=1e-8, maxfev=5000)
    }
    if (m == "SANN") {
      lower_arg <- c(delta_KL=1e-6, delta_VAE=1e-6)
      upper_arg <- c(delta_KL=1-1e-6, delta_VAE=1-1e-6)
      control_arg <- list(maxit=5000, temp=10, tmax=50)
    }
    if (m == "DE") {
      lower_arg <- c(gamma=1, lambda=0.001, delta_KL=1e-6, delta_VAE=1e-6, nu=1)
      upper_arg <- c(gamma=10, lambda=0.2, delta_KL=1-1e-6, delta_VAE=1-1e-6, nu=10)
      control_arg <- list(itermax=100)
    }
    
    
    args_list <- list(
      yName = "Ys",
      xNames = c("Ks","Ls","Es"),
      tName = "t",
      data = d_num,
      vrs = TRUE,
      multErr = TRUE,
      method = m,
      rho1 = rhoGrid_KL,
      rho = rhoGrid_VAE,
      returnGridAll = TRUE
    )
    if (!is.null(start_arg))   args_list$start   <- start_arg
    if (!is.null(lower_arg))   args_list$lower   <- lower_arg
    if (!is.null(upper_arg))   args_list$upper   <- upper_arg
    if (!is.null(control_arg)) args_list$control <- control_arg
    
    fit_try <- try(suppressWarnings(do.call(cesEst, args_list)), silent=TRUE)
    
    runtime <- as.numeric(difftime(Sys.time(), t0, units="secs"))
    times_all[m] <- runtime
    
    success   <- inherits(fit_try,"cesEst")
    conv_flag <- if (success && !is.null(fit_try$convergence)) fit_try$convergence else FALSE
    msg_flag  <- if (inherits(fit_try,"try-error")) as.character(fit_try)[1]
    else if (success && !is.null(fit_try$message)) fit_try$message else ""
    
    fit_all[[m]] <- if (success) fit_try else NULL
    conv_all[m]  <- conv_flag
    msg_all[m]   <- msg_flag
    
    message(sprintf("  %s %-10s in %.1fs%s",
                    if (conv_flag) "✓" else "✗", m, runtime,
                    if (!conv_flag && nzchar(msg_flag)) paste0(" msg: ", msg_flag) else ""))
  }
  
  list(fits=fit_all, conv=conv_all, msg=msg_all, times=times_all, data=d_num)
}

# ---- EXTRACT RESULTS ----
extract_region <- function(region_name, region_fits) {
  if (is.null(region_fits)) return(NULL)
  
  diag_tbl     <- tibble()
  timevary_tbl <- tibble()
  grid_tbl     <- tibble()
  coef_long    <- tibble()
  
  for (m in names(region_fits$fits)) {
    fit_obj   <- region_fits$fits[[m]]
    conv_flag <- region_fits$conv[[m]] %||% FALSE
    msg_flag  <- region_fits$msg[[m]] %||% NA_character_
    runtime   <- region_fits$times[[m]] %||% NA_real_
    
    gamma_est <- lambda_est <- delta_KL_est <- delta_VAE_est <- nu_est <- NA_real_
    rho1 <- rhoT <- sigma_KL <- sigma_VAE <- NA_real_
    rss_val <- R2_val <- adjR2_val <- NA_real_
    iter_val <- NA_integer_
    aic_naive <- bic_naive <- aic_plusRho <- bic_plusRho <- NA_real_
    p_gamma <- p_lambda <- p_delta_KL <- p_delta_VAE <- p_nu <- NA_real_
    
    if (inherits(fit_obj,"cesEst")) {
      s        <- try(summary(fit_obj), silent = TRUE)
      coef_mat <- if (!inherits(s, "try-error")) {
        # coef(summary()) already is a numeric matrix with rows like "gamma", "lambda", "delta_1", "delta", "nu"
        suppressWarnings( tryCatch(as.matrix(coef(s)), error = function(e) NULL) )
      } else NULL
      
      get_coef <- function(mat, par, col = "Estimate") {
        if (!is.null(mat) && par %in% rownames(mat) && col %in% colnames(mat)) mat[par, col] else NA_real_
      }
      
      gamma_est     <- get_coef(coef_mat, "gamma",   "Estimate")
      lambda_est    <- get_coef(coef_mat, "lambda",  "Estimate")
      delta_KL_est  <- get_coef(coef_mat, "delta_1", "Estimate")
      delta_VAE_est <- get_coef(coef_mat, "delta",   "Estimate")
      nu_est        <- get_coef(coef_mat, "nu",      "Estimate")
      
      # p-values directly from the same table
      p_gamma     <- get_coef(coef_mat, "gamma",   "Pr(>|t|)")
      p_lambda    <- get_coef(coef_mat, "lambda",  "Pr(>|t|)")
      p_delta_KL  <- get_coef(coef_mat, "delta_1", "Pr(>|t|)")
      p_delta_VAE <- get_coef(coef_mat, "delta",   "Pr(>|t|)")
      p_nu        <- get_coef(coef_mat, "nu",      "Pr(>|t|)")
      
      # AIC/BIC
      aic_val <- suppressWarnings(tryCatch(AIC(fit_obj), error = function(e) NA_real_))
      bic_val <- suppressWarnings(tryCatch(BIC(fit_obj), error = function(e) NA_real_))
      
      # Smaller version of coefficients for export
      if (!is.null(coef_mat)) {
        coef_long <- bind_rows(
          coef_long,
          as_tibble(coef_mat, rownames = "parameter") |>
            mutate(r = region_name, method = m, .before = 1)
        )
      }
      
      # Best rho from grid (if present), and elasticities
      if (!is.null(fit_obj$allRhoSum) && nrow(fit_obj$allRhoSum) > 0) {
        best <- fit_obj$allRhoSum |>
          tidyr::drop_na(rss) |>
          dplyr::slice_min(rss, n = 1)
        rho1 <- best$rho1 %||% NA_real_
        rhoT <- best$rho  %||% NA_real_
        grid_tbl <- bind_rows(grid_tbl, dplyr::mutate(fit_obj$allRhoSum, r = region_name, method = m))
      } else {
        cf   <- tryCatch(coef(fit_obj), error = function(e) NULL)
        rho1 <- if (!is.null(cf) && "rho1" %in% names(cf)) cf[["rho1"]] else NA_real_
        rhoT <- if (!is.null(cf) && "rho"  %in% names(cf)) cf[["rho"]]  else NA_real_
      }
      sigma_KL  <- ifelse(is.finite(1/(1+rho1)), 1/(1+rho1), NA_real_)
      sigma_VAE <- ifelse(is.finite(1/(1+rhoT)), 1/(1+rhoT), NA_real_)
      
      # Goodness of fit & information criteria on the log scale
      iter_val <- fit_obj$iter %||% NA_integer_
      obs_log  <- log(region_fits$data$Ys + 1e-12)
      fit_log  <- log(as.numeric(fit_obj$fitted.values) + 1e-12)
      resid_log <- obs_log - fit_log
      R2_val <- 1 - sum(resid_log^2) / sum((obs_log - mean(obs_log))^2)
      
      k_hat <- length(tryCatch(coef(fit_obj), error=function(e) numeric(0)))
      n_obs <- length(obs_log)
      adjR2_val <- 1 - (1 - R2_val) * (n_obs - 1) / (n_obs - k_hat - 1)
      
      ic0 <- aic_bic_from_rss(resid_log, k = k_hat, add_k_rho = 0L)
      icR <- aic_bic_from_rss(resid_log, k = k_hat, add_k_rho = rho_penalty)
      
      rss_val      <- ic0$RSS
      aic_naive    <- ic0$AIC
      bic_naive    <- ic0$BIC
      aicc_naive   <- ic0$AICc
      aic_plusRho  <- icR$AIC
      bic_plusRho  <- icR$BIC
      aicc_plusRho <- icR$AICc

      # Time-varying TFP series for IAM table
      timevary_tbl <- bind_rows(timevary_tbl, tibble(
        r        = region_name,
        method   = m,
        t        = region_fits$data$t,
        fitted   = exp(fit_log),                  # back on original scale
        residual = obs_log - fit_log,             # log residual
        TFP      = gamma_est * exp(lambda_est * region_fits$data$t)
      ))
    }
    
    # Building the regional table of parameters and statistics
    diag_tbl <- bind_rows(diag_tbl, tibble(
      r           = region_name,
      method      = m,
      rss         = rss_val,
      R2          = R2_val,
      adjR2       = adjR2_val,
      iter        = iter_val,
      conv        = conv_flag,
      message     = msg_flag,
      gamma       = gamma_est,
      p_gamma     = p_gamma,
      lambda      = lambda_est,
      p_lambda    = p_lambda,
      delta_KL    = delta_KL_est,
      p_delta_KL  = p_delta_KL,
      delta_VAE   = delta_VAE_est,
      p_delta_VAE = p_delta_VAE,
      nu          = nu_est,
      p_nu        = p_nu,
      rho_KL      = rho1,
      rho_VAE     = rhoT,
      sigma_KL    = sigma_KL,
      sigma_VAE   = sigma_VAE,
      AIC_naive   = aic_naive,
      BIC_naive   = bic_naive,
      AICc_naive  = aicc_naive,
      AIC_plusRho = aic_plusRho,
      BIC_plusRho = bic_plusRho,
      AICc_plusRho= aicc_plusRho,
      runtime     = runtime
    ))
  }
  
  list(diag=diag_tbl, timevary=timevary_tbl, grid=grid_tbl, coef_long=coef_long)
}


# Region splits
splits <- dfS %>% group_split(r, .keep=TRUE)
region_names <- dfS %>% distinct(r) %>% pull(r)

# ---- RUN ALL ----
fits_all <- future_map2(
  splits, region_names,
  ~ tryCatch(
    estimate_region(.x, .y),
    error = function(e) {
      message("Completely failed region: ", .y, " → ", e$message)
      list(fits=list(), conv=list(), msg=list(error=e$message), times=list(), data=.x)
    }
  ),
  .progress=TRUE,
  .options=furrr_options(
    packages=c("micEconCES","dplyr","tidyr","purrr","readr"),
    globals=c("rhoGrid_KL","rhoGrid_VAE","estimate_region","extract_region"),
    seed=TRUE
  ),
  .env_globals=globalenv()
)




# --- Extract diagnostics ---
results            <- map2(region_names, fits_all, extract_region)
results_table      <- bind_rows(map(results, "diag"))
results_table_time <- bind_rows(map(results, "timevary"))
results_grid       <- bind_rows(compact(map(results, "grid")))
results_coef_long  <- bind_rows(compact(map(results, "coef_long")))


# ---- EXPORT ----
write_csv(results_table, "CES_region_method.csv")
write_csv(results_table_time, "CES_region_method_year.csv")
write_csv(results_grid, "CES_gridsearch.csv")
write_csv(results_coef_long,  "CES_coefficients_long.csv") 

# Save all results
saveRDS(results, file = "results_run1.rds")
}
results_table      <- read_csv("CES_region_method.csv",      show_col_types = FALSE)
results_table_time <- read_csv("CES_region_method_year.csv", show_col_types = FALSE)
results_grid       <- if (file.exists("CES_gridsearch.csv")) read_csv("CES_gridsearch.csv", show_col_types = FALSE) else tibble()
results_coef_long  <- if (file.exists("CES_coefficients_long.csv")) read_csv("CES_coefficients_long.csv", show_col_types = FALSE) else tibble()

# ---- Convergence summary ----
convergence_summary <- results_table %>%
  select(r, method, conv) %>%
  distinct() %>%
  count(method, conv, name = "count") %>%
  mutate(status = ifelse(conv, "Converged", "Failed"))

# ---- Validity flags (computed here only) ----
results_table_valid <- results_table %>%
  mutate(
    on_edge_KL  = on_grid_edge(rho_KL,  rhoGrid_KL),
    on_edge_VAE = on_grid_edge(rho_VAE, rhoGrid_VAE),
    valid = (conv == TRUE) &
      is.finite(delta_KL) & delta_KL > 0 & delta_KL < 1 &
      is.finite(delta_VAE) & delta_VAE > 0 & delta_VAE < 1 &
      is.finite(gamma) & gamma > 0 & gamma < 1e6 &
      is.finite(sigma_KL)  & sigma_KL  > 0 & sigma_KL  < 100 &
      is.finite(sigma_VAE) & sigma_VAE > 0 & sigma_VAE < 100
  ) %>%
  mutate(
    invalid_reason = if_else(
      valid,
      NA_character_,
      paste(
        if_else(conv != TRUE, "no_convergence", NA_character_),
        if_else(!is.finite(delta_KL) | delta_KL <= 0 | delta_KL >= 1, "bad_delta_KL", NA_character_),
        if_else(!is.finite(delta_VAE) | delta_VAE <= 0 | delta_VAE >= 1, "bad_delta_VAE", NA_character_),
        if_else(!is.finite(gamma)     | gamma     <= 0 | gamma     >= 1e6, "bad_gamma",   NA_character_),
        if_else(!is.finite(sigma_KL)  | sigma_KL  <= 0 | sigma_KL  >= 100,  "bad_sigma_KL",  NA_character_),
        if_else(!is.finite(sigma_VAE) | sigma_VAE <= 0 | sigma_VAE >= 100,  "bad_sigma_VAE", NA_character_),
        if_else(on_edge_KL,  "rho_KL_on_edge",  NA_character_),
        if_else(on_edge_VAE, "rho_VAE_on_edge", NA_character_),
        sep = "|"
      )
    )
  )

# ---- Best method per region (AICc+rho → RSS → runtime) ----
order_cols <- c("AICc_plusRho", "rss", "runtime")

best_valid <- results_table_valid %>%
  filter(valid) %>%
  group_by(r) %>%
  arrange(across(all_of(order_cols))) %>%
  slice_head(n = 1) %>%
  ungroup()

still_missing <- setdiff(unique(results_table_valid$r), best_valid$r)

fallback <- results_table_valid %>%
  filter(r %in% still_missing) %>%
  filter(conv == TRUE) %>%
  filter(
    is.finite(delta_KL), delta_KL > 0, delta_KL < 1,
    is.finite(delta_VAE), delta_VAE > 0, delta_VAE < 1,
    is.finite(gamma), gamma > 0, gamma < 1e6,   # looser gamma upper bound
    is.finite(sigma_KL),  sigma_KL > 0, sigma_KL < 50, # looser elasticities
    is.finite(sigma_VAE), sigma_VAE > 0, sigma_VAE < 50
  ) %>%
  arrange(rss, runtime) %>%
  group_by(r) %>%
  slice_head(n = 1) %>%
  mutate(valid = FALSE, fallback_tier = "relaxed") %>%
  ungroup()

best_methods <- bind_rows(best_valid, fallback)

# ---- Invalid runs table ----
invalid_runs <- results_table_valid %>%
  filter(!valid) %>%
  select(r, method, invalid_reason, everything())

# ---- AICc weights (among valid only) ----
aic_weights <- results_table_valid %>%
  filter(valid) %>%
  group_by(r) %>%
  mutate(
    dAICc = AICc_plusRho - min(AICc_plusRho, na.rm = TRUE),
    wAICc = exp(-0.5 * dAICc) / sum(exp(-0.5 * dAICc), na.rm = TRUE)
  ) %>%
  ungroup()

# ---- Share of grid converged ----
grid_conv_share <- results_grid %>%
  mutate(convergence = as.logical(convergence %||% !is.na(rss))) %>%
  group_by(r, method) %>%
  summarise(
    grid_points = n(),
    n_converged = sum(convergence, na.rm = TRUE),
    share_converged = n_converged / grid_points,
    .groups = "drop"
  )

# ---- IAM parameter table (join best methods) ----
iam_table <- results_table_time %>%
  inner_join(
    best_methods %>% 
      select(r, method, gamma, nu, delta_KL, delta_VAE, rho_KL, rho_VAE, sigma_KL, sigma_VAE, conv),
    by = c("r","method")
  ) %>%
  distinct(r, t, .keep_all = TRUE) %>%
  transmute(
    year   = t,
    region = r,
    total_factor_productivity   = TFP,
    share_capital_valueadded    = delta_KL,
    share_valueadded_output     = delta_VAE,
    elasticity_substitution_capital_labour = sigma_KL,
    elasticity_substitution_valueadded_energy = sigma_VAE,
    valueadded_CES_exponent     = rho_KL,
    output_CES_exponent         = rho_VAE,
    gamma_parameter             = gamma,
    nu_parameter                = nu,
    converged                   = conv
  )

# ---- EXPORT derived artifacts ----
write_csv(results_table_valid,  "CES_region_method_valid.csv")
write_csv(convergence_summary,  "CES_convergence_summary.csv")
write_csv(best_methods,         "CES_best_methods.csv")
write_csv(invalid_runs,         "CES_invalid_runs.csv")
write_csv(aic_weights,          "CES_AICc_weights.csv")
write_csv(grid_conv_share,      "CES_grid_convergence_share.csv")
write_csv(iam_table,            "IAM_params.csv")