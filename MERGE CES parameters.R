options(scipen = 999) # avoids scientific notation unless necessary
setTimeLimit(cpu = Inf, elapsed = Inf, transient = TRUE)

# ---- PACKAGES ----
#install.packages(c("micEconCES","dplyr","readr","purrr","parallel","tibble","progressr"))

library(micEconCES)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(furrr)
library(parallel)
library(tibble)

# ---- SETTINGS ----
options(future.rng.onMisuse = "ignore")
options(future.wait.timeout = 0)   # disables waiting timeout


setwd("C:/Users/escami_g/OneDrive - Paul Scherrer Institut/05.Models/MERGE updates/CES-parametrisation/stage1")
infile <- "MERGE macro.csv"

# TRUE to run a fresh estimation or FALSE to reuse saved results (.rds file)
RUN_ESTIMATION <- FALSE

# Stage 1 grid (broad scope)
rhoGrid_KL  <- c(seq(-0.9, -0.15, by = 0.15),
                 seq(-0.1, 0.5, by = 0.1),
                 seq(0.7, 2, by = 0.3),
                 3, 5, 10, 20, 50)
rhoGrid_VAE <- c(seq(-0.8, -0.2, by = 0.2),
                 seq(-0.1, 0.5, by = 0.1),
                 seq(0.7, 3, by = 0.3),
                 4, 6, 10, 20)

# Stage 2 grid (refined)
# rhoGrid_KL  <- c(seq(-0.8, -0.5, by = 0.01),
#                  seq(-0.495, 0.5, by = 0.005),
#                  seq(0.51, 2, by = 0.01),
#                  seq(2.05, 5, by = 0.05),
#                  seq(5.1, 7, by = 0.1),
#                  8, 9, 10, 15, 20, 50, 100)
# rhoGrid_VAE <- c(seq(-0.8, 0, by = 0.01),
#                  seq(0.05, 6, by = 0.05),
#                  seq(6.2, 8, by = 0.2),
#                  9, 10, 15, 20, 50, 100)

# Test grid
# rhoGrid_KL <- c(seq(-0.5, 1, by = 0.25))
# rhoGrid_VAE <- c(seq(-0.5, 1, by = 0.25))

# Helper functions
# Log scale modification of information criteria: AIC, BIC, AICc from residual sums of squares
aic_bic_from_rss <- function(resid_log, k, add_k_rho = 0L) {
  n   <- length(resid_log)
  RSS <- sum(resid_log^2)
  k0  <- k + add_k_rho
  AIC <- n*log(RSS/n) + 2*k0
  BIC <- n*log(RSS/n) + k0*log(n)
  AICc <- if (n - k0 - 1 > 0) AIC + (2*k0*(k0 + 1))/(n - k0 - 1) else NA_real_
  list(AIC = AIC, BIC = BIC, AICc = AICc, RSS = RSS, n = n, k = k0)
}

# Counts number of rho treated as parameters for conservative AIC penalties
rho_penalty <- as.integer(length(rhoGrid_KL) > 1) +
  as.integer(length(rhoGrid_VAE) > 1)

# Custom operator function to return a fallback value if the left side is NULL or NA
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || isTRUE(all(is.na(x)))) {
    y
  } else {
    # collapse to scalar if needed
    if (length(x) > 1) x[[1]] else x
  }
}

# Extractor of boolean TRUE/FALSE defaults
safe_bool <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) {
    default
  } else {
    out <- suppressWarnings(as.logical(x[[1]]))
    if (is.na(out)) default else out
  }
}

# Extractor of text defaults
safe_chr <- function(x, default = NA_character_) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) {
    default
  } else {
    as.character(x[[1]])
  }
}

# Extractor of number defaults
safe_num <- function(x, default = NA_real_) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) {
    default
  } else {
    suppressWarnings(as.numeric(x[[1]]))
  }
}

# Detects if the rho values are at the edge of their respective grids
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

# Standardises the coefficients across different summary formats
coef_table_safe <- function(fit_obj) {
  # Try several common locations/dispatch paths
  try_list <- list(
    function(x) coef(summary(x)),
    function(x) summary(x)$coefficients,
    function(x) summary(x)$coef,
    function(x) summary(x)$coefTable
  )
  for (f in try_list) {
    cm <- try(f(fit_obj), silent = TRUE)
    if (!inherits(cm, "try-error") && !is.null(cm)) {
      cm <- try(as.matrix(cm), silent = TRUE)
      if (!inherits(cm, "try-error") && is.matrix(cm) && nrow(cm) > 0)
        return(cm)
    }
  }
  return(NULL)
}

# Extract a single value from a coef table if present
coef_get <- function(cm, par, col = "Estimate") {
  if (!is.null(cm) && par %in% rownames(cm) && col %in% colnames(cm)) {
    val <- suppressWarnings(as.numeric(cm[par, col]))
    if (is.finite(val)) return(val)
  }
  NA_real_
}

# Tries different approaches to extract the number of iterations from the solvers
iter_safe <- function(fit_obj) {
  candidates <- list(
    tryCatch(as.numeric(fit_obj$iter), error = function(e) NA_real_),
    tryCatch(as.numeric(fit_obj$iterations), error = function(e) NA_real_),
    tryCatch(as.numeric(fit_obj$niter), error = function(e) NA_real_),
    tryCatch(as.numeric(fit_obj$counts[["function"]]), error = function(e) NA_real_),
    tryCatch(as.numeric(fit_obj$optim$counts[["function"]]), error = function(e) NA_real_)
  )
  it <- NA_real_
  for (v in candidates) {
    if (length(v) == 1 && is.finite(v)) { it <- v; break }
  }
  it
}

# Helper: grid-point stats extractor
extract_grid_coeffs <- function(fit_sub) {
  # fit_sub is a single cesEst object from allRhoFull
  
  if (!inherits(fit_sub, "cesEst")) return(NULL)
  
  cm <- coef_table_safe(fit_sub)
  if (is.null(cm) || !is.matrix(cm) || nrow(cm) == 0L) return(NULL)
  
  rn <- rownames(cm)
  cn <- colnames(cm)
  
  # flexible column pickers
  pick_col <- function(patterns) {
    ix <- which(vapply(
      cn,
      function(z) any(grepl(patterns, z, ignore.case = TRUE)),
      logical(1)
    ))
    if (length(ix) == 0) NA_integer_ else ix[1]
  }
  col_est <- pick_col("^(estimate|coef|value)$|^estimate$|^coef$|^coeff")
  col_se  <- pick_col("(std\\.? ?error|se)")
  col_t   <- pick_col("^(t.?value|z|t)$")
  col_p   <- pick_col("^(pr\\(|p.?value|p$)")
  
  # flexible row pickers for parameter names
  pick_row <- function(patterns, exclude = NULL) {
    ok <- which(vapply(
      rn,
      function(z) any(grepl(patterns, z, ignore.case = TRUE)),
      logical(1)
    ))
    if (!is.null(exclude) && length(ok) > 0L) {
      ok <- ok[!vapply(
        rn[ok],
        function(z) any(grepl(exclude, z, ignore.case = TRUE)),
        logical(1)
      )]
    }
    if (length(ok) == 0L) NA_integer_ else ok[1]
  }
  r_gamma     <- pick_row("^gamma$")
  r_lambda    <- pick_row("^lambda|^lam$")
  r_delta1    <- pick_row("^delta[_ ]?1$|^delta-?1$|delta[_]?kl")
  r_deltaMain <- pick_row("^delta$", exclude = "1")
  r_nu        <- pick_row("^nu$")
  
  get_val <- function(ri, ci) {
    if (is.na(ri) || is.na(ci)) return(NA_real_)
    v <- suppressWarnings(as.numeric(cm[ri, ci, drop = TRUE]))
    if (is.finite(v)) v else NA_real_
  }
  
  est <- list(
    gamma     = get_val(r_gamma,     col_est),
    lambda    = get_val(r_lambda,    col_est),
    delta_KVA = get_val(r_delta1,    col_est),
    delta_VAY = get_val(r_deltaMain, col_est),
    nu        = get_val(r_nu,        col_est)
  )
  
  se  <- list(
    gamma     = get_val(r_gamma,     col_se),
    lambda    = get_val(r_lambda,    col_se),
    delta_KVA = get_val(r_delta1,    col_se),
    delta_VAY = get_val(r_deltaMain, col_se),
    nu        = get_val(r_nu,        col_se)
  )
  
  tval <- list(
    gamma     = get_val(r_gamma,     col_t),
    lambda    = get_val(r_lambda,    col_t),
    delta_KVA = get_val(r_delta1,    col_t),
    delta_VAY = get_val(r_deltaMain, col_t),
    nu        = get_val(r_nu,        col_t)
  )
  pval <- list(
    gamma     = get_val(r_gamma,     col_p),
    lambda    = get_val(r_lambda,    col_p),
    delta_KVA = get_val(r_delta1,    col_p),
    delta_VAY = get_val(r_deltaMain, col_p),
    nu        = get_val(r_nu,        col_p)
  )
  
  # fallback SEs from vcov if SEs are missing
  if (any(!is.finite(unlist(se)))) {
    vc <- try(vcov(fit_sub), silent = TRUE)
    if (!inherits(vc, "try-error") && is.matrix(vc)) {
      se_v <- try(sqrt(diag(vc)), silent = TRUE)
      if (!inherits(se_v, "try-error")) {
        nms <- names(se_v); if (is.null(nms)) nms <- rownames(vc)
        if (length(nms)) {
          get_se <- function(pattern) {
            ix <- which(grepl(pattern, nms, ignore.case = TRUE))
            if (length(ix) == 0) NA_real_ else {
              v <- unname(se_v[ix[1]])
              if (is.finite(v) && v > 0) v else NA_real_
            }
          }
          if (!is.finite(se$gamma))     se$gamma     <- get_se("^gamma$")
          if (!is.finite(se$lambda))    se$lambda    <- get_se("^lambda|^lam$")
          if (!is.finite(se$delta_KVA)) se$delta_KVA <- get_se("^delta[_ ]?1$|^delta-?1$|delta[_]?kl")
          if (!is.finite(se$delta_VAY)) se$delta_VAY <- get_se("^delta$")
          if (!is.finite(se$nu))        se$nu        <- get_se("^nu$")
        }
      }
    }
  }
  
  # fill in t if missing
  for (nm in names(est)) {
    if (!is.finite(tval[[nm]]) &&
        is.finite(est[[nm]]) &&
        is.finite(se[[nm]]) && se[[nm]] > 0) {
      tval[[nm]] <- est[[nm]] / se[[nm]]
    }
  }
  # fill in p if missing (normal approx)
  for (nm in names(tval)) {
    if (!is.finite(pval[[nm]]) && is.finite(tval[[nm]])) {
      pval[[nm]] <- 2 * (1 - pnorm(abs(tval[[nm]])))
    }
  }
  
  # rhos for joining
  cf <- try(stats::coef(fit_sub), silent = TRUE)
  rho_KL  <- NA_real_
  rho_VAE <- NA_real_
  if (!inherits(cf, "try-error") && length(cf) > 0L) {
    nms <- names(cf)
    if (any(grepl("^rho[_]?1$", nms))) {
      rho_KL <- as.numeric(cf[grep("^rho[_]?1$", nms)[1]])
    }
    if ("rho" %in% nms) {
      rho_VAE <- as.numeric(cf[["rho"]])
    } else if (any(grepl("^rho2$", nms))) {
      rho_VAE <- as.numeric(cf[grep("^rho2$", nms)[1]])
    }
  }
  
  tibble(
    rho_KL  = rho_KL,
    rho_VAE = rho_VAE,
    gamma     = est$gamma,
    lambda    = est$lambda,
    delta_KVA = est$delta_KVA,
    delta_VAY = est$delta_VAY,
    nu        = est$nu,
    se_gamma     = se$gamma,
    se_lambda    = se$lambda,
    se_delta_KVA = se$delta_KVA,
    se_delta_VAY = se$delta_VAY,
    se_nu        = se$nu,
    t_gamma      = tval$gamma,
    t_lambda     = tval$lambda,
    t_delta_KVA  = tval$delta_KVA,
    t_delta_VAY  = tval$delta_VAY,
    t_nu         = tval$nu,
    p_gamma      = pval$gamma,
    p_lambda     = pval$lambda,
    p_delta_KVA  = pval$delta_KVA,
    p_delta_VAY  = pval$delta_VAY,
    p_nu         = pval$nu,
    ci_lo_gamma     = ifelse(is.finite(est$gamma)     & is.finite(se$gamma),     est$gamma     - 1.96 * se$gamma,     NA_real_),
    ci_hi_gamma     = ifelse(is.finite(est$gamma)     & is.finite(se$gamma),     est$gamma     + 1.96 * se$gamma,     NA_real_),
    ci_lo_lambda    = ifelse(is.finite(est$lambda)    & is.finite(se$lambda),    est$lambda    - 1.96 * se$lambda,    NA_real_),
    ci_hi_lambda    = ifelse(is.finite(est$lambda)    & is.finite(se$lambda),    est$lambda    + 1.96 * se$lambda,    NA_real_),
    ci_lo_delta_KVA = ifelse(is.finite(est$delta_KVA) & is.finite(se$delta_KVA), est$delta_KVA - 1.96 * se$delta_KVA, NA_real_),
    ci_hi_delta_KVA = ifelse(is.finite(est$delta_KVA) & is.finite(se$delta_KVA), est$delta_KVA + 1.96 * se$delta_KVA, NA_real_),
    ci_lo_delta_VAY = ifelse(is.finite(est$delta_VAY) & is.finite(se$delta_VAY), est$delta_VAY - 1.96 * se$delta_VAY, NA_real_),
    ci_hi_delta_VAY = ifelse(is.finite(est$delta_VAY) & is.finite(se$delta_VAY), est$delta_VAY + 1.96 * se$delta_VAY, NA_real_),
    ci_lo_nu        = ifelse(is.finite(est$nu)        & is.finite(se$nu),        est$nu        - 1.96 * se$nu,        NA_real_),
    ci_hi_nu        = ifelse(is.finite(est$nu)        & is.finite(se$nu),        est$nu        + 1.96 * se$nu,        NA_real_)
  )
}




# ---- ESTIMATION ----
if (isTRUE(RUN_ESTIMATION)) {
  # Run in multiple sessions for faster solve times. Includes a safe fallback to sequential solve if it fails
  suppressWarnings({
    tryCatch({
      future::plan(future::multisession, workers = max(1, parallel::detectCores() - 1))
    }, error = function(e) {
      warning("multisession failed (", conditionMessage(e), "); falling back to sequential.")
      future::plan(future::sequential)
    })
  })
  
  # Data loading and normalisation
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

  # Region estimation loop
  estimate_region <- function(d, region_name) {
    message("\nEstimating region: ", region_name)
    d_num <- d %>% transmute(t, Ys, Ks, Ls, Es)
    # Setting a deterministic seed for replicability
    seed_val <- sum(utf8ToInt(region_name))
    seed_val <- abs(seed_val) %% .Machine$integer.max
    set.seed(seed_val)
    # Fast methods: NM, Nelder-Mead, Newton, L-BFGS-B
    # Slow methods: LM, PORT, BFGS
    # Terribly slow: SANN, DE, CG
    # opt_methods <- unique(c("Newton", "L-BFGS-B", "BFGS"))
    opt_methods <- unique(c("LM", "NM", "Nelder-Mead", "BFGS", "PORT", "Newton", "CG", "L-BFGS-B", "SANN", "DE"))
    
    fit_all <- setNames(vector("list", length(opt_methods)), opt_methods)
    conv_all <- setNames(rep(FALSE, length(opt_methods)), opt_methods)
    msg_all <- setNames(rep(NA_character_, length(opt_methods)), opt_methods)
    times_all <- setNames(rep(NA_real_, length(opt_methods)), opt_methods)
    
    for (m in opt_methods) {
      t0 <- Sys.time()
      
      start_arg <- NULL; lower_arg <- NULL; upper_arg <- NULL; control_arg <- NULL
      
      if (m == "Newton") {
        start_arg <- c(
          gamma = runif(1, 0.9, 1.1),
          lambda = runif(1, -0.005, 0.005),
          delta_KVA = runif(1, 0.4, 0.6),
          delta_VAY = runif(1, 0.4, 0.6),
          nu = runif(1, 0.9, 1.1)
        )
        # no lower/upper/control
      }
      if (m == "L-BFGS-B") {
        start_arg <- c(gamma = 1, lambda = 0.001, delta_KVA = 0.5, delta_VAY = 0.5, nu = 1)
        lower_arg <- c(gamma = 0.1, lambda = -0.3, delta_KVA = 0.1, delta_VAY = 0.1, nu = 0.3)
        upper_arg <- c(gamma = 10, lambda = 0.3, delta_KVA = 0.9, delta_VAY = 0.9, nu = 5)
        control_arg <- list(maxit = 10000, factr = 1e9)
      }
      if (m == "PORT") control_arg <- list(eval.max=1e5, iter.max=1e5, reltol=1e-8)
      if (m == "BFGS") {
        start_arg <- c(gamma = 1, lambda = 0.001, delta_KVA = 0.5, delta_VAY = 0.5, nu = 1)
        lower_arg <- c(delta_KVA = 0.1, delta_VAY = 0.1)
        upper_arg <- c(delta_KVA = 0.9, delta_VAY = 0.9)
        control_arg <- list(maxit = 10000, reltol = 1e-8)
      }
      if (m == "CG") {
        start_arg <- c(gamma = 1, lambda = 0.001, delta_KVA = 0.5, delta_VAY = 0.5, nu = 1)
        lower_arg <- c(delta_KVA = 0.1, delta_VAY = 0.1)
        upper_arg <- c(delta_KVA = 0.9, delta_VAY = 0.9)
        control_arg <- list(maxit = 1000, reltol = 1e-8)
      }    
      if (m %in% c("NM","Nelder-Mead")) {
        lower_arg <- c(delta_KVA = 0.1, delta_VAY = 0.1)
        upper_arg <- c(delta_KVA = 0.9, delta_VAY = 0.9)
        control_arg <- list(maxit = 10000, reltol = 1e-8)
      }
      if (m == "LM") {
        lower_arg <- c(delta_KVA = 0.1, delta_VAY = 0.1)
        upper_arg <- c(delta_KVA = 0.9, delta_VAY = 0.9)
        control_arg <- list(maxiter = 10000, ftol = 1e-8, maxfev = 5000)
      }
      if (m == "SANN") {
        lower_arg <- c(delta_KVA = 0.1, delta_VAY = 0.1)
        upper_arg <- c(delta_KVA = 0.9, delta_VAY = 0.9)
        control_arg <- list(maxit = 10000, temp = 10, tmax = 50)
      }
      if (m == "DE") {
        lower_arg <- c(gamma = 0.1, lambda = -0.3, delta_KVA = 0.1, delta_VAY = 0.1, nu = 0.3)
        upper_arg <- c(gamma = 10, lambda = 0.3, delta_KVA = 0.9, delta_VAY = 0.9, nu = 5)
        control_arg <- list(itermax = 500)
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
      if (!is.null(start_arg)) args_list$start <- start_arg
      if (!is.null(lower_arg)) args_list$lower <- lower_arg
      if (!is.null(upper_arg)) args_list$upper <- upper_arg
      if (!is.null(control_arg)) args_list$control <- control_arg
      
      fit_try <- try(
        suppressWarnings(
          do.call(cesEst, args_list)
          ),
        silent = TRUE
      )
      
      runtime <- as.numeric(difftime(Sys.time(), t0, units="secs"))
      times_all[[m]] <- runtime
      
      success <- inherits(fit_try,"cesEst")
      conv_flag <- if (success && !is.null(fit_try$convergence)) fit_try$convergence else FALSE
      msg_flag <- if (inherits(fit_try,"try-error")) as.character(fit_try)[1]
      else if (success && !is.null(fit_try$message)) fit_try$message else ""
      
      fit_all[[m]] <- if (success) fit_try else NULL
      conv_all[[m]] <- conv_flag
      msg_all[[m]] <- msg_flag
      
      base::message(sprintf("  %s %-10s in %.1fs%s",
                            if (conv_flag) "✓" else "✗", m, runtime,
                            if (!conv_flag && nzchar(msg_flag)) paste0(" msg: ", msg_flag) else ""))
    }
    
    list(fits = fit_all, conv = conv_all, msg = msg_all, times = times_all, data = d_num)
}

  # ---- EXTRACT RESULTS ----
  # This function builds a parameter grid from the rhos tested, attaches coefficient estimates, flags run validity and produces a table per region/method/grid point
  extract_region <- function(region_name, region_fits) {
    # empty tibble helper 
    if (is.null(region_fits) || is.null(region_fits$data) || length(region_fits$fits) == 0L) {
      return(tibble(
        r = character(), 
        method = character(), 
        rho_KL = numeric(), 
        rho_VAE = numeric(),
        conv = logical(), 
        msg = character(), 
        rss = numeric(),
        gamma = numeric(), 
        lambda = numeric(), 
        delta_KVA = numeric(), 
        delta_VAY = numeric(), 
        nu = numeric(),
        se_gamma = numeric(), 
        se_lambda = numeric(), 
        se_delta_KVA = numeric(), 
        se_delta_VAY = numeric(), 
        se_nu = numeric(),
        t_gamma = numeric(), 
        t_lambda = numeric(), 
        t_delta_KVA = numeric(), 
        t_delta_VAY = numeric(), 
        t_nu = numeric(),
        p_gamma = numeric(), 
        p_lambda = numeric(), 
        p_delta_KVA = numeric(), 
        p_delta_VAY = numeric(), 
        p_nu = numeric(),
        ci_lo_gamma = numeric(), 
        ci_hi_gamma = numeric(),
        ci_lo_lambda = numeric(), 
        ci_hi_lambda = numeric(),
        ci_lo_delta_KVA = numeric(), 
        ci_hi_delta_KVA = numeric(),
        ci_lo_delta_VAY = numeric(), 
        ci_hi_delta_VAY = numeric(),
        ci_lo_nu = numeric(), 
        ci_hi_nu = numeric(),
        R2 = numeric(), 
        adjR2 = numeric(),
        AIC_naive = numeric(), 
        AICc_naive = numeric(), 
        AIC_plusRho = numeric(), 
        AICc_plusRho = numeric(),
        iter = numeric(),
        sigma_KL = numeric(), 
        sigma_VAE = numeric(),
        on_edge_KL = logical(), 
        on_edge_VAE = logical(),
        n_grid = integer(), 
        runtime_total = numeric(), 
        runtime_per_grid = numeric(),
        valid = logical()
      ))
    }
    
    d_num <- region_fits$data
    grid_tbl <- tibble()
    
    for (m in names(region_fits$fits)) {
      fit_obj   <- region_fits$fits[[m]]
      conv_flag <- safe_bool(region_fits$conv[[m]], FALSE)
      msg_flag  <- safe_chr(region_fits$msg[[m]], NA_character_)
      runtime   <- safe_num(region_fits$times[[m]], NA_real_)
      
      # base grid for this region × method
      full_grid <- expand_grid(rho_KL = rhoGrid_KL, rho_VAE = rhoGrid_VAE) %>%
        mutate(
          r = region_name,
          method = m,
          n_grid = n(),
          runtime_total = runtime,
          runtime_per_grid = runtime / n_grid,
          msg   = msg_flag,
          conv  = conv_flag,  # will be overridden by grid-level convergence if available
          sigma_KL  = ifelse(is.finite(1/(1+rho_KL)), 1/(1+rho_KL), NA_real_),
          sigma_VAE = ifelse(is.finite(1/(1+rho_VAE)), 1/(1+rho_VAE), NA_real_),
          on_edge_KL  = on_grid_edge(rho_KL,  rhoGrid_KL),
          on_edge_VAE = on_grid_edge(rho_VAE, rhoGrid_VAE),
          rss = NA_real_
        )
      
      # pre-fill parameter/stat columns as NA
      full_grid <- full_grid %>%
        mutate(
          gamma = NA_real_, lambda = NA_real_,
          delta_KVA = NA_real_, delta_VAY = NA_real_, nu = NA_real_,
          se_gamma = NA_real_, se_lambda = NA_real_,
          se_delta_KVA = NA_real_, se_delta_VAY = NA_real_, se_nu = NA_real_,
          t_gamma = NA_real_, t_lambda = NA_real_,
          t_delta_KVA = NA_real_, t_delta_VAY = NA_real_, t_nu = NA_real_,
          p_gamma = NA_real_, p_lambda = NA_real_,
          p_delta_KVA = NA_real_, p_delta_VAY = NA_real_, p_nu = NA_real_,
          ci_lo_gamma = NA_real_, ci_hi_gamma = NA_real_,
          ci_lo_lambda = NA_real_, ci_hi_lambda = NA_real_,
          ci_lo_delta_KVA = NA_real_, ci_hi_delta_KVA = NA_real_,
          ci_lo_delta_VAY = NA_real_, ci_hi_delta_VAY = NA_real_,
          ci_lo_nu = NA_real_, ci_hi_nu = NA_real_,
          R2 = NA_real_, adjR2 = NA_real_,
          AIC_naive = NA_real_, AICc_naive = NA_real_,
          AIC_plusRho = NA_real_, AICc_plusRho = NA_real_,
          iter = NA_real_
        )
      
      if (inherits(fit_obj, "cesEst")) {
        # --- 1. RSS and convergence per grid from allRhoSum ---
        if (!is.null(fit_obj$allRhoSum) && nrow(fit_obj$allRhoSum) > 0L) {
          sum_tbl <- fit_obj$allRhoSum
          # assume columns rho1 & rho & convergence exist for nested case
          sum_tbl <- sum_tbl %>%
            mutate(
              rho1 = as.numeric(.data[["rho1"]]),
              rho  = as.numeric(.data[["rho"]]),
              rss  = as.numeric(.data[["rss"]]),
              conv_grid = safe_bool(.data[["convergence"]], TRUE)
            ) %>%
            select(rho1, rho, rss, conv_grid) %>%
            rename(rho_KL = rho1, rho_VAE = rho)
          
          full_grid <- full_grid %>%
            left_join(sum_tbl, by = c("rho_KL", "rho_VAE"), suffix = c("", ".sum")) %>%
            mutate(
              rss  = coalesce(rss.sum, rss),
              conv = ifelse(is.na(conv_grid), conv, conv_grid)
            ) %>%
            select(-rss.sum, -conv_grid)
        }
        
        # --- 2. Per-grid coefficients from allRhoFull ---
        coef_df <- NULL
        if (!is.null(fit_obj$allRhoFull) && length(fit_obj$allRhoFull) > 0L) {
          coef_df <- purrr::map_dfr(fit_obj$allRhoFull, extract_grid_coeffs)
          
          full_grid <- full_grid %>%
            left_join(coef_df, by = c("rho_KL","rho_VAE"), suffix = c("", ".sub")) %>%
            mutate(
              gamma     = coalesce(gamma.sub,     gamma),
              lambda    = coalesce(lambda.sub,    lambda),
              delta_KVA = coalesce(delta_KVA.sub, delta_KVA),
              delta_VAY = coalesce(delta_VAY.sub, delta_VAY),
              nu        = coalesce(nu.sub,        nu),
              se_gamma     = coalesce(se_gamma.sub,     se_gamma),
              se_lambda    = coalesce(se_lambda.sub,    se_lambda),
              se_delta_KVA = coalesce(se_delta_KVA.sub, se_delta_KVA),
              se_delta_VAY = coalesce(se_delta_VAY.sub, se_delta_VAY),
              se_nu        = coalesce(se_nu.sub,        se_nu),
              t_gamma      = coalesce(t_gamma.sub,      t_gamma),
              t_lambda     = coalesce(t_lambda.sub,     t_lambda),
              t_delta_KVA  = coalesce(t_delta_KVA.sub,  t_delta_KVA),
              t_delta_VAY  = coalesce(t_delta_VAY.sub,  t_delta_VAY),
              t_nu         = coalesce(t_nu.sub,         t_nu),
              p_gamma      = coalesce(p_gamma.sub,      p_gamma),
              p_lambda     = coalesce(p_lambda.sub,     p_lambda),
              p_delta_KVA  = coalesce(p_delta_KVA.sub,  p_delta_KVA),
              p_delta_VAY  = coalesce(p_delta_VAY.sub,  p_delta_VAY),
              p_nu         = coalesce(p_nu.sub,         p_nu),
              ci_lo_gamma     = coalesce(ci_lo_gamma.sub,     ci_lo_gamma),
              ci_hi_gamma     = coalesce(ci_hi_gamma.sub,     ci_hi_gamma),
              ci_lo_lambda    = coalesce(ci_lo_lambda.sub,    ci_lo_lambda),
              ci_hi_lambda    = coalesce(ci_hi_lambda.sub,    ci_hi_lambda),
              ci_lo_delta_KVA = coalesce(ci_lo_delta_KVA.sub, ci_lo_delta_KVA),
              ci_hi_delta_KVA = coalesce(ci_hi_delta_KVA.sub, ci_hi_delta_KVA),
              ci_lo_delta_VAY = coalesce(ci_lo_delta_VAY.sub, ci_lo_delta_VAY),
              ci_hi_delta_VAY = coalesce(ci_hi_delta_VAY.sub, ci_hi_delta_VAY),
              ci_lo_nu        = coalesce(ci_lo_nu.sub,        ci_lo_nu),
              ci_hi_nu        = coalesce(ci_hi_nu.sub,        ci_hi_nu)
            ) %>%
            select(-ends_with(".sub"))
        }
        
        # --- 3. R² and AIC per grid cell (using rss) ---
        obs_log <- try(log(d_num$Ys + 1e-12), silent = TRUE)
        if (!inherits(obs_log, "try-error")) {
          n_obs   <- length(obs_log)
          TSS_log <- sum((obs_log - mean(obs_log))^2)
          k_hat   <- length(tryCatch(stats::coef(fit_obj), error = function(e) numeric(0)))
          k0      <- k_hat + rho_penalty
          
          full_grid <- full_grid %>%
            mutate(
              R2 = ifelse(is.finite(rss) & TSS_log > 0,
                          1 - rss / TSS_log,
                          NA_real_),
              adjR2 = ifelse(
                is.finite(R2) & (n_obs - k_hat - 1) > 0,
                1 - (1 - R2) * (n_obs - 1) / (n_obs - k_hat - 1),
                NA_real_
              ),
              AIC_naive = ifelse(
                is.finite(rss),
                n_obs * log(rss / n_obs) + 2 * k_hat,
                NA_real_
              ),
              AICc_naive = ifelse(
                is.finite(AIC_naive) & (n_obs - k_hat - 1) > 0,
                AIC_naive + (2 * k_hat * (k_hat + 1)) / (n_obs - k_hat - 1),
                NA_real_
              ),
              AIC_plusRho = ifelse(
                is.finite(rss),
                n_obs * log(rss / n_obs) + 2 * k0,
                NA_real_
              ),
              AICc_plusRho = ifelse(
                is.finite(AIC_plusRho) & (n_obs - k0 - 1) > 0,
                AIC_plusRho + (2 * k0 * (k0 + 1)) / (n_obs - k0 - 1),
                NA_real_
              )
            )
        }
        
        # --- 4. iterations (method-level) ---
        full_grid$iter <- iter_safe(fit_obj)
      }
      
      grid_tbl <- bind_rows(grid_tbl, full_grid)
    }
    
    # cast numerics
    grid_tbl %>%
      mutate(
        across(c(
          delta_KVA, delta_VAY, gamma, nu, lambda, sigma_KL, sigma_VAE,
          se_gamma, se_lambda, se_delta_KVA, se_delta_VAY, se_nu,
          t_gamma, t_lambda, t_delta_KVA, t_delta_VAY, t_nu,
          p_gamma, p_lambda, p_delta_KVA, p_delta_VAY, p_nu,
          ci_lo_gamma, ci_hi_gamma, ci_lo_lambda, ci_hi_lambda,
          ci_lo_delta_KVA, ci_hi_delta_KVA, ci_lo_delta_VAY, ci_hi_delta_VAY,
          ci_lo_nu, ci_hi_nu,
          rss, R2, adjR2, AIC_naive, AICc_naive, AIC_plusRho, AICc_plusRho,
          iter, rho_KL, rho_VAE, runtime_total, runtime_per_grid, n_grid
        ), ~ suppressWarnings(as.numeric(.)))
      )
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
  .progress = TRUE,
  .options = furrr_options(
    packages = c("micEconCES","dplyr","tidyr","purrr","readr"),
    globals = c("rhoGrid_KL","rhoGrid_VAE","estimate_region","extract_region"),
    seed = TRUE
  )
)

# Extract grid results
results_grid <- map2_dfr(region_names, fits_all, extract_region)

# ---- SAVE MASTER TABLE ----
write_csv(results_grid, "CES_results_grid.csv")
saveRDS(results_grid, "results_run1.rds")
} else {
  if (file.exists("results_run1.rds")) {
    message("Loading results from results_run1.rds")
    results_grid <- readRDS("results_run1.rds")
  } else if (file.exists("CES_results_grid.csv")) {
    message("Loading results from CES_results_grid.csv")
    results_grid <- read_csv("CES_results_grid.csv", show_col_types = FALSE)
  } else {
    stop("No saved results found. Run with RUN_ESTIMATION = TRUE first.")
  }
  
  # Load input data so dfS is available for IAM table
  df <- read_csv(infile, show_col_types = FALSE)
  
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
}

# Convergence summary by region
convergence_summary <- results_grid %>%
  select(r, method, conv) %>%
  distinct() %>%
  count(method, conv, name = "count") %>%
  mutate(status = ifelse(conv, "Converged", "Failed"))

# Adding validity to the runs, economically feasible and converged from the solver
# Adding validity to the runs, with detailed status diagnostics
add_validity <- function(df) {
  df %>%
    mutate(
      across(
        any_of(c("delta_KVA","delta_VAY","gamma","nu","lambda")),
        ~ suppressWarnings(as.numeric(.))
      ),
      
      valid_reason = pmap_chr(
        list(conv, rss, R2, delta_KVA, delta_VAY, gamma, nu, lambda),
        function(conv, rss, R2, dK, dV, g, n, lam) {
          reasons <- c()
          
          if (!isTRUE(conv))                     reasons <- c(reasons, "Solver did not converge")
          if (!is.finite(rss) || rss <= 0)       reasons <- c(reasons, "RSS invalid")
          if (!is.finite(R2)  || R2 <= 0 || R2 > 1) reasons <- c(reasons, "R2 invalid")
          if (!is.finite(dK)  || dK < 0 || dK > 1) reasons <- c(reasons, "dK-VA out of [0,1]")
          if (!is.finite(dV)  || dV < 0 || dV > 1) reasons <- c(reasons, "dVA-Y out of [0,1]")
          if (!is.finite(g)   || g < 0.5 || g > 3) reasons <- c(reasons, "gamma out of [0.5,3]")
          if (!is.finite(n)   || n < 0.7 || n > 1.3) reasons <- c(reasons, "v out of [0.7,1.3]")
          if (!is.finite(lam) || lam < -0.05 || lam > 0.05) reasons <- c(reasons, "lambda out of [-0.05,0.05]")
          
          if (length(reasons) == 0) "OK" else paste(reasons, collapse = "; ")
        }
      ),
      
      valid = valid_reason == "OK",
      
      solver_reason = msg_group_map(msg),
      
      status = case_when(
        valid                 ~ "Valid",
        !conv                 ~ solver_reason,
        valid_reason != "OK"  ~ valid_reason,
        TRUE                  ~ "Unspecified"
      ),
      
      status = factor(
        status,
        levels = c(
          "False convergence","Max iterations","Reduction criterion",
          "Bounds/tolerance","Unspecified",
          "RSS invalid","R2 invalid",
          "dK-VA out of [0,1]","dVA-Y out of [0,1]",
          "gamma out of [0.5,3]","v out of [0.7,1.3]","lambda out of [-0.05,0.05]",
          "Valid"
        )
      )
    )
}

  

results_grid <- results_grid %>% select(-any_of("valid")) %>% add_validity()

# Warning messages
message("Finite gamma:      ", sum(is.finite(results_grid$gamma)))
message("Finite delta_KVA:  ", sum(is.finite(results_grid$delta_KVA)))
message("Finite delta_VAY:  ", sum(is.finite(results_grid$delta_VAY)))
message("Valid grid points: ", sum(results_grid$valid))

if (all(is.na(results_grid$gamma))) {
  warning("All gamma are NA – something is still wrong with coefficient extraction.")
}
if (sum(results_grid$valid) == 0) {
  warning("No grid point passed validity – consider relaxing thresholds.")
}


# Valid runs table
results_grid_valid <- results_grid %>% filter(valid)

# Invalid runs table
results_grid_invalid <- results_grid %>% filter(!valid)

# Best method per region
# Strict criteria
best_methods <- results_grid %>%
  mutate(validity = ifelse(valid, "valid", "invalid")) %>% 
  filter(valid,
         is.finite(R2), R2 > 0.7,
         is.finite(delta_KVA) & delta_KVA >= 0.2 & delta_KVA <= 0.8 &
         is.finite(delta_VAY) & delta_VAY >= 0.2 & delta_VAY <= 0.8 &
         rho_KL < 50 & rho_VAE < 50
         ) %>%
  group_by(r) %>%
  arrange(AICc_plusRho, rss, desc(R2)) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(best_tier = "strict")

# Relaxed criteria, valid runs regardless of goodness of fit
still_missing <- setdiff(unique(results_grid$r), unique(best_methods$r))
if (length(still_missing) > 0) {
  best_relaxed <- results_grid %>%
    mutate(validity = ifelse(valid, "valid", "invalid")) %>%
    filter(r %in% still_missing, valid) %>%
    group_by(r) %>%
    filter(!on_edge_KL & !on_edge_VAE) %>%
    arrange(AICc_plusRho, rss, desc(R2)) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    mutate(best_tier = "relaxed")
  
  best_methods <- bind_rows(best_methods, best_relaxed) %>%
    distinct(r, .keep_all = TRUE)
}

# Forced criteria, first near-valid ranges and if not, arranged by goodness of fit regardless of range
still_missing <- setdiff(unique(results_grid$r), unique(best_methods$r))
if (length(still_missing) > 0) {
  best_forced <- results_grid %>%
    mutate(validity = ifelse(valid, "valid", "invalid")) %>%
    filter(r %in% still_missing) %>%
    group_by(r) %>%
    # Step 1: prioritise "near-valid" (parameters just outside thresholds)
    mutate(
      near_valid = (
        between(delta_KVA, 0.02, 0.98) &
          between(delta_VAY, 0.02, 0.98) &
          between(gamma, 0.1, 100) &
          between(nu, 0.3, 10) &
          between(lambda, -0.5, 0.5)
      )
    ) %>%
    arrange(desc(near_valid), AICc_plusRho, rss, desc(R2), runtime_total) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    mutate(best_tier = "forced")

  best_methods <- bind_rows(best_methods, best_forced) %>%
    distinct(r, .keep_all = TRUE)
}


# Median valid parameters
results_grid_valid_summary <- results_grid %>%
  mutate(flag = case_when(
    paste(r, method) %in% paste(best_methods$r, best_methods$method) ~ "best",
    valid ~ "valid",
    TRUE ~ "invalid"
  )) %>%
  group_by(r, method, flag) %>%
  summarise(
    n_runs   = n(),
    # RSS / R² ranges
    min_RSS  = min(rss, na.rm = TRUE),
    max_RSS  = max(rss, na.rm = TRUE),
    med_RSS  = median(rss, na.rm = TRUE),
    min_R2   = min(R2, na.rm = TRUE),
    max_R2   = max(R2, na.rm = TRUE),
    med_R2   = median(R2, na.rm = TRUE),
    # Parameters (medians across runs)
    gamma_med     = median(gamma, na.rm = TRUE),
    lambda_med    = median(lambda, na.rm = TRUE),
    delta_KVA_med  = median(delta_KVA, na.rm = TRUE),
    delta_VAY_med = median(delta_VAY, na.rm = TRUE),
    nu_med        = median(nu, na.rm = TRUE),
    sigma_KL_med  = median(sigma_KL, na.rm = TRUE),
    sigma_VAE_med = median(sigma_VAE, na.rm = TRUE),
    # Iterations and runtime
    med_iter      = median(iter, na.rm = TRUE),
    med_runtime   = median(runtime_total, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(r, method, flag)



# Robustness summary by region
robustness_summary <- results_grid %>%
  group_by(r) %>%
  summarise(
    n_grid_total = n(),
    n_valid = sum(valid, na.rm = TRUE),
    share_valid = n_valid / n_grid_total,
    gamma_min = if (any(valid & is.finite(gamma))) min(gamma[valid], na.rm = TRUE) else NA_real_,
    gamma_max = if (any(valid & is.finite(gamma))) max(gamma[valid], na.rm = TRUE) else NA_real_,
    deltaKL_min = if (any(valid & is.finite(delta_KVA))) min(delta_KVA[valid], na.rm = TRUE) else NA_real_,
    deltaKL_max = if (any(valid & is.finite(delta_KVA))) max(delta_KVA[valid], na.rm = TRUE) else NA_real_,
    deltaVAE_min = if (any(valid & is.finite(delta_VAY))) min(delta_VAY[valid], na.rm = TRUE) else NA_real_,
    deltaVAE_max = if (any(valid & is.finite(delta_VAY))) max(delta_VAY[valid], na.rm = TRUE) else NA_real_,
    nu_min = if (any(valid & is.finite(nu))) min(nu[valid], na.rm = TRUE) else NA_real_,
    nu_max = if (any(valid & is.finite(nu))) max(nu[valid], na.rm = TRUE) else NA_real_,
    sigmaKL_min = if (any(valid & is.finite(sigma_KL))) min(sigma_KL[valid], na.rm = TRUE) else NA_real_,
    sigmaKL_max = if (any(valid & is.finite(sigma_KL))) max(sigma_KL[valid], na.rm = TRUE) else NA_real_,
    sigmaVAE_min = if (any(valid & is.finite(sigma_VAE))) min(sigma_VAE[valid], na.rm = TRUE) else NA_real_,
    sigmaVAE_max = if (any(valid & is.finite(sigma_VAE))) max(sigma_VAE[valid], na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>%
  left_join(
    best_methods %>%
      select(r, method, best_tier, gamma, delta_KVA, delta_VAY, nu, sigma_KL, sigma_VAE),
    by = "r"
  )



# AICc weights (among valid only)
aic_weights <- results_grid %>%
  filter(valid, is.finite(AICc_plusRho)) %>%
  group_by(r) %>%
  mutate(
    dAICc = AICc_plusRho - min(AICc_plusRho, na.rm = TRUE),
    wAICc = exp(-0.5 * dAICc) / sum(exp(-0.5 * dAICc), na.rm = TRUE)
  ) %>%
  ungroup()


# Share of grid converged
grid_conv_share <- results_grid %>%
  group_by(r, method) %>%
  summarise(
    grid_points = n(),
    n_converged = sum(conv, na.rm = TRUE),
    share_converged = n_converged/grid_points,
    .groups = "drop"
  )


# IAM parameter table (join best methods)
years_by_region <- dfS %>% distinct(r, t)
iam_table <- best_methods %>%
  inner_join(years_by_region, by = "r") %>%
  transmute(
    year = t,
    region = r,
    total_factor_productivity = gamma * exp(lambda * t),
    share_capital_valueadded = delta_KVA,
    share_valueadded_output = delta_VAY,
    elasticity_substitution_KL = sigma_KL,
    elasticity_substitution_VAE = sigma_VAE,
    exponent_rho_KL = rho_KL,
    exponent_rho_VAE = rho_VAE,
    TFP_scale_factor = gamma,
    returns_scale = nu,
    TFP_growth_rate = lambda,
    p_gamma, p_lambda, p_delta_KVA, p_delta_VAY, p_nu
  )



# ---- EXPORT derived artifacts ----
#write_csv(convergence_summary, "CES_convergence_summary.csv")
write_csv(best_methods, "CES_best_methods.csv")
write_csv(results_grid, "CES_results_grid.csv")
#write_csv(results_grid_valid, "CES_results_grid_valid.csv")
#write_csv(results_grid_invalid, "CES_results_grid_invalid.csv")
write_csv(robustness_summary, "CES_robustness_summary.csv")
#write_csv(aic_weights, "CES_AICc_weights.csv")
write_csv(grid_conv_share, "CES_grid_convergence_share.csv")
write_csv(iam_table, "IAM_params.csv")
