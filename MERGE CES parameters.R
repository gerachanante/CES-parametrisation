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


setwd("C:/Users/escami_g/OneDrive - Paul Scherrer Institut/05.Models/MERGE updates/CES-parametrisation/-.9 to 1 by 0.15, 1.5 to 5 by 0.5, 10, 50, 100")
infile <- "MERGE macro.csv"

# TRUE to run a fresh estimation or FALSE to reuse saved results (.rds file)
RUN_ESTIMATION <- TRUE

# Full grid
#rhoGrid_KL <- c(seq(-0.86, 0.4, by = 0.02), seq(0.5, 3, by = 0.1), 4, 5, 10, 15, 20, 30, 50, 100)
#rhoGrid_VAE <- c(seq(-0.8, 0.4, by = 0.02), seq(0.5, 3, by = 0.1), 4, 5, 10, 15, 20, 30, 50, 100)

# Test grid
rhoGrid_KL <- c(seq(-0.9, 1, by = 0.15), seq(1.5, 5, by = 0.5), 10, 50, 100)
rhoGrid_VAE <- c(seq(-0.9, 1, by = 0.15), seq(1.5, 5, by = 0.5), 10, 50, 100)

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
    # Fast methods: NM, Nelder-Mead, Newton, L-BFGS-B
    # Slow methods: LM, PORT, BFGS
    # Terribly slow: SANN, DE, CG
    #opt_methods <- unique(c("LM", "NM", "Nelder-Mead", "PORT", "Newton", "L-BFGS-B"))
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
          delta_KL = runif(1, 0.4, 0.6),
          delta_VAE = runif(1, 0.4, 0.6),
          nu = runif(1, 0.9, 1.1)
        )
        # no lower/upper/control
      }
      if (m == "L-BFGS-B") {
        start_arg <- c(gamma = 1, lambda = 0.001, delta_KL = 0.5, delta_VAE = 0.5, nu = 1)
        lower_arg <- c(gamma = 0.1, lambda = -0.3, delta_KL = 0.1, delta_VAE = 0.1, nu = 0.3)
        upper_arg <- c(gamma = 10, lambda = 0.3, delta_KL = 0.9, delta_VAE = 0.9, nu = 5)
        control_arg <- list(maxit = 1000, factr = 1e9)
      }
      if (m == "PORT") control_arg <- list(eval.max=1e5, iter.max=2e4, reltol=1e-8)
      if (m == "BFGS") {
        start_arg <- c(gamma = 1, lambda = 0.001, delta_KL = 0.5, delta_VAE = 0.5, nu = 1)
        lower_arg <- c(delta_KL = 0.1, delta_VAE = 0.1)
        upper_arg <- c(delta_KL = 0.9, delta_VAE = 0.9)
        control_arg <- list(maxit = 1000, reltol = 1e-8)
      }
      if (m == "CG") {
        start_arg <- c(gamma = 1, lambda = 0.001, delta_KL = 0.5, delta_VAE = 0.5, nu = 1)
        lower_arg <- c(delta_KL = 0.1, delta_VAE = 0.1)
        upper_arg <- c(delta_KL = 0.9, delta_VAE = 0.9)
        control_arg <- list(maxit = 500, reltol = 1e-8)
      }    
      if (m %in% c("NM","Nelder-Mead")) {
        lower_arg <- c(delta_KL = 0.1, delta_VAE = 0.1)
        upper_arg <- c(delta_KL = 0.9, delta_VAE = 0.9)
        control_arg <- list(maxit = 5000, reltol = 1e-8)
      }
      if (m == "LM") {
        lower_arg <- c(delta_KL = 0.1, delta_VAE = 0.1)
        upper_arg <- c(delta_KL = 0.9, delta_VAE = 0.9)
        control_arg <- list(maxiter = 10000, ftol = 1e-8, maxfev = 5000)
      }
      if (m == "SANN") {
        lower_arg <- c(delta_KL = 0.1, delta_VAE = 0.1)
        upper_arg <- c(delta_KL = 0.9, delta_VAE = 0.9)
        control_arg <- list(maxit = 5000, temp = 10, tmax = 50)
      }
      if (m == "DE") {
        lower_arg <- c(gamma = 0.1, lambda = -0.3, delta_KL = 0.1, delta_VAE = 0.1, nu = 0.3)
        upper_arg <- c(gamma = 10, lambda = 0.3, delta_KL = 0.9, delta_VAE = 0.9, nu = 5)
        control_arg <- list(itermax = 100)
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
    if (is.null(region_fits) || is.null(region_fits$data) || length(region_fits$fits) == 0L) {
      return(tibble(
        r = character(), method = character(), rho_KL = numeric(), rho_VAE = numeric(),
        conv = logical(), msg = character(), rss = numeric(),
        gamma = numeric(), lambda = numeric(), delta_KL = numeric(), delta_VAE = numeric(), nu = numeric(),
        se_gamma = numeric(), se_lambda = numeric(), se_delta_KL = numeric(), se_delta_VAE = numeric(), se_nu = numeric(),
        t_gamma = numeric(), t_lambda = numeric(), t_delta_KL = numeric(), t_delta_VAE = numeric(), t_nu = numeric(),
        p_gamma = numeric(), p_lambda = numeric(), p_delta_KL = numeric(), p_delta_VAE = numeric(), p_nu = numeric(),
        ci_lo_gamma = numeric(), ci_hi_gamma = numeric(),
        ci_lo_lambda = numeric(), ci_hi_lambda = numeric(),
        ci_lo_delta_KL = numeric(), ci_hi_delta_KL = numeric(),
        ci_lo_delta_VAE = numeric(), ci_hi_delta_VAE = numeric(),
        ci_lo_nu = numeric(), ci_hi_nu = numeric(),
        R2 = numeric(), adjR2 = numeric(),
        AIC_naive = numeric(), AICc_naive = numeric(), AIC_plusRho = numeric(), AICc_plusRho = numeric(),
        iter = numeric(),
        sigma_KL = numeric(), sigma_VAE = numeric(),
        on_edge_KL = logical(), on_edge_VAE = logical(),
        n_grid = integer(), runtime_total = numeric(), runtime_per_grid = numeric(),
        valid = logical()
      ))
    }
    
    d_num <- region_fits$data
    grid_tbl <- tibble()
    
    for (m in names(region_fits$fits)) {
      fit_obj <- region_fits$fits[[m]]
      conv_flag <- safe_bool(region_fits$conv[[m]], FALSE)
      msg_flag <- safe_chr(region_fits$msg[[m]], NA_character_)
      runtime <- safe_num(region_fits$times[[m]], NA_real_)

      # Build the full grid
      full_grid <- expand_grid(rho_KL = rhoGrid_KL, rho_VAE = rhoGrid_VAE) %>%
        mutate(
          r = region_name,
          method = m,
          conv = conv_flag,
          msg = msg_flag,
          n_grid = nrow(expand_grid(rho_KL = rhoGrid_KL, rho_VAE = rhoGrid_VAE)),
          runtime_total = runtime,
          runtime_per_grid = runtime / n_grid,
          sigma_KL = ifelse(is.finite(1/(1+rho_KL)), 1/(1+rho_KL), NA_real_),
          sigma_VAE = ifelse(is.finite(1/(1+rho_VAE)),  1/(1+rho_VAE),  NA_real_),
          on_edge_KL = on_grid_edge(rho_KL, rhoGrid_KL),
          on_edge_VAE = on_grid_edge(rho_VAE,  rhoGrid_VAE),
          
          rss = NA_real_,
          gamma = NA_real_, lambda = NA_real_, delta_KL = NA_real_, delta_VAE = NA_real_, nu = NA_real_,
          se_gamma = NA_real_, se_lambda = NA_real_, se_delta_KL = NA_real_, se_delta_VAE = NA_real_, se_nu = NA_real_,
          t_gamma = NA_real_, t_lambda = NA_real_, t_delta_KL = NA_real_, t_delta_VAE = NA_real_, t_nu = NA_real_,
          p_gamma = NA_real_, p_lambda = NA_real_, p_delta_KL = NA_real_, p_delta_VAE = NA_real_, p_nu = NA_real_,
          ci_lo_gamma = NA_real_, ci_hi_gamma = NA_real_,
          ci_lo_lambda = NA_real_, ci_hi_lambda = NA_real_,
          ci_lo_delta_KL = NA_real_, ci_hi_delta_KL = NA_real_,
          ci_lo_delta_VAE = NA_real_, ci_hi_delta_VAE = NA_real_,
          ci_lo_nu = NA_real_, ci_hi_nu = NA_real_,
          R2 = NA_real_, adjR2 = NA_real_,
          AIC_naive = NA_real_, AICc_naive = NA_real_,
          AIC_plusRho = NA_real_, AICc_plusRho = NA_real_,
          iter = NA_real_
        )
      
      # If the fit succeeded, merge results
      if (inherits(fit_obj, "cesEst")) {
        # RSS per grid combination
        if (!is.null(fit_obj$allRhoSum) && nrow(fit_obj$allRhoSum) > 0) {
          allRho_tbl <- fit_obj$allRhoSum %>%
            select(rho1, rho, rss) %>%
            rename(rho_KL = rho1, rho_VAE = rho) %>%
            suppressWarnings(drop_na(rss))
          if (nrow(allRho_tbl) > 0) {
            full_grid <- full_grid %>%
              left_join(allRho_tbl, by = c("rho_KL","rho_VAE"), suffix = c("", ".fit")) %>%
              mutate(rss = coalesce(rss.fit, rss)) %>%
              select(-rss.fit)
          }
        }
        
        # Coefficients, SE,t,p,R2,AIC,iter attached to grid combinations from each method.
        cm <- NULL
        s1 <- try(summary(fit_obj), silent = TRUE)
        s2 <- try(summary(fit_obj, ela = TRUE), silent = TRUE)
        
        # grab the first available coef-like table
        for (s in list(s1, s2)) {
          if (inherits(s, "try-error") || is.null(s)) next
          for (slot in c("coefficients", "coef", "coefTable")) {
            if (!is.null(s[[slot]])) {
              tmp <- try(as.data.frame(s[[slot]]), silent = TRUE)
              if (!inherits(tmp, "try-error") && nrow(tmp) > 0) { cm <- tmp; break }
            }
          }
          if (!is.null(cm)) break
        }
        # sometimes summary returns a matrix directly
        if (is.null(cm) && !inherits(s1, "try-error") && is.matrix(s1)) cm <- as.data.frame(s1)
        if (is.null(cm) && !inherits(s2, "try-error") && is.matrix(s2)) cm <- as.data.frame(s2)
        
        est <- se <- tval <- pval <- list(gamma=NA_real_, lambda=NA_real_, delta_KL=NA_real_, delta_VAE=NA_real_, nu=NA_real_)
        
        if (!is.null(cm) && nrow(cm) > 0) {
          rn <- rownames(cm); cn <- colnames(cm)
          
          # flexible column pickers (case-insensitive)
          pick_col <- function(patterns) {
            ix <- which(vapply(cn, function(z) any(grepl(patterns, z, ignore.case = TRUE)), logical(1)))
            if (length(ix) == 0) NA_integer_ else ix[1]
          }
          col_est <- pick_col("^(estimate|coef|value)$|^estimate$|^coef$|^coeff")
          col_se <- pick_col("(std\\.? ?error|se)")
          col_t <- pick_col("^(t.?value|z|t)$")
          col_p <- pick_col("^(pr\\(|p.?value|p$)")
          
          # flexible row pickers
          pick_row <- function(patterns, exclude=NULL) {
            ok <- which(vapply(rn, function(z) any(grepl(patterns, z, ignore.case = TRUE)), logical(1)))
            if (!is.null(exclude)) ok <- ok[!vapply(rn[ok], function(z) any(grepl(exclude, z, ignore.case = TRUE)), logical(1))]
            if (length(ok) == 0) NA_integer_ else ok[1]
          }
          r_gamma <- pick_row("^gamma$")
          r_lambda <- pick_row("^lambda|^lam$")
          r_delta1 <- pick_row("^delta[_ ]?1$|^delta-?1$|delta[_]?kl")
          r_deltaMain <- pick_row("^delta$", exclude="1")  # avoid delta_1
          r_nu <- pick_row("^nu$")
          
          get_val <- function(ri, ci) {
            if (is.na(ri) || is.na(ci)) return(NA_real_)
            v <- suppressWarnings(as.numeric(cm[ri, ci, drop=TRUE]))
            if (is.finite(v)) v else NA_real_
          }
          
          # estimates
          est$gamma <- get_val(r_gamma, col_est)
          est$lambda <- get_val(r_lambda, col_est)
          est$delta_KL <- get_val(r_delta1, col_est)
          est$delta_VAE <- get_val(r_deltaMain, col_est)
          est$nu <- get_val(r_nu, col_est)
          
          # Standard errors
          se$gamma <- get_val(r_gamma, col_se)
          se$lambda <- get_val(r_lambda, col_se)
          se$delta_KL <- get_val(r_delta1, col_se)
          se$delta_VAE <- get_val(r_deltaMain, col_se)
          se$nu <- get_val(r_nu, col_se)
          
          # t-tests
          tval$gamma <- get_val(r_gamma, col_t)
          tval$lambda <- get_val(r_lambda, col_t)
          tval$delta_KL <- get_val(r_delta1, col_t)
          tval$delta_VAE <- get_val(r_deltaMain, col_t)
          tval$nu <- get_val(r_nu, col_t)
          
          # p-values
          pval$gamma <- get_val(r_gamma, col_p)
          pval$lambda <- get_val(r_lambda, col_p)
          pval$delta_KL <- get_val(r_delta1, col_p)
          pval$delta_VAE <- get_val(r_deltaMain, col_p)
          pval$nu <- get_val(r_nu, col_p)
        }
        
        # Fallback from vcov if needed (compute t, p)
        if (any(!is.finite(unlist(se)))) {
          vc <- try(vcov(fit_obj), silent = TRUE)
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
                if (!is.finite(se$gamma)) se$gamma <- get_se("^gamma$")
                if (!is.finite(se$lambda)) se$lambda <- get_se("^lambda|^lam$")
                if (!is.finite(se$delta_KL)) se$delta_KL <- get_se("^delta[_ ]?1$|^delta-?1$|delta[_]?kl")
                if (!is.finite(se$delta_VAE)) se$delta_VAE <- get_se("^delta$")
                if (!is.finite(se$nu)) se$nu <- get_se("^nu$")
              }
            }
          }
        }
        
        # t from est/SE when missing
        for (nm in names(est)) {
          if (!is.finite(tval[[nm]]) && is.finite(est[[nm]]) && is.finite(se[[nm]]) && se[[nm]] > 0)
            tval[[nm]] <- est[[nm]]/se[[nm]]
        }
        # p from t (normal)
        for (nm in names(tval)) {
          if (!is.finite(pval[[nm]]) && is.finite(tval[[nm]]))
            pval[[nm]] <- 2*(1 - pnorm(abs(tval[[nm]])))
        }
        
        # residual diagnostics in log space
        obs_log  <- try(log(d_num$Ys + 1e-12), silent = TRUE)
        fit_logv <- try(log(as.numeric(fit_obj$fitted.values) + 1e-12), silent = TRUE)
        R2_val <- adjR2_val <- AIC0 <- AICc0 <- AICR <- AICcR <- NA_real_
        if (!inherits(obs_log,"try-error") && !inherits(fit_logv,"try-error") && length(obs_log) == length(fit_logv)) {
          resid_log <- obs_log - fit_logv
          R2_val <- 1 - sum(resid_log^2) / sum((obs_log - mean(obs_log))^2)
          k_hat <- length(tryCatch(stats::coef(fit_obj), error = function(e) numeric(0)))
          n_obs <- length(obs_log)
          adjR2_val <- if (n_obs - k_hat - 1 > 0) 1 - (1 - R2_val)*(n_obs - 1)/(n_obs - k_hat - 1) else NA_real_
          ic0 <- aic_bic_from_rss(resid_log, k = k_hat, add_k_rho = 0L)
          icR <- aic_bic_from_rss(resid_log, k = k_hat, add_k_rho = rho_penalty)
          AIC0  <- ic0$AIC;  AICc0 <- ic0$AICc
          AICR  <- icR$AIC;  AICcR <- icR$AICc
        }
        
        # Iterations
        iter_val <- iter_safe(fit_obj)
        
        # Attach to all rows in the grid combo
        full_grid <- full_grid %>%
          mutate(
            # Estimated values
            gamma = est$gamma, 
            lambda = est$lambda, 
            delta_KL = est$delta_KL, 
            delta_VAE = est$delta_VAE, 
            nu = est$nu,
            # Standard errors
            se_gamma = se$gamma, 
            se_lambda = se$lambda, 
            se_delta_KL = se$delta_KL, 
            se_delta_VAE = se$delta_VAE, 
            se_nu = se$nu,
            # t-tests
            t_gamma = tval$gamma, 
            t_lambda = tval$lambda, 
            t_delta_KL = tval$delta_KL, 
            t_delta_VAE = tval$delta_VAE, 
            t_nu = tval$nu,
            # P-values
            p_gamma = pval$gamma, 
            p_lambda = pval$lambda, 
            p_delta_KL = pval$delta_KL,
            p_delta_VAE = pval$delta_VAE, 
            p_nu = pval$nu,
            # Confidence intervals
            ci_lo_gamma = ifelse(is.finite(gamma) & is.finite(se_gamma), gamma - 1.96*se_gamma, NA_real_),
            ci_hi_gamma = ifelse(is.finite(gamma) & is.finite(se_gamma), gamma + 1.96*se_gamma, NA_real_),
            ci_lo_lambda = ifelse(is.finite(lambda) & is.finite(se_lambda), lambda - 1.96*se_lambda, NA_real_),
            ci_hi_lambda = ifelse(is.finite(lambda) & is.finite(se_lambda), lambda + 1.96*se_lambda, NA_real_),
            ci_lo_delta_KL = ifelse(is.finite(delta_KL) & is.finite(se_delta_KL), delta_KL - 1.96*se_delta_KL, NA_real_),
            ci_hi_delta_KL = ifelse(is.finite(delta_KL) & is.finite(se_delta_KL), delta_KL + 1.96*se_delta_KL, NA_real_),
            ci_lo_delta_VAE = ifelse(is.finite(delta_VAE) & is.finite(se_delta_VAE), delta_VAE - 1.96*se_delta_VAE, NA_real_),
            ci_hi_delta_VAE = ifelse(is.finite(delta_VAE) & is.finite(se_delta_VAE), delta_VAE + 1.96*se_delta_VAE, NA_real_),
            ci_lo_nu = ifelse(is.finite(nu) & is.finite(se_nu), nu - 1.96*se_nu, NA_real_),
            ci_hi_nu = ifelse(is.finite(nu) & is.finite(se_nu), nu + 1.96*se_nu, NA_real_),
            # Goodness of fit
            R2 = R2_val, 
            adjR2 = adjR2_val,
            # Information criteria
            AIC_naive = AIC0, 
            AICc_naive = AICc0, 
            AIC_plusRho = AICR, 
            AICc_plusRho = AICcR,
            # Number of iterations
            iter = iter_val
          )
      }
      
      grid_tbl <- bind_rows(grid_tbl, full_grid)
    }
    
    
    # Final validity flags. These determine if a run is "valid", allowing it to be a candidate for best method run
    grid_tbl %>%
      mutate(
        across(c(
          delta_KL, delta_VAE, gamma, nu, lambda, sigma_KL, sigma_VAE,
          se_gamma, se_lambda, se_delta_KL, se_delta_VAE, se_nu,
          t_gamma, t_lambda, t_delta_KL, t_delta_VAE, t_nu,
          p_gamma, p_lambda, p_delta_KL, p_delta_VAE, p_nu,
          ci_lo_gamma, ci_hi_gamma, ci_lo_lambda, ci_hi_lambda,
          ci_lo_delta_KL, ci_hi_delta_KL, ci_lo_delta_VAE, ci_hi_delta_VAE,
          ci_lo_nu, ci_hi_nu,
          rss, R2, adjR2, AIC_naive, AICc_naive, AIC_plusRho, AICc_plusRho,
          iter, rho_KL, rho_VAE, runtime_total, runtime_per_grid, n_grid
        ), ~ suppressWarnings(as.numeric(.)))
      ) %>%
      rowwise() %>%
      mutate(
        valid = isTRUE(conv) && # only runs that converged
          !is.na(delta_KL)  && delta_KL  > 0.1 && delta_KL  < 0.9 && # deltas in the range
          !is.na(delta_VAE) && delta_VAE > 0.1 && delta_VAE < 0.9 && # deltas in the range
          !is.na(sigma_KL)  && sigma_KL  > 0.1 && sigma_KL  < 5   && # sigmas in the range
          !is.na(sigma_VAE) && sigma_VAE > 0.1 && sigma_VAE < 5   && # sigmas in the range
          !is.na(gamma)     && gamma     > 0.5 && gamma     < 5   && # gamma in the range
          !is.na(nu)        && nu        > 0.5 && nu        < 5   && # nu in the range
          !is.na(lambda)    && lambda    > -0.4 && lambda   < 0.4 # lambda in the range
      ) %>%
      ungroup()
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
  results_grid <- readRDS("results_run1.rds")
}

# Read the results when running the code without estimation.
results_grid <- if(file.exists("CES_results_grid.csv")) read_csv("CES_results_grid.csv", show_col_types = FALSE) else tibble()

# Convergence summary by region
convergence_summary <- results_grid %>%
  select(r, method, conv) %>%
  distinct() %>%
  count(method, conv, name = "count") %>%
  mutate(status = ifelse(conv, "Converged", "Failed"))

# Best method per region
# Strict criteria
best_methods <- results_grid %>%
  filter(valid) %>%
  group_by(r) %>%
  filter(!on_edge_KL | !on_edge_VAE | n() == 1) %>%
  arrange(AICc_plusRho, rss) %>%
  slice_head(n = 1) %>%
  mutate(best_tier = "strict valid") %>%
  ungroup()

# Regions not yet covered by strict criteria
still_missing <- setdiff(unique(results_grid$r), unique(best_methods$r))
if (length(still_missing)) {
  best_relaxed <- results_grid %>%
    filter(r %in% still_missing, valid) %>%
    arrange(AICc_plusRho, rss) %>%
    group_by(r) %>%
    slice_head(n = 1) %>%
    mutate(best_tier = "relaxed fallback") %>%
    ungroup()
  best_methods <- bind_rows(best_methods, best_relaxed) %>%
    distinct(r, .keep_all = TRUE)
}

# Valid runs table
results_grid_valid <- results_grid %>% filter(valid)

# Invalid runs table
results_grid_invalid <- results_grid %>% filter(!valid)


# Robustness summary by region
robustness_summary <- results_grid %>%
  filter(valid) %>%
  group_by(r) %>%
  summarise(
    n_grid_total = n(),
    n_valid = n(),
    share_valid = 1,
    gamma_min = min(gamma, na.rm = TRUE),  gamma_max = max(gamma, na.rm = TRUE),
    deltaKL_min = min(delta_KL, na.rm = TRUE), deltaKL_max = max(delta_KL, na.rm = TRUE),
    deltaVAE_min = min(delta_VAE, na.rm = TRUE), deltaVAE_max = max(delta_VAE, na.rm = TRUE),
    nu_min = min(nu, na.rm = TRUE), nu_max = max(nu, na.rm = TRUE),
    sigmaKL_min = min(sigma_KL, na.rm = TRUE), sigmaKL_max = max(sigma_KL, na.rm = TRUE),
    sigmaVAE_min = min(sigma_VAE, na.rm = TRUE), sigmaVAE_max = max(sigma_VAE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    best_methods %>% select(r, method, best_tier, gamma, delta_KL, delta_VAE, nu, sigma_KL, sigma_VAE),
    by = "r"
  )


# AICc weights (among valid only)
aic_weights <- results_grid %>%
  filter(valid) %>%
  group_by(r) %>%
  mutate(
    dAICc = AICc_plusRho - min(AICc_plusRho, na.rm = TRUE),
    wAICc = exp(-0.5*dAICc)/sum(exp(-0.5*dAICc), na.rm = TRUE)
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
    share_capital_valueadded = delta_KL,
    share_valueadded_output = delta_VAE,
    elasticity_substitution_KL = sigma_KL,
    elasticity_substitution_VAE = sigma_VAE,
    exponent_rho_KL = rho_KL,
    exponent_rho_VAE = rho_VAE,
    gamma_parameter = gamma,
    nu_parameter = nu,
    lambda_parameter = lambda,
    p_gamma, p_lambda, p_delta_KL, p_delta_VAE, p_nu
  )



# ---- EXPORT derived artifacts ----
write_csv(convergence_summary, "CES_convergence_summary.csv")
write_csv(best_methods, "CES_best_methods.csv")
write_csv(results_grid_valid, "CES_results_grid_valid.csv")
write_csv(results_grid_invalid, "CES_results_grid_invalid.csv")
write_csv(robustness_summary, "CES_robustness_summary.csv")
write_csv(aic_weights, "CES_AICc_weights.csv")
write_csv(grid_conv_share, "CES_grid_convergence_share.csv")
write_csv(iam_table, "IAM_params.csv")
