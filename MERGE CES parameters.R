options(scipen = 999) # avoids scientific notation unless necessary
setTimeLimit(cpu = Inf, elapsed = Inf, transient = TRUE) # allows long runs without R aborting due to time limits (high resolution grids can take days to run)

#### PACKAGES ####
#install.packages(c("micEconCES","dplyr","readr","purrr","parallel","tibble","progressr"))

library(micEconCES)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(furrr)
library(parallel)
library(tibble)
library(stringr)

#### SETTINGS ####
options(future.rng.onMisuse = "ignore") # suppresses warnings
options(future.wait.timeout = 0)   # disables waiting timeout for parallel workers


# setwd("C:/Users/escami_g/OneDrive - Paul Scherrer Institut/05.Models/MERGE updates/Macroeconomic submodel/CES-parametrisation/stage1")
setwd("C:/Users/escami_g/Desktop/CES data/stage1")
infile <- "MERGE macro.csv"

# TRUE to run a fresh estimation or FALSE to reuse saved results (.rds and CSV files from previous estimation)
RUN_ESTIMATION <- TRUE

# Stage 1 grid (broad scope)
# Each sequence (seq) defines a range of rho values; together they comprise the full search grid for each rho
rhoGrid_KL  <- c(seq(-0.9, -0.15, by = 0.15),
                 seq(-0.1, 0.5, by = 0.1),
                 seq(0.7, 2, by = 0.3),
                 3, 5, 10, 20, 50)
rhoGrid_VAE <- c(seq(-0.8, -0.2, by = 0.2),
                 seq(-0.1, 0.5, by = 0.1),
                 seq(0.7, 3, by = 0.3),
                 4, 6, 10, 20)

# Stage 2 grid (refined)
# rhoGrid_KL  <- c(seq(-0.90, -0.60, by = 0.10),
#                  seq(-0.60, 0.20, by = 0.03),
#                  seq(0.20, 1.90, by = 0.05),
#                  seq(1.90, 5.00, by = 0.20),
#                  seq(2.00, 21.00, by = 0.50))
# 
# rhoGrid_VAE <- c(seq(-0.72, -0.40, by = 0.08),
#                  seq(-0.40, 0.80, by = 0.03),
#                  seq(0.80, 2.50, by = 0.06),
#                  seq(2.50, 6.50, by = 0.20),
#                  seq(6.50, 10.00, by = 0.50))

# Test grid. Very coarse for debugging and testing, not for estimation
# rhoGrid_KL <- c(seq(-0.5, 1, by = 0.25))
# rhoGrid_VAE <- c(seq(-0.5, 1, by = 0.25))

# Grid sorting and rounding
rhoGrid_KL <- sort(unique(round(rhoGrid_KL, 4)))
rhoGrid_VAE <- sort(unique(round(rhoGrid_VAE, 4)))

#### FUNCTIONS ####

# Counts 1 rho_KL and 1 rho_VAE if there is more than one in each search grid. 
rho_penalty <- as.integer(length(rhoGrid_KL) > 1) + as.integer(length(rhoGrid_VAE) > 1)

# Function to return a fallback value if the left side is NULL or NA. Used to extract values from missing/empty objects
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || isTRUE(all(is.na(x)))) {
    y
  } else {
    # collapse to scalar if needed taking the first element
    if (length(x) > 1) x[[1]] else x
  }
}

# Function to extract boolean TRUE/FALSE from NULL, empty, or non-logical data
safe_bool <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) {
    default
  } else {
    out <- suppressWarnings(as.logical(x[[1]]))
    if (is.na(out)) default else out
  }
}

# Function to extract text from NULL or empty data
safe_chr <- function(x, default = NA_character_) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) {
    default
  } else {
    as.character(x[[1]])
  }
}

# Function to extract numbers from NULL or empty data
safe_num <- function(x, default = NA_real_) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) {
    default
  } else {
    suppressWarnings(as.numeric(x[[1]]))
  }
}

# Function to detect rho values at the edge of their grids. Used to reduce the priority of edge rho values in best estimates.
on_grid_edge <- function(val, grid, tol = 1e-12) {
  if (length(grid) == 0L) {
    return(rep(FALSE, length(val)))
  }
  g <- sort(unique(grid)) # removes rho duplicates and sorts in ascending order
  gmin <- g[1] # smallest rho value
  gmax <- g[length(g)] # largest rho value
  # check if rho is within tolerance of either edge. 
  # Near is from dplyr and checks if two numbers are equal within a tolerance
  near(val, gmin, tol = tol) | near(val, gmax, tol = tol)
}

# Function to standardise coefficients across formats from each solver, so the extract_grid_coeffs works regardless of the solver
coef_table_safe <- function(fit_obj) {
  # Try several places where coefficients may be stored by different solvers
  try_list <- list(
    function(x) coef(summary(x)),
    function(x) summary(x)$coefficients,
    function(x) summary(x)$coef,
    function(x) summary(x)$coefTable
  )
  for (f in try_list) {
    cm <- try(f(fit_obj), silent = TRUE)
    if (!inherits(cm, "try-error") && !is.null(cm)) {
      cm <- try(as.matrix(cm), silent = TRUE) # add to matrix if not empty
      if (!inherits(cm, "try-error") && is.matrix(cm) && nrow(cm) > 0)
        return(cm) # returns contents in matrix
    }
  }
  return(NULL) # if nothing works
}


# Function to extract the number of iterations from the solvers
# Each solver saves iteration counts in different ways, this function looks for them
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

# Function to extract statistics from each grid combination
# Extracts standard errors, t-statistics, p-values, confidence intervals fromthe cesEST objects created by the package
extract_grid_coeffs <- function(fit_sub) {
  # fit_sub is a single cesEst object from allRhoFull
  if (!inherits(fit_sub, "cesEst")) return(NULL)
  
  cm <- coef_table_safe(fit_sub) # harmonised coefficient table
  if (is.null(cm) || !is.matrix(cm) || nrow(cm) == 0L) return(NULL)
  
  rn <- rownames(cm)
  cn <- colnames(cm)
  
  # flexible columns, for different naming conventions for the statistics across solvers
  pick_col <- function(patterns) {
    ix <- which(vapply( # which returns indices of columns that match. vapply with logical makes the result a logical boolean value
      cn, # column names from the coefficient table
      function(z) any(grepl(patterns, z, ignore.case = TRUE)), # for each column name z in cn, check if it matches a regular expression
      logical(1)
    ))
    if (length(ix) == 0) NA_integer_ else ix[1] # return first matching index or NA if empty
  }
  col_est <- pick_col("^(estimate|coef|value)$|^estimate$|^coef$|^coeff")
  col_se <- pick_col("(std\\.? ?error|se)")
  col_t <- pick_col("^(t.?value|z|t)$")
  col_p <- pick_col("^(pr\\(|p.?value|p$)")
  
  # flexible rows for various names of the parameters across solvers
  pick_row <- function(patterns, exclude = NULL) {
    ok <- which(vapply(
      rn, # row names (parameters)
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
  r_gamma <- pick_row("^gamma$")
  r_lambda <- pick_row("^lambda|^lam$")
  r_delta1 <- pick_row("^delta[_ ]?1$|^delta-?1$|delta[_]?kl")
  r_deltaMain <- pick_row("^delta$", exclude = "1")
  r_nu <- pick_row("^nu$")
  
  get_val <- function(ri, ci) {
    if (is.na(ri) || is.na(ci)) return(NA_real_)
    v <- suppressWarnings(as.numeric(cm[ri, ci, drop = TRUE]))
    if (is.finite(v)) v else NA_real_
  }
  
  # Builds the estimates of each parameters from the list of statistics
  est <- list(
    gamma = get_val(r_gamma, col_est),
    lambda = get_val(r_lambda, col_est),
    delta_KVA = get_val(r_delta1, col_est),
    delta_VAY = get_val(r_deltaMain, col_est),
    nu = get_val(r_nu, col_est)
  )
  
  # Builds the standard errors for each parameter from the cleaned standard error column
  se  <- list(
    gamma = get_val(r_gamma, col_se),
    lambda = get_val(r_lambda, col_se),
    delta_KVA = get_val(r_delta1, col_se),
    delta_VAY = get_val(r_deltaMain, col_se),
    nu = get_val(r_nu, col_se)
  )
  
  # Builds the t-statistic for each parameter from the cleaned t-statistic column
  tval <- list(
    gamma = get_val(r_gamma, col_t),
    lambda = get_val(r_lambda, col_t),
    delta_KVA = get_val(r_delta1, col_t),
    delta_VAY = get_val(r_deltaMain, col_t),
    nu = get_val(r_nu, col_t)
  )
  
  # Builds the p-value for each parameter from the cleaned p-value column
  pval <- list(
    gamma = get_val(r_gamma, col_p),
    lambda = get_val(r_lambda, col_p),
    delta_KVA = get_val(r_delta1, col_p),
    delta_VAY = get_val(r_deltaMain, col_p),
    nu = get_val(r_nu, col_p)
  )
  
  # fallback standard errors if missing from summary tables from the variance-covariance matrix
  if (any(!is.finite(unlist(se)))) { # when extracted SE are missing or non-finite
    vc <- try(vcov(fit_sub), silent = TRUE) # we call the variance-covariance matrix from the solvers
    if (!inherits(vc, "try-error") && is.matrix(vc)) {
      se_v <- try(sqrt(diag(vc)), silent = TRUE) # square root of diagonals to get standard errors of parameters
      if (!inherits(se_v, "try-error")) {
        nms <- names(se_v); if (is.null(nms)) nms <- rownames(vc)
        if (length(nms)) {
          get_se <- function(pattern) {
            ix <- which(grepl(pattern, nms, ignore.case = TRUE))
            if (length(ix) == 0) NA_real_ else {
              v <- unname(se_v[ix[1]])
              if (is.finite(v) && v > 0) v else NA_real_
            }
          } # various naming conventions would be gathered
          if (!is.finite(se$gamma)) se$gamma <- get_se("^gamma$")
          if (!is.finite(se$lambda)) se$lambda <- get_se("^lambda|^lam$")
          if (!is.finite(se$delta_KVA)) se$delta_KVA <- get_se("^delta[_ ]?1$|^delta-?1$|delta[_]?kl")
          if (!is.finite(se$delta_VAY)) se$delta_VAY <- get_se("^delta$")
          if (!is.finite(se$nu)) se$nu <- get_se("^nu$")
        }
      }
    }
  }
  
  # Build t-statistic if missing: we have 
  for (nm in names(est)) {
    if (!is.finite(tval[[nm]]) && # if solver or summary doesn't give t-satistic
        is.finite(est[[nm]]) &&
        is.finite(se[[nm]]) && se[[nm]] > 0) { # but there are standard errors and estimates
      tval[[nm]] <- est[[nm]]/se[[nm]] # we calculate t-statistic
    }
  }
  # fill in p-values if missing. Just an approximation
  for (nm in names(tval)) {
    if (!is.finite(pval[[nm]]) && is.finite(tval[[nm]])) { # p-values are missing
      pval[[nm]] <- 2*(1 - pnorm(abs(tval[[nm]]))) # calculate with normal distribution pnorm multiplied by 2 to get a two-sided p-value
    }
  }
  
  # Extract the rho values to be merged back to the grid
  cf <- try(stats::coef(fit_sub), silent = TRUE)
  rho_KL  <- NA_real_
  rho_VAE <- NA_real_
  if (!inherits(cf, "try-error") && length(cf) > 0L) {
    nms <- names(cf)
    if (any(grepl("^rho[_]?1$", nms))) { # if any naming convention about rho1 exists
      rho_KL <- as.numeric(cf[grep("^rho[_]?1$", nms)[1]]) # map it to rho_KL
    }
    if ("rho" %in% nms) { # if rho or rho2 is in the names, it maps it to rho_VAE
      rho_VAE <- as.numeric(cf[["rho"]])
    } else if (any(grepl("^rho2$", nms))) {
      rho_VAE <- as.numeric(cf[grep("^rho2$", nms)[1]])
    }
  }
  
  # Build a table with all the extracted statistics and parameters. Also calculates the confidence intervals
  tibble(
    rho_KL = rho_KL,
    rho_VAE = rho_VAE,
    gamma = est$gamma,
    lambda = est$lambda,
    delta_KVA = est$delta_KVA,
    delta_VAY = est$delta_VAY,
    nu = est$nu,
    se_gamma = se$gamma,
    se_lambda = se$lambda,
    se_delta_KVA = se$delta_KVA,
    se_delta_VAY = se$delta_VAY,
    se_nu = se$nu,
    t_gamma = tval$gamma,
    t_lambda = tval$lambda,
    t_delta_KVA  = tval$delta_KVA,
    t_delta_VAY  = tval$delta_VAY,
    t_nu = tval$nu,
    p_gamma = pval$gamma,
    p_lambda = pval$lambda,
    p_delta_KVA = pval$delta_KVA,
    p_delta_VAY = pval$delta_VAY,
    p_nu = pval$nu,
    # Newly calculated confidence intervals once the data has been placed in the table
    ci_lo_gamma = ifelse(is.finite(est$gamma) & is.finite(se$gamma), est$gamma - 1.96*se$gamma, NA_real_),
    ci_hi_gamma = ifelse(is.finite(est$gamma) & is.finite(se$gamma), est$gamma + 1.96*se$gamma, NA_real_),
    ci_lo_lambda = ifelse(is.finite(est$lambda) & is.finite(se$lambda), est$lambda - 1.96*se$lambda, NA_real_),
    ci_hi_lambda = ifelse(is.finite(est$lambda) & is.finite(se$lambda), est$lambda + 1.96*se$lambda, NA_real_),
    ci_lo_delta_KVA = ifelse(is.finite(est$delta_KVA) & is.finite(se$delta_KVA), est$delta_KVA - 1.96*se$delta_KVA, NA_real_),
    ci_hi_delta_KVA = ifelse(is.finite(est$delta_KVA) & is.finite(se$delta_KVA), est$delta_KVA + 1.96*se$delta_KVA, NA_real_),
    ci_lo_delta_VAY = ifelse(is.finite(est$delta_VAY) & is.finite(se$delta_VAY), est$delta_VAY - 1.96*se$delta_VAY, NA_real_),
    ci_hi_delta_VAY = ifelse(is.finite(est$delta_VAY) & is.finite(se$delta_VAY), est$delta_VAY + 1.96*se$delta_VAY, NA_real_),
    ci_lo_nu = ifelse(is.finite(est$nu) & is.finite(se$nu), est$nu - 1.96*se$nu, NA_real_),
    ci_hi_nu = ifelse(is.finite(est$nu) & is.finite(se$nu), est$nu + 1.96*se$nu, NA_real_)
  )
}


#### ESTIMATION ####
if (isTRUE(RUN_ESTIMATION)) {
  # Run in multiple sessions for faster solve time but more computing power
  suppressWarnings({
    tryCatch({
      future::plan(future::multisession, workers = max(1, parallel::detectCores() - 1)) # Runs with all cores - 1 to avoid freezing the computer
    }, error = function(e) {
      warning("multisession failed (", conditionMessage(e), "); falling back to sequential.")
      future::plan(future::sequential) # fallback if the multisession fails, it runs sequentially
    })
  })
  
  # Data loading and normalisation
  df <- read_csv(infile, show_col_types = TRUE) # read the dataset
  
  dfS <- df %>%
    group_by(r) %>% # normalising the data to 2022. X'(t,r) = X(t,r)/X(2022,r)
    mutate(
      Ybase = Y[t == 2022][1], 
      Kbase = K[t == 2022][1],
      Lbase = L[t == 2022][1],
      Ebase = E[t == 2022][1]
    ) %>%
    mutate(
      Ys = Y/Ybase, # Economic output
      Ks = K/Kbase, # Capital
      Ls = L/Lbase, # Labour
      Es = E/Ebase  # Energy
    ) %>%
    ungroup()

  # Region estimation loop. Estimates the parameters for each region for each grid point and with each solver
  # The code in this estimate_region runs for each region
  estimate_region <- function(d, region_name) {
    message("\nEstimating region: ", region_name)
    d_num <- d %>% transmute(t, Ys, Ks, Ls, Es) # grabbing numeric values from the data
    # Setting a deterministic seed for replicability. Different seeds per region
    seed_val <- sum(utf8ToInt(region_name)) # converts region text to numbers and sums them to create an ID
    seed_val <- abs(seed_val) %% .Machine$integer.max # limiting seed to a valid range
    set.seed(seed_val) # sets the seed to the region numeric ID
    
    # Fast methods: NM, Nelder-Mead, Newton, L-BFGS-B
    # Slow methods: LM, PORT, BFGS
    # Terribly slow: SANN, DE, CG
    
    # Choose the optimisation methods to use. Two options here.
    # Sage 1
    opt_methods <- unique(c("LM", "NM", "Nelder-Mead", "BFGS", "PORT", "Newton", "CG", "L-BFGS-B", "SANN", "DE"))
    # Stage 2
    # opt_methods <- unique(c("Newton", "L-BFGS-B"))
    
    # Setting up objects to store the estimations by method
    fit_all <- setNames(vector("list", length(opt_methods)), opt_methods)
    conv_all <- setNames(rep(FALSE, length(opt_methods)), opt_methods)
    msg_all <- setNames(rep(NA_character_, length(opt_methods)), opt_methods)
    times_all <- setNames(rep(NA_real_, length(opt_methods)), opt_methods)
    
    for (m in opt_methods) {
      t0 <- Sys.time() # time start for each method
      
      # Emptying any arguments we will give to the solvers for starting points, min/max and other control settings
      start_arg <- NULL; lower_arg <- NULL; upper_arg <- NULL; control_arg <- NULL
      
      # Tuned (via trial and error) starting values, solver bounds, and control settings
      if (m == "Newton") {
        start_arg <- c(
          gamma = runif(1, 0.9, 1.1),
          lambda = runif(1, -0.005, 0.005),
          delta_KVA = runif(1, 0.4, 0.6),
          delta_VAY = runif(1, 0.4, 0.6),
          nu = runif(1, 0.9, 1.1)
        )
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
      
      # Arguments passed to the MicEconCES::cesEst() package
      args_list <- list(
        yName = "Ys",
        xNames = c("Ks","Ls","Es"),
        tName = "t",
        data = d_num,
        vrs = TRUE, # nu is estimated (variable returns to scale)
        multErr = TRUE, # multiplicative error term in CES
        method = m,
        rho1 = rhoGrid_KL, # assigning the rhoGrid_KL to rho1 (package name)
        rho = rhoGrid_VAE, # assigning the rhoGrid_VAE to rho (package name)
        returnGridAll = TRUE # Returns estimates for each point of the grid, rather than only summary
        )
      # These ifs pass the method-specific starting, min/max and control arguments. If empty, it uses the package standards
      if (!is.null(start_arg)) args_list$start <- start_arg
      if (!is.null(lower_arg)) args_list$lower <- lower_arg
      if (!is.null(upper_arg)) args_list$upper <- upper_arg
      if (!is.null(control_arg)) args_list$control <- control_arg
      
      # Fitting the estimation and silencing warnings so failed ones don't stop the process
      fit_try <- try(
        suppressWarnings(
          do.call(cesEst, args_list)
          ),
        silent = TRUE
      )
      
      runtime <- as.numeric(difftime(Sys.time(), t0, units="secs")) # defining runtime as seconds since start
      times_all[[m]] <- runtime # building a runtime per method object
      
      success <- inherits(fit_try,"cesEst") # if the estimation returns a valid object
      conv_flag <- if (success && !is.null(fit_try$convergence)) fit_try$convergence else FALSE # flag for convergence
      msg_flag <- if (inherits(fit_try,"try-error")) as.character(fit_try)[1] # captures the solver's error message
      else if (success && !is.null(fit_try$message)) fit_try$message else ""
      
      fit_all[[m]] <- if (success) fit_try else NULL # gets the successful estimations per method
      conv_all[[m]] <- conv_flag # gets all convergence flags per method
      msg_all[[m]] <- msg_flag # gets all failure messages per method
      
      # Progress status message per region-method
      base::message(sprintf("  %s %-10s in %.1fs%s",
                            if (conv_flag) "✓" else "✗", m, runtime,
                            if (!conv_flag && nzchar(msg_flag)) paste0(" msg: ", msg_flag) else ""))
    }
    
    # Returns all estimations, convergence, messages and runtimes for the region
    list(fits = fit_all, conv = conv_all, msg = msg_all, times = times_all, data = d_num)
}

  #### EXTRACT RESULTS ####
  # This function builds a parameter grid from the rhos tested, attaches coefficient estimates, flags run validity and produces a table per region/method/grid point
  extract_region <- function(region_name, region_fits) {
    # Create an empty tibble to receive the data
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
    
    # Placing the normalised data into d_num
    d_num <- region_fits$data
    # Creating an empty table to store all the results later
    grid_tbl <- tibble()
    
    for (m in names(region_fits$fits)) {
      fit_obj <- region_fits$fits[[m]]
      conv_flag <- safe_bool(region_fits$conv[[m]], FALSE)
      msg_flag <- safe_chr(region_fits$msg[[m]], NA_character_)
      runtime <- safe_num(region_fits$times[[m]], NA_real_)
      
      # Creating a full_grid table for each region × method
      # testing the full grid of all combinations of rho_KL and rho_VAE
      full_grid <- expand_grid(rho_KL = rhoGrid_KL, rho_VAE = rhoGrid_VAE) %>%
        mutate(
          r = region_name,
          method = m,
          n_grid = n(), # grid points per region-method
          runtime_total = runtime, # runtime for the region-method run
          runtime_per_grid = runtime / n_grid, # averge runtime per grid point
          msg = msg_flag, # solver message per method
          conv = conv_flag,  # will be overridden by grid-level convergence if available
          # Calculating constant elasticity of substitution
          sigma_KL = ifelse(is.finite(1/(1+rho_KL)), 1/(1+rho_KL), NA_real_),
          sigma_VAE = ifelse(is.finite(1/(1+rho_VAE)), 1/(1+rho_VAE), NA_real_),
          # Flags for rho at the boundary of the grids
          on_edge_KL = on_grid_edge(rho_KL,  rhoGrid_KL),
          on_edge_VAE = on_grid_edge(rho_VAE, rhoGrid_VAE),
          rss = NA_real_
        )
      
      # pre-fill parameter/statistics as NA
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
        # RSS and convergence per grid from allRhoSum
        # allRhoSum is the grid point summary table produced by the package
        if (!is.null(fit_obj$allRhoSum) && nrow(fit_obj$allRhoSum) > 0L) {
          sum_tbl <- fit_obj$allRhoSum
          # assume columns rho1 & rho & convergence exist for nested case
          sum_tbl <- sum_tbl %>%
            mutate(
              rho1 = as.numeric(.data[["rho1"]]),
              rho = as.numeric(.data[["rho"]]),
              rss = as.numeric(.data[["rss"]]),
              conv_grid = safe_bool(.data[["convergence"]], TRUE)
            ) %>%
            select(rho1, rho, rss, conv_grid) %>%
            rename(rho_KL = rho1, rho_VAE = rho)
          
          full_grid <- full_grid %>%
            left_join(sum_tbl, by = c("rho_KL", "rho_VAE"), suffix = c("", ".sum")) %>%
            mutate(
              rss = coalesce(rss.sum, rss), # override RSS with the grid point value
              conv = ifelse(is.na(conv_grid), conv, conv_grid) # use convergence at grid point if available
            ) %>%
            select(-rss.sum, -conv_grid)
        }
        
        # Per-grid coefficients from allRhoFull
        # allRhoFull has all the solver objects per grid point, so we extract parameter values and statistics
        coef_df <- NULL
        if (!is.null(fit_obj$allRhoFull) && length(fit_obj$allRhoFull) > 0L) {
          coef_df <- purrr::map_dfr(fit_obj$allRhoFull, extract_grid_coeffs)
          
          full_grid <- full_grid %>%
            left_join(coef_df, by = c("rho_KL","rho_VAE"), suffix = c("", ".sub")) %>%
            mutate(
              gamma = coalesce(gamma.sub, gamma),
              lambda = coalesce(lambda.sub, lambda),
              delta_KVA = coalesce(delta_KVA.sub, delta_KVA),
              delta_VAY = coalesce(delta_VAY.sub, delta_VAY),
              nu = coalesce(nu.sub, nu),
              se_gamma = coalesce(se_gamma.sub, se_gamma),
              se_lambda = coalesce(se_lambda.sub, se_lambda),
              se_delta_KVA = coalesce(se_delta_KVA.sub, se_delta_KVA),
              se_delta_VAY = coalesce(se_delta_VAY.sub, se_delta_VAY),
              se_nu = coalesce(se_nu.sub, se_nu),
              t_gamma = coalesce(t_gamma.sub, t_gamma),
              t_lambda = coalesce(t_lambda.sub, t_lambda),
              t_delta_KVA = coalesce(t_delta_KVA.sub, t_delta_KVA),
              t_delta_VAY = coalesce(t_delta_VAY.sub, t_delta_VAY),
              t_nu = coalesce(t_nu.sub, t_nu),
              p_gamma = coalesce(p_gamma.sub, p_gamma),
              p_lambda = coalesce(p_lambda.sub, p_lambda),
              p_delta_KVA = coalesce(p_delta_KVA.sub, p_delta_KVA),
              p_delta_VAY = coalesce(p_delta_VAY.sub, p_delta_VAY),
              p_nu = coalesce(p_nu.sub, p_nu),
              ci_lo_gamma = coalesce(ci_lo_gamma.sub, ci_lo_gamma),
              ci_hi_gamma = coalesce(ci_hi_gamma.sub, ci_hi_gamma),
              ci_lo_lambda = coalesce(ci_lo_lambda.sub, ci_lo_lambda),
              ci_hi_lambda = coalesce(ci_hi_lambda.sub, ci_hi_lambda),
              ci_lo_delta_KVA = coalesce(ci_lo_delta_KVA.sub, ci_lo_delta_KVA),
              ci_hi_delta_KVA = coalesce(ci_hi_delta_KVA.sub, ci_hi_delta_KVA),
              ci_lo_delta_VAY = coalesce(ci_lo_delta_VAY.sub, ci_lo_delta_VAY),
              ci_hi_delta_VAY = coalesce(ci_hi_delta_VAY.sub, ci_hi_delta_VAY),
              ci_lo_nu = coalesce(ci_lo_nu.sub, ci_lo_nu),
              ci_hi_nu = coalesce(ci_hi_nu.sub, ci_hi_nu)
            ) %>%
            select(-ends_with(".sub"))
        }
        
        # R² and AIC per grid cell (using rss)
        # Calculating goodness of fit and information criteria in log for Ys
        obs_log <- try(log(d_num$Ys + 1e-12), silent = TRUE) # adding 1e-12 to avoid log(0)
        if (!inherits(obs_log, "try-error")) {
          n_obs <- length(obs_log)
          TSS_log <- sum((obs_log - mean(obs_log))^2) # total sum of squares in log
          k_hat <- length(tryCatch(stats::coef(fit_obj), error = function(e) numeric(0))) # number of estimated parameters (should be 5)
          k0 <- k_hat + rho_penalty # parameters + rho penalty (should be 7 in most cases)
          
          full_grid <- full_grid %>%
            mutate(
              # R² calculated in log
              R2 = ifelse(is.finite(rss) & TSS_log > 0,
                          1 - rss/TSS_log,
                          NA_real_),
              # Adjusted R² penalising models with more estimated parameters k_hat
              adjR2 = ifelse(
                is.finite(R2) & (n_obs - k_hat - 1) > 0,
                1 - (1 - R2)*(n_obs - 1)/(n_obs - k_hat - 1),
                NA_real_
              ),
              # Akaike Information Criterion with number of estimated parameter k_hat
              AIC_naive = ifelse(
                is.finite(rss),
                n_obs * log(rss/n_obs) + 2*k_hat,
                NA_real_
              ),
              # Small sample correction AIC
              AICc_naive = ifelse(
                is.finite(AIC_naive) & (n_obs - k_hat - 1) > 0,
                AIC_naive + (2*k_hat*(k_hat + 1))/(n_obs - k_hat - 1),
                NA_real_
              ),
              # Same as AIC plus a rho penalty inside k0 = khat + rho_penalty
              AIC_plusRho = ifelse(
                is.finite(rss),
                n_obs*log(rss/n_obs) + 2*k0,
                NA_real_
              ),
              AICc_plusRho = ifelse(
                is.finite(AIC_plusRho) & (n_obs - k0 - 1) > 0,
                AIC_plusRho + (2*k0*(k0 + 1))/(n_obs - k0 - 1),
                NA_real_
              )
            )
        }
        
        # Iterations (method-level)
        full_grid$iter <- iter_safe(fit_obj) # adding iteration count of the method to all grid points
      }
      
      grid_tbl <- bind_rows(grid_tbl, full_grid)
    }
    
    # Filling in the grid table with the previously calculated or extracted parameters/statistics
    grid_tbl %>%
      mutate(
        across(c(
          delta_KVA, 
          delta_VAY, 
          gamma, 
          nu, 
          lambda, 
          sigma_KL, 
          sigma_VAE,
          se_gamma, 
          se_lambda, 
          se_delta_KVA, 
          se_delta_VAY, 
          se_nu,
          t_gamma, 
          t_lambda, 
          t_delta_KVA, 
          t_delta_VAY, 
          t_nu,
          p_gamma, 
          p_lambda, 
          p_delta_KVA, 
          p_delta_VAY, 
          p_nu,
          ci_lo_gamma, 
          ci_hi_gamma, 
          ci_lo_lambda, 
          ci_hi_lambda,
          ci_lo_delta_KVA, 
          ci_hi_delta_KVA, 
          ci_lo_delta_VAY, 
          ci_hi_delta_VAY,
          ci_lo_nu, 
          ci_hi_nu,
          rss, 
          R2, 
          adjR2, 
          AIC_naive, 
          AICc_naive, 
          AIC_plusRho, 
          AICc_plusRho,
          iter, 
          rho_KL, 
          rho_VAE, 
          runtime_total, 
          runtime_per_grid, 
          n_grid
        ), ~ suppressWarnings(as.numeric(.)))
      )
  }
  
  
# Region splits of the normalised data (dFS)
splits <- dfS %>% group_split(r, .keep=TRUE) 
region_names <- dfS %>% distinct(r) %>% pull(r) # extract the region names from dFS



#### RUN ALL ####
# Running estimations in parallel across regions
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
    packages = c("micEconCES","dplyr","tidyr","purrr","readr"), # making sure these are loaded for the parallel processes
    globals = c("rhoGrid_KL","rhoGrid_VAE","estimate_region","extract_region"), # exporting objects to parallel processes
    seed = TRUE
  )
)

# Extract grid results from all region-method combinations into a single table
results_grid <- map2_dfr(region_names, fits_all, extract_region)

#### SAVE MASTER TABLE ####
write_csv(results_grid, "CES_results_grid.csv") # main table with all grid results and statistics
saveRDS(results_grid, "results_run1.rds") # R object with all results
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
  
  # Load input data so df is available for IAM table
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

#### BEST RESULTS AND OUTPUTS ####
# Convergence summary by region
convergence_summary <- results_grid %>%
  select(r, method, conv) %>%
  distinct() %>%
  count(method, conv, name = "count") %>%
  mutate(status = ifelse(conv, "Converged", "Failed"))

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

  
# REbuild with validity columns for the whole grid
results_grid <- results_grid %>% add_validity()

# Valid runs table. Making a subset with only valid estimations and adding model difference in AICc
results_grid_valid <- results_grid %>% 
  filter(valid) %>%
  group_by(r) %>%
  mutate(
    valid_estimations = n(),
    dAICc = AICc_plusRho - min(AICc_plusRho, na.rm = TRUE),
    wAICc_raw = exp(-0.5 * dAICc)
  ) %>%
  ungroup()


# Model-average best methods
dAICcmin <- 4

best_methods_average <- results_grid_valid %>%
  group_by(r) %>%
  mutate(
    any_valid_strict = any(valid_strict, na.rm = TRUE),
    any_strict_support = any(valid_strict & dAICc <= dAICcmin, na.rm = TRUE),
    any_support = any(dAICc <= dAICcmin, na.rm = TRUE),
    support_mode = case_when(
      any_support & any_strict_support ~ "strong AICc support", # dAICc <= 4 AND at least one strict
      any_support & !any_strict_support ~ "some AICc support", # dAICc <= 4 but no strict inside support
      !any_support ~ "AICc min support", # no dAICc <= 4 → fallback to min
      TRUE ~ "AICc min support" # safety catch
    )
  ) %>%
  # 1) Keep only rows belonging to the chosen support set for that region
  filter(
    (support_mode %in% c("strong AICc support", "some AICc support") & dAICc <= dAICcmin) |
      (support_mode == "AICc min support" & dAICc == min(dAICc, na.rm = TRUE))
  ) %>%
  # 2) If support_mode is strong AICc support, keep only strict points
  filter(
    support_mode != "strong AICc support" | valid_strict
  ) %>%
  group_by(r, support_mode) %>%
  # 3) Renormalise AICc weights within the final averaging set
  mutate(
    wAICc = wAICc_raw/sum(wAICc_raw, na.rm = TRUE),
    avg_estimations = n(), # grid points actually used in averaging
    valid_estimations = first(valid_estimations),  # total valid grid points in region
    share_valid = avg_estimations/valid_estimations
  ) %>%
  summarise(
    # model-averaged parameters (AICc weights)
    gamma = sum(wAICc*gamma, na.rm = TRUE),
    delta_KVA = sum(wAICc*delta_KVA, na.rm = TRUE),
    delta_VAY = sum(wAICc*delta_VAY, na.rm = TRUE),
    nu = sum(wAICc*nu, na.rm = TRUE),
    lambda = sum(wAICc*lambda, na.rm = TRUE),
    rho_KL = sum(wAICc*rho_KL, na.rm = TRUE),
    rho_VAE = sum(wAICc*rho_VAE, na.rm = TRUE),
    sigma_KL = 1/(1 + rho_KL),
    sigma_VAE = 1/(1 + rho_VAE),
    avg_estimations = first(avg_estimations),
    valid_estimations = first(valid_estimations),
    share_valid = first(share_valid),
    .groups = "drop"
  )



best_methods <- results_grid %>%
  mutate(
    tier = case_when(
      valid_strict ~ "strict-valid",
      valid ~ "valid",
      converged & is.finite(R2) & R2 > 0 ~ "conv-only",
      TRUE ~ "any"
    ),
    tier_rank = factor(tier, levels = c("strict-valid","valid","conv-only","any")),
    aicc_rank = if_else(is.finite(AICc_plusRho), AICc_plusRho, Inf),
    rss_rank = if_else(is.finite(rss), rss, Inf),
    R2_rank = if_else(is.finite(R2), R2, -Inf)
  ) %>%
  group_by(r) %>%
  arrange(
    tier_rank,
    aicc_rank,
    rss_rank,
    desc(R2_rank)
    ) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(best_tier = as.character(tier_rank))


best_methods_AICc <- results_grid_valid %>%
  group_by(r) %>%
  filter(dAICc == min(dAICc, na.rm = TRUE)) %>%
  slice_head(n = 1) %>%
  ungroup()





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
    min_RSS = min(rss, na.rm = TRUE),
    max_RSS = max(rss, na.rm = TRUE),
    med_RSS = median(rss, na.rm = TRUE),
    min_R2 = min(R2, na.rm = TRUE),
    max_R2 = max(R2, na.rm = TRUE),
    med_R2 = median(R2, na.rm = TRUE),
    # Parameters (medians across runs)
    gamma_med = median(gamma, na.rm = TRUE),
    lambda_med = median(lambda, na.rm = TRUE),
    delta_KVA_med = median(delta_KVA, na.rm = TRUE),
    delta_VAY_med = median(delta_VAY, na.rm = TRUE),
    nu_med = median(nu, na.rm = TRUE),
    sigma_KL_med = median(sigma_KL, na.rm = TRUE),
    sigma_VAE_med = median(sigma_VAE, na.rm = TRUE),
    # Iterations and runtime
    med_iter = median(iter, na.rm = TRUE),
    med_runtime = median(runtime_total, na.rm = TRUE),
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




# ---- EXPORT derived datasets ----
#write_csv(convergence_summary, "CES_convergence_summary.csv")
write_csv(best_methods, "CES_best_methods.csv")
write_csv(best_methods_average, "CES_best_methods_average.csv")
write_csv(results_grid, "CES_results_grid.csv")
#write_csv(results_grid_valid, "CES_results_grid_valid.csv")
write_csv(robustness_summary, "CES_robustness_summary.csv")
#write_csv(aic_weights, "CES_AICc_weights.csv")
write_csv(grid_conv_share, "CES_grid_convergence_share.csv")
write_csv(iam_table, "IAM_params.csv")
