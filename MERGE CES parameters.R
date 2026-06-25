options(scipen = 999) # avoids scientific notation unless necessary
setTimeLimit(cpu = Inf, elapsed = Inf, transient = TRUE) # allows long runs without R aborting due to time limits (high resolution grids can take days to run)

#### PURPOSE ####
# Estimate the nested CES production-function parameters used by MERGE-ETL.

# The economic structure is:
#   Y = gamma*exp(lambda*time_from_base_year)*[delta_VAY*VA^rho_VAE + (1 - delta_VAY)*E^rho_VAE]^(nu/rho_VAE)

#   VA = [delta_KVA*K^rho_KL + (1 - delta_KVA)*L^rho_KL]^(1/rho_KL)

# micEconCES uses rho as the curvature parameter and sigma = 1/(1 + rho)
# as the substitution elasticity. The script builds sigma grids because sigma
# is easier to interpret economically, then converts them to rho for estimation.

# Main outputs:
# - stage*_results.rds keeps the full estimation record and diagnostics.
# - stage2_merge_iam_parameter_table.csv is the compact regional table for MERGE.
# - gamma is the base-year total factor productivity scale.
# - lambda is the continuous annual TFP growth rate relative to the base year.
# - delta_KVA and delta_VAY are nest-specific CES distribution weights.
# - rho_KL and rho_VAE are CES curvature terms; sigma_KL and sigma_VAE are the
#   corresponding substitution elasticities.
# - nu is the returns-to-scale exponent.

#### PACKAGES ####
#install.packages(c("micEconCES","dplyr","tidyr","readr","purrr","furrr","parallel","tibble","stringr"))

pkgs <- c("micEconCES", "dplyr", "tidyr", "readr", "purrr", "furrr", "parallel", "tibble", "stringr")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) {
  stop("Install required packages before running this script: ", paste(missing_pkgs, collapse = ", "))
}
invisible(lapply(pkgs, function(pkg) {
  suppressPackageStartupMessages(suppressWarnings(library(pkg, character.only = TRUE)))
}))

#### SETTINGS ####
options(future.rng.onMisuse = "ignore") # suppresses warnings
options(future.wait.timeout = 0)   # disables waiting timeout for parallel workers

# Quick map:
# - Run controls: choose Stage 1, Stage 2, or a small test run.
# - Data paths: point the workflow to the macro input and output folder.
# - Search design: set Stage 1/Stage 2 sigma grids and solver counts.
# - Economic validity: edit plausibility and final-selection thresholds here,
#   not inside the validity code below.
# - Solver bounds: keep optimizer bounds aligned with the economic validity
#   thresholds unless there is a documented reason to loosen them.

# RStudio use: edit this block, then click Source on this file. Do not paste a
# source("MERGE CES parameters.R") line inside this file; source() is only for
# running the script from the console or from another driver script.
# Stage 1 searches broadly and learns where each region has statistical support.
# Stage 2 uses those Stage 1 support windows to produce the final MERGE table.
# The separate Stage 1 sequence script is optional; it only automates repeated
# targeted Stage 1 runs when several rounds are needed.
setting <- function(name, default) {
  if (exists(name, envir = .GlobalEnv, inherits = FALSE)) get(name, envir = .GlobalEnv) else default
}

previous_stage1_file_given <- exists("previous_stage1_file", envir = .GlobalEnv, inherits = FALSE)

run_stage <- setting("run_stage", "stage1") # "stage1", "stage2", or "test"
stage1 <- setting("stage1", "stage1.1") # stage1.1 is broad; stage1.2+ are targeted audits
stage2 <- setting("stage2", "stage2.10") # latest final MERGE table; increase only when rerunning final refinement
test <- setting("test", "stage1") # in test mode, choose whether the test mimics Stage 1 or Stage 2

estimate <- setting("estimate", TRUE) # TRUE fits models; FALSE rebuilds diagnostics from a saved RDS
saved_results <- setting("saved_results", NA_character_) # optional RDS file to reanalyse when estimate is FALSE
data_folder <- setting("data_folder", "C:/Users/escami_g/Desktop/CES data") # input and output root folder
previous_stage1_file <- setting("previous_stage1_file", NA_character_) # optional override; stage1.2+ infer this automatically
stage1_file_for_stage2 <- setting("stage1_file_for_stage2", file.path(data_folder, "stage1", "stage1_results.rds")) # Stage 2 reads this
test_n <- setting("test_n", 2L) # number of regions in quick test runs
estimate_flagged_regions <- setting("estimate_flagged_regions", FALSE) # TRUE estimates regions that fail the Stage 0 input gate

stage1_expand_points <- setting("stage1_expand_points", 11L) # points per nest in targeted Stage 1 expansion grids
stage1_sigma_min <- setting("stage1_sigma_min", 0.01) # stress-test lower elasticity bound; very close to Leontief
stage1_sigma_max <- setting("stage1_sigma_max", 12) # stress-test upper elasticity bound; very high substitution
stage1_max_solvers <- setting("stage1_max_solvers", 4L) # solvers kept per region after Stage 1.1 evidence

stage2_grid_points <- setting("stage2_grid_points", 17L) # points per nest in final Stage 2 grids before extra focus points
stage2_window_padding <- setting("stage2_window_padding", 0.25) # widens the Stage 1 support window before Stage 2
stage2_max_solvers <- setting("stage2_max_solvers", 3L) # solvers kept per region in Stage 2
stage2_sigma_min <- setting("stage2_sigma_min", 0.05) # final MERGE lower elasticity bound
stage2_sigma_max <- setting("stage2_sigma_max", 5) # final MERGE upper elasticity bound
stage2_boundary_sigma_min <- setting("stage2_boundary_sigma_min", stage1_sigma_min) # boundary-supported Stage 2 can retain Stage 1 stress-test elasticity limits
stage2_boundary_sigma_max <- setting("stage2_boundary_sigma_max", stage1_sigma_max)
stage2_solver_pool <- setting("stage2_solver_pool", c("PORT", "L-BFGS-B")) # micEconCES source supports bounds for these methods
aicc_support_delta <- setting("aicc_support_delta", 4) # delta AICc support window for Stage 1 support and final Stage 2 selection

# Economic validity and final-selection thresholds.
# Edit these lists when the scientific review changes the admissible parameter
# ranges. The add_validity() function below reads these objects directly.
limits <- setting("limits", list(
  soft = list(
    gamma = c(0.5, 3),
    lambda = c(-0.05, 0.05),
    nu = c(0.7, 1.3)
  ),
  final = list(
    delta_KVA = c(0.02, 0.98),
    delta_VAY = c(0.02, 0.98),
    sigma_KL = c(stage2_sigma_min, stage2_sigma_max),
    sigma_VAE = c(stage2_sigma_min, stage2_sigma_max),
    lambda = c(-0.03, 0.03),
    nu = c(0.5, 2)
  ),
  strict = list(
    delta_KVA = c(0.2, 0.8),
    delta_VAY = c(0.2, 0.8)
  )
))

# Solver bounds and starting values. Stage 2 uses the tighter final bounds
# because it is producing the MERGE-facing table; Stage 1 stays wider so the
# broad audit can identify whether a region is forcing boundary behaviour.
start_vals <- setting("start_vals", c(
  gamma = 1,
  lambda = 0.001,
  delta_KVA = 0.5,
  delta_VAY = 0.5,
  nu = 1
))
stage2_lower <- setting("stage2_lower", c(
  gamma = 0.1,
  lambda = limits$soft$lambda[[1]],
  delta_KVA = limits$final$delta_KVA[[1]],
  delta_VAY = limits$final$delta_VAY[[1]],
  nu = limits$final$nu[[1]]
))
stage2_upper <- setting("stage2_upper", c(
  gamma = 10,
  lambda = limits$final$lambda[[2]],
  delta_KVA = limits$final$delta_KVA[[2]],
  delta_VAY = limits$final$delta_VAY[[2]],
  nu = 2.5
))
stage1_lower <- setting("stage1_lower", c(
  gamma = 0.1,
  lambda = -0.3,
  delta_KVA = 0.1,
  delta_VAY = 0.1,
  nu = 0.3
))
stage1_upper <- setting("stage1_upper", c(
  gamma = 10,
  lambda = 0.3,
  delta_KVA = 0.9,
  delta_VAY = 0.9,
  nu = 5
))
share_lower <- setting("share_lower", c(delta_KVA = 0.1, delta_VAY = 0.1))
share_upper <- setting("share_upper", c(delta_KVA = 0.9, delta_VAY = 0.9))
stage1_solvers <- setting("stage1_solvers", c("LM", "NM", "BFGS", "PORT", "Newton", "CG", "L-BFGS-B", "SANN", "DE"))
test_solvers <- setting("test_solvers", c("L-BFGS-B"))

backend <- setting("backend", "multisession") # fastest; use "sequential" only for one-region checks
workers <- setting("workers", min(6L, max(1L, parallel::detectCores() - 1L))) # stable long-run pool; raise only after a clean full pass
resume_regions <- setting("resume_regions", TRUE) # reuse completed region runs only when their settings match this run
save_run_history <- setting("save_run_history", TRUE) # keeps a dated scientific record under stage*/history/

# Optional command-line overrides for isolated checks. These are not needed for normal RStudio use.
command_line_args <- commandArgs(trailingOnly = TRUE)

get_command_line_value <- function(prefix, default = NA_character_) {
  matching_argument <- command_line_args[grepl(paste0("^", prefix, "="), command_line_args)]
  if (length(matching_argument) == 0L) {
    return(default)
  }
  sub(paste0("^", prefix, "="), "", matching_argument[[1]])
}

previous_stage1_file_from_stage <- function(stage1_label, data_folder) {
  match <- regexec("^stage1\\.([0-9]+)$", stage1_label)
  parts <- regmatches(stage1_label, match)[[1]]
  if (length(parts) < 2L) return(NA_character_)
  stage1_number <- suppressWarnings(as.integer(parts[[2]]))
  if (!is.finite(stage1_number) || stage1_number <= 1L) return(NA_character_)

  previous_stage1 <- paste0("stage1.", stage1_number - 1L)
  file.path(
    data_folder,
    "stage1",
    "history",
    previous_stage1,
    paste0(previous_stage1, "_results.rds")
  )
}

if ("--test" %in% command_line_args) {
  run_stage <- "test"
}
if ("--test-stage2" %in% command_line_args) {
  run_stage <- "test"
  test <- "stage2"
}
if ("--multisession" %in% command_line_args) {
  backend <- "multisession"
}
if ("--no-checkpoint-resume" %in% command_line_args) {
  resume_regions <- FALSE
}
previous_stage1_override <- get_command_line_value("--previous-stage1-file")
if (!is.na(previous_stage1_override) && nzchar(previous_stage1_override)) {
  previous_stage1_file <- previous_stage1_override
  previous_stage1_file_given <- TRUE
}
run_label_override <- get_command_line_value("--run-label")
if (!is.na(run_label_override) && nzchar(run_label_override)) {
  if (identical(run_stage, "stage2")) {
    stage2 <- run_label_override
  } else {
    stage1 <- run_label_override
  }
}
test_region_arg <- command_line_args[grepl("^--test-regions=", command_line_args)]
if (length(test_region_arg) > 0L) {
  test_n <- suppressWarnings(
    as.integer(sub("^--test-regions=", "", test_region_arg[[1]]))
  )
  if (!is.finite(test_n) || test_n < 1L) {
    stop("--test-regions must be a positive integer.")
  }
}
if (!run_stage %in% c("stage1", "stage2", "test")) {
  stop('run_stage must be "stage1", "stage2", or "test".')
}
if (!test %in% c("stage1", "stage2")) {
  stop('test must be "stage1" or "stage2".')
}
if (!isTRUE(previous_stage1_file_given) && identical(run_stage, "stage1")) {
  previous_stage1_file <- previous_stage1_file_from_stage(stage1, data_folder)
}

# Derived settings. These translate the human run choice above into concrete
# folders, file labels, checkpoint locations, and test-mode behaviour.
estimation_stage <- if (identical(run_stage, "test")) test else run_stage
base_year <- 2022
stage_output_folder <- if (identical(run_stage, "stage2")) "stage2" else "stage1"
output_folder <- file.path(data_folder, stage_output_folder)
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
stage_output_folder_absolute <- normalizePath(output_folder, winslash = "/", mustWork = FALSE)
input_file <- file.path(data_folder, "stage1", "MERGE macro.csv")
setwd(output_folder)
output_folder <- "."
test_mode <- identical(run_stage, "test") # non-scientific integration test only; writes test-labelled outputs
run_label <- if (isTRUE(test_mode)) {
  "test"
} else if (identical(estimation_stage, "stage1")) {
  stage1
} else {
  stage2
}
run_file_label <- gsub("[^A-Za-z0-9.]+", "_", run_label)

output_tag <- if (isTRUE(test_mode)) {
  "test"
} else {
  estimation_stage
}
output_file <- function(file_stem, extension) {
  file.path(output_folder, paste0(file_stem, "_", output_tag, extension))
}
main_results_file <- file.path(output_folder, paste0(output_tag, "_results.rds"))
input_data_summary_file <- file.path(output_folder, paste0(output_tag, "_input_data_summary.csv"))
run_history_folder <- file.path(output_folder, "history", run_file_label)
run_history_results_file <- file.path(
  run_history_folder,
  paste0(run_file_label, "_results.rds")
)
run_history_index_file <- file.path(output_folder, paste0(output_tag, "_run_history.csv"))
region_checkpoint_folder <- file.path(
  stage_output_folder_absolute,
  paste0(output_tag, "_", run_file_label, "_region_checkpoints")
)

safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) "unnamed" else x
}

region_fit_checkpoint_file <- function(region_name) {
  file.path(
    region_checkpoint_folder,
    paste0(output_tag, "_", safe_filename(region_name), "_fit.rds")
  )
}

file_fingerprint <- function(file) {
  if (is.na(file) || !nzchar(file)) {
    return(list(path = NA_character_, exists = FALSE))
  }
  normalized_file <- normalizePath(file, winslash = "/", mustWork = FALSE)
  if (!file.exists(file)) {
    return(list(path = normalized_file, exists = FALSE))
  }
  info <- file.info(file)
  list(
    path = normalized_file,
    exists = TRUE,
    size = unname(info$size),
    modified = as.character(unname(info$mtime))
  )
}

checkpoint_run_key <- function(region_name, d) {
  list(
    stage = estimation_stage,
    run_label = run_label,
    region = region_name,
    test_mode = test_mode,
    input_file = file_fingerprint(input_file),
    previous_stage1_file = file_fingerprint(previous_stage1_file),
    stage1_file_for_stage2 = file_fingerprint(stage1_file_for_stage2),
    sigma_grid_KL = sigma_grid_KL,
    sigma_grid_VAE = sigma_grid_VAE,
    rho_grid_KL = rho_grid_KL,
    rho_grid_VAE = rho_grid_VAE,
    stage1_expand_points = stage1_expand_points,
    stage1_sigma_min = stage1_sigma_min,
    stage1_sigma_max = stage1_sigma_max,
    stage1_max_solvers = stage1_max_solvers,
    stage1_solvers = stage1_solvers,
    stage2_grid_points = stage2_grid_points,
    stage2_window_padding = stage2_window_padding,
    stage2_max_solvers = stage2_max_solvers,
    stage2_sigma_min = stage2_sigma_min,
    stage2_sigma_max = stage2_sigma_max,
    stage2_solver_pool = stage2_solver_pool,
    aicc_support_delta = aicc_support_delta,
    limits = limits,
    start_vals = start_vals,
    stage1_lower = stage1_lower,
    stage1_upper = stage1_upper,
    stage2_lower = stage2_lower,
    stage2_upper = stage2_upper,
    share_lower = share_lower,
    share_upper = share_upper,
    region_data = d %>% arrange(t) %>% select(any_of(c("r", "t", "Y", "K", "L", "E")))
  )
}

checkpoint_reuse_status <- function(checkpoint, current_key) {
  if (!is.list(checkpoint) || length(checkpoint$fits %||% list()) == 0L) {
    return("incomplete")
  }
  if (!"checkpoint_key" %in% names(checkpoint)) {
    return("missing settings key")
  }
  if (!identical(checkpoint$checkpoint_key, current_key)) {
    return("settings changed")
  }
  "ok"
}

safe_write_csv <- function(x, file) {
  tryCatch(
    {
      readr::write_csv(x, file)
      TRUE
    },
    error = function(e) {
      warning(
        "Could not write CSV '", file, "'. ",
        "The RDS results file will still be written. Close the file in Excel/RStudio if you need this CSV refreshed. ",
        "Original error: ", conditionMessage(e),
        call. = FALSE
      )
      FALSE
    }
  )
}

safe_save_rds <- function(object, file) {
  tryCatch(
    {
      saveRDS(object, file)
      file
    },
    error = function(e) {
      dated_copy_file <- file.path(
        dirname(file),
        paste0(
          tools::file_path_sans_ext(basename(file)),
          "_",
          format(Sys.time(), "%Y%m%d_%H%M%S"),
          ".rds"
        )
      )
      tryCatch(
        {
          saveRDS(object, dated_copy_file)
          warning(
            "Could not write RDS '", file, "'. ",
            "Wrote dated copy instead: ", dated_copy_file, ". ",
            "Close the locked file before rerunning if you want the standard filename refreshed. ",
            "Original error: ", conditionMessage(e),
            call. = FALSE
          )
          dated_copy_file
        },
        error = function(e2) {
          warning(
            "Could not write RDS '", file, "' or dated fallback '", dated_copy_file, "'. ",
            "Original error: ", conditionMessage(e), ". ",
            "Fallback error: ", conditionMessage(e2),
            call. = FALSE
          )
          NA_character_
        }
      )
    }
  )
}

save_region_checkpoint <- function(region_fit_result, checkpoint_file) {
  tryCatch(
    {
      dir.create(dirname(checkpoint_file), recursive = TRUE, showWarnings = FALSE)
      written_file <- safe_save_rds(region_fit_result, checkpoint_file)
      if (is.na(written_file)) {
        warning(
          "Region fit finished, but its checkpoint could not be written: ",
          checkpoint_file,
          call. = FALSE
        )
      }
      invisible(written_file)
    },
    error = function(e) {
      warning(
        "Region fit finished, but checkpoint handling failed for '",
        checkpoint_file,
        "': ",
        conditionMessage(e),
        call. = FALSE
      )
      invisible(NA_character_)
    }
  )
}

collapse_for_history <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(NA_character_)
  }
  paste(x, collapse = "; ")
}

write_run_history <- function(results_object, main_file, history_file, index_file) {
  dir.create(dirname(history_file), recursive = TRUE, showWarnings = FALSE)
  written_history_file <- safe_save_rds(results_object, history_file)
  metadata <- results_object$metadata
  pick <- function(x, y) {
    if (is.null(x) || length(x) == 0L || isTRUE(all(is.na(x)))) {
      y
    } else {
      x[[1]]
    }
  }
  history_row <- tibble(
    run_label = pick(metadata$run_label, NA_character_),
    estimation_stage = pick(metadata$estimation_stage, NA_character_),
    created_at = as.character(pick(metadata$created_at, NA_character_)),
    run_mode = pick(metadata$run_mode, NA_character_),
    estimate = pick(metadata$estimate, NA),
    reused_saved_results = pick(metadata$reused_saved_results, NA),
    reused_results_source_file = pick(metadata$reused_results_source_file, NA_character_),
    input_regions = pick(metadata$input_regions, NA_integer_),
    grid_rows = nrow(results_object$grid_results),
    valid_grid_rows = nrow(results_object$valid_grid_results),
    sigma_grid_KL = collapse_for_history(metadata$obs_sigma_KL),
    sigma_grid_VAE = collapse_for_history(metadata$obs_sigma_VAE),
    rho_grid_KL = collapse_for_history(metadata$obs_rho_KL),
    rho_grid_VAE = collapse_for_history(metadata$obs_rho_VAE),
    stage1_solvers = collapse_for_history(metadata$stage1_solvers),
    stage2_solver_source = pick(metadata$stage2_solver_source, NA_character_),
    main_results_file = normalizePath(main_file, winslash = "/", mustWork = FALSE),
    run_history_results_file = normalizePath(written_history_file, winslash = "/", mustWork = FALSE)
  ) %>%
    mutate(across(everything(), as.character))
  history_index <- if (file.exists(index_file)) {
    suppressWarnings(readr::read_csv(index_file, col_types = readr::cols(.default = readr::col_character())))
  } else {
    tibble()
  }
  history_index <- bind_rows(history_index, history_row) %>%
    arrange(estimation_stage, run_label, created_at) %>%
    group_by(estimation_stage, run_label) %>%
    slice_tail(n = 1) %>%
    ungroup()
  safe_write_csv(history_index, index_file)
  written_history_file
}

# Stage 1 grid (broad scope). Sigma is the elasticity of substitution: values
# below 1 imply more complementarity, values around 1 are Cobb-Douglas-like, and
# values above 1 imply greater substitution. micEconCES documents
# sigma = 1 / (1 + rho), so the broad Stage 1 grid is specified in sigma space
# and converted to rho.
sigma_grid_KL_stage1 <- c(
  0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80,
  0.90, 1.00, 1.10, 1.25, 1.50, 2.00, 3.00, 5.00
)
sigma_grid_VAE_stage1 <- c(
  0.05, 0.075, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50,
  0.60, 0.75, 1.00, 1.25, 1.50, 2.00, 3.00
)
rho_from_sigma <- function(sigma) {
  (1 / sigma) - 1
}
rho_grid_KL_stage1 <- rho_from_sigma(sigma_grid_KL_stage1)
rho_grid_VAE_stage1 <- rho_from_sigma(sigma_grid_VAE_stage1)

# Test grid. Very coarse for quick checks, not for estimation.
sigma_grid_KL_test <- c(0.50, 1.00, 2.00)
sigma_grid_VAE_test <- c(0.50, 1.00, 2.00)
rho_grid_KL_test <- rho_from_sigma(sigma_grid_KL_test)
rho_grid_VAE_test <- rho_from_sigma(sigma_grid_VAE_test)

if (isTRUE(test_mode)) {
  rhoGrid_KL <- rho_grid_KL_test
  rhoGrid_VAE <- rho_grid_VAE_test
} else if (identical(estimation_stage, "stage1")) {
  rhoGrid_KL <- rho_grid_KL_stage1
  rhoGrid_VAE <- rho_grid_VAE_stage1
} else if (identical(estimation_stage, "stage2")) {
  # Stage 2 has no global rho grid. Region-specific sigma/rho grids are built
  # from the completed Stage 1 support set inside build_stage2_plan().
  rhoGrid_KL <- numeric(0)
  rhoGrid_VAE <- numeric(0)
} else {
  stop("Unknown estimation_stage: ", estimation_stage)
}

# Grid sorting and rounding
rhoGrid_KL <- sort(unique(round(rhoGrid_KL, 4)))
rhoGrid_VAE <- sort(unique(round(rhoGrid_VAE, 4)))
rho_grid_KL <- rhoGrid_KL
rho_grid_VAE <- rhoGrid_VAE
sigma_grid_KL <- sort(unique(round(1 / (1 + rho_grid_KL), 4)))
sigma_grid_VAE <- sort(unique(round(1 / (1 + rho_grid_VAE), 4)))

#### FUNCTIONS ####

# Counts 1 rho_KL and 1 rho_VAE if there is more than one in each search grid. 
rho_penalty <- as.integer(length(rhoGrid_KL) > 1) + as.integer(length(rhoGrid_VAE) > 1)

# Function to return a default value if the left side is NULL or NA. Used to extract values from missing/empty objects
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || isTRUE(all(is.na(x)))) {
    y
  } else {
    # collapse to scalar if needed taking the first element
    if (length(x) > 1) x[[1]] else x
  }
}

safe_first <- function(x, default) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) {
    return(default)
  }
  x[[1]]
}

# Function to extract boolean TRUE/FALSE from NULL, empty, or non-logical data
safe_bool <- function(x, default = FALSE) {
  out <- suppressWarnings(as.logical(safe_first(x, default)))
  if (is.na(out)) default else out
}

# Function to extract text from NULL or empty data
safe_chr <- function(x, default = NA_character_) {
  as.character(safe_first(x, default))
}

# Function to extract numbers from NULL or empty data
safe_num <- function(x, default = NA_real_) {
  suppressWarnings(as.numeric(safe_first(x, default)))
}

finite_values <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.finite(x)]
}

finite_stat <- function(x, fn) {
  x <- finite_values(x)
  if (length(x) == 0L) NA_real_ else fn(x)
}

finite_min <- function(x) finite_stat(x, min)
finite_max <- function(x) finite_stat(x, max)
finite_median <- function(x) finite_stat(x, median)

finite_first_at_min <- function(value, rank_by) {
  rank_by <- suppressWarnings(as.numeric(rank_by))
  ok <- is.finite(rank_by)
  if (!any(ok)) {
    return(NA_real_)
  }
  value[which(ok)[which.min(rank_by[ok])]][1]
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

make_refined_sigma_grid <- function(sigma_min, sigma_max, sigma_best,
                                    points = 17L, padding_share = 0.25,
                                    lower_bound = 0.03, upper_bound = 8) {
  values <- suppressWarnings(as.numeric(c(sigma_min, sigma_max, sigma_best)))
  values <- values[is.finite(values) & values > 0]
  if (length(values) == 0L) {
    return(numeric(0))
  }
  values <- pmin(pmax(values, lower_bound), upper_bound)

  grid_min <- min(values)
  grid_max <- max(values)
  if (near(grid_min, grid_max, tol = 1e-10)) {
    grid_min <- grid_min * (1 - padding_share)
    grid_max <- grid_max * (1 + padding_share)
  } else {
    width <- grid_max - grid_min
    grid_min <- grid_min - padding_share * width
    grid_max <- grid_max + padding_share * width
  }

  grid_min <- max(lower_bound, grid_min)
  grid_max <- min(upper_bound, grid_max)
  if (grid_min <= 0 || grid_max <= 0 || grid_min > grid_max) {
    return(sort(unique(round(values[values >= lower_bound & values <= upper_bound], 5))))
  }

  # Sigma is positive and often spans orders of magnitude, so log spacing gives
  # useful resolution at both low- and high-substitution ends. Extra points near
  # sigma_best sharpen Stage 2 around the Stage 1-supported optimum.
  global_grid <- exp(seq(log(grid_min), log(grid_max), length.out = points))
  sigma_best <- suppressWarnings(as.numeric(sigma_best))
  best <- if (is.finite(sigma_best) && sigma_best > 0) {
    pmin(pmax(sigma_best, lower_bound), upper_bound)
  } else {
    exp(stats::median(log(values)))
  }
  focus_ratio <- max(1.25, min(2.5, sqrt(grid_max / grid_min)))
  focus_min <- max(grid_min, best / focus_ratio)
  focus_max <- min(grid_max, best * focus_ratio)
  focus_grid <- exp(seq(log(focus_min), log(focus_max), length.out = max(7L, ceiling(points / 2))))
  final_grid <- sort(unique(round(c(global_grid, focus_grid, values), 5)))
  final_grid[final_grid >= lower_bound & final_grid <= upper_bound]
}

make_stage2_boundary_sigma_grid <- function(sigma_min, sigma_max, sigma_best,
                                            points = 17L,
                                            lower_bound = 0.01,
                                            upper_bound = 12) {
  values <- suppressWarnings(as.numeric(c(sigma_min, sigma_max, sigma_best)))
  values <- values[is.finite(values) & values > 0]
  if (length(values) == 0L) {
    return(numeric(0))
  }
  grid_min <- max(lower_bound, min(values))
  grid_max <- min(upper_bound, max(values))
  if (near(grid_min, grid_max, tol = 1e-10)) {
    grid_min <- max(lower_bound, grid_min / 2)
    grid_max <- min(upper_bound, grid_max * 2)
  }
  make_refined_sigma_grid(
    sigma_min = grid_min,
    sigma_max = grid_max,
    sigma_best = sigma_best,
    points = points,
    padding_share = 0,
    lower_bound = lower_bound,
    upper_bound = upper_bound
  )
}

split_optimizer_list <- function(x, max_methods = 3L) {
  if (is.null(x) || length(x) == 0L || all(is.na(x)) || !nzchar(as.character(x[[1]]))) {
    return(character(0))
  }
  methods <- trimws(unlist(strsplit(as.character(x[[1]]), ",")))
  methods <- methods[nzchar(methods)]
  methods <- unique(methods)
  head(methods, max_methods)
}

make_expanded_sigma_grid <- function(best_sigma, support_min, support_max,
                                             previous_grid, points = 9L,
                                             lower_bound = 0.03, upper_bound = 8) {
  previous_grid <- sort(unique(suppressWarnings(as.numeric(previous_grid))))
  previous_grid <- previous_grid[is.finite(previous_grid) & previous_grid > 0]
  best_sigma <- suppressWarnings(as.numeric(best_sigma))
  support_values <- suppressWarnings(as.numeric(c(best_sigma, support_min, support_max)))
  support_values <- support_values[is.finite(support_values) & support_values > 0]
  if (length(previous_grid) == 0L || !is.finite(best_sigma) || best_sigma <= 0) {
    return(sort(unique(round(support_values, 5))))
  }
  grid_min <- previous_grid[[1]]
  grid_max <- previous_grid[[length(previous_grid)]]
  if (near(best_sigma, grid_min, tol = 1e-8)) {
    upper <- max(c(best_sigma * 3, support_values, previous_grid[previous_grid <= 3 * grid_min]), na.rm = TRUE)
    upper <- min(upper_bound, max(upper, lower_bound * 3))
    values <- seq(lower_bound, upper, length.out = points)
  } else if (near(best_sigma, grid_max, tol = 1e-8)) {
    lower <- min(c(best_sigma / 3, support_values, previous_grid[previous_grid >= grid_max / 3]), na.rm = TRUE)
    lower <- max(lower_bound, min(lower, upper_bound / 3))
    values <- seq(lower, upper_bound, length.out = points)
  } else {
    values <- make_refined_sigma_grid(
      sigma_min = support_min,
      sigma_max = support_max,
      sigma_best = best_sigma,
      points = points,
      padding_share = 0.25,
      lower_bound = lower_bound,
      upper_bound = upper_bound
    )
  }
  sort(unique(round(c(values, support_values, best_sigma), 5)))
}

combine_stage1_iteration_grid <- function(current_grid, previous_stage1_results, current_run_label) {
  if (!is.list(previous_stage1_results) || !"grid_results" %in% names(previous_stage1_results)) {
    return(list(grid_results = current_grid, carried_forward_regions = character(0)))
  }
  previous_grid <- previous_stage1_results$grid_results
  if (!all(c("r", "method", "rho_KL", "rho_VAE") %in% names(current_grid)) ||
      !"r" %in% names(previous_grid)) {
    return(list(grid_results = current_grid, carried_forward_regions = character(0)))
  }
  current_regions <- sort(unique(current_grid$r))
  carried_forward_regions <- sort(setdiff(unique(previous_grid$r), current_regions))
  if (length(carried_forward_regions) == 0L) {
    return(list(grid_results = current_grid, carried_forward_regions = character(0)))
  }
  previous_run_label <- previous_stage1_results$metadata$run_label %||% "previous_stage1"
  carried_forward_grid <- previous_grid %>%
    filter(r %in% carried_forward_regions) %>%
    mutate(stage1_source = paste0("carried_forward_from_", previous_run_label))
  current_grid <- current_grid %>%
    mutate(stage1_source = paste0("estimated_in_", current_run_label))
  list(
    grid_results = bind_rows(carried_forward_grid, current_grid),
    carried_forward_regions = carried_forward_regions
  )
}

build_next_stage1_plan <- function(stage1_results,
                                             grid_points = 9L,
                                             sigma_lower_bound = 0.03,
                                             sigma_upper_bound = 8,
                                             max_solvers = 3L) {
  required_objects <- c("metadata", "stage2_design", "method_choice")
  missing_objects <- setdiff(required_objects, names(stage1_results))
  if (length(missing_objects) > 0L) {
    stop("Stage 1 results object is missing: ", paste(missing_objects, collapse = ", "))
  }
  previous_sigma_KL <- stage1_results$metadata$obs_sigma_KL
  previous_sigma_VAE <- stage1_results$metadata$obs_sigma_VAE
  previous_region_grids <- if ("grid_results" %in% names(stage1_results)) {
    stage1_results$grid_results %>%
      group_by(r) %>%
      summarise(
        prev_sigma_KL = list(sort(unique(sigma_KL[is.finite(sigma_KL) & sigma_KL > 0]))),
        prev_sigma_VAE = list(sort(unique(sigma_VAE[is.finite(sigma_VAE) & sigma_VAE > 0]))),
        .groups = "drop"
      )
  } else {
    stage1_results$stage2_design %>%
      distinct(r) %>%
      mutate(
        prev_sigma_KL = list(previous_sigma_KL),
        prev_sigma_VAE = list(previous_sigma_VAE)
      )
  }
  optimizer_candidates <- stage1_results$method_choice %>%
    filter(valid_share > 0) %>%
    group_by(r) %>%
    arrange(
      desc(convergence_share),
      desc(valid_share),
      replace_na(edge_solution_share, 1),
      median_valid_AICc,
      median_runtime_seconds,
      .by_group = TRUE
    ) %>%
    summarise(
      selected_stage1_solvers = paste(head(method, max_solvers), collapse = ", "),
      .groups = "drop"
    )
  region_plan <- stage1_results$stage2_design %>%
    left_join(previous_region_grids, by = "r") %>%
    left_join(optimizer_candidates, by = "r") %>%
    rowwise() %>%
    mutate(
      run_next_stage1 = rho_refinement_action == "expand_stage1_grid_before_stage2" &&
        !is.na(selected_stage1_solvers) && nzchar(selected_stage1_solvers),
      sigma_grid_KL = list(if (isTRUE(run_next_stage1)) {
        make_expanded_sigma_grid(
          best_sigma = best_stage1_sigma_KL,
          support_min = sigma_KL_supported_min,
          support_max = sigma_KL_supported_max,
          previous_grid = prev_sigma_KL,
          points = grid_points,
          lower_bound = sigma_lower_bound,
          upper_bound = sigma_upper_bound
        )
      } else {
        numeric(0)
      }),
      sigma_grid_VAE = list(if (isTRUE(run_next_stage1)) {
        make_expanded_sigma_grid(
          best_sigma = best_stage1_sigma_VAE,
          support_min = sigma_VAE_supported_min,
          support_max = sigma_VAE_supported_max,
          previous_grid = prev_sigma_VAE,
          points = grid_points,
          lower_bound = sigma_lower_bound,
          upper_bound = sigma_upper_bound
        )
      } else {
        numeric(0)
      }),
      rho_grid_KL = list(rho_from_sigma(sigma_grid_KL)),
      rho_grid_VAE = list(rho_from_sigma(sigma_grid_VAE)),
      selected_stage1_solvers_list = list(split_optimizer_list(selected_stage1_solvers, max_solvers)),
      stage1_next_grid_cell_count = length(rho_grid_KL) *
        length(rho_grid_VAE) *
        length(selected_stage1_solvers_list),
      stage1_next_plan_note = case_when(
        run_next_stage1 ~
          "Run targeted Stage 1 expansion: broader boundary-aware sigma grid with Stage 1-supported optimizers.",
        rho_refinement_action == "ready_for_stage2_refinement" ~
          "Do not rerun in Stage 1 expansion; this region is ready for Stage 2 refinement.",
        rho_refinement_action == "boundary_supported_review_before_stage2" ~
          "Do not rerun automatically; targeted Stage 1 still supports a boundary solution.",
        TRUE ~ "Do not rerun automatically; review data, validity, or optimizer support."
      )
    ) %>%
    ungroup()
  list(
    metadata = list(
      created_at = Sys.time(),
    source_stage1 = stage1_results$metadata$run_label %||% NA_character_,
      grid_points = grid_points,
      sigma_lower_bound = sigma_lower_bound,
      sigma_upper_bound = sigma_upper_bound,
      max_solvers = max_solvers
    ),
    region_plan = region_plan,
    regions_to_estimate = region_plan %>%
      filter(run_next_stage1) %>%
      pull(r),
    regions_ready_for_stage2 = region_plan %>%
      filter(rho_refinement_action == "ready_for_stage2_refinement") %>%
      pull(r)
  )
}

build_stage2_plan <- function(stage1_results,
                                  grid_points = 9L,
                                  window_padding = 0.25,
                                  max_solvers = 3L,
                                  sigma_lower_bound = 0.05,
                                  sigma_upper_bound = 5,
                                  boundary_sigma_lower_bound = sigma_lower_bound,
                                  boundary_sigma_upper_bound = sigma_upper_bound,
                                  solver_pool = c("PORT", "L-BFGS-B")) {
  required_objects <- c("stage2_design", "rho_windows", "method_choice")
  missing_objects <- setdiff(required_objects, names(stage1_results))
  if (length(missing_objects) > 0L) {
    stop("Stage 1 results object is missing: ", paste(missing_objects, collapse = ", "))
  }

  design <- stage1_results$stage2_design
  region_plan <- design %>%
    rowwise() %>%
    mutate(
      selected_stage2_solvers_list = list({
        ranked_solvers <- split_optimizer_list(selected_stage2_solvers, max_methods = max_solvers)
        bounded_solvers <- intersect(ranked_solvers, solver_pool)
        if (length(bounded_solvers) > 0L) bounded_solvers else solver_pool
      }),
      has_stage2_optimizer_support = length(selected_stage2_solvers_list) > 0L,
      stage2_estimation_class = case_when(
        rho_refinement_action == "ready_for_stage2_refinement" ~ "interior-supported",
        rho_refinement_action == "boundary_supported_review_before_stage2" ~ "boundary-supported",
        TRUE ~ "not-estimated"
      ),
      estimate_in_stage2 =
        rho_refinement_action %in% c(
          "ready_for_stage2_refinement",
          "boundary_supported_review_before_stage2"
        ) &
        isTRUE(has_stage2_optimizer_support),
      sigma_grid_KL = list(
        if (isTRUE(estimate_in_stage2)) {
          if (stage2_estimation_class == "boundary-supported") {
            make_stage2_boundary_sigma_grid(
              sigma_KL_supported_min, sigma_KL_supported_max, best_stage1_sigma_KL,
              points = grid_points,
              lower_bound = boundary_sigma_lower_bound,
              upper_bound = boundary_sigma_upper_bound
            )
          } else {
            make_refined_sigma_grid(
              sigma_KL_supported_min, sigma_KL_supported_max, best_stage1_sigma_KL,
              points = grid_points,
              padding_share = window_padding,
              lower_bound = sigma_lower_bound,
              upper_bound = sigma_upper_bound
            )
          }
        } else {
          numeric(0)
        }
      ),
      sigma_grid_VAE = list(
        if (isTRUE(estimate_in_stage2)) {
          if (stage2_estimation_class == "boundary-supported") {
            make_stage2_boundary_sigma_grid(
              sigma_VAE_supported_min, sigma_VAE_supported_max, best_stage1_sigma_VAE,
              points = grid_points,
              lower_bound = boundary_sigma_lower_bound,
              upper_bound = boundary_sigma_upper_bound
            )
          } else {
            make_refined_sigma_grid(
              sigma_VAE_supported_min, sigma_VAE_supported_max, best_stage1_sigma_VAE,
              points = grid_points,
              padding_share = window_padding,
              lower_bound = sigma_lower_bound,
              upper_bound = sigma_upper_bound
            )
          }
        } else {
          numeric(0)
        }
      ),
      rho_grid_KL = list(rho_from_sigma(sigma_grid_KL)),
      rho_grid_VAE = list(rho_from_sigma(sigma_grid_VAE)),
      stage2_grid_cell_count =
        length(rho_grid_KL) *
        length(rho_grid_VAE) *
        length(selected_stage2_solvers_list),
      stage2_plan_note = case_when(
        estimate_in_stage2 & stage2_estimation_class == "interior-supported" ~
          "Estimate final interior-supported Stage 2 grid with region-specific sigma windows and Stage 1-ranked optimizers.",
        estimate_in_stage2 & stage2_estimation_class == "boundary-supported" ~
          "Estimate final boundary-supported Stage 2 grid; keep the boundary status in the MERGE parameter review.",
        rho_refinement_action == "ready_for_stage2_refinement" & !has_stage2_optimizer_support ~
          "Do not run Stage 2 yet; Stage 1 did not identify a reliable optimizer for this region.",
        rho_refinement_action == "expand_stage1_grid_before_stage2" ~
          "Do not run Stage 2 yet; expand Stage 1 grid because the supported optimum is on a boundary.",
        rho_refinement_action == "boundary_supported_review_before_stage2" ~
          "Do not run Stage 2 yet; Stage 1 did not identify a reliable optimizer for this boundary-supported region.",
        TRUE ~ "Do not estimate in Stage 2 without review."
      )
    ) %>%
    ungroup()

  list(
    metadata = list(
      created_at = Sys.time(),
      source_stage1_created_at = stage1_results$metadata$created_at %||% NA,
      grid_points = grid_points,
      window_padding = window_padding,
      max_solvers = max_solvers,
      sigma_lower_bound = sigma_lower_bound,
      sigma_upper_bound = sigma_upper_bound,
      boundary_sigma_lower_bound = boundary_sigma_lower_bound,
      boundary_sigma_upper_bound = boundary_sigma_upper_bound,
      solver_pool = paste(solver_pool, collapse = ", ")
    ),
    region_plan = region_plan,
    regions_to_estimate = region_plan %>% filter(estimate_in_stage2) %>% pull(r),
    regions_requiring_stage1_expansion = region_plan %>%
      filter(rho_refinement_action == "expand_stage1_grid_before_stage2") %>%
      pull(r),
    regions_requiring_manual_review = region_plan %>%
      filter(!estimate_in_stage2, rho_refinement_action != "expand_stage1_grid_before_stage2") %>%
      pull(r)
  )
}

prepare_grid_results_for_postprocessing <- function(grid_results_raw, normalised_input_data,
                                                    obs_rho_KL,
                                                    obs_rho_VAE,
                                                    obs_rho_penalty) {
  required_columns <- c("r", "method", "conv", "rho_KL", "rho_VAE", "rss")
  missing_required_columns <- setdiff(required_columns, names(grid_results_raw))
  if (length(missing_required_columns) > 0L) {
    stop(
      "Saved grid results are missing required columns for post-processing: ",
      paste(missing_required_columns, collapse = ", ")
    )
  }

  observation_diagnostics <- normalised_input_data %>%
    group_by(r) %>%
    summarise(
      n_obs_for_fit = n(),
      TSS_log_for_fit = {
        log_output <- log(Ys)
        sum((log_output - mean(log_output, na.rm = TRUE))^2, na.rm = TRUE)
      },
      .groups = "drop"
    )

  grid_results_prepared <- grid_results_raw %>%
    left_join(observation_diagnostics, by = "r") %>%
    mutate(
      across(
        any_of(c(
          "rho_KL", "rho_VAE", "rss", "gamma", "lambda", "delta_KVA", "delta_VAY",
          "nu", "runtime_total", "runtime_per_grid", "iter"
        )),
        ~ suppressWarnings(as.numeric(.))
      ),
      conv = replace_na(suppressWarnings(as.logical(conv)), FALSE),
      sigma_KL = if_else(is.finite(rho_KL) & rho_KL > -1, 1 / (1 + rho_KL), NA_real_),
      sigma_VAE = if_else(is.finite(rho_VAE) & rho_VAE > -1, 1 / (1 + rho_VAE), NA_real_),
      k_hat_postprocessed = 5L,
      k_plus_rho_postprocessed = k_hat_postprocessed + obs_rho_penalty,
      R2 = if_else(
        is.finite(rss) & is.finite(TSS_log_for_fit) & TSS_log_for_fit > 0,
        1 - rss / TSS_log_for_fit,
        NA_real_
      ),
      adjR2 = if_else(
        is.finite(R2) & (n_obs_for_fit - k_hat_postprocessed - 1) > 0,
        1 - (1 - R2) * (n_obs_for_fit - 1) / (n_obs_for_fit - k_hat_postprocessed - 1),
        NA_real_
      ),
      AIC_naive = if_else(
        is.finite(rss) & rss > 0 & n_obs_for_fit > 0,
        n_obs_for_fit * log(rss / n_obs_for_fit) + 2 * k_hat_postprocessed,
        NA_real_
      ),
      AICc_naive = if_else(
        is.finite(AIC_naive) & (n_obs_for_fit - k_hat_postprocessed - 1) > 0,
        AIC_naive + (2 * k_hat_postprocessed * (k_hat_postprocessed + 1)) /
          (n_obs_for_fit - k_hat_postprocessed - 1),
        NA_real_
      ),
      AIC_plusRho = if_else(
        is.finite(rss) & rss > 0 & n_obs_for_fit > 0,
        n_obs_for_fit * log(rss / n_obs_for_fit) + 2 * k_plus_rho_postprocessed,
        NA_real_
      ),
      AICc_plusRho = if_else(
        is.finite(AIC_plusRho) & (n_obs_for_fit - k_plus_rho_postprocessed - 1) > 0,
        AIC_plusRho + (2 * k_plus_rho_postprocessed * (k_plus_rho_postprocessed + 1)) /
          (n_obs_for_fit - k_plus_rho_postprocessed - 1),
        NA_real_
      )
    ) %>%
    group_by(r) %>%
    mutate(
      on_edge_KL = on_grid_edge(rho_KL, unique(rho_KL[is.finite(rho_KL)])),
      on_edge_VAE = on_grid_edge(rho_VAE, unique(rho_VAE[is.finite(rho_VAE)]))
    ) %>%
    ungroup() %>%
    select(-n_obs_for_fit, -TSS_log_for_fit, -k_hat_postprocessed, -k_plus_rho_postprocessed)

  for (missing_column in c(
    "msg", "gamma", "lambda", "delta_KVA", "delta_VAY", "nu",
    "se_gamma", "se_lambda", "se_delta_KVA", "se_delta_VAY", "se_nu",
    "t_gamma", "t_lambda", "t_delta_KVA", "t_delta_VAY", "t_nu",
    "p_gamma", "p_lambda", "p_delta_KVA", "p_delta_VAY", "p_nu",
    "iter", "runtime_total", "runtime_per_grid", "n_grid"
  )) {
    if (!missing_column %in% names(grid_results_prepared)) {
      grid_results_prepared[[missing_column]] <- if (missing_column == "msg") NA_character_ else NA_real_
    }
  }

  grid_results_prepared
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

# Function to count free non-rho parameters for information criteria in grid-search estimates.
# micEconCES documents rho1/rho2/rho grid values as fixed during each cell's least-squares fit; this
# keeps AIC_naive tied to the fitted free parameters and AIC_plusRho tied to the grid-selected rhos.
count_free_params <- function(fit_obj, rho_penalty) {
  cf <- tryCatch(stats::coef(fit_obj), error = function(e) numeric(0))
  if (length(cf) == 0L) {
    return(NA_integer_)
  }
  nms <- names(cf)
  if (!is.null(nms) && length(nms) == length(cf)) {
    free <- !grepl("^rho(_?1|_?2)?$|^rho[12]$", nms, ignore.case = TRUE)
    return(sum(free, na.rm = TRUE))
  }
  max(1L, length(cf) - rho_penalty)
}

ces_solver_args <- function(method, stage) {
  stage_lower <- if (identical(stage, "stage2")) stage2_lower else stage1_lower
  stage_upper <- if (identical(stage, "stage2")) stage2_upper else stage1_upper

  switch(
    method,
    Newton = list(
      start = c(
        gamma = runif(1, 0.9, 1.1),
        lambda = runif(1, -0.005, 0.005),
        delta_KVA = runif(1, 0.4, 0.6),
        delta_VAY = runif(1, 0.4, 0.6),
        nu = runif(1, 0.9, 1.1)
      )
    ),
    `L-BFGS-B` = list(start = start_vals, lower = stage_lower, upper = stage_upper, control = list(maxit = 10000, factr = 1e9)),
    PORT = list(
      start = start_vals,
      lower = if (identical(stage, "stage2")) stage2_lower else NULL,
      upper = if (identical(stage, "stage2")) stage2_upper else NULL,
      control = list(eval.max = 1e5, iter.max = 1e5, reltol = 1e-8)
    ),
    BFGS = list(start = start_vals, lower = share_lower, upper = share_upper, control = list(maxit = 10000, reltol = 1e-8)),
    CG = list(start = start_vals, lower = share_lower, upper = share_upper, control = list(maxit = 1000, reltol = 1e-8)),
    NM = list(lower = share_lower, upper = share_upper, control = list(maxit = 10000, reltol = 1e-8)),
    LM = list(lower = share_lower, upper = share_upper, control = list(maxiter = 10000, ftol = 1e-8, maxfev = 5000)),
    SANN = list(lower = share_lower, upper = share_upper, control = list(maxit = 10000, temp = 10, tmax = 50)),
    DE = list(lower = stage_lower, upper = stage_upper, control = list(itermax = 500)),
    list()
  )
}

# Function to inspect source macro data before estimation. This is a research-methods gate:
# regions with missing/non-positive inputs or unusable base-year anchors should be repaired before
# Stage 1 unless estimate_flagged_regions is explicitly set to TRUE.
build_input_data_summary <- function(raw_input_data, base_year) {
  required_columns <- c("r", "t", "Y", "K", "L", "E")
  missing_columns <- setdiff(required_columns, names(raw_input_data))
  if (length(missing_columns) > 0L) {
    stop("Input file is missing required columns: ", paste(missing_columns, collapse = ", "))
  }

  duplicate_rows <- raw_input_data %>%
    count(r, t, name = "duplicate_count") %>%
    filter(duplicate_count > 1L)
  if (nrow(duplicate_rows) > 0L) {
    stop("Input file has duplicated r,t rows. Repair before estimation.")
  }

  diagnostic_data <- raw_input_data %>%
    arrange(r, t) %>%
    group_by(r) %>%
    mutate(
      K_Y = K / Y,
      E_Y = E / Y,
      L_Y = L / Y,
      Y_growth = log(Y) - lag(log(Y)),
      K_growth = log(K) - lag(log(K)),
      L_growth = log(L) - lag(log(L)),
      E_growth = log(E) - lag(log(E))
    ) %>%
    ungroup()

  base_year_values <- raw_input_data %>%
    filter(t == base_year) %>%
    transmute(
      r,
      Ybase = Y,
      Kbase = K,
      Lbase = L,
      Ebase = E
    )

  normalized_ranges <- raw_input_data %>%
    left_join(base_year_values, by = "r") %>%
    mutate(
      Y_norm = Y / Ybase,
      K_norm = K / Kbase,
      L_norm = L / Lbase,
      E_norm = E / Ebase
    ) %>%
    group_by(r) %>%
    summarise(
      Y_norm_min = min(Y_norm, na.rm = TRUE),
      Y_norm_max = max(Y_norm, na.rm = TRUE),
      K_norm_min = min(K_norm, na.rm = TRUE),
      K_norm_max = max(K_norm, na.rm = TRUE),
      L_norm_min = min(L_norm, na.rm = TRUE),
      L_norm_max = max(L_norm, na.rm = TRUE),
      E_norm_min = min(E_norm, na.rm = TRUE),
      E_norm_max = max(E_norm, na.rm = TRUE),
      .groups = "drop"
    )

  diagnostic_data %>%
    group_by(r) %>%
    summarise(
      n_observations = n(),
      first_year = min(t, na.rm = TRUE),
      last_year = max(t, na.rm = TRUE),
      has_base_year = any(t == base_year),
      missing_input_count = sum(is.na(Y) | is.na(K) | is.na(L) | is.na(E)),
      nonpositive_input_count = sum(Y <= 0 | K <= 0 | L <= 0 | E <= 0, na.rm = TRUE),
      K_Y_base = K_Y[t == base_year][1],
      E_Y_base = E_Y[t == base_year][1],
      L_Y_base = L_Y[t == base_year][1],
      K_Y_min = min(K_Y, na.rm = TRUE),
      K_Y_max = max(K_Y, na.rm = TRUE),
      E_Y_min = min(E_Y, na.rm = TRUE),
      E_Y_max = max(E_Y, na.rm = TRUE),
      Y_growth_min = min(Y_growth, na.rm = TRUE),
      Y_growth_max = max(Y_growth, na.rm = TRUE),
      K_growth_min = min(K_growth, na.rm = TRUE),
      K_growth_max = max(K_growth, na.rm = TRUE),
      E_growth_min = min(E_growth, na.rm = TRUE),
      E_growth_max = max(E_growth, na.rm = TRUE),
      L_growth_min = min(L_growth, na.rm = TRUE),
      L_growth_max = max(L_growth, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(base_year_values, by = "r") %>%
    left_join(normalized_ranges, by = "r") %>%
    mutate(
      base_year_anchor_invalid =
        !has_base_year |
        !is.finite(Ybase) | !is.finite(Kbase) |
        !is.finite(Lbase) | !is.finite(Ebase) |
        Ybase <= 0 | Kbase <= 0 | Lbase <= 0 | Ebase <= 0,
      very_large_annual_jump =
        pmax(abs(Y_growth_min), abs(Y_growth_max),
             abs(K_growth_min), abs(K_growth_max),
             abs(L_growth_min), abs(L_growth_max),
             abs(E_growth_min), abs(E_growth_max),
             na.rm = TRUE) > 0.5,
      flat_series_warning =
        (Y_norm_max - Y_norm_min) < 0.01 |
        (K_norm_max - K_norm_min) < 0.01 |
        (L_norm_max - L_norm_min) < 0.01 |
        (E_norm_max - E_norm_min) < 0.01,
      norm_explosion_warning =
        Y_norm_max > 10 | K_norm_max > 10 |
        L_norm_max > 10 | E_norm_max > 10,
      K_Y_warning =
        K_Y_base < 1 | K_Y_base > 12 |
        K_Y_min < 1 | K_Y_max > 12,
      short_sample_warning = n_observations < 20,
      input_data_class = case_when(
        missing_input_count > 0 | nonpositive_input_count > 0 | base_year_anchor_invalid ~
          "exclude_or_repair_before_stage1",
        very_large_annual_jump | flat_series_warning | norm_explosion_warning |
          K_Y_warning | short_sample_warning ~
          "usable_with_warning",
        TRUE ~ "usable"
      ),
      input_data_note = pmap_chr(
        list(
          missing_input_count, nonpositive_input_count, base_year_anchor_invalid,
          very_large_annual_jump, flat_series_warning, norm_explosion_warning,
          K_Y_warning, short_sample_warning
        ),
        function(missing_count, nonpositive_count, anchor_invalid, large_jump,
                 flat_series, norm_explosion, K_Y_warning, short_sample) {
          notes <- c()
          if (missing_count > 0) notes <- c(notes, "missing Y/K/L/E")
          if (nonpositive_count > 0) notes <- c(notes, "non-positive Y/K/L/E")
          if (isTRUE(anchor_invalid)) notes <- c(notes, "invalid or missing base-year anchor")
          if (isTRUE(large_jump)) notes <- c(notes, "large annual log jump")
          if (isTRUE(flat_series)) notes <- c(notes, "near-flat normalised series")
          if (isTRUE(norm_explosion)) notes <- c(notes, "normalised variable exceeds 10")
          if (isTRUE(K_Y_warning)) notes <- c(notes, "K/Y outside [1,12]")
          if (isTRUE(short_sample)) notes <- c(notes, "short sample")
          if (length(notes) == 0L) "OK" else paste(notes, collapse = "; ")
        }
      )
    )
}

# Function to normalise macro data using the checked base-year anchors from input_data_summary.
# This prevents Stage 1 from silently estimating regions where the normalisation denominator is
# missing, zero, or non-positive.
build_normalised_input_data <- function(raw_input_data, input_data_summary, base_year,
                                        estimate_flagged_regions = FALSE) {
  regions_to_exclude <- input_data_summary %>%
    filter(input_data_class == "exclude_or_repair_before_stage1") %>%
    pull(r)

  if (length(regions_to_exclude) > 0L && !isTRUE(estimate_flagged_regions)) {
    message(
      "Excluding regions that require repair before Stage 1: ",
      paste(regions_to_exclude, collapse = ", ")
    )
  }

  regions_to_estimate <- input_data_summary %>%
    filter(input_data_class != "exclude_or_repair_before_stage1" | isTRUE(estimate_flagged_regions)) %>%
    pull(r)

  if (length(regions_to_estimate) == 0L) {
    stop("No regions passed the input-data gate. Repair input data or set estimate_flagged_regions = TRUE explicitly.")
  }

  raw_input_data %>%
    filter(r %in% regions_to_estimate) %>%
    left_join(
      input_data_summary %>%
        select(r, Ybase, Kbase, Lbase, Ebase),
      by = "r"
    ) %>%
    mutate(
      Ys = Y / Ybase,
      Ks = K / Kbase,
      Ls = L / Lbase,
      Es = E / Ebase,
      time_from_base_year = t - base_year
    )
}

# Function to extract statistics from each grid combination
# Extracts standard errors, t-statistics, p-values, and confidence intervals from the cesEst objects created by the package
# These inference statistics are diagnostics only. Regional macro time-series estimates are exposed
# to trending variables, serial correlation, heteroskedasticity, endogeneity, and small-sample limits.
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
  
  # standard errors from the variance-covariance matrix if summary tables do not provide them
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
saved_stage_results <- NULL
previous_stage1_results <- NULL
next_stage1_plan <- NULL
stage2_plan <- NULL
carried_forward_stage1_regions <- character(0)
if (isTRUE(estimate)) {
  if (identical(backend, "multisession")) {
    message(
      "Parallel backend: multisession with ", workers,
      " worker(s). Use backend <- \"sequential\" only when checking or when a single quiet R session is required."
    )
    suppressWarnings({
      tryCatch({
        future::plan(future::multisession, workers = workers)
      }, error = function(e) {
        warning("multisession failed (", conditionMessage(e), "); using sequential backend.")
        backend <<- "sequential"
        future::plan(future::sequential)
      })
    })
  } else if (identical(backend, "sequential")) {
    message("Parallel backend: sequential. This avoids Windows worker popups and is easiest to reproduce.")
    future::plan(future::sequential)
  } else {
    stop('backend must be "sequential" or "multisession".')
  }
  
  # Data loading and normalisation
  df <- read_csv(input_file, show_col_types = FALSE) # read the dataset
  input_data_summary <- build_input_data_summary(df, base_year)
  safe_write_csv(input_data_summary, input_data_summary_file)

  normalised_input_data <- build_normalised_input_data(
    raw_input_data = df,
    input_data_summary = input_data_summary,
    base_year = base_year,
    estimate_flagged_regions = estimate_flagged_regions
  )
  if (isTRUE(test_mode)) {
    test_r <- normalised_input_data %>%
      distinct(r) %>%
      arrange(r) %>%
      slice_head(n = test_n) %>%
      pull(r)
    message(
      "Test run: estimating ",
      length(test_r),
      " region(s): ",
      paste(test_r, collapse = ", ")
    )
    normalised_input_data <- normalised_input_data %>%
      filter(r %in% test_r)
  }
    dfS <- normalised_input_data # short name used by the estimation code below

  if (
    identical(estimation_stage, "stage1") &&
    !is.na(previous_stage1_file) &&
    nzchar(previous_stage1_file)
  ) {
    if (!file.exists(previous_stage1_file)) {
      stop("Stage 1 targeted iteration requires a previous Stage 1 results file. Not found: ", previous_stage1_file)
    }
    previous_stage1_results <- readRDS(previous_stage1_file)
    next_stage1_plan <- build_next_stage1_plan(
      stage1_results = previous_stage1_results,
      grid_points = stage1_expand_points,
      sigma_lower_bound = stage1_sigma_min,
      sigma_upper_bound = stage1_sigma_max,
      max_solvers = stage1_max_solvers
    )
    next_stage1_regions <- next_stage1_plan$regions_to_estimate
    if (length(next_stage1_regions) == 0L) {
      stop("No regions require targeted Stage 1 expansion according to the previous Stage 1 results.")
    }
    message("Targeted Stage 1 iteration plan: estimating ", length(next_stage1_regions), " region(s).")
    dfS <- dfS %>% filter(r %in% next_stage1_regions)
  }

  if (identical(estimation_stage, "stage2")) {
    if (!file.exists(stage1_file_for_stage2)) {
      stop("Stage 2 requires a Stage 1 results file. Not found: ", stage1_file_for_stage2)
    }
    stage1_results_for_stage2 <- readRDS(stage1_file_for_stage2)
    stage2_plan <- build_stage2_plan(
      stage1_results = stage1_results_for_stage2,
      grid_points = stage2_grid_points,
      window_padding = stage2_window_padding,
      max_solvers = stage2_max_solvers,
      sigma_lower_bound = stage2_sigma_min,
      sigma_upper_bound = stage2_sigma_max,
      boundary_sigma_lower_bound = stage2_boundary_sigma_min,
      boundary_sigma_upper_bound = stage2_boundary_sigma_max,
      solver_pool = stage2_solver_pool
    )
    if (length(stage2_plan$regions_requiring_stage1_expansion) > 0L) {
      stop(
        "Stage 2 is not ready. Run the next targeted Stage 1 iteration for ",
        length(stage2_plan$regions_requiring_stage1_expansion),
        " region(s) first. Use MERGE CES stage1 sequence.R to do this automatically."
      )
    }
    regions_ready_for_stage2 <- stage2_plan$regions_to_estimate
    if (length(regions_ready_for_stage2) == 0L) {
      stop("No regions are ready for Stage 2 estimation according to the Stage 1 run plan.")
    }
    message("Stage 2 run plan: estimating ", length(regions_ready_for_stage2), " region(s).")
    dfS <- dfS %>% filter(r %in% regions_ready_for_stage2)
  }

  # Region estimation loop. Estimates the parameters for each region for each grid point and with each solver
  # The code in this estimate_region runs for each region
  quiet_worker_packages <- function() {
    invisible(lapply(c("micEconCES", "dplyr", "tidyr", "purrr", "readr"), function(pkg) {
      suppressPackageStartupMessages(suppressWarnings(library(pkg, character.only = TRUE)))
    }))
  }

  estimate_region <- function(d, region_name) {
    quiet_worker_packages()
    message("\nEstimating region: ", region_name)
    d_num <- d %>% transmute(t, time_from_base_year, Ys, Ks, Ls, Es) # grabbing numeric values from the data
    local_rhoGrid_KL <- rhoGrid_KL
    local_rhoGrid_VAE <- rhoGrid_VAE
    local_stage1_optimizers <- NULL
    local_stage2_optimizers <- NULL
    if (
      identical(estimation_stage, "stage1") &&
      !is.null(next_stage1_plan)
    ) {
      region_plan <- next_stage1_plan$region_plan %>% filter(r == region_name) %>% slice_head(n = 1)
      if (nrow(region_plan) == 1L && isTRUE(region_plan$run_next_stage1[[1]])) {
        local_rhoGrid_KL <- sort(unique(round(unlist(region_plan$rho_grid_KL[[1]]), 5)))
        local_rhoGrid_VAE <- sort(unique(round(unlist(region_plan$rho_grid_VAE[[1]]), 5)))
        local_stage1_optimizers <- unlist(region_plan$selected_stage1_solvers_list[[1]])
      }
    }
    if (identical(estimation_stage, "stage2") && !is.null(stage2_plan)) {
      region_plan <- stage2_plan$region_plan %>% filter(r == region_name) %>% slice_head(n = 1)
      if (nrow(region_plan) == 1L && isTRUE(region_plan$estimate_in_stage2[[1]])) {
        local_rhoGrid_KL <- sort(unique(round(unlist(region_plan$rho_grid_KL[[1]]), 5)))
        local_rhoGrid_VAE <- sort(unique(round(unlist(region_plan$rho_grid_VAE[[1]]), 5)))
        local_stage2_optimizers <- unlist(region_plan$selected_stage2_solvers_list[[1]])
      }
    }
    # Setting a deterministic seed for replicability. Different seeds per region
    seed_val <- sum(utf8ToInt(region_name)) # converts region text to numbers and sums them to create an ID
    seed_val <- abs(seed_val) %% .Machine$integer.max # limiting seed to a valid range
    set.seed(seed_val) # sets the seed to the region numeric ID
    
    # Choose the optimisation methods to use. micEconCES documents analytical gradients for LM,
    # BFGS, CG, L-BFGS-B, Newton, and PORT. In micEconCES::cesEst source, method == "NM" is
    # immediately recoded to "Nelder-Mead", so only NM is retained here to avoid duplicate work.
    opt_methods <- if (isTRUE(test_mode)) {
      test_solvers
    } else if (identical(estimation_stage, "stage1")) {
      if (!is.null(local_stage1_optimizers) && length(local_stage1_optimizers) > 0L) {
        local_stage1_optimizers
      } else {
        stage1_solvers
      }
    } else if (identical(estimation_stage, "stage2")) {
      if (!is.null(local_stage2_optimizers) && length(local_stage2_optimizers) > 0L) {
        local_stage2_optimizers
      } else {
        stop("Stage 2 requires region-specific optimizers selected from Stage 1 diagnostics for region: ", region_name)
      }
    } else {
      stop("Unknown estimation_stage: ", estimation_stage)
    }
    
    # Setting up objects to store the estimations by method
    fit_all <- setNames(vector("list", length(opt_methods)), opt_methods)
    conv_all <- setNames(rep(FALSE, length(opt_methods)), opt_methods)
    msg_all <- setNames(rep(NA_character_, length(opt_methods)), opt_methods)
    times_all <- setNames(rep(NA_real_, length(opt_methods)), opt_methods)
    
    for (m in opt_methods) {
      t0 <- Sys.time() # time start for each method
      solver_args <- ces_solver_args(m, estimation_stage)
      
      # Arguments passed to the MicEconCES::cesEst() package
      args_list <- list(
        yName = "Ys",
        xNames = c("Ks","Ls","Es"),
        tName = "time_from_base_year",
        data = d_num,
        vrs = TRUE, # nu is estimated (variable returns to scale)
        multErr = TRUE, # multiplicative error term in CES
        method = m,
        rho1 = local_rhoGrid_KL, # assigning the K-L rho grid to rho1 (package name)
        rho = local_rhoGrid_VAE, # assigning the VA-E rho grid to rho (package name)
        returnGridAll = TRUE
        )
      # These ifs pass the method-specific starting, min/max and control arguments. If empty, it uses the package standards
      args_list <- c(args_list, solver_args[!vapply(solver_args, is.null, logical(1))])
      
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
                            if (conv_flag) "OK" else "NO", m, runtime,
                            if (!conv_flag && nzchar(msg_flag)) paste0(" msg: ", msg_flag) else ""))
    }
    
    # Returns all estimations, convergence, messages and runtimes for the region
    region_fit_result <- list(
      checkpoint_key = checkpoint_run_key(region_name, d),
      checkpoint_created_at = Sys.time(),
      fits = fit_all,
      conv = conv_all,
      msg = msg_all,
      times = times_all,
      data = d_num,
      rhoGrid_KL = local_rhoGrid_KL,
      rhoGrid_VAE = local_rhoGrid_VAE
    )

    checkpoint_file <- region_fit_checkpoint_file(region_name)
    save_region_checkpoint(region_fit_result, checkpoint_file)

    region_fit_result
}

estimate_region_or_checkpoint <- function(d, region_name) {
  quiet_worker_packages()
  checkpoint_file <- region_fit_checkpoint_file(region_name)
  current_key <- checkpoint_run_key(region_name, d)
  if (isTRUE(resume_regions) && file.exists(checkpoint_file)) {
    checkpoint <- readRDS(checkpoint_file)
    reuse_status <- checkpoint_reuse_status(checkpoint, current_key)
    if (identical(reuse_status, "ok")) {
      message("\nUsing checkpoint for region: ", region_name)
      return(checkpoint)
    }
    message("\nIgnoring checkpoint for region: ", region_name, " (", reuse_status, ")")
  }
  estimate_region(d, region_name)
}

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
    region_rhoGrid_KL <- if (!is.null(region_fits$rhoGrid_KL)) {
      sort(unique(round(as.numeric(region_fits$rhoGrid_KL), 5)))
    } else {
      rhoGrid_KL
    }
    region_rhoGrid_VAE <- if (!is.null(region_fits$rhoGrid_VAE)) {
      sort(unique(round(as.numeric(region_fits$rhoGrid_VAE), 5)))
    } else {
      rhoGrid_VAE
    }
    region_rho_penalty <- as.integer(length(region_rhoGrid_KL) > 1) +
      as.integer(length(region_rhoGrid_VAE) > 1)
    # Creating an empty table to store all the results later
    grid_tbl <- tibble()
    
    for (m in names(region_fits$fits)) {
      fit_obj <- region_fits$fits[[m]]
      conv_flag <- safe_bool(region_fits$conv[[m]], FALSE)
      msg_flag <- safe_chr(region_fits$msg[[m]], NA_character_)
      runtime <- safe_num(region_fits$times[[m]], NA_real_)
      
      # Creating a full_grid table for each region × method
      # testing the full grid of all combinations of rho_KL and rho_VAE
      full_grid <- expand_grid(rho_KL = region_rhoGrid_KL, rho_VAE = region_rhoGrid_VAE) %>%
        mutate(
          r = region_name,
          method = m,
          n_grid = n(), # grid points per region-method
          runtime_total = runtime, # runtime for the region-method run
          runtime_per_grid = runtime / n_grid, # average runtime per grid point
          msg = msg_flag, # solver message per method
          conv = conv_flag,  # will be overridden by grid-level convergence if available
          # Calculating constant elasticity of substitution
          sigma_KL = ifelse(is.finite(1/(1+rho_KL)), 1/(1+rho_KL), NA_real_),
          sigma_VAE = ifelse(is.finite(1/(1+rho_VAE)), 1/(1+rho_VAE), NA_real_),
          # Flags for rho at the boundary of the grids
          on_edge_KL = on_grid_edge(rho_KL,  region_rhoGrid_KL),
          on_edge_VAE = on_grid_edge(rho_VAE, region_rhoGrid_VAE),
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
          if (!"convergence" %in% names(sum_tbl)) {
            sum_tbl$convergence <- NA
          }
          # assume columns rho1 & rho & convergence exist for nested case
          sum_tbl <- sum_tbl %>%
            mutate(
              rho1 = as.numeric(.data[["rho1"]]),
              rho = as.numeric(.data[["rho"]]),
              rss = as.numeric(.data[["rss"]]),
              conv_grid = as.logical(.data[["convergence"]])
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
        
        # Per-grid coefficients from allRhoFull. This is useful for refined/final analysis, but expensive
        # in broad Stage 1 because it keeps every grid-cell fit object for every optimiser.
        coef_df <- NULL
        if (!is.null(fit_obj$allRhoFull) && length(fit_obj$allRhoFull) > 0L) {
          coef_df <- purrr::map_dfr(fit_obj$allRhoFull, extract_grid_coeffs)
        } else {
          # If returnGridAll = FALSE, micEconCES still returns the best grid-cell estimate. Attach
          # those coefficients to the corresponding rho cell so Stage 1 can rank optimisers without
          # materialising every grid-cell fit.
          coef_df <- extract_grid_coeffs(fit_obj)
        }

        if (!is.null(coef_df) && nrow(coef_df) > 0L) {
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
        # AIC/AICc are model-selection heuristics here. The likelihood is not separately estimated;
        # selected fixed-grid rho values are counted in AIC_plusRho because they are chosen over the grid.
        obs_log <- try(log(d_num$Ys + 1e-12), silent = TRUE) # adding 1e-12 to avoid log(0)
        if (!inherits(obs_log, "try-error")) {
          n_obs <- length(obs_log)
          TSS_log <- sum((obs_log - mean(obs_log))^2) # total sum of squares in log
          k_hat <- count_free_params(fit_obj, region_rho_penalty) # free non-rho parameters, usually gamma/lambda/delta_1/delta/nu
          k0 <- k_hat + region_rho_penalty # adds selected grid-search rho parameters, usually rho_KL and rho_VAE
          
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
  
  
if (isTRUE(estimate)) {

# Region splits of the normalised data (dfS)
splits <- dfS %>% group_split(r, .keep=TRUE) 
region_names <- dfS %>% distinct(r) %>% pull(r) # extract the region names from dfS



#### RUN ALL ####
# Running estimations in parallel across regions
fits_all <- future_map2(
  splits, region_names,
  ~ tryCatch(
    estimate_region_or_checkpoint(.x, .y),
    error = function(e) {
      message("Completely failed region: ", .y, " -> ", e$message)
      list(
        fits = list(),
        conv = list(),
        msg = list(error = e$message),
        times = list(),
        data = .x,
        rhoGrid_KL = rhoGrid_KL,
        rhoGrid_VAE = rhoGrid_VAE
      )
    }
  ),
  .progress = TRUE,
  .options = furrr_options(
    packages = character(0), # worker packages are loaded quietly by quiet_worker_packages()
    conditions = structure("condition", exclude = "warning"), # avoids repeated package-version warnings from workers
    globals = c(
      "rhoGrid_KL", "rhoGrid_VAE", "estimation_stage", "test_mode",
      "next_stage1_plan", "stage2_plan", "stage1_solvers", "test_solvers",
      "ces_solver_args", "start_vals", "stage1_lower", "stage1_upper",
      "stage2_lower", "stage2_upper", "share_lower", "share_upper",
      "run_label", "input_file", "previous_stage1_file", "stage1_file_for_stage2",
      "sigma_grid_KL", "sigma_grid_VAE", "rho_grid_KL", "rho_grid_VAE",
      "stage1_expand_points", "stage1_sigma_min", "stage1_sigma_max",
      "stage1_max_solvers", "stage2_grid_points", "stage2_window_padding",
      "stage2_max_solvers", "stage2_sigma_min", "stage2_sigma_max",
      "stage2_boundary_sigma_min", "stage2_boundary_sigma_max",
      "stage2_solver_pool", "aicc_support_delta", "limits",
      "estimate_region", "estimate_region_or_checkpoint", "extract_region",
      "region_checkpoint_folder", "output_tag", "safe_filename",
      "safe_save_rds", "region_fit_checkpoint_file", "resume_regions",
      "save_region_checkpoint",
      "file_fingerprint", "checkpoint_run_key", "checkpoint_reuse_status",
      "quiet_worker_packages", "%||%"
    ), # exporting objects to parallel processes
    seed = TRUE
  )
)

# Extract grid results from all region-method combinations into a single table
results_grid <- map2_dfr(region_names, fits_all, extract_region)
if (nrow(results_grid) == 0L) {
  stop("No grid results were produced. Check worker errors above before saving this run.")
}

#### MASTER TABLE ####
} else {
  stage_rds_file <- if (!is.na(saved_results) && nzchar(saved_results)) {
    saved_results
  } else {
    main_results_file
  }
  if (file.exists(stage_rds_file)) {
    message("Loading results from ", stage_rds_file)
    saved_stage_results <- readRDS(stage_rds_file)
    if (!is.list(saved_stage_results) || !"grid_results" %in% names(saved_stage_results)) {
      stop("Saved results must be a current-format RDS with a grid_results table.")
    }
    results_grid <- saved_stage_results$grid_results
    if (nrow(results_grid) == 0L) {
      stop("Saved results contain an empty grid_results table.")
    }
    if (is.list(saved_stage_results) && "stage2_plan" %in% names(saved_stage_results)) {
      stage2_plan <- saved_stage_results$stage2_plan
    }
  } else {
    stop("No saved results found. Run with estimate = TRUE first.")
  }
  
  # Load input data so df is available for IAM table
  df <- read_csv(input_file, show_col_types = FALSE)
  input_data_summary <- build_input_data_summary(df, base_year)
  safe_write_csv(input_data_summary, input_data_summary_file)

  normalised_input_data <- build_normalised_input_data(
    raw_input_data = df,
    input_data_summary = input_data_summary,
    base_year = base_year,
    estimate_flagged_regions = estimate_flagged_regions
  )
  dfS <- normalised_input_data # short name used by the analysis code below
}

#### BEST RESULTS AND OUTPUTS ####
# This section turns the full grid of region/solver/rho candidates into
# interpretable outputs. The first tables preserve diagnostics for review; the
# final Stage 2 block builds the compact MERGE parameter table.
if (
  identical(estimation_stage, "stage1") &&
  is.null(previous_stage1_results) &&
  !is.na(previous_stage1_file) &&
  nzchar(previous_stage1_file) &&
  file.exists(previous_stage1_file)
) {
  previous_stage1_results <- readRDS(previous_stage1_file)
}

if (identical(estimation_stage, "stage1") && !is.null(previous_stage1_results)) {
  combined_stage1_grid <- combine_stage1_iteration_grid(
    current_grid = results_grid,
    previous_stage1_results = previous_stage1_results,
    current_run_label = run_label
  )
  results_grid <- combined_stage1_grid$grid_results
  carried_forward_stage1_regions <- combined_stage1_grid$carried_forward_regions
  if (length(carried_forward_stage1_regions) > 0L) {
    message(
      "Carried forward ", length(carried_forward_stage1_regions),
      " region(s) from the previous Stage 1 run: ",
      paste(carried_forward_stage1_regions, collapse = ", ")
    )
  }
}

obs_rho_KL <- results_grid %>%
  transmute(rho_KL = suppressWarnings(as.numeric(rho_KL))) %>%
  filter(is.finite(rho_KL)) %>%
  distinct(rho_KL) %>%
  arrange(rho_KL) %>%
  pull(rho_KL)

obs_rho_VAE <- results_grid %>%
  transmute(rho_VAE = suppressWarnings(as.numeric(rho_VAE))) %>%
  filter(is.finite(rho_VAE)) %>%
  distinct(rho_VAE) %>%
  arrange(rho_VAE) %>%
  pull(rho_VAE)

obs_rho_penalty <- as.integer(length(obs_rho_KL) > 1) +
  as.integer(length(obs_rho_VAE) > 1)

obs_sigma_KL <- sort(unique(round(1 / (1 + obs_rho_KL), 4)))
obs_sigma_VAE <- sort(unique(round(1 / (1 + obs_rho_VAE), 4)))

results_grid <- prepare_grid_results_for_postprocessing(
  grid_results_raw = results_grid,
  normalised_input_data = normalised_input_data,
  obs_rho_KL = obs_rho_KL,
  obs_rho_VAE = obs_rho_VAE,
  obs_rho_penalty = obs_rho_penalty
)

# Convergence summary by region
convergence_summary <- results_grid %>%
  select(r, method, conv) %>%
  distinct() %>%
  count(method, conv, name = "count") %>%
  mutate(status = ifelse(conv, "Converged", "Failed"))

# Function that adds validity to the runs, economically feasible and converged from the solver
add_validity <- function(df) {
  soft <- limits$soft
  final <- limits$final
  strict <- limits$strict

  df %>%
    mutate(
      across(
        any_of(c("delta_KVA","delta_VAY","gamma","nu","lambda")),
        ~ suppressWarnings(as.numeric(.)) # Force to be numeric
      ),
      
      # Solver level convergence
      converged = replace_na(suppressWarnings(as.logical(conv)), FALSE),
      
      # Hard numerical and CES-domain checks. These follow the fitted CES form:
      # positive scale, positive returns curvature, within-nest shares in [0,1],
      # finite substitution elasticities, and finite log-residual RSS.
      valid_numeric =
        is.finite(rss) & rss >= 0 &
        is.finite(AICc_plusRho) &
        is.finite(gamma) & gamma > 0 &
        is.finite(nu) & nu > 0 &
        is.finite(lambda) &
        is.finite(delta_KVA) & between(delta_KVA, 0, 1) &
        is.finite(delta_VAY) & between(delta_VAY, 0, 1) &
        is.finite(sigma_KL) & sigma_KL > 0 &
        is.finite(sigma_VAE) & sigma_VAE > 0,

      # Fit-quality diagnostics. R2 can be negative for a poor nonlinear fit; that is poor fit,
      # not a mathematical impossibility. Keep the flag separate from core CES validity.
      valid_stat = valid_numeric & is.finite(R2) & R2 <= 1,
      positive_R2_fit = valid_stat & R2 > 0,

      # Soft economic plausibility screen. These thresholds are research judgement, not
      # micEconCES package constraints, so they should guide ranking/review rather than define
      # whether a grid cell is mathematically admissible.
      plausible_econ =
        is.finite(gamma) & between(gamma, soft$gamma[[1]], soft$gamma[[2]]) &
        is.finite(nu) & between(nu, soft$nu[[1]], soft$nu[[2]]) &
        is.finite(lambda) & between(lambda, soft$lambda[[1]], soft$lambda[[2]]),

      plausible_final =
        plausible_econ &
        between(delta_KVA, final$delta_KVA[[1]], final$delta_KVA[[2]]) &
        between(delta_VAY, final$delta_VAY[[1]], final$delta_VAY[[2]]) &
        between(sigma_KL, final$sigma_KL[[1]], final$sigma_KL[[2]]) &
        between(sigma_VAE, final$sigma_VAE[[1]], final$sigma_VAE[[2]]) &
        between(nu, final$nu[[1]], final$nu[[2]]) &
        between(lambda, final$lambda[[1]], final$lambda[[2]]),

      valid_econ = valid_numeric & plausible_econ,
      
      # Core Stage 1 validity: enough to use a grid cell in the broad audit.
      # Plausibility and R2 are kept as separate diagnostics.
      valid = converged & valid_stat,
      
      # Strict support used for final ranking and Stage 2 preference, not for censoring Stage 1.
      valid_strict =
        valid &
        positive_R2_fit &
        plausible_final &
        between(delta_KVA, strict$delta_KVA[[1]], strict$delta_KVA[[2]]) &
        between(delta_VAY, strict$delta_VAY[[1]], strict$delta_VAY[[2]]) &
        !on_edge_KL &
        !on_edge_VAE,
      
      # high-level “which filter failed”
      solver_issue = !converged,
      stat_issue = converged & !valid_stat,
      econ_issue = converged & valid_stat & !plausible_econ,
      
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
        list(converged, rss, AICc_plusRho, R2, delta_KVA, delta_VAY, gamma, nu, lambda, sigma_KL, sigma_VAE),
        function(conv, rss, AICc_plusRho, R2, delta_KVA, delta_VAY, gamma, nu, lambda, sigma_KL, sigma_VAE) {
          reasons <- c()
          
          if (!isTRUE(conv)) reasons <- c(reasons, "Solver did not converge")
          if (!is.finite(rss) || rss < 0) reasons <- c(reasons, "log-RSS invalid")
          if (!is.finite(AICc_plusRho)) reasons <- c(reasons, "AICc unavailable")
          if (!is.finite(R2) || R2 > 1) reasons <- c(reasons, "R2 invalid")
          if (!is.finite(delta_KVA) || delta_KVA < 0 || delta_KVA > 1) reasons <- c(reasons, "dK-VA out of [0,1]")
          if (!is.finite(delta_VAY) || delta_VAY < 0 || delta_VAY > 1) reasons <- c(reasons, "dVA-Y out of [0,1]")
          if (!is.finite(gamma) || gamma <= 0) reasons <- c(reasons, "gamma not positive")
          if (!is.finite(nu) || nu <= 0) reasons <- c(reasons, "nu not positive")
          if (!is.finite(lambda)) reasons <- c(reasons, "lambda not finite")
          if (!is.finite(sigma_KL) || sigma_KL <= 0) reasons <- c(reasons, "sigma K-L not positive")
          if (!is.finite(sigma_VAE) || sigma_VAE <= 0) reasons <- c(reasons, "sigma VA-E not positive")
          if (length(reasons) == 0 && (!between(gamma, soft$gamma[[1]], soft$gamma[[2]]) ||
              !between(nu, soft$nu[[1]], soft$nu[[2]]) ||
              !between(lambda, soft$lambda[[1]], soft$lambda[[2]]))) {
            reasons <- c(reasons, "outside soft plausibility band")
          }
          
          if (length(reasons) == 0) "OK" else paste(reasons, collapse = "; ") 
        }
      ),
      
      # Status combines solver and validity info into one
      status = case_when(
        valid ~ "Valid",
        # solver failed, why?
        !converged & solver_reason != "Valid" ~ solver_reason,
        # solver ok, statistical failure, why?
        grepl("log-RSS invalid", valid_reason) ~ "RSS invalid",
        grepl("AICc unavailable", valid_reason) ~ "AICc unavailable",
        grepl("R2 invalid", valid_reason) ~ "R2 invalid",
        # solver & statistics ok, economic failure, why?
        grepl("dK-VA", valid_reason) ~ "dK-VA out of [0,1]",
        grepl("dVA-Y", valid_reason) ~ "dVA-Y out of [0,1]",
        grepl("gamma not positive", valid_reason) ~ "gamma not positive",
        grepl("nu not positive", valid_reason) ~ "nu not positive",
        grepl("lambda not finite", valid_reason) ~ "lambda not finite",
        grepl("sigma", valid_reason) ~ "sigma not positive",
        grepl("soft plausibility", valid_reason) ~ "Outside plausibility band",
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
          "AICc unavailable",
          "R2 invalid",
          # economic failures
          "dK-VA out of [0,1]",
          "dVA-Y out of [0,1]", 
          "gamma not positive",
          "nu not positive",
          "lambda not finite",
          "sigma not positive",
          "Outside plausibility band",
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
    dAICc = AICc_plusRho - finite_min(AICc_plusRho),
    wAICc_raw = exp(-0.5 * dAICc)
  ) %>%
  ungroup()


# Model-average best methods
dAICcmin <- aicc_support_delta

best_methods_average <- results_grid_valid %>%
  group_by(r) %>%
  mutate(
    any_valid_strict = any(valid_strict, na.rm = TRUE),
    any_strict_support = any(valid_strict & dAICc <= dAICcmin, na.rm = TRUE),
    any_support = any(dAICc <= dAICcmin, na.rm = TRUE),
    support_mode = case_when(
      any_support & any_strict_support ~ "strong AICc support", # dAICc <= 4 AND at least one strict
      any_support & !any_strict_support ~ "some AICc support", # dAICc <= 4 but no strict inside support
      !any_support ~ "AICc min support", # no dAICc <= 4 -> use the minimum AICc cell
      TRUE ~ "AICc min support" # safety catch
    )
  ) %>%
  # 1) Keep only rows belonging to the chosen support set for that region
  filter(
      (support_mode %in% c("strong AICc support", "some AICc support") & dAICc <= dAICcmin) |
      (support_mode == "AICc min support" & dAICc == finite_min(dAICc))
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
  filter(dAICc == finite_min(dAICc)) %>%
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
    min_RSS = finite_min(rss),
    max_RSS = finite_max(rss),
    med_RSS = finite_median(rss),
    min_R2 = finite_min(R2),
    max_R2 = finite_max(R2),
    med_R2 = finite_median(R2),
    # Parameters (medians across runs)
    gamma_med = finite_median(gamma),
    lambda_med = finite_median(lambda),
    delta_KVA_med = finite_median(delta_KVA),
    delta_VAY_med = finite_median(delta_VAY),
    nu_med = finite_median(nu),
    sigma_KL_med = finite_median(sigma_KL),
    sigma_VAE_med = finite_median(sigma_VAE),
    # Iterations and runtime
    med_iter = finite_median(iter),
    med_runtime = finite_median(runtime_total),
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
    # micEconCES defines the fitted trend as gamma * exp(lambda * tName). The estimation uses
    # time_from_base_year = year - base_year, so gamma is the base-year scale and the index below
    # is aligned with the fitted time convention.
    TFP_index = exp(lambda * (t - base_year)),
    delta_KVA,
    delta_VAY,
    sigma_KL,
    sigma_VAE,
    rho_KL,
    rho_VAE,
    gamma,
    lambda,
    nu,
    p_gamma, p_lambda, p_delta_KVA, p_delta_VAY, p_nu
  )

# ---- Stage 1 audit and Stage 2 design tables ----
grid_results <- results_grid
valid_grid_results <- results_grid_valid
best_method_by_region <- best_methods
iam_parameter_table <- iam_table
estimated_regions <- dfS %>%
  distinct(r) %>%
  pull(r)

method_runtime <- grid_results %>%
  distinct(r, method, runtime_total)

optimizer_diagnostic <- grid_results %>%
  mutate(
    edge_solution = valid & (on_edge_KL | on_edge_VAE),
    gamma_at_optimizer_bound = is.finite(gamma) &
      (
        near(gamma, stage1_lower[["gamma"]], tol = 1e-6) |
          near(gamma, stage1_upper[["gamma"]], tol = 1e-6)
      )
  ) %>%
  group_by(method) %>%
  summarise(
    n_grid_cells = n(),
    n_region_method_pairs = n_distinct(paste(r, method, sep = " | ")),
    convergence_share = mean(converged, na.rm = TRUE),
    valid_share = mean(valid, na.rm = TRUE),
    valid_non_edge_share = mean(valid & !on_edge_KL & !on_edge_VAE, na.rm = TRUE),
    edge_solution_share = if_else(sum(valid, na.rm = TRUE) > 0,
                                  sum(edge_solution, na.rm = TRUE) / sum(valid, na.rm = TRUE),
                                  NA_real_),
    median_valid_RSS = finite_median(rss[valid]),
    median_valid_AICc = finite_median(AICc_plusRho[valid]),
    regions_with_valid_results = n_distinct(r[valid]),
    gamma_optimizer_bound_hits = sum(gamma_at_optimizer_bound, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    method_runtime %>%
      group_by(method) %>%
      summarise(
        median_runtime_seconds = finite_median(runtime_total),
        max_runtime_seconds = finite_max(runtime_total),
        .groups = "drop"
      ),
    by = "method"
  ) %>%
  arrange(desc(convergence_share), desc(valid_share), edge_solution_share, median_valid_AICc, median_runtime_seconds)

method_choice <- grid_results %>%
  mutate(edge_solution = valid & (on_edge_KL | on_edge_VAE)) %>%
  group_by(r, method) %>%
  summarise(
    n_grid_cells = n(),
    convergence_share = mean(converged, na.rm = TRUE),
    valid_share = mean(valid, na.rm = TRUE),
    valid_non_edge_share = mean(valid & !on_edge_KL & !on_edge_VAE, na.rm = TRUE),
    edge_solution_share = if_else(sum(valid, na.rm = TRUE) > 0,
                                  sum(edge_solution, na.rm = TRUE) / sum(valid, na.rm = TRUE),
                                  NA_real_),
    median_valid_RSS = finite_median(rss[valid]),
    median_valid_AICc = finite_median(AICc_plusRho[valid]),
    median_runtime_seconds = finite_median(runtime_total),
    .groups = "drop"
  ) %>%
  group_by(r) %>%
  arrange(
    desc(convergence_share),
    desc(valid_share),
    replace_na(edge_solution_share, 1),
    median_valid_AICc,
    median_runtime_seconds,
    .by_group = TRUE
  ) %>%
  mutate(
    method_rank = row_number(),
    optimizer_recommendation = case_when(
      valid_share == 0 ~ "do_not_use_for_stage2_without_review",
      edge_solution_share >= 0.75 ~ "review_before_stage2_edge_heavy",
      method_rank <= 3 ~ "candidate_for_stage2",
      TRUE ~ "audit_only"
    )
  ) %>%
  ungroup()

gamma_time_origin_audit <- grid_results %>%
  mutate(
    gamma_base = gamma * exp(lambda * base_year),
    raw_gamma_valid = is.finite(gamma) & between(gamma, limits$soft$gamma[[1]], limits$soft$gamma[[2]]),
    base_year_gamma_valid = is.finite(gamma_base) & between(gamma_base, limits$soft$gamma[[1]], limits$soft$gamma[[2]]),
    gamma_validity_changed = raw_gamma_valid != base_year_gamma_valid,
    gamma_at_optimizer_bound = is.finite(gamma) &
      (
        near(gamma, stage1_lower[["gamma"]], tol = 1e-6) |
          near(gamma, stage1_upper[["gamma"]], tol = 1e-6)
      )
  ) %>%
  group_by(r, method) %>%
  summarise(
    n_grid_cells = n(),
    raw_gamma_valid_count = sum(raw_gamma_valid, na.rm = TRUE),
    base_year_gamma_valid_count = sum(base_year_gamma_valid, na.rm = TRUE),
    gamma_validity_changed_count = sum(gamma_validity_changed, na.rm = TRUE),
    gamma_validity_changed_share = gamma_validity_changed_count / n_grid_cells,
    gamma_optimizer_bound_hits = sum(gamma_at_optimizer_bound, na.rm = TRUE),
    .groups = "drop"
  )

stage1_support_set <- valid_grid_results %>%
  group_by(r) %>%
  mutate(
    region_best_AICc = finite_min(AICc_plusRho),
    stage1_best_on_edge_KL = any(near(AICc_plusRho, region_best_AICc, tol = 1e-8) & on_edge_KL, na.rm = TRUE),
    stage1_best_on_edge_VAE = any(near(AICc_plusRho, region_best_AICc, tol = 1e-8) & on_edge_VAE, na.rm = TRUE),
    stage1_best_on_edge = stage1_best_on_edge_KL | stage1_best_on_edge_VAE
  ) %>%
  ungroup()

stage1_supported_non_edge <- stage1_support_set %>%
  filter(dAICc <= dAICcmin, !on_edge_KL, !on_edge_VAE)

stage1_region_grid_limits <- grid_results %>%
  group_by(r) %>%
  summarise(
    sigma_KL_grid_min = finite_min(sigma_KL),
    sigma_KL_grid_max = finite_max(sigma_KL),
    sigma_VAE_grid_min = finite_min(sigma_VAE),
    sigma_VAE_grid_max = finite_max(sigma_VAE),
    .groups = "drop"
  )

rho_windows <- input_data_summary %>%
  select(r, input_data_class, input_data_note) %>%
  left_join(stage1_region_grid_limits, by = "r") %>%
  left_join(
    stage1_support_set %>%
      group_by(r) %>%
      summarise(
        n_valid_stage1_grid_cells = n(),
        stage1_best_on_edge = first(stage1_best_on_edge),
        stage1_best_on_edge_KL = first(stage1_best_on_edge_KL),
        stage1_best_on_edge_VAE = first(stage1_best_on_edge_VAE),
        best_stage1_AICc = finite_min(AICc_plusRho),
        best_stage1_rho_KL = finite_first_at_min(rho_KL, AICc_plusRho),
        best_stage1_rho_VAE = finite_first_at_min(rho_VAE, AICc_plusRho),
        best_stage1_sigma_KL = finite_first_at_min(sigma_KL, AICc_plusRho),
        best_stage1_sigma_VAE = finite_first_at_min(sigma_VAE, AICc_plusRho),
        .groups = "drop"
      ),
    by = "r"
  ) %>%
  left_join(
    stage1_supported_non_edge %>%
      group_by(r) %>%
      summarise(
        supported_non_edge_grid_cells = n(),
        rho_KL_supported_min = finite_min(rho_KL),
        rho_KL_supported_max = finite_max(rho_KL),
        rho_VAE_supported_min = finite_min(rho_VAE),
        rho_VAE_supported_max = finite_max(rho_VAE),
        sigma_KL_supported_min = finite_min(sigma_KL),
        sigma_KL_supported_max = finite_max(sigma_KL),
        sigma_VAE_supported_min = finite_min(sigma_VAE),
        sigma_VAE_supported_max = finite_max(sigma_VAE),
        supporting_methods = paste(sort(unique(method)), collapse = ", "),
        .groups = "drop"
      ),
    by = "r"
  ) %>%
  mutate(
    stage1_can_expand_KL = replace_na(stage1_best_on_edge_KL, FALSE) &
      (
        (near(best_stage1_sigma_KL, sigma_KL_grid_min, tol = 1e-8) & stage1_sigma_min < sigma_KL_grid_min) |
          (near(best_stage1_sigma_KL, sigma_KL_grid_max, tol = 1e-8) & stage1_sigma_max > sigma_KL_grid_max)
      ),
    stage1_can_expand_VAE = replace_na(stage1_best_on_edge_VAE, FALSE) &
      (
        (near(best_stage1_sigma_VAE, sigma_VAE_grid_min, tol = 1e-8) & stage1_sigma_min < sigma_VAE_grid_min) |
          (near(best_stage1_sigma_VAE, sigma_VAE_grid_max, tol = 1e-8) & stage1_sigma_max > sigma_VAE_grid_max)
      ),
    stage1_can_expand_grid = replace_na(stage1_can_expand_KL | stage1_can_expand_VAE, FALSE),
    rho_refinement_action = case_when(
      isTRUE(test_mode) & !r %in% estimated_regions ~
        "not_estimated_in_test_run",
      input_data_class == "exclude_or_repair_before_stage1" & !isTRUE(estimate_flagged_regions) ~
        "repair_input_data_before_stage1",
      is.na(n_valid_stage1_grid_cells) | n_valid_stage1_grid_cells == 0 ~
        "no_valid_stage1_result_review_data_or_optimizer_settings",
      stage1_best_on_edge & stage1_can_expand_grid ~
        "expand_stage1_grid_before_stage2",
      stage1_best_on_edge ~
        "boundary_supported_review_before_stage2",
      is.na(supported_non_edge_grid_cells) | supported_non_edge_grid_cells == 0 ~
        "review_edge_only_or_weak_stage1_support",
      TRUE ~ "ready_for_stage2_refinement"
    )
  )

residuals <- grid_results %>%
  left_join(input_data_summary %>% select(r, n_observations), by = "r") %>%
  mutate(log_residual_rmse = sqrt(rss / n_observations)) %>%
  group_by(r, method) %>%
  summarise(
    n_grid_cells = n(),
    n_valid_grid_cells = sum(valid, na.rm = TRUE),
    median_valid_log_residual_rmse = if (any(valid & is.finite(log_residual_rmse))) {
      finite_median(log_residual_rmse[valid])
    } else {
      NA_real_
    },
    best_valid_log_residual_rmse = if (any(valid & is.finite(log_residual_rmse))) {
      finite_min(log_residual_rmse[valid])
    } else {
      NA_real_
    },
    median_valid_R2_log = finite_median(R2[valid]),
    .groups = "drop"
  )

stability <- stage1_supported_non_edge %>%
  group_by(r) %>%
  summarise(
    support_grid_cells = n(),
    support_methods = paste(sort(unique(method)), collapse = ", "),
    gamma_min = finite_min(gamma),
    gamma_max = finite_max(gamma),
    lambda_min = finite_min(lambda),
    lambda_max = finite_max(lambda),
    delta_KVA_min = finite_min(delta_KVA),
    delta_KVA_max = finite_max(delta_KVA),
    delta_VAY_min = finite_min(delta_VAY),
    delta_VAY_max = finite_max(delta_VAY),
    nu_min = finite_min(nu),
    nu_max = finite_max(nu),
    sigma_KL_min = finite_min(sigma_KL),
    sigma_KL_max = finite_max(sigma_KL),
    sigma_VAE_min = finite_min(sigma_VAE),
    sigma_VAE_max = finite_max(sigma_VAE),
    .groups = "drop"
  )

stage2_optimizer_candidates <- method_choice %>%
  filter(optimizer_recommendation == "candidate_for_stage2") %>%
  group_by(r) %>%
  summarise(
    selected_stage2_solvers = paste(method, collapse = ", "),
    .groups = "drop"
  )

stage2_design <- rho_windows %>%
  left_join(stage2_optimizer_candidates, by = "r") %>%
  mutate(
    selected_stage2_solvers = replace_na(selected_stage2_solvers, ""),
    stage2_design_note = case_when(
      rho_refinement_action == "ready_for_stage2_refinement" ~
        "Use valid non-edge Stage 1 support set with delta AICc <= 4; keep AICc as heuristic.",
      rho_refinement_action == "expand_stage1_grid_before_stage2" ~
        "Best supported result is on a Stage 1 grid edge; expand Stage 1 before Stage 2.",
      rho_refinement_action == "boundary_supported_review_before_stage2" ~
        "Targeted Stage 1 still selects a boundary solution; review as a boundary-supported case before MERGE use.",
      rho_refinement_action == "repair_input_data_before_stage1" ~
        "Input data failed the Stage 1 gate and should not be estimated without explicit override.",
      rho_refinement_action == "not_estimated_in_test_run" ~
        "Region was intentionally skipped by the non-scientific test run.",
      TRUE ~ "Review Stage 1 support before finalising Stage 2."
    )
  )

if (identical(estimation_stage, "stage1")) {
  stage2_plan <- build_stage2_plan(
    stage1_results = list(
      metadata = list(
        created_at = Sys.time(),
        obs_sigma_KL = obs_sigma_KL,
        obs_sigma_VAE = obs_sigma_VAE
      ),
      stage2_design = stage2_design,
      rho_windows = rho_windows,
      method_choice = method_choice
    ),
      grid_points = stage2_grid_points,
      window_padding = stage2_window_padding,
      max_solvers = stage2_max_solvers,
      sigma_lower_bound = stage2_sigma_min,
      sigma_upper_bound = stage2_sigma_max,
      boundary_sigma_lower_bound = stage2_boundary_sigma_min,
      boundary_sigma_upper_bound = stage2_boundary_sigma_max,
      solver_pool = stage2_solver_pool
    )
}

if (identical(estimation_stage, "stage1") && !isTRUE(test_mode)) {
  next_stage1_plan <- build_next_stage1_plan(
    stage1_results = list(
      metadata = list(
        run_label = run_label,
        obs_sigma_KL = obs_sigma_KL,
        obs_sigma_VAE = obs_sigma_VAE
      ),
      stage2_design = stage2_design,
      method_choice = method_choice,
      grid_results = grid_results
    ),
    grid_points = stage1_expand_points,
    sigma_lower_bound = stage1_sigma_min,
    sigma_upper_bound = stage1_sigma_max,
    max_solvers = stage1_max_solvers
  )
}

stage2_final_support_set <- tibble()
stage2_review <- tibble()
stage2_fit <- tibble()
stage2_fit_summary <- tibble()
merge_regional_parameter_table <- tibble()
merge_tfp_time_series <- tibble()
merge_annual_diagnostic_table <- tibble()
merge_long_parameter_table <- tibble()
merge_iam_parameter_table <- tibble()

if (identical(estimation_stage, "stage2")) {
  stage2_final_support_set <- valid_grid_results %>%
    left_join(
      stage2_design %>% select(r, rho_refinement_action),
      by = "r"
    ) %>%
    mutate(non_edge = !on_edge_KL & !on_edge_VAE) %>%
    group_by(r) %>%
    mutate(
      has_strict_supported = any(valid_strict & dAICc <= dAICcmin, na.rm = TRUE),
      has_plausible_non_edge_supported = any(valid & plausible_final & non_edge & dAICc <= dAICcmin, na.rm = TRUE),
      has_plausible_non_edge_any = any(valid & plausible_final & non_edge, na.rm = TRUE),
      best_plausible_AICc = if (any(valid & plausible_final & non_edge, na.rm = TRUE)) {
        finite_min(AICc_plusRho[valid & plausible_final & non_edge])
      } else {
        NA_real_
      },
      has_valid_non_edge_supported = any(valid & non_edge & dAICc <= dAICcmin, na.rm = TRUE),
      has_boundary_supported = any(
        valid & !non_edge & dAICc <= dAICcmin,
        na.rm = TRUE
      ),
      final_support_tier = case_when(
        has_strict_supported & valid_strict & dAICc <= dAICcmin ~ "strict_supported",
        !has_strict_supported & has_plausible_non_edge_supported &
          valid & plausible_final & non_edge & dAICc <= dAICcmin ~
          "plausible_non_edge_supported",
        !has_strict_supported & !has_plausible_non_edge_supported & has_plausible_non_edge_any &
          valid & plausible_final & non_edge &
          AICc_plusRho <= best_plausible_AICc + dAICcmin ~
          "plausible_fit_tradeoff",
        !has_strict_supported & !has_plausible_non_edge_supported & !has_plausible_non_edge_any &
          has_valid_non_edge_supported & valid & non_edge & dAICc <= dAICcmin ~
          "valid_non_edge_supported",
        !has_strict_supported & !has_plausible_non_edge_supported & !has_plausible_non_edge_any &
          !has_valid_non_edge_supported &
          has_boundary_supported & valid & !non_edge & dAICc <= dAICcmin ~
          "boundary_supported",
        TRUE ~ NA_character_
      )
    ) %>%
    ungroup() %>%
    filter(!is.na(final_support_tier)) %>%
    group_by(r) %>%
    mutate(
      final_dAICc = AICc_plusRho - finite_min(AICc_plusRho),
      final_weight_raw = exp(-0.5 * final_dAICc),
      final_AICc_weight = final_weight_raw / sum(final_weight_raw, na.rm = TRUE)
    ) %>%
    ungroup()

  stage2_review <- stage2_final_support_set %>%
    left_join(input_data_summary %>% select(r, n_observations), by = "r") %>%
    group_by(r) %>%
    arrange(AICc_plusRho, .by_group = TRUE) %>%
    summarise(
      final_selection_rule = first(final_support_tier),
      final_support_grid_cells = n(),
      final_support_methods = paste(sort(unique(method)), collapse = ", "),
      final_selected_solver = first(method),
      final_best_AICc = finite_min(AICc_plusRho),
      final_log_residual_RMSE = sqrt(first(rss) / first(n_observations)),
      gamma = first(gamma),
      lambda = first(lambda),
      delta_KVA = first(delta_KVA),
      delta_VAY = first(delta_VAY),
      nu = first(nu),
      rho_KL = first(rho_KL),
      rho_VAE = first(rho_VAE),
      sigma_KL = 1 / (1 + rho_KL),
      sigma_VAE = 1 / (1 + rho_VAE),
      final_any_edge_solution = any(on_edge_KL | on_edge_VAE, na.rm = TRUE),
      final_any_gamma_bound_hit = any(
        is.finite(gamma) &
          (
            near(gamma, stage2_lower[["gamma"]], tol = 1e-6) |
              near(gamma, stage2_upper[["gamma"]], tol = 1e-6)
          ),
        na.rm = TRUE
      ),
      final_review_note = case_when(
        first(final_support_tier) == "strict_supported" & !final_any_edge_solution ~
          "MERGE-ready subject to substantive regional review.",
        first(final_support_tier) == "plausible_non_edge_supported" ~
          "MERGE-ready using plausible non-edge support; document that strict support was not available.",
        first(final_support_tier) == "plausible_fit_tradeoff" ~
          "MERGE-ready as a plausibility-weighted trade-off; best statistical support was economically implausible.",
        first(final_support_tier) == "valid_non_edge_supported" ~
          "Use as fallback for complete MERGE coverage; no plausible support set was available.",
        first(final_support_tier) == "boundary_supported" ~
          "MERGE input is boundary-supported; use because MERGE needs complete coverage, but report the boundary status.",
        TRUE ~ "Manual review required before MERGE use."
      ),
      .groups = "drop"
    )

  stage2_fit <- stage2_review %>%
    select(
      r, final_selection_rule, final_selected_solver,
      gamma, lambda, delta_KVA, delta_VAY, rho_KL, rho_VAE, nu
    ) %>%
    inner_join(
      dfS %>% select(r, t, time_from_base_year, Ys, Ks, Ls, Es, Ybase),
      by = "r"
    ) %>%
    group_by(r) %>%
    group_modify(function(region_data, key) {
      selected <- region_data[1, ]
      ces_coefficients <- c(
        selected$gamma, selected$lambda, selected$delta_KVA, selected$delta_VAY,
        selected$rho_KL, selected$rho_VAE, selected$nu
      )
      predicted_Y_norm <- micEconCES::cesCalc(
        xNames = c("Ks", "Ls", "Es"),
        data = region_data,
        coef = ces_coefficients,
        tName = "time_from_base_year",
        nested = TRUE
      )
      region_data %>%
        mutate(
          predicted_Y_norm = predicted_Y_norm,
          predicted_output = predicted_Y_norm * Ybase,
          log_observed_output = log(Ys),
          log_predicted_output = log(predicted_Y_norm),
          log_residual = log_observed_output - log_predicted_output,
          percent_residual = 100 * (Ys / predicted_Y_norm - 1)
        )
    }) %>%
    ungroup()

  stage2_fit_summary <- stage2_fit %>%
    group_by(r) %>%
    summarise(
      selected_fit_log_RMSE = sqrt(mean(log_residual^2, na.rm = TRUE)),
      selected_fit_log_bias = mean(log_residual, na.rm = TRUE),
      selected_fit_mean_absolute_percent_error = mean(abs(percent_residual), na.rm = TRUE),
      selected_fit_log_R2 = {
        rss_selected <- sum(log_residual^2, na.rm = TRUE)
        tss_selected <- sum((log_observed_output - mean(log_observed_output, na.rm = TRUE))^2, na.rm = TRUE)
        if (is.finite(tss_selected) && tss_selected > 0) 1 - rss_selected / tss_selected else NA_real_
      },
      .groups = "drop"
    )

  stage2_review <- stage2_review %>%
    left_join(stage2_fit_summary, by = "r")

  # Final regional MERGE parameter table. Each row is one MERGE region and the
  # values come from one coherent selected CES fit, not from averaging unrelated
  # nonlinear parameter vectors.
  #
  # Parameter interpretation:
  # - delta_KVA and delta_VAY are CES distribution weights in the K-L and VA-E nests.
  # - rho_KL and rho_VAE are the CES curvature terms used by the MERGE equation.
  # - sigma_KL and sigma_VAE are the substitution elasticities, sigma = 1 / (1 + rho).
  # - gamma is base-year TFP scale because time_from_base_year equals zero in base_year.
  # - lambda is continuous annual TFP growth; 0.01 is roughly 1 percent per year before compounding.
  # - nu is the fitted returns-to-scale exponent.
  merge_regional_parameter_table <- stage2_review %>%
    inner_join(dfS %>% distinct(r, Ybase, Kbase, Lbase, Ebase), by = "r") %>%
    transmute(
      region = r,
      delta_KVA,
      delta_VAY,
      rho_KL,
      rho_VAE,
      nu,
      Ybase,
      Kbase,
      Lbase,
      Ebase,
      sigma_KL,
      sigma_VAE,
      gamma,
      lambda,
      final_selected_solver,
      final_selection_rule,
      final_review_note
    )

  # TFP diagnostic path implied by gamma and lambda. MERGE currently consumes
  # only the regional gamma and lambda fields, but the annual path is retained
  # to make the growth-rate interpretation visible during review.
  merge_tfp_time_series <- stage2_review %>%
    inner_join(dfS %>% distinct(r, t), by = "r") %>%
    transmute(
      region = r,
      year = t,
      TFP = gamma * exp(lambda * (t - base_year)),
      TFP_index = exp(lambda * (t - base_year))
    )

  merge_annual_diagnostic_table <- merge_tfp_time_series %>%
    left_join(merge_regional_parameter_table, by = "region")

  merge_long_parameter_table <- bind_rows(
    merge_regional_parameter_table %>%
      select(
        region,
        delta_KVA, delta_VAY,
        rho_KL, rho_VAE,
        nu, Ybase, Kbase, Lbase, Ebase,
        sigma_KL, sigma_VAE,
        gamma, lambda
      ) %>%
      pivot_longer(
        cols = c(
          delta_KVA, delta_VAY,
          rho_KL, rho_VAE,
          nu, Ybase, Kbase, Lbase, Ebase,
          sigma_KL, sigma_VAE,
          gamma, lambda
        ),
        names_to = "parameter",
        values_to = "value"
      ) %>%
      mutate(year = NA_integer_, dimension = "region"),
    merge_tfp_time_series %>%
      pivot_longer(c(TFP, TFP_index), names_to = "parameter", values_to = "value") %>%
      mutate(dimension = "region-year")
  ) %>%
    left_join(
      merge_regional_parameter_table %>%
        select(region, final_selected_solver, final_selection_rule, final_review_note),
      by = "region"
    ) %>%
    select(parameter, region, year, value, dimension, final_selected_solver, final_selection_rule, final_review_note) %>%
    arrange(parameter, region, year)

  merge_iam_parameter_table <- merge_regional_parameter_table
}

metadata <- list(
  estimation_stage = estimation_stage,
  run_label = run_label,
  base_year = base_year,
  created_at = Sys.time(),
  run_mode = if (isTRUE(estimate)) "estimation_plus_postprocessing" else "postprocess_saved_grid_results",
  estimate = estimate,
  resume_regions = resume_regions,
  reused_saved_results = !isTRUE(estimate),
  reused_results_source_file = if (!isTRUE(estimate)) normalizePath(stage_rds_file, winslash = "/", mustWork = FALSE) else NA_character_,
  previous_metadata = if (is.list(saved_stage_results) && "metadata" %in% names(saved_stage_results)) saved_stage_results$metadata else NULL,
  input_file = input_file,
  input_rows = nrow(df),
  input_regions = n_distinct(df$r),
  input_year_min = min(df$t, na.rm = TRUE),
  input_year_max = max(df$t, na.rm = TRUE),
  sigma_grid_KL = sigma_grid_KL,
  sigma_grid_VAE = sigma_grid_VAE,
  rho_grid_KL = rho_grid_KL,
  rho_grid_VAE = rho_grid_VAE,
  obs_sigma_KL = obs_sigma_KL,
  obs_sigma_VAE = obs_sigma_VAE,
  obs_rho_KL = obs_rho_KL,
  obs_rho_VAE = obs_rho_VAE,
  obs_rho_penalty = obs_rho_penalty,
  stage1_solvers = stage1_solvers,
  stage2_solver_source = "region-specific Stage 1 optimizer diagnostics",
  stage1_file_for_stage2 = stage1_file_for_stage2,
  previous_stage1_file = if (!is.na(previous_stage1_file) && nzchar(previous_stage1_file)) {
    normalizePath(previous_stage1_file, winslash = "/", mustWork = FALSE)
  } else {
    NA_character_
  },
  carried_forward_stage1_regions = carried_forward_stage1_regions,
  carried_forward_stage1_region_count = length(carried_forward_stage1_regions),
  stage1_expand_points = stage1_expand_points,
  stage1_sigma_min = stage1_sigma_min,
  stage1_sigma_max = stage1_sigma_max,
  stage1_max_solvers = stage1_max_solvers,
  run_history_results_file = normalizePath(run_history_results_file, winslash = "/", mustWork = FALSE),
  run_history_index_file = normalizePath(run_history_index_file, winslash = "/", mustWork = FALSE),
  stage2_grid_points = stage2_grid_points,
  stage2_window_padding = stage2_window_padding,
  stage2_max_solvers = stage2_max_solvers,
  stage2_sigma_min = stage2_sigma_min,
  stage2_sigma_max = stage2_sigma_max,
  stage2_boundary_sigma_min = stage2_boundary_sigma_min,
  stage2_boundary_sigma_max = stage2_boundary_sigma_max,
  stage2_solver_pool = stage2_solver_pool,
  aicc_support_delta = aicc_support_delta,
  limits = limits,
  start_vals = start_vals,
  stage1_lower = stage1_lower,
  stage1_upper = stage1_upper,
  stage2_lower = stage2_lower,
  stage2_upper = stage2_upper,
  share_lower = share_lower,
  share_upper = share_upper,
  test_mode = test_mode,
  test = if (isTRUE(test_mode)) test else NA_character_,
  test_n = if (isTRUE(test_mode)) test_n else NA_integer_,
  test_r = if (isTRUE(test_mode)) test_r else character(0),
  estimate_flagged_regions = estimate_flagged_regions,
  aicc_label = "heuristic: RSS-based log-residual information criterion with grid-selected rho penalty",
  micEconCES_version = as.character(packageVersion("micEconCES"))
)

results <- list(
  metadata = metadata,
  input_data_summary = input_data_summary,
  grid_results = grid_results,
  optimizer_diagnostic = optimizer_diagnostic,
  method_choice = method_choice,
  rho_windows = rho_windows,
  residuals = residuals,
  stability = stability,
  next_stage1_plan = next_stage1_plan,
  stage2_design = stage2_design,
  stage2_plan = stage2_plan,
  stage2_final_support_set = stage2_final_support_set,
  stage2_review = stage2_review,
  stage2_fit = stage2_fit,
  stage2_fit_summary = stage2_fit_summary,
  merge_regional_parameter_table = merge_regional_parameter_table,
  merge_tfp_time_series = merge_tfp_time_series,
  merge_annual_diagnostic_table = merge_annual_diagnostic_table,
  merge_long_parameter_table = merge_long_parameter_table,
  merge_iam_parameter_table = merge_iam_parameter_table,
  gamma_time_origin_audit = gamma_time_origin_audit,
  valid_grid_results = valid_grid_results,
  best_method_by_region = best_method_by_region,
  best_average_by_region = best_methods_average,
  iam_parameter_table = iam_parameter_table,
  convergence_summary = convergence_summary,
  robustness_summary = robustness_summary
)

written_main_results_file <- safe_save_rds(results, main_results_file)
if (is.na(written_main_results_file)) {
  stop(
    "The run finished but the main results RDS could not be written. ",
    "Close any program using the output folder and rerun this stage."
  )
}

if (isTRUE(save_run_history) && !isTRUE(test_mode)) {
  written_run_history_results_file <- write_run_history(
    results_object = results,
    main_file = written_main_results_file,
    history_file = run_history_results_file,
    index_file = run_history_index_file
  )
  message("Run history written to: ", normalizePath(written_run_history_results_file, winslash = "/"))
}

if (identical(estimation_stage, "stage2")) {
  safe_write_csv(merge_iam_parameter_table, output_file("merge_iam_parameter_table", ".csv"))
}
