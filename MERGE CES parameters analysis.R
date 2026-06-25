options(scipen = 999)

#### PURPOSE ####
# Analyse a completed MERGE-ETL nested-CES run.
# This script does not estimate anything. It reads the main RDS object from
# "MERGE CES parameters.R", updates compact review tables, and draws paper
# figures from the saved estimation results.
#
# The analysis outputs are meant for human review: they explain whether Stage 1
# found interior support, which solvers behaved well, how Stage 2 selected final
# parameter vectors, and what the final MERGE table contains.

#### PACKAGES ####
packages <- c("dplyr", "tidyr", "readr", "purrr", "tibble", "ggplot2", "scales", "stringr")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Install required packages before running this analysis: ", paste(missing_packages, collapse = ", "))
}
invisible(lapply(packages, function(pkg) {
  suppressPackageStartupMessages(suppressWarnings(library(pkg, character.only = TRUE)))
}))

#### SETTINGS ####
setting <- function(name, default) {
  if (exists(name, envir = .GlobalEnv, inherits = FALSE)) get(name, envir = .GlobalEnv) else default
}

data_folder <- setting("data_folder", setting("data_root_folder", "C:/Users/escami_g/Desktop/CES data"))
analysis_stage <- setting("analysis_stage", "stage2") # "stage1", "stage2", or "test"
write_csv_tables <- setting("write_csv_tables", TRUE)
write_full_plan_tables <- setting("write_full_plan_tables", FALSE)
write_individual_figures <- setting("write_individual_figures", TRUE)
show_plots <- setting("show_plots", interactive())
save_history <- setting("save_history", TRUE)

args <- commandArgs(trailingOnly = TRUE)
if ("--test" %in% args) analysis_stage <- "test"
if ("--stage2" %in% args) analysis_stage <- "stage2"
if (!analysis_stage %in% c("stage1", "stage2", "test")) {
  stop('analysis_stage must be "stage1", "stage2", or "test".')
}

stage_folder <- if (analysis_stage == "stage2") "stage2" else "stage1"
stage_path <- file.path(data_folder, stage_folder)
results_file <- file.path(stage_path, paste0(analysis_stage, "_results.rds"))
analysis_file <- file.path(stage_path, paste0(analysis_stage, "_analysis_results.rds"))
report_file <- file.path(stage_path, paste0(analysis_stage, "_scientific_report.md"))

#### HELPERS ####
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || isTRUE(all(is.na(x)))) y else x
}

safe_write_csv <- function(x, file) {
  tryCatch(
    readr::write_csv(x, file),
    error = function(e) warning("Could not write CSV '", file, "'. Close it if it is open. ", conditionMessage(e), call. = FALSE)
  )
}

safe_save_rds <- function(x, file) {
  tryCatch(
    saveRDS(x, file),
    error = function(e) {
      dated_copy <- file.path(dirname(file), paste0(tools::file_path_sans_ext(basename(file)), "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds"))
      saveRDS(x, dated_copy)
      warning("Could not write RDS '", file, "'. Wrote dated copy: ", dated_copy, ". ", conditionMessage(e), call. = FALSE)
    }
  )
}

as_yes_no <- function(x) if_else(replace_na(as.logical(x), FALSE), "Yes", "No")

font_family <- function() {
  windows_fonts <- tryCatch(names(grDevices::windowsFonts()), error = function(e) character())
  if ("Aptos" %in% windows_fonts) "Aptos" else "sans"
}

paper_theme <- function(base_size = 10.5) {
  theme_minimal(base_size = base_size, base_family = font_family()) +
    theme(
      text = element_text(colour = "#1B1B1B"),
      plot.title = element_text(face = "bold", size = base_size + 3, margin = margin(b = 4)),
      plot.subtitle = element_text(size = base_size, colour = "#3A3A3A", margin = margin(b = 10)),
      plot.caption = element_text(size = base_size - 2, colour = "#555555", hjust = 0, margin = margin(t = 8)),
      panel.grid.major = element_line(linewidth = 0.22, colour = "#E6E6E6"),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "#333333"),
      strip.text = element_text(face = "bold", hjust = 0),
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(9, 12, 9, 12)
    )
}

colour <- list(
  ready = "#26734D",
  expand = "#D66A00",
  review = "#B3222A",
  blue = "#1F5A93",
  teal = "#178F8F",
  purple = "#6F4EA2",
  grey = "#656565",
  pale = "#F2F2F2",
  line = "#B8B8B8"
)

stage_colours <- c(
  ready_for_stage2_refinement = colour$ready,
  expand_stage1_grid_before_stage2 = colour$expand,
  boundary_supported_review_before_stage2 = colour$review,
  review_edge_only_or_weak_stage1_support = colour$review,
  no_valid_stage1_result_review_data_or_optimizer_settings = colour$review,
  repair_input_data_before_stage1 = colour$review,
  not_estimated_in_test_run = colour$grey
)

selection_colours <- c(
  strict_supported = colour$ready,
  plausible_non_edge_supported = colour$teal,
  plausible_fit_tradeoff = colour$purple,
  valid_non_edge_supported = colour$expand,
  boundary_supported = colour$review,
  `strict-valid` = colour$ready,
  valid = colour$expand,
  `conv-only` = colour$review,
  any = colour$grey
)

action_label <- function(x) {
  recode(
    x,
    ready_for_stage2_refinement = "Ready for Stage 2",
    expand_stage1_grid_before_stage2 = "Expand Stage 1 grid",
    boundary_supported_review_before_stage2 = "Boundary-supported review",
    review_edge_only_or_weak_stage1_support = "Review weak support",
    no_valid_stage1_result_review_data_or_optimizer_settings = "Review optimizer/data",
    repair_input_data_before_stage1 = "Repair input data",
    not_estimated_in_test_run = "Skipped in test",
    .default = str_to_sentence(str_replace_all(x, "_", " "))
  )
}

selection_label <- function(x) {
  recode(
    x,
    strict_supported = "Strict supported estimate",
    plausible_non_edge_supported = "Plausible non-edge supported estimate",
    plausible_fit_tradeoff = "Plausible fit trade-off estimate",
    valid_non_edge_supported = "Valid non-edge supported estimate",
    boundary_supported = "Boundary-supported estimate",
    `strict-valid` = "Strict valid estimate",
    valid = "Valid estimate",
    `conv-only` = "Converged estimate only",
    any = "Best available estimate",
    .default = str_to_sentence(str_replace_all(x %||% "", "_", " "))
  )
}

safe_write_lines <- function(x, file) {
  tryCatch(
    writeLines(x, file),
    error = function(e) warning("Could not write report '", file, "'. Close it if it is open. ", conditionMessage(e), call. = FALSE)
  )
}

grid_text <- function(x) {
  values <- suppressWarnings(as.numeric(unlist(strsplit(x %||% "", ";", fixed = TRUE))))
  values <- values[is.finite(values)]
  if (length(values) == 0L) return("")
  paste0(length(values), " values: ", number(min(values), accuracy = 0.001), " to ", number(max(values), accuracy = 0.001))
}

collapse_grid <- function(x) paste(round(x, 5), collapse = "; ")

clean_plan <- function(plan, stage) {
  if (is.null(plan) || !"region_plan" %in% names(plan)) return(tibble())
  df <- plan$region_plan
  if (stage == "stage1") {
    df %>%
      mutate(
        across(c(sigma_grid_KL, sigma_grid_VAE, rho_grid_KL, rho_grid_VAE), \(x) map_chr(x, collapse_grid)),
        selected_stage1_solvers = map_chr(selected_stage1_solvers_list, paste, collapse = ", ")
      ) %>%
      select(-selected_stage1_solvers_list)
  } else {
    df %>%
      mutate(
        across(c(sigma_grid_KL, sigma_grid_VAE, rho_grid_KL, rho_grid_VAE), \(x) map_chr(x, collapse_grid)),
        selected_stage2_solvers = map_chr(selected_stage2_solvers_list, paste, collapse = ", ")
      ) %>%
      select(-selected_stage2_solvers_list)
  }
}

write_named_tables <- function(tables, folder, prefix) {
  if (!write_csv_tables) return(invisible(NULL))
  tables <- keep(tables, ~ nrow(.x) > 0L)
  iwalk(tables, ~ safe_write_csv(.x, file.path(folder, paste0(prefix, "_", .y, ".csv"))))
}

save_figure_book <- function(plots, file) {
  grDevices::pdf(file, width = 9.4, height = 7.0, onefile = TRUE, useDingbats = FALSE)
  walk(plots, print)
  grDevices::dev.off()
}

save_individual_figures <- function(plots, folder) {
  if (!write_individual_figures) return(invisible(NULL))
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  iwalk(plots, function(plot, name) {
    ggsave(file.path(folder, paste0(name, ".pdf")), plot, width = 9.4, height = 7.0, device = cairo_pdf)
    ggsave(file.path(folder, paste0(name, ".png")), plot, width = 9.4, height = 7.0, dpi = 400, bg = "white")
  })
}

#### LOAD ####
if (!file.exists(results_file)) stop("Results file not found: ", results_file)
ces <- readRDS(results_file)

required_objects <- c(
  "metadata", "input_data_summary", "grid_results", "optimizer_diagnostic",
  "method_choice", "rho_windows", "residuals",
  "stability", "stage2_design",
  "valid_grid_results", "best_method_by_region", "best_average_by_region",
  "iam_parameter_table"
)
missing_objects <- setdiff(required_objects, names(ces))
if (length(missing_objects) > 0L) stop("Results object is missing: ", paste(missing_objects, collapse = ", "))
if (!"stage2_plan" %in% names(ces)) {
  stop("Results object is missing: stage2_plan")
}

metadata <- ces$metadata
run_label <- metadata$run_label %||% analysis_stage
run_file_label <- gsub("[^A-Za-z0-9.]+", "_", run_label)
stage_tag <- if (analysis_stage == "test") "Test" else if ((metadata$estimation_stage %||% analysis_stage) == "stage2") "Stage2" else "Stage1"
figure_prefix <- paste(stage_tag, run_file_label, sep = ".")
figure_book_file <- file.path(stage_path, paste0(figure_prefix, ".figures.pdf"))
figure_folder <- file.path(stage_path, paste0(figure_prefix, ".figures"))
history_folder <- file.path(stage_path, "history", run_file_label)
history_analysis_file <- file.path(history_folder, paste0(run_file_label, "_analysis_results.rds"))
history_figure_file <- file.path(history_folder, paste0(figure_prefix, ".figures.pdf"))
history_report_file <- file.path(history_folder, paste0(run_file_label, "_scientific_report.md"))

input_data <- ces$input_data_summary
grid <- ces$grid_results
valid_grid <- ces$valid_grid_results
optimizer <- ces$optimizer_diagnostic
windows <- ces$rho_windows
residuals <- ces$residuals
stage2_design <- ces$stage2_design
stage1_plan_object <- ces$next_stage1_plan %||% NULL
stage2_plan_object <- ces$stage2_plan
best_region <- ces$best_method_by_region
stage2_final <- ces$stage2_review %||% tibble()
stage2_fit_summary <- ces$stage2_fit_summary %||% tibble()
merge_table <- ces$merge_iam_parameter_table %||% tibble()

#### TABLES ####
run_summary <- tibble(
  analysis_stage,
  estimation_stage = metadata$estimation_stage %||% NA_character_,
  run_label,
  run_mode = metadata$run_mode %||% NA_character_,
  created_at = as.character(metadata$created_at %||% NA),
  input_regions = metadata$input_regions %||% n_distinct(input_data$r),
  grid_rows = nrow(grid),
  valid_grid_rows = nrow(valid_grid),
  valid_grid_share = if_else(nrow(grid) > 0, nrow(valid_grid) / nrow(grid), NA_real_),
  stage2_ready_regions = length(stage2_plan_object$regions_to_estimate %||% character()),
  next_stage1_regions = length(stage1_plan_object$regions_to_estimate %||% character()),
  stage1_expansion_regions = length(stage2_plan_object$regions_requiring_stage1_expansion %||% character()),
  final_parameter_regions = nrow(stage2_final),
  merge_parameter_rows = nrow(merge_table)
)

stage1_plan_full <- clean_plan(stage1_plan_object, "stage1")
stage1_plan <- if (nrow(stage1_plan_full) > 0L) {
  stage1_plan_full %>%
    transmute(
      Region = r,
      `Next action` = action_label(rho_refinement_action),
      `Run in next Stage 1` = as_yes_no(run_next_stage1),
      `Selected solvers` = selected_stage1_solvers,
      `Best sigma K-L` = round(best_stage1_sigma_KL, 4),
      `Best sigma VA-E` = round(best_stage1_sigma_VAE, 4),
      `Best at grid edge` = as_yes_no(stage1_best_on_edge),
      `Supported non-edge cells` = supported_non_edge_grid_cells,
      `Next sigma grid K-L` = map_chr(sigma_grid_KL, grid_text),
      `Next sigma grid VA-E` = map_chr(sigma_grid_VAE, grid_text),
      `Model evaluations` = stage1_next_grid_cell_count,
      `Scientific reason` = stage1_next_plan_note
    )
} else tibble()

stage2_plan_full <- clean_plan(stage2_plan_object, "stage2")
stage2_plan <- if (nrow(stage2_plan_full) > 0L) {
  stage2_plan_full %>%
    transmute(
      Region = r,
      `Next action` = action_label(rho_refinement_action),
      `Run in Stage 2` = as_yes_no(estimate_in_stage2),
      `Estimate class` = str_to_sentence(stage2_estimation_class),
      `Selected solvers` = selected_stage2_solvers,
      `Supported sigma K-L` = if_else(is.finite(sigma_KL_supported_min), paste0(number(sigma_KL_supported_min, accuracy = 0.001), " to ", number(sigma_KL_supported_max, accuracy = 0.001)), ""),
      `Supported sigma VA-E` = if_else(is.finite(sigma_VAE_supported_min), paste0(number(sigma_VAE_supported_min, accuracy = 0.001), " to ", number(sigma_VAE_supported_max, accuracy = 0.001)), ""),
      `Stage 2 sigma grid K-L` = map_chr(sigma_grid_KL, grid_text),
      `Stage 2 sigma grid VA-E` = map_chr(sigma_grid_VAE, grid_text),
      `Model evaluations` = stage2_grid_cell_count,
      `Scientific reason` = stage2_plan_note
    )
} else tibble()

readiness <- stage2_design %>%
  count(rho_refinement_action, name = "regions") %>%
  mutate(decision = action_label(rho_refinement_action)) %>%
  arrange(desc(regions))

solver_summary <- optimizer %>%
  transmute(
    method,
    convergence_share,
    valid_share,
    valid_non_edge_share,
    edge_solution_share,
    median_valid_AICc,
    median_runtime_seconds,
    regions_with_valid_results
  ) %>%
  arrange(desc(valid_non_edge_share), median_runtime_seconds)

parameter_review <- if (nrow(stage2_final) > 0L) {
  stage2_final %>% arrange(r)
} else {
  best_region %>%
    transmute(
      r,
      final_selection_rule = best_tier,
      final_support_grid_cells = NA_integer_,
      final_support_methods = method,
      final_selected_solver = method,
      final_best_AICc = AICc_plusRho,
      final_log_residual_RMSE = NA_real_,
      gamma, lambda, delta_KVA, delta_VAY, nu, rho_KL, rho_VAE, sigma_KL, sigma_VAE,
      final_any_edge_solution = on_edge_KL | on_edge_VAE,
      final_any_gamma_bound_hit = is.finite(gamma) & (near(gamma, 0.1, tol = 1e-6) | near(gamma, 10, tol = 1e-6)),
      final_review_note = "Stage 1 best method only; do not treat as final MERGE input."
    ) %>%
    arrange(r)
}

# Stage 1 and test runs do not always carry the Stage 2-only review fields.
# Add blank placeholders so the human review tables keep the same columns.
review_defaults <- list(
  final_selected_solver = NA_character_,
  final_support_grid_cells = NA_integer_,
  final_support_methods = NA_character_,
  final_log_residual_RMSE = NA_real_,
  final_any_edge_solution = NA,
  final_any_gamma_bound_hit = NA,
  final_review_note = NA_character_
)
for (review_column in names(review_defaults)) {
  if (!review_column %in% names(parameter_review)) {
    parameter_review[[review_column]] <- review_defaults[[review_column]]
  }
}

run_summary_table <- run_summary %>%
  transmute(
    `Analysis stage` = analysis_stage,
    `Estimation stage` = estimation_stage,
    `Run` = run_label,
    `Created at` = created_at,
    `Input regions` = input_regions,
    `Grid rows` = grid_rows,
    `Valid grid rows` = valid_grid_rows,
    `Valid grid share` = percent(valid_grid_share, accuracy = 0.1),
    `Regions scheduled for Stage 2` = stage2_ready_regions,
    `Regions needing another Stage 1` = next_stage1_regions,
    `Regions needing Stage 1 expansion` = stage1_expansion_regions,
    `Final parameter regions` = final_parameter_regions,
    `MERGE parameter rows` = merge_parameter_rows
  )

solver_summary_table <- solver_summary %>%
  transmute(
    Solver = method,
    `Convergence share` = percent(convergence_share, accuracy = 0.1),
    `Valid share` = percent(valid_share, accuracy = 0.1),
    `Valid non-edge share` = percent(valid_non_edge_share, accuracy = 0.1),
    `Edge share among valid results` = percent(edge_solution_share, accuracy = 0.1),
    `Median valid AICc` = round(median_valid_AICc, 2),
    `Median runtime seconds` = round(median_runtime_seconds, 2),
    `Regions with valid results` = regions_with_valid_results
  )

parameter_review_table <- parameter_review %>%
  transmute(
    Region = r,
    `Selection rule` = selection_label(final_selection_rule),
    `Support grid cells` = final_support_grid_cells,
    `Supporting solvers` = final_support_methods,
    `Selected solver` = final_selected_solver,
    `Best AICc` = round(final_best_AICc, 2),
    `Log residual RMSE` = round(final_log_residual_RMSE, 5),
    Gamma = round(gamma, 6),
    Lambda = round(lambda, 6),
    `Delta K-VA` = round(delta_KVA, 6),
    `Delta VA-Y` = round(delta_VAY, 6),
    Nu = round(nu, 6),
    `Rho K-L` = round(rho_KL, 6),
    `Rho VA-E` = round(rho_VAE, 6),
    `Sigma K-L` = round(sigma_KL, 6),
    `Sigma VA-E` = round(sigma_VAE, 6),
    `Any grid-edge solution` = as_yes_no(final_any_edge_solution),
    `Gamma bound hit` = as_yes_no(final_any_gamma_bound_hit),
    `Review note` = final_review_note
  )

if (analysis_stage == "stage2" && nrow(stage2_final) > 0L) {
  readiness_table <- stage2_final %>%
    count(final_selection_rule, name = "Regions") %>%
    transmute(`Final selection status` = selection_label(final_selection_rule), Regions)
} else {
  readiness_table <- readiness %>%
    transmute(`Stage 1 decision` = decision, Regions = regions)
}

decision_lines <- if (nrow(readiness_table) > 0L) {
  paste0("- ", readiness_table[[1]], ": ", readiness_table$Regions, " regions")
} else {
  "- No decisions were available."
}

next_step <- if (analysis_stage == "stage2" && nrow(stage2_final) > 0L) {
  "Review the final parameter table, residual diagnostics, and publication figures before copying the MERGE input table into the model workflow."
} else if ((run_summary$next_stage1_regions %||% 0) > 0L) {
  paste0(
    "Run the next targeted Stage 1 round for ",
    run_summary$next_stage1_regions,
    " regions. These regions still have boundary-driven Stage 1 evidence."
  )
} else if (any(readiness$rho_refinement_action == "boundary_supported_review_before_stage2")) {
  "Run Stage 2 for complete MERGE coverage. Interior-supported regions and boundary-supported regions must be labelled separately in the final parameter table."
} else {
  "Stage 1 has produced a Stage 2 plan for the full set of regions."
}

what_this_run_does <- if (analysis_stage == "stage2") {
  c(
    "Stage 2 is the final estimation stage for MERGE inputs.",
    "",
    "It uses the region-specific sigma grids and solver subsets learned in Stage 1. The final table covers every MERGE region available in the estimation data and labels the selection rule behind each estimate."
  )
} else {
  c(
    "Stage 1 is a scientific audit of substitution-elasticity grids and solver behaviour. It is not the final MERGE parameter-selection stage.",
    "",
    "Stage 2 produces the complete parameter table required by MERGE. Regions with valid non-edge Stage 1 support are treated as interior-supported. Regions that still prefer a grid boundary after targeted Stage 1 are treated as boundary-supported: they are estimated for MERGE coverage, but their boundary status is retained in the final review table."
  )
}

stage2_boundary_note <- if (
  analysis_stage == "stage2" &&
  nrow(stage2_final) > 0L &&
  !any(stage2_final$final_selection_rule == "boundary_supported", na.rm = TRUE)
) {
  "- Stage 2 found valid non-edge supported final estimates for all regions; Stage 1 boundary history remains documented in the saved RDS and run-plan tables."
} else {
  "- Stage 2 refines the Stage 1-supported sigma space. Boundary-supported regions remain boundary-supported in the final MERGE table rather than being re-labelled as interior estimates."
}

report_lines <- c(
  paste0("# MERGE CES ", toupper(analysis_stage), " scientific report"),
  "",
  paste0("Run analysed: **", run_label, "**"),
  paste0("Source file: `", normalizePath(results_file, winslash = "/", mustWork = FALSE), "`"),
  "",
  "## What this run does",
  "",
  what_this_run_does,
  "",
  "## Decision summary",
  "",
  decision_lines,
  "",
  "## Next step",
  "",
  next_step,
  "",
  "## Statistical interpretation",
  "",
  "- Residual diagnostics use log residuals because the CES estimation uses multiplicative errors.",
  "- AICc is used as an RSS-based log-residual heuristic for comparing grid cells and solvers; it is not treated as a fully specified likelihood result.",
  "- The two share parameters, Delta K-VA and Delta VA-Y, belong to different CES nests and are not summed.",
  stage2_boundary_note
)

#### FIGURE DATA ####
caption <- paste0("Source: ", basename(results_file), ". AICc is an RSS-based log-residual heuristic; p-values and standard errors are diagnostics.")

region_status <- stage2_design %>%
  transmute(r, decision = action_label(rho_refinement_action), rho_refinement_action)

solver_plot_data <- optimizer %>%
  mutate(runtime_minutes = median_runtime_seconds / 60)

support_intervals <- windows %>%
  mutate(decision = action_label(rho_refinement_action)) %>%
  select(r, decision, rho_refinement_action, sigma_KL_supported_min, sigma_KL_supported_max, sigma_VAE_supported_min, sigma_VAE_supported_max) %>%
  pivot_longer(
    starts_with("sigma_"),
    names_to = c("nest", ".value"),
    names_pattern = "sigma_(KL|VAE)_supported_(min|max)"
  ) %>%
  filter(is.finite(min), is.finite(max)) %>%
  mutate(
    nest = recode(nest, KL = "K-L", VAE = "VA-E"),
    width = max - min,
    r = reorder(r, width)
  )

elasticity_data <- parameter_review %>%
  left_join(region_status, by = c("r")) %>%
  filter(is.finite(sigma_KL), is.finite(sigma_VAE))

residual_data <- if (nrow(stage2_fit_summary) > 0L) {
  stage2_fit_summary %>%
    transmute(r, fit_log_rmse = selected_fit_log_RMSE) %>%
    left_join(region_status, by = "r")
} else {
  residuals %>%
    group_by(r) %>%
    summarise(
      fit_log_rmse = {
        values <- best_valid_log_residual_rmse[is.finite(best_valid_log_residual_rmse)]
        if (length(values) > 0L) min(values) else NA_real_
      },
      .groups = "drop"
    ) %>%
    left_join(region_status, by = "r")
} %>%
  filter(is.finite(fit_log_rmse))

parameter_data <- parameter_review %>%
  select(r, lambda, nu, delta_KVA, delta_VAY) %>%
  pivot_longer(-r, names_to = "parameter", values_to = "estimate") %>%
  mutate(parameter = recode(
    parameter,
    # lambda is the continuous annual TFP growth rate from
    # gamma * exp(lambda * (year - base_year)).
    lambda = "TFP growth rate",
    nu = "Returns-to-scale exponent",
    delta_KVA = "K share in VA nest",
    delta_VAY = "Value-added share in output nest"
  )) %>%
  filter(is.finite(estimate))

tfp_data <- parameter_review %>%
  filter(is.finite(lambda), is.finite(gamma)) %>%
  transmute(r, gamma, lambda, annual_tfp_growth_pct = 100 * (exp(lambda) - 1))

#### FIGURES ####
fig_readiness <- readiness %>%
  ggplot(aes(x = reorder(decision, regions), y = regions, fill = rho_refinement_action)) +
  geom_col(width = 0.62, colour = "white") +
  geom_text(aes(label = regions), hjust = -0.25, family = font_family(), size = 3.4) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = stage_colours, guide = "none", drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Stage 1 decision status by region",
    subtitle = "Regions move to Stage 2 only when the supported optimum is valid and not boundary-driven",
    x = NULL, y = "Regions", caption = caption
  ) +
  paper_theme()

fig_solver_frontier <- ggplot(solver_plot_data, aes(x = runtime_minutes, y = valid_non_edge_share)) +
  geom_point(aes(size = regions_with_valid_results, colour = edge_solution_share), alpha = 0.9) +
  geom_text(aes(label = method), nudge_y = 0.025, family = font_family(), size = 3.1, check_overlap = TRUE) +
  scale_x_continuous(labels = label_number(accuracy = 0.1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_colour_gradient(low = colour$ready, high = colour$review, labels = percent_format(accuracy = 1)) +
  scale_size_continuous(range = c(3, 9), breaks = pretty_breaks(4)) +
  labs(
    title = "Solver reliability and speed",
    subtitle = "Best Stage 2 candidates combine many valid non-edge solutions with low runtime and few edge hits",
    x = "Median runtime per region-method run (minutes)",
    y = "Valid non-edge share",
    colour = "Edge share\namong valid",
    size = "Regions with\nvalid results",
    caption = caption
  ) +
  paper_theme()

fig_support_intervals <- ggplot(support_intervals, aes(y = r, colour = nest)) +
  geom_segment(aes(x = min, xend = max, yend = r), linewidth = 1.0, alpha = 0.85) +
  geom_point(aes(x = min), size = 1.7) +
  geom_point(aes(x = max), size = 1.7) +
  facet_wrap(~ nest, scales = "free_x") +
  scale_x_continuous(trans = "log10", breaks = c(0.03, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 8), labels = label_number(accuracy = 0.01)) +
  scale_colour_manual(values = c("K-L" = colour$blue, "VA-E" = colour$purple), guide = "none") +
  labs(
    title = "Supported substitution-elasticity windows",
    subtitle = "Intervals show valid non-edge Stage 1 support within the delta AICc rule",
    x = expression("Substitution elasticity, " * sigma),
    y = NULL,
    caption = caption
  ) +
  paper_theme(base_size = 9)

fig_elasticity_map <- ggplot(elasticity_data, aes(x = sigma_KL, y = sigma_VAE)) +
  annotate("rect", xmin = 0.03, xmax = 0.5, ymin = 0.03, ymax = 8, fill = "#FFF3E6", alpha = 0.45) +
  annotate("rect", xmin = 0.03, xmax = 8, ymin = 0.03, ymax = 0.5, fill = "#FFF3E6", alpha = 0.45) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = colour$line, linewidth = 0.45) +
  geom_point(aes(colour = rho_refinement_action, shape = final_any_edge_solution), size = 3.0, alpha = 0.9) +
  scale_x_continuous(trans = "log10", breaks = c(0.03, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 8), labels = label_number(accuracy = 0.01)) +
  scale_y_continuous(trans = "log10", breaks = c(0.03, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 8), labels = label_number(accuracy = 0.01)) +
  scale_colour_manual(values = stage_colours, labels = action_label, drop = FALSE) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17), labels = c("Interior", "Grid edge")) +
  labs(
    title = "Selected substitution elasticities",
    subtitle = "Shaded zones mark very low substitution; the dashed line marks equal elasticity across nests",
    x = expression(sigma["K-L"]),
    y = expression(sigma["VA-E"]),
    colour = "Stage decision",
    shape = "Grid status",
    caption = caption
  ) +
  paper_theme()

fig_residuals <- ggplot(residual_data, aes(x = reorder(r, fit_log_rmse), y = fit_log_rmse, colour = rho_refinement_action)) +
  geom_segment(aes(xend = r, y = 0, yend = fit_log_rmse), linewidth = 0.45, colour = colour$line) +
  geom_point(size = 2.2) +
  coord_flip() +
  scale_colour_manual(values = stage_colours, labels = action_label, drop = FALSE) +
  labs(
    title = "Selected CES fit by region",
    subtitle = "Observed output is compared with output reconstructed from the final selected CES parameter vector",
    x = NULL,
    y = "Selected fit log-residual RMSE",
    colour = "Stage decision",
    caption = caption
  ) +
  paper_theme(base_size = 9)

fig_parameters <- ggplot(parameter_data, aes(x = estimate)) +
  geom_histogram(aes(y = after_stat(density)), bins = 18, fill = colour$pale, colour = "white") +
  geom_density(linewidth = 0.8, colour = colour$blue, na.rm = TRUE) +
  facet_wrap(~ parameter, scales = "free", ncol = 2) +
  labs(
    title = "Distribution of selected CES parameters",
    subtitle = "Share parameters are nest-specific and are not summed across the two nests",
    x = "Estimate",
    y = "Density",
    caption = caption
  ) +
  paper_theme()

fig_final_status <- if (nrow(stage2_final) > 0L) {
  stage2_final %>%
    count(final_selection_rule, name = "regions") %>%
    ggplot(aes(x = reorder(final_selection_rule, regions), y = regions, fill = final_selection_rule)) +
    geom_col(width = 0.62, colour = "white") +
    geom_text(aes(label = regions), hjust = -0.25, family = font_family(), size = 3.4) +
    coord_flip(clip = "off") +
    scale_fill_manual(values = selection_colours, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title = "Final parameter-selection status",
      subtitle = "Final MERGE inputs require non-edge support and substantive review",
      x = NULL, y = "Regions", caption = caption
    ) +
    paper_theme()
} else {
  fig_readiness +
    labs(title = "Final selection not yet available", subtitle = "Run Stage 2 after all required Stage 1 grid expansions are complete")
}

fig_tfp_growth <- ggplot(tfp_data, aes(x = reorder(r, annual_tfp_growth_pct), y = annual_tfp_growth_pct)) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = colour$line) +
  geom_col(fill = colour$teal, width = 0.62) +
  coord_flip() +
  labs(
    title = "Implied annual TFP growth",
    subtitle = "Computed as 100 * (exp(lambda) - 1), so values are comparable annual percentage rates",
    x = NULL,
    y = "Annual TFP growth (%)",
    caption = caption
  ) +
  paper_theme(base_size = 9)

all_plots <- list(
  readiness = fig_readiness,
  solver_frontier = fig_solver_frontier,
  support_windows = fig_support_intervals,
  residual_fit = fig_residuals,
  final_selection = fig_final_status,
  elasticities = fig_elasticity_map,
  parameters = fig_parameters,
  tfp_growth = fig_tfp_growth
)

stage_plot_order <- if (analysis_stage == "stage2") {
  c("final_selection", "elasticities", "parameters", "tfp_growth", "residual_fit", "solver_frontier", "support_windows")
} else {
  c("readiness", "solver_frontier", "support_windows", "residual_fit", "elasticities", "parameters")
}
stage_plot_order <- intersect(stage_plot_order, names(all_plots))
plots <- all_plots[stage_plot_order]
names(plots) <- paste0(figure_prefix, ".", sprintf("%02d", seq_along(plots)), ".", names(plots))
assign(paste0(stage_tag, "_plots"), plots, envir = .GlobalEnv)

#### SAVE ####
dir.create(stage_path, recursive = TRUE, showWarnings = FALSE)
save_figure_book(plots, figure_book_file)
save_individual_figures(plots, figure_folder)
if (isTRUE(show_plots)) walk(plots, print)
safe_write_lines(report_lines, report_file)

if (save_history && analysis_stage != "test") {
  dir.create(history_folder, recursive = TRUE, showWarnings = FALSE)
  invisible(file.copy(figure_book_file, history_figure_file, overwrite = TRUE))
  invisible(file.copy(report_file, history_report_file, overwrite = TRUE))
}

analysis_results <- list(
  metadata = list(
    analysis_stage = analysis_stage,
    run_label = run_label,
    created_at = Sys.time(),
    source_results_file = normalizePath(results_file, winslash = "/", mustWork = FALSE),
    figure_book_file = normalizePath(figure_book_file, winslash = "/", mustWork = FALSE),
    figure_folder = normalizePath(figure_folder, winslash = "/", mustWork = FALSE),
    report_file = normalizePath(report_file, winslash = "/", mustWork = FALSE)
  ),
  run_summary = run_summary,
  run_summary_table = run_summary_table,
  solver_summary = solver_summary,
  solver_summary_table = solver_summary_table,
  stage2_readiness_summary = readiness,
  stage2_readiness_table = readiness_table,
  next_stage1_plan_full = stage1_plan_full,
  next_stage1_plan = stage1_plan,
  stage2_plan_full = stage2_plan_full,
  stage2_plan = stage2_plan,
  final_parameter_table_for_review = parameter_review,
  final_parameter_review_table = parameter_review_table,
  plots = plots
)
assign(paste0(stage_tag, "_analysis"), analysis_results, envir = .GlobalEnv)

safe_save_rds(analysis_results, analysis_file)
if (save_history && analysis_stage != "test") safe_save_rds(analysis_results, history_analysis_file)

tables <- list(
  run_summary = run_summary_table,
  solver_summary = solver_summary_table,
  stage2_readiness_summary = readiness_table,
  next_stage1_plan = stage1_plan,
  stage2_plan = stage2_plan,
  candidate_parameter_review = parameter_review_table
)
if (analysis_stage == "stage2") names(tables)[names(tables) == "candidate_parameter_review"] <- "final_parameter_review"
if (nrow(merge_table) > 0L) tables$merge_iam_parameter_table <- merge_table
write_named_tables(tables, stage_path, analysis_stage)

if (write_full_plan_tables) {
  write_named_tables(list(next_stage1_plan_full = stage1_plan_full, stage2_plan_full = stage2_plan_full), stage_path, analysis_stage)
}

message("Analysis results written to: ", normalizePath(analysis_file, winslash = "/"))
message("Publication figure book written to: ", normalizePath(figure_book_file, winslash = "/"))
if (isTRUE(write_individual_figures)) {
  message("Individual figures written to: ", normalizePath(figure_folder, winslash = "/"))
}
message("Scientific report written to: ", normalizePath(report_file, winslash = "/"))
