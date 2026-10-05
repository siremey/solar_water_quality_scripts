# Otto and Markum mixed-effects model analysis
# Fits the primary models, sensitivity models, tables, and SI figures.

# 1. Load packages -----------------------------------------------------------
required_packages <- c(
  "tidyverse", "lme4", "lmerTest", "MuMIn", "performance",
  "writexl", "patchwork"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install these packages before running the script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(MuMIn)
  library(performance)
  library(writexl)
  library(patchwork)
})

options(na.action = "na.fail")

# 2. Set file paths ----------------------------------------------------------
PROJECT_DIR <- "C:/Users/User/Box/My Box Notes/manuscript_revision/WQ_revision"
OTTO_INPUT_FILE <- file.path(
  PROJECT_DIR, "data", "05_model_components", "otto_model_components_v1.csv"
)
MARKUM_INPUT_FILE <- file.path(
  PROJECT_DIR, "data", "05_model_components", "markum_model_components_v1.csv"
)
OUTPUT_PARENT_DIR <- file.path(PROJECT_DIR, "outputs")

if (!dir.exists(PROJECT_DIR)) {
  stop("Project folder not found: ", PROJECT_DIR)
}
if (!file.exists(OTTO_INPUT_FILE)) {
  stop("Otto input file not found: ", OTTO_INPUT_FILE)
}
if (!file.exists(MARKUM_INPUT_FILE)) {
  stop("Markum input file not found: ", MARKUM_INPUT_FILE)
}

project_dir <- normalizePath(PROJECT_DIR, winslash = "/", mustWork = TRUE)
otto_input <- normalizePath(OTTO_INPUT_FILE, winslash = "/", mustWork = TRUE)
markum_input <- normalizePath(MARKUM_INPUT_FILE, winslash = "/", mustWork = TRUE)
message("Otto input: ", otto_input)
message("Markum input: ", markum_input)
message("Output parent: ", OUTPUT_PARENT_DIR)

make_new_output_dir <- function(parent_dir, stem) {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  candidate <- file.path(parent_dir, paste0(stem, "_", timestamp))
  suffix <- 0L

  while (dir.exists(candidate) || file.exists(candidate)) {
    suffix <- suffix + 1L
    candidate <- file.path(
      parent_dir,
      paste0(stem, "_", timestamp, "_", sprintf("%02d", suffix))
    )
  }

  dir.create(candidate, recursive = TRUE, showWarnings = FALSE)
  normalizePath(candidate, winslash = "/", mustWork = TRUE)
}

output_parent <- OUTPUT_PARENT_DIR
dir.create(output_parent, recursive = TRUE, showWarnings = FALSE)

output_dir <- make_new_output_dir(
  output_parent,
  "combined_Otto_Markum_final_models"
)

table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")
model_dir <- file.path(output_dir, "models")
audit_dir <- file.path(output_dir, "audit")
log_dir <- file.path(output_dir, "logs")

purrr::walk(
  c(table_dir, figure_dir, model_dir, audit_dir, log_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
)

# 3. Define data-checking functions -----------------------------------------
read_component_file <- function(path) {
  readr::read_csv(
    path,
    col_types = cols(
      Date = col_character(),
      Site = col_character(),
      station_uid = col_character(),
      sampling_event_id = col_character(),
      observation_id = col_character(),
      .default = col_guess()
    ),
    show_col_types = FALSE
  )
}

require_columns <- function(data, required, site_name) {
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns) > 0) {
    stop(
      site_name, " is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
}

z_scale_checked <- function(x, variable_name, site_name) {
  if (any(!is.finite(x))) {
    stop(site_name, ": ", variable_name, " contains non-finite values.")
  }
  x_sd <- stats::sd(x)
  if (!is.finite(x_sd) || x_sd <= 0) {
    stop(site_name, ": ", variable_name, " has zero or undefined SD.")
  }
  as.numeric(scale(x))
}

validate_station_mapping <- function(data, site_name) {
  crosswalk <- data %>%
    distinct(station_uid, Site) %>%
    arrange(station_uid, Site)

  check <- crosswalk %>%
    summarise(
      n_station_uid = n_distinct(station_uid),
      n_Site = n_distinct(Site),
      station_uid_maps_to_one_Site = all(table(station_uid) == 1),
      Site_maps_to_one_station_uid = all(table(Site) == 1)
    )

  if (
    !check$station_uid_maps_to_one_Site ||
    !check$Site_maps_to_one_station_uid
  ) {
    stop(site_name, ": Site and station_uid are not one-to-one.")
  }

  list(crosswalk = crosswalk, check = check)
}

check_missingness <- function(data, variables, site_name) {
  out <- tibble(
    site = site_name,
    variable = variables,
    missing_n = vapply(
      data[variables],
      function(x) sum(is.na(x)),
      integer(1)
    )
  )

  if (any(out$missing_n > 0)) {
    stop(site_name, ": missing values exist in required model variables.")
  }
  out
}

check_solar_components <- function(data, site_name) {
  by_station <- data %>%
    group_by(station_uid) %>%
    summarise(
      n = n(),
      mean_solar_within = mean(solar_within),
      solar_between_unique_n = n_distinct(solar_between),
      .groups = "drop"
    )

  if (any(abs(by_station$mean_solar_within) > 1e-8)) {
    warning(site_name, ": solar_within is not exactly centered within station.")
  }
  if (any(by_station$solar_between_unique_n != 1)) {
    stop(site_name, ": solar_between is not constant within station.")
  }

  overall <- tibble(
    site = site_name,
    component = c("solar_within", "solar_between"),
    mean = c(mean(data$solar_within), mean(data$solar_between)),
    sd = c(sd(data$solar_within), sd(data$solar_between))
  )

  if (any(abs(overall$mean) > 1e-7) || any(abs(overall$sd - 1) > 1e-7)) {
    warning(
      site_name,
      ": within/between solar components are not globally mean 0 and SD 1. ",
      "Confirm that the frozen component variables were standardized separately."
    )
  }

  list(by_station = by_station, overall = overall)
}

# 4. Read and prepare data ---------------------------------------------------
otto_required <- c(
  "Turbidity", "Turbidity_sqrt", "Site", "station_uid",
  "sampling_event_id", "observation_id", "drainage_area_km2",
  "max_rain_event", "NDVI", "solar_panel_pct", "solar_between",
  "solar_within", "A_mean_t_ha_yr", "treeVegetation_pct", "max_air_temp",
  "cultivated_bareland_pct", "builtup_pct", "water_pct", "season"
)

markum_required <- c(
  "Turbidity", "Turbidity_log10", "Site", "station_uid",
  "sampling_event_id", "observation_id", "drainage_area_km2",
  "max_rain_event", "NDVI", "solar_panel_pct", "solar_between",
  "solar_within", "A_mean_t_ha_yr", "treeVegetation_pct", "max_air_temp"
)

otto_raw <- read_component_file(otto_input)
markum_raw <- read_component_file(markum_input)

require_columns(otto_raw, otto_required, "Otto")
require_columns(markum_raw, markum_required, "Markum")

if ("included_in_candidate_models" %in% names(otto_raw)) {
  otto_raw <- otto_raw %>% filter(as.logical(included_in_candidate_models))
}
if ("included_in_candidate_models" %in% names(markum_raw)) {
  markum_raw <- markum_raw %>% filter(as.logical(included_in_candidate_models))
}

otto <- otto_raw %>%
  mutate(
    station_uid = factor(station_uid),
    sampling_event_id = factor(sampling_event_id),
    Site = factor(Site),
    season = factor(season)
  )

markum <- markum_raw %>%
  mutate(
    station_uid = factor(station_uid),
    sampling_event_id = factor(sampling_event_id),
    Site = factor(Site)
  )

if (nrow(otto) != 129) stop("Expected 129 eligible Otto observations.")
if (n_distinct(otto$station_uid) != 13) stop("Expected 13 Otto stations.")
if (n_distinct(otto$sampling_event_id) != 15) {
  stop("Expected 15 Otto sampling-event levels.")
}

if (nrow(markum) != 56) stop("Expected 56 eligible Markum observations.")
if (n_distinct(markum$station_uid) != 4) stop("Expected 4 Markum stations.")
if (n_distinct(markum$sampling_event_id) != 27) {
  stop("Expected 27 Markum sampling-event levels.")
}

if (anyDuplicated(otto$observation_id) > 0) {
  stop("Otto observation_id is not unique.")
}
if (anyDuplicated(markum$observation_id) > 0) {
  stop("Markum observation_id is not unique.")
}

if (any(!is.finite(otto$Turbidity)) || any(otto$Turbidity <= 0)) {
  stop("All Otto turbidity observations must be finite and greater than zero.")
}
if (any(!is.finite(markum$Turbidity)) || any(markum$Turbidity <= 0)) {
  stop("All Markum turbidity observations must be finite and greater than zero.")
}

sqrt_difference <- max(
  abs(otto$Turbidity_sqrt - sqrt(otto$Turbidity)),
  na.rm = TRUE
)
log10_difference <- max(
  abs(markum$Turbidity_log10 - log10(markum$Turbidity)),
  na.rm = TRUE
)

if (!is.finite(sqrt_difference) || sqrt_difference > 1e-8) {
  stop("Otto Turbidity_sqrt does not equal sqrt(Turbidity).")
}
if (!is.finite(log10_difference) || log10_difference > 1e-10) {
  stop("Markum Turbidity_log10 does not equal log10(Turbidity).")
}

otto_missingness <- check_missingness(otto, otto_required, "Otto")
markum_missingness <- check_missingness(markum, markum_required, "Markum")

otto_mapping <- validate_station_mapping(otto, "Otto")
markum_mapping <- validate_station_mapping(markum, "Markum")

otto_solar_component_check <- check_solar_components(otto, "Otto")
markum_solar_component_check <- check_solar_components(markum, "Markum")

otto_scale_variables <- c(
  "drainage_area_km2", "max_rain_event", "NDVI", "solar_panel_pct",
  "A_mean_t_ha_yr", "treeVegetation_pct", "max_air_temp",
  "cultivated_bareland_pct", "builtup_pct", "water_pct"
)

markum_scale_variables <- c(
  "drainage_area_km2", "max_rain_event", "NDVI", "solar_panel_pct",
  "A_mean_t_ha_yr", "treeVegetation_pct", "max_air_temp"
)

otto <- otto %>%
  mutate(
    drainage_area_z = z_scale_checked(drainage_area_km2, "drainage_area_km2", "Otto"),
    max_rain_event_z = z_scale_checked(max_rain_event, "max_rain_event", "Otto"),
    NDVI_z = z_scale_checked(NDVI, "NDVI", "Otto"),
    solar_panel_pct_z = z_scale_checked(solar_panel_pct, "solar_panel_pct", "Otto"),
    A_mean_t_ha_yr_z = z_scale_checked(A_mean_t_ha_yr, "A_mean_t_ha_yr", "Otto"),
    treeVegetation_pct_z = z_scale_checked(treeVegetation_pct, "treeVegetation_pct", "Otto"),
    max_air_temp_z = z_scale_checked(max_air_temp, "max_air_temp", "Otto"),
    cultivated_bareland_pct_z = z_scale_checked(
      cultivated_bareland_pct,
      "cultivated_bareland_pct",
      "Otto"
    ),
    builtup_pct_z = z_scale_checked(builtup_pct, "builtup_pct", "Otto"),
    water_pct_z = z_scale_checked(water_pct, "water_pct", "Otto")
  )

markum <- markum %>%
  mutate(
    drainage_area_z = z_scale_checked(drainage_area_km2, "drainage_area_km2", "Markum"),
    max_rain_event_z = z_scale_checked(max_rain_event, "max_rain_event", "Markum"),
    NDVI_z = z_scale_checked(NDVI, "NDVI", "Markum"),
    solar_panel_pct_z = z_scale_checked(solar_panel_pct, "solar_panel_pct", "Markum"),
    A_mean_t_ha_yr_z = z_scale_checked(A_mean_t_ha_yr, "A_mean_t_ha_yr", "Markum"),
    treeVegetation_pct_z = z_scale_checked(treeVegetation_pct, "treeVegetation_pct", "Markum"),
    max_air_temp_z = z_scale_checked(max_air_temp, "max_air_temp", "Markum")
  )

scaling_parameters <- bind_rows(
  tibble(site = "Otto", variable = otto_scale_variables) %>%
    mutate(
      center = map_dbl(variable, ~ mean(otto[[.x]])),
      scale = map_dbl(variable, ~ sd(otto[[.x]]))
    ),
  tibble(site = "Markum", variable = markum_scale_variables) %>%
    mutate(
      center = map_dbl(variable, ~ mean(markum[[.x]])),
      scale = map_dbl(variable, ~ sd(markum[[.x]]))
    )
)

data_audit <- tibble(
  site = c("Otto", "Markum"),
  input_file = c(otto_input, markum_input),
  observations = c(nrow(otto), nrow(markum)),
  station_levels = c(n_distinct(otto$station_uid), n_distinct(markum$station_uid)),
  sampling_event_levels = c(
    n_distinct(otto$sampling_event_id),
    n_distinct(markum$sampling_event_id)
  ),
  response = c("Turbidity_sqrt", "Turbidity_log10"),
  response_scale = c("sqrt(NTU)", "log10(NTU)"),
  minimum_turbidity_NTU = c(min(otto$Turbidity), min(markum$Turbidity)),
  maximum_turbidity_NTU = c(max(otto$Turbidity), max(markum$Turbidity)),
  transformation_check_max_abs_difference = c(sqrt_difference, log10_difference)
)

readr::write_csv(otto, file.path(audit_dir, "otto_model_cohort_used.csv"))
readr::write_csv(markum, file.path(audit_dir, "markum_model_cohort_used.csv"))
readr::write_csv(data_audit, file.path(audit_dir, "data_audit.csv"))
readr::write_csv(
  bind_rows(otto_missingness, markum_missingness),
  file.path(audit_dir, "missingness.csv")
)
readr::write_csv(scaling_parameters, file.path(audit_dir, "scaling_parameters.csv"))
readr::write_csv(
  bind_rows(
    otto_solar_component_check$overall,
    markum_solar_component_check$overall
  ),
  file.path(audit_dir, "solar_component_scaling_check.csv")
)
readr::write_csv(
  bind_rows(
    otto_solar_component_check$by_station %>% mutate(site = "Otto", .before = 1),
    markum_solar_component_check$by_station %>% mutate(site = "Markum", .before = 1)
  ),
  file.path(audit_dir, "solar_component_station_check.csv")
)
readr::write_csv(
  bind_rows(
    otto_mapping$crosswalk %>% mutate(site = "Otto", .before = 1),
    markum_mapping$crosswalk %>% mutate(site = "Markum", .before = 1)
  ),
  file.path(audit_dir, "station_crosswalk.csv")
)

# 5. Define model and output functions --------------------------------------
fit_control <- lme4::lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 200000)
)

random_terms <- c(
  "(1 | station_uid)",
  "(1 | sampling_event_id)"
)

make_formula <- function(response, fixed_terms) {
  stats::as.formula(
    paste(
      response,
      "~",
      paste(c(fixed_terms, random_terms), collapse = " + ")
    )
  )
}

fit_ml <- function(data, response, fixed_terms) {
  lmerTest::lmer(
    make_formula(response, fixed_terms),
    data = data,
    REML = FALSE,
    control = fit_control
  )
}

# Refit one ML model with REML.
refit_reml <- function(model, model_data) {
  lmerTest::lmer(
    stats::formula(model),
    data = model_data,
    REML = TRUE,
    control = fit_control
  )
}

fit_pair <- function(data, response, base_terms) {
  list(
    no_solar = fit_ml(data, response, base_terms),
    solar = fit_ml(data, response, c(base_terms, "solar_panel_pct_z"))
  )
}

fit_solar_formulations <- function(data, response, base_terms, primary_pair) {
  list(
    no_solar = primary_pair$no_solar,
    overall_solar = primary_pair$solar,
    within_solar = fit_ml(data, response, c(base_terms, "solar_within")),
    between_solar = fit_ml(data, response, c(base_terms, "solar_between")),
    within_between = fit_ml(
      data,
      response,
      c(base_terms, "solar_within", "solar_between")
    )
  )
}

assert_comparable <- function(models, comparison_name) {
  if (!all(vapply(models, inherits, logical(1), what = "merMod"))) {
    stop(comparison_name, ": one or more objects are not fitted merMod models.")
  }
  if (any(vapply(models, lme4::isREML, logical(1)))) {
    stop(comparison_name, ": model comparisons must use REML = FALSE.")
  }

  responses <- vapply(
    models,
    function(x) all.vars(stats::formula(x))[[1]],
    character(1)
  )
  if (length(unique(responses)) != 1) {
    stop(comparison_name, ": response variables differ.")
  }

  row_sets <- lapply(models, function(x) rownames(stats::model.frame(x)))
  if (!all(vapply(row_sets[-1], identical, logical(1), y = row_sets[[1]]))) {
    stop(comparison_name, ": models do not use identical observations/order.")
  }

  random_sets <- lapply(models, function(x) names(lme4::getME(x, "flist")))
  if (!all(vapply(random_sets[-1], identical, logical(1), y = random_sets[[1]]))) {
    stop(comparison_name, ": random-effects grouping structures differ.")
  }

  invisible(TRUE)
}

convergence_details <- function(model) {
  message <- model@optinfo$conv$lme4$messages
  optimizer_code <- model@optinfo$conv$opt
  gradient <- model@optinfo$derivs$gradient

  message_text <- if (is.null(message) || length(message) == 0) NA_character_ else {
    paste(message, collapse = "; ")
  }
  code_value <- if (is.null(optimizer_code) || length(optimizer_code) == 0) NA_integer_ else {
    as.integer(optimizer_code[[1]])
  }

  tibble(
    convergence_status = ifelse(
      (is.na(code_value) || code_value == 0L) && is.na(message_text),
      "OK",
      "Check"
    ),
    optimizer_code = code_value,
    convergence_message = message_text,
    maximum_absolute_gradient = if (is.null(gradient) || length(gradient) == 0) {
      NA_real_
    } else {
      max(abs(gradient))
    }
  )
}

variance_components <- function(model) {
  vc <- as.data.frame(lme4::VarCorr(model))

  get_variance <- function(group_name) {
    value <- vc$vcov[vc$grp == group_name]
    if (length(value) == 0) NA_real_ else value[[1]]
  }

  tibble(
    station_variance = get_variance("station_uid"),
    sampling_event_variance = get_variance("sampling_event_id"),
    residual_variance = get_variance("Residual")
  )
}

safe_r2 <- function(model) {
  tryCatch(
    {
      result <- suppressWarnings(performance::r2_nakagawa(model))
      tibble(
        marginal_R2 = unname(result$R2_marginal),
        conditional_R2 = unname(result$R2_conditional)
      )
    },
    error = function(e) {
      tibble(marginal_R2 = NA_real_, conditional_R2 = NA_real_)
    }
  )
}

fixed_formula_text <- function(model) {
  paste(deparse(lme4::nobars(stats::formula(model))), collapse = " ")
}

model_metrics_row <- function(
    model,
    panel,
    site,
    analysis_role,
    family,
    model_id,
    solar_included
) {
  model_frame <- stats::model.frame(model)
  log_likelihood <- logLik(model)
  variance <- variance_components(model)
  convergence <- convergence_details(model)

  tibble(
    panel = panel,
    site = site,
    analysis_role = analysis_role,
    sensitivity_family = family,
    model_id = model_id,
    fixed_effect_specification = fixed_formula_text(model),
    solar_included = solar_included,
    drainage_included = "drainage_area_z" %in% all.vars(stats::formula(model)),
    n = nrow(model_frame),
    station_levels = n_distinct(model_frame$station_uid),
    sampling_event_levels = n_distinct(model_frame$sampling_event_id),
    estimated_parameters = attr(log_likelihood, "df"),
    log_likelihood = as.numeric(log_likelihood),
    AICc = MuMIn::AICc(model),
    station_variance = variance$station_variance,
    sampling_event_variance = variance$sampling_event_variance,
    residual_variance = variance$residual_variance,
    singular_fit = lme4::isSingular(model, tol = 1e-4),
    convergence_status = convergence$convergence_status,
    optimizer_code = convergence$optimizer_code,
    convergence_message = convergence$convergence_message,
    maximum_absolute_gradient = convergence$maximum_absolute_gradient
  )
}

lrt_values <- function(smaller_model, larger_model, comparison_name) {
  assert_comparable(
    list(smaller = smaller_model, larger = larger_model),
    comparison_name
  )

  small_ll <- logLik(smaller_model)
  large_ll <- logLik(larger_model)
  df_difference <- attr(large_ll, "df") - attr(small_ll, "df")
  if (df_difference <= 0) {
    stop(comparison_name, ": larger model does not have more parameters.")
  }

  chi_square <- 2 * (as.numeric(large_ll) - as.numeric(small_ll))

  tibble(
    LRT_comparison = comparison_name,
    LRT_chi_square = chi_square,
    LRT_df = df_difference,
    LRT_p_value = stats::pchisq(
      chi_square,
      df = df_difference,
      lower.tail = FALSE
    )
  )
}

pair_panel_table <- function(
    pair,
    panel,
    site,
    analysis_role,
    family,
    model_ids,
    comparison_name
) {
  assert_comparable(pair, comparison_name)

  table <- bind_rows(
    model_metrics_row(
      pair$no_solar,
      panel,
      site,
      analysis_role,
      family,
      model_ids[[1]],
      FALSE
    ),
    model_metrics_row(
      pair$solar,
      panel,
      site,
      analysis_role,
      family,
      model_ids[[2]],
      TRUE
    )
  ) %>%
    mutate(
      delta_AICc = AICc - min(AICc),
      relative_likelihood = exp(-0.5 * delta_AICc),
      Akaike_weight = relative_likelihood / sum(relative_likelihood)
    )

  lrt <- lrt_values(pair$no_solar, pair$solar, comparison_name)

  table %>%
    mutate(
      LRT_comparison = if_else(solar_included, lrt$LRT_comparison, NA_character_),
      LRT_chi_square = if_else(solar_included, lrt$LRT_chi_square, NA_real_),
      LRT_df = if_else(solar_included, as.numeric(lrt$LRT_df), NA_real_),
      LRT_p_value = if_else(solar_included, lrt$LRT_p_value, NA_real_)
    ) %>%
    select(
      panel, site, analysis_role, sensitivity_family, model_id,
      fixed_effect_specification, solar_included, drainage_included,
      n, station_levels, sampling_event_levels, estimated_parameters,
      log_likelihood, AICc, delta_AICc, Akaike_weight,
      LRT_comparison, LRT_chi_square, LRT_df, LRT_p_value,
      station_variance, sampling_event_variance, residual_variance,
      convergence_status, singular_fit, optimizer_code,
      maximum_absolute_gradient, convergence_message
    )
}

five_model_table <- function(models, panel, site, analysis_role, model_ids) {
  assert_comparable(models, paste(site, "five-model solar comparison"))

  labels <- c(
    no_solar = "No solar",
    overall_solar = "Overall solar",
    within_solar = "Within-station solar",
    between_solar = "Between-station solar",
    within_between = "Within + between solar"
  )

  result <- bind_rows(lapply(seq_along(models), function(i) {
    model_name <- names(models)[[i]]
    model_metrics_row(
      models[[i]],
      panel,
      site,
      analysis_role,
      "Five solar formulations",
      model_ids[[i]],
      model_name != "no_solar"
    ) %>%
      mutate(solar_formulation = labels[[model_name]])
  })) %>%
    mutate(
      delta_AICc = AICc - min(AICc),
      relative_likelihood = exp(-0.5 * delta_AICc),
      Akaike_weight = relative_likelihood / sum(relative_likelihood),
      AICc_rank = rank(AICc, ties.method = "first")
    ) %>%
    select(
      panel, site, analysis_role, model_id, solar_formulation,
      fixed_effect_specification, n, station_levels, sampling_event_levels,
      estimated_parameters, log_likelihood, AICc, delta_AICc,
      Akaike_weight, AICc_rank, station_variance,
      sampling_event_variance, residual_variance, convergence_status,
      singular_fit, optimizer_code, maximum_absolute_gradient,
      convergence_message
    )

  result
}

readable_predictor <- function(term) {
  labels <- c(
    "(Intercept)" = "Intercept",
    "drainage_area_z" = "Drainage area",
    "max_rain_event_z" = "Maximum rainfall event",
    "NDVI_z" = "NDVI",
    "solar_panel_pct_z" = "Solar-panel cover",
    "solar_within" = "Within-station solar cover",
    "solar_between" = "Between-station solar cover",
    "cultivated_bareland_pct_z" = "Cultivated/bare land",
    "treeVegetation_pct_z" = "Tree vegetation",
    "builtup_pct_z" = "Built-up land",
    "water_pct_z" = "Open water",
    "A_mean_t_ha_yr_z" = "RUSLE mean soil loss",
    "max_air_temp_z" = "Maximum air temperature"
  )

  mapped <- unname(labels[term])
  fallback <- stringr::str_replace_all(term, "_", " ")
  ifelse(is.na(mapped), fallback, mapped)
}

vif_values <- function(model) {
  result <- tryCatch(
    as.data.frame(suppressWarnings(performance::check_collinearity(model))),
    error = function(e) data.frame()
  )

  if (nrow(result) == 0) {
    return(tibble(model_term = character(), VIF = numeric()))
  }

  term_column <- intersect(
    c("Parameter", "Term", "parameter", "term"),
    names(result)
  )
  vif_column <- intersect(c("VIF", "vif"), names(result))

  if (length(term_column) == 0 || length(vif_column) == 0) {
    return(tibble(model_term = character(), VIF = numeric()))
  }

  tibble(
    model_term = as.character(result[[term_column[[1]]]]),
    VIF = as.numeric(result[[vif_column[[1]]]])
  )
}

coefficient_table <- function(
    model,
    panel,
    site,
    analysis_role,
    model_id,
    response_scale
) {
  coefficients <- as.data.frame(coef(summary(model))) %>%
    rownames_to_column("model_term")

  if (!"df" %in% names(coefficients)) {
    coefficients$df <- Inf
  }

  p_column <- grep("^Pr\\(", names(coefficients), value = TRUE)[[1]]

  out <- coefficients %>%
    transmute(
      panel = panel,
      site = site,
      analysis_role = analysis_role,
      model_id = model_id,
      model_term = model_term,
      predictor = readable_predictor(model_term),
      estimate = Estimate,
      standard_error = `Std. Error`,
      Satterthwaite_df = df,
      t_value = `t value`,
      p_value = .data[[p_column]],
      CI_95_low = estimate - qt(0.975, df = Satterthwaite_df) * standard_error,
      CI_95_high = estimate + qt(0.975, df = Satterthwaite_df) * standard_error,
      response_scale = response_scale,
      predictor_increment = if_else(
        model_term == "(Intercept)",
        NA_character_,
        "1 SD"
      )
    )

  out %>%
    left_join(vif_values(model), by = "model_term")
}

model_summary <- function(
    model,
    panel,
    site,
    analysis_role,
    model_id,
    response_scale
) {
  model_frame <- model.frame(model)
  variance <- variance_components(model)
  r2 <- safe_r2(model)
  convergence <- convergence_details(model)

  tibble(
    panel = panel,
    site = site,
    analysis_role = analysis_role,
    model_id = model_id,
    fixed_effect_specification = fixed_formula_text(model),
    response_scale = response_scale,
    observations = nrow(model_frame),
    station_levels = n_distinct(model_frame$station_uid),
    sampling_event_levels = n_distinct(model_frame$sampling_event_id),
    station_variance = variance$station_variance,
    sampling_event_variance = variance$sampling_event_variance,
    residual_variance = variance$residual_variance,
    marginal_R2 = r2$marginal_R2,
    conditional_R2 = r2$conditional_R2,
    singular_fit = lme4::isSingular(model, tol = 1e-4),
    convergence_status = convergence$convergence_status,
    optimizer_code = convergence$optimizer_code,
    maximum_absolute_gradient = convergence$maximum_absolute_gradient,
    convergence_message = convergence$convergence_message
  )
}

write_table_package <- function(filename_stem, sheets) {
  workbook_path <- file.path(table_dir, paste0(filename_stem, ".xlsx"))
  writexl::write_xlsx(sheets, workbook_path)

  csv_dir <- file.path(table_dir, filename_stem)
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

  purrr::iwalk(sheets, function(data, sheet_name) {
    safe_name <- stringr::str_replace_all(sheet_name, "[^A-Za-z0-9]+", "_")
    readr::write_csv(data, file.path(csv_dir, paste0(safe_name, ".csv")), na = "")
  })

  workbook_path
}

# 6. Fit primary and within-between models ----------------------------------
otto_primary_base <- c(
  "max_rain_event_z", "NDVI_z", "cultivated_bareland_pct_z",
  "treeVegetation_pct_z", "builtup_pct_z", "water_pct_z"
)

otto_drainage_base <- c("drainage_area_z", otto_primary_base)

markum_primary_base <- c(
  "drainage_area_z", "max_rain_event_z", "NDVI_z"
)

otto_primary_pair <- fit_pair(otto, "Turbidity_sqrt", otto_primary_base)
otto_drainage_pair <- fit_pair(otto, "Turbidity_sqrt", otto_drainage_base)
markum_primary_pair <- fit_pair(markum, "Turbidity_log10", markum_primary_base)

otto_solar_models <- fit_solar_formulations(
  otto,
  "Turbidity_sqrt",
  otto_primary_base,
  otto_primary_pair
)

markum_solar_models <- fit_solar_formulations(
  markum,
  "Turbidity_log10",
  markum_primary_base,
  markum_primary_pair
)

# Refit final models with REML.
otto_primary_REML <- refit_reml(otto_primary_pair$solar, otto)
markum_primary_REML <- refit_reml(markum_primary_pair$solar, markum)
otto_drainage_REML <- refit_reml(otto_drainage_pair$solar, otto)

# Refit within-between models with REML.
otto_within_between_REML <- refit_reml(
  otto_solar_models$within_between,
  otto
)
markum_within_between_REML <- refit_reml(
  markum_solar_models$within_between,
  markum
)

# 7. Fit sensitivity models -------------------------------------------------
# Exclude NDVI from the Otto RUSLE model because RUSLE C is NDVI-derived.
otto_rusle_pair <- fit_pair(
  otto,
  "Turbidity_sqrt",
  c("drainage_area_z", "max_rain_event_z", "A_mean_t_ha_yr_z")
)

otto_season_pair <- fit_pair(
  otto,
  "Turbidity_sqrt",
  c("drainage_area_z", "max_rain_event_z", "NDVI_z", "season")
)

# Retain the earlier Markum RUSLE and NDVI model for comparison.
markum_rusle_ndvi_pair <- fit_pair(
  markum,
  "Turbidity_log10",
  c("A_mean_t_ha_yr_z", "NDVI_z")
)

markum_rusle_rain_pair <- fit_pair(
  markum,
  "Turbidity_log10",
  c("max_rain_event_z", "A_mean_t_ha_yr_z")
)

# Fit the remaining sensitivity models.
otto_core_pair <- fit_pair(
  otto,
  "Turbidity_sqrt",
  c("drainage_area_z", "max_rain_event_z", "NDVI_z")
)

otto_compact_land_pair <- fit_pair(
  otto,
  "Turbidity_sqrt",
  c(
    "drainage_area_z", "max_rain_event_z", "NDVI_z",
    "cultivated_bareland_pct_z", "treeVegetation_pct_z"
  )
)

otto_temperature_pair <- fit_pair(
  otto,
  "Turbidity_sqrt",
  c("drainage_area_z", "max_rain_event_z", "NDVI_z", "max_air_temp_z")
)

otto_reduced_vegetation_pair <- fit_pair(
  otto,
  "Turbidity_sqrt",
  c(
    "drainage_area_z", "max_rain_event_z", "NDVI_z",
    "treeVegetation_pct_z"
  )
)

markum_tree_pair <- fit_pair(
  markum,
  "Turbidity_log10",
  c("drainage_area_z", "max_rain_event_z", "treeVegetation_pct_z")
)

markum_temperature_pair <- fit_pair(
  markum,
  "Turbidity_log10",
  c(
    "drainage_area_z", "max_rain_event_z", "NDVI_z",
    "max_air_temp_z"
  )
)

# 8. Build matched model-comparison tables ---------------------------------
comparison_panel_a <- pair_panel_table(
  otto_primary_pair,
  "A",
  "Otto",
  "Primary no-drainage analysis",
  "Primary",
  c("O-P0", "O-P1"),
  "Otto primary no-drainage: add overall solar"
)

comparison_panel_b <- pair_panel_table(
  otto_drainage_pair,
  "B",
  "Otto",
  "Drainage-adjusted sensitivity",
  "Drainage sensitivity",
  c("O-D0", "O-D1"),
  "Otto drainage-adjusted sensitivity: add overall solar"
)

comparison_panel_c <- pair_panel_table(
  markum_primary_pair,
  "C",
  "Markum",
  "Primary drainage-adjusted analysis",
  "Primary",
  c("M-P0", "M-P1"),
  "Markum primary drainage-adjusted: add overall solar"
)

comparison_panel_d <- bind_rows(
  pair_panel_table(
    otto_rusle_pair,
    "D",
    "Otto",
    "Other matched sensitivity",
    "Drainage + rainfall + RUSLE A",
    c("O-RUSLE0", "O-RUSLE1"),
    "Otto RUSLE sensitivity: add overall solar"
  ),
  pair_panel_table(
    otto_season_pair,
    "D",
    "Otto",
    "Other matched sensitivity",
    "Drainage + rainfall + NDVI + season",
    c("O-SEAS0", "O-SEAS1"),
    "Otto seasonal sensitivity: add overall solar"
  ),
  pair_panel_table(
    markum_rusle_ndvi_pair,
    "D",
    "Markum",
    "Other matched sensitivity",
    "RUSLE A + NDVI sensitivity",
    c("M-RN0", "M-RN1"),
    "Markum RUSLE + NDVI sensitivity: add overall solar"
  ),
  pair_panel_table(
    markum_rusle_rain_pair,
    "D",
    "Markum",
    "Other matched sensitivity",
    "Rainfall + RUSLE A",
    c("M-RR0", "M-RR1"),
    "Markum rainfall + RUSLE sensitivity: add overall solar"
  )
)

comparison_all_panels <- bind_rows(
  comparison_panel_a,
  comparison_panel_b,
  comparison_panel_c,
  comparison_panel_d
)

comparison_notes <- tibble(
  note_id = 1:4,
  note = c(
    "Models were fitted by maximum likelihood.",
    "AICc differences and weights were calculated within each matched pair.",
    "Likelihood-ratio results appear on the solar-model row.",
    "The Markum RUSLE and NDVI model is included as a sensitivity comparison."
  )
)

stopifnot(
  nrow(comparison_panel_a) == 2,
  nrow(comparison_panel_b) == 2,
  nrow(comparison_panel_c) == 2,
  nrow(comparison_panel_d) == 8
)

comparison_path <- write_table_package(
  "model_comparisons",
  list(
    Panel_A_Otto_primary = comparison_panel_a,
    Panel_B_Otto_drainage = comparison_panel_b,
    Panel_C_Markum_primary = comparison_panel_c,
    Panel_D_other_sensitivity = comparison_panel_d,
    All_panels = comparison_all_panels,
    Notes = comparison_notes
  )
)

# 9. Build final coefficient tables -----------------------------------------
otto_primary_coefficients <- coefficient_table(
  otto_primary_REML,
  "A",
  "Otto",
  "Primary no-drainage REML model",
  "O-P1-REML",
  "sqrt(NTU)"
)

markum_primary_coefficients <- coefficient_table(
  markum_primary_REML,
  "B",
  "Markum",
  "Primary drainage-adjusted REML model",
  "M-P1-REML",
  "log10(NTU)"
)

otto_drainage_coefficients <- coefficient_table(
  otto_drainage_REML,
  "C",
  "Otto",
  "Drainage-adjusted REML sensitivity",
  "O-D1-REML",
  "sqrt(NTU)"
)

final_coefficients <- bind_rows(otto_primary_coefficients, markum_primary_coefficients, otto_drainage_coefficients)

final_model_summary <- bind_rows(
  model_summary(
    otto_primary_REML,
    "A",
    "Otto",
    "Primary no-drainage REML model",
    "O-P1-REML",
    "sqrt(NTU)"
  ),
  model_summary(
    markum_primary_REML,
    "B",
    "Markum",
    "Primary drainage-adjusted REML model",
    "M-P1-REML",
    "log10(NTU)"
  ),
  model_summary(
    otto_drainage_REML,
    "C",
    "Otto",
    "Drainage-adjusted REML sensitivity",
    "O-D1-REML",
    "sqrt(NTU)"
  )
)

coefficient_notes <- tibble(
  note_id = 1:4,
  note = c(
    "Coefficients and variance components were estimated by REML.",
    "Continuous coefficients represent a 1-SD predictor increase.",
    "Responses were transformed but not standardized.",
    "Confidence intervals use Satterthwaite degrees of freedom."
  )
)

coefficient_path <- write_table_package(
  "final_model_coefficients_and_summaries",
  list(
    Panel_A_Otto_primary = otto_primary_coefficients,
    Panel_B_Markum_primary = markum_primary_coefficients,
    Panel_C_Otto_drainage = otto_drainage_coefficients,
    Model_summary = final_model_summary,
    All_coefficients = final_coefficients,
    Notes = coefficient_notes
  )
)

# 10. Build within-between tables -------------------------------------------
within_between_model_comparison <- bind_rows(
  five_model_table(
    otto_solar_models,
    "A",
    "Otto",
    "Primary no-drainage adjustment",
    c("O-WB0", "O-WB1", "O-WBW", "O-WBB", "O-WBWB")
  ),
  five_model_table(
    markum_solar_models,
    "A",
    "Markum",
    "Primary drainage-adjusted adjustment",
    c("M-WB0", "M-WB1", "M-WBW", "M-WBB", "M-WBWB")
  )
)

within_between_lrt_rows <- function(site, models) {
  bind_rows(
    lrt_values(
      models$no_solar,
      models$within_solar,
      paste(site, "add within-station solar to no-solar model")
    ),
    lrt_values(
      models$no_solar,
      models$between_solar,
      paste(site, "add between-station solar to no-solar model")
    ),
    lrt_values(
      models$within_solar,
      models$within_between,
      paste(site, "add between-station solar to within model")
    ),
    lrt_values(
      models$between_solar,
      models$within_between,
      paste(site, "add within-station solar to between model")
    )
  ) %>%
    mutate(panel = "B", site = site, .before = 1)
}

within_between_lrt <- bind_rows(
  within_between_lrt_rows("Otto", otto_solar_models),
  within_between_lrt_rows("Markum", markum_solar_models)
)

component_coefficient_rows <- function(
    model,
    site,
    analysis_role,
    model_id,
    response_scale
) {
  model_frame <- model.frame(model)

  coefficient_table(
    model,
    "C",
    site,
    analysis_role,
    model_id,
    response_scale
  ) %>%
    filter(model_term %in% c("solar_within", "solar_between")) %>%
    mutate(
      predictor_increment = "1 SD; component standardized separately",
      observations = nrow(model_frame),
      station_levels = n_distinct(model_frame$station_uid),
      sampling_event_levels = n_distinct(model_frame$sampling_event_id)
    )
}

within_between_coefficients <- bind_rows(
  component_coefficient_rows(
    otto_within_between_REML,
    "Otto",
    "No-drainage within-between REML model",
    "O-WBWB-REML",
    "sqrt(NTU)"
  ),
  component_coefficient_rows(
    markum_within_between_REML,
    "Markum",
    "Drainage-adjusted within-between REML model",
    "M-WBWB-REML",
    "log10(NTU)"
  )
)

within_between_model_summary <- bind_rows(
  model_summary(
    otto_within_between_REML,
    "C",
    "Otto",
    "No-drainage within-between REML model",
    "O-WBWB-REML",
    "sqrt(NTU)"
  ),
  model_summary(
    markum_within_between_REML,
    "C",
    "Markum",
    "Drainage-adjusted within-between REML model",
    "M-WBWB-REML",
    "log10(NTU)"
  )
)

within_between_notes <- tibble(
  note_id = 1:4,
  note = c(
    "Model comparisons used ML.",
    "Coefficient estimates used the joint within-between REML model.",
    "Within-station solar represents change around each station mean.",
    "Between-station solar represents differences among station means."
  )
)

stopifnot(
  sum(within_between_model_comparison$site == "Otto") == 5,
  sum(within_between_model_comparison$site == "Markum") == 5,
  sum(within_between_lrt$site == "Otto") == 4,
  sum(within_between_lrt$site == "Markum") == 4,
  sum(within_between_coefficients$site == "Otto") == 2,
  sum(within_between_coefficients$site == "Markum") == 2
)

within_between_path <- write_table_package(
  "within_between_analysis",
  list(
    Panel_A_five_models = within_between_model_comparison,
    Panel_B_nested_LRT = within_between_lrt,
    Panel_C_REML_components = within_between_coefficients,
    Panel_C_model_summary = within_between_model_summary,
    Notes = within_between_notes
  )
)

# 11. Build sensitivity tables ----------------------------------------------
additional_model_comparisons <- bind_rows(
  pair_panel_table(
    otto_core_pair,
    "A",
    "Otto",
    "Additional matched sensitivity",
    "Drainage + rainfall + NDVI",
    c("O-CORE0", "O-CORE1"),
    "Otto core sensitivity: add overall solar"
  ),
  pair_panel_table(
    otto_compact_land_pair,
    "A",
    "Otto",
    "Additional matched sensitivity",
    "Compact land-management adjustment",
    c("O-LAND0", "O-LAND1"),
    "Otto compact land-management sensitivity: add overall solar"
  ),
  pair_panel_table(
    otto_temperature_pair,
    "A",
    "Otto",
    "Additional matched sensitivity",
    "Core + maximum temperature",
    c("O-TEMP0", "O-TEMP1"),
    "Otto temperature sensitivity: add overall solar"
  ),
  pair_panel_table(
    otto_reduced_vegetation_pair,
    "A",
    "Otto",
    "Additional matched sensitivity",
    "Reduced vegetation adjustment",
    c("O-VEG0", "O-VEG1"),
    "Otto reduced vegetation sensitivity: add overall solar"
  ),
  pair_panel_table(
    markum_tree_pair,
    "A",
    "Markum",
    "Additional matched sensitivity",
    "Drainage + rainfall + tree cover",
    c("M-TREE0", "M-TREE1"),
    "Markum tree-cover sensitivity: add overall solar"
  ),
  pair_panel_table(
    markum_temperature_pair,
    "A",
    "Markum",
    "Additional matched sensitivity",
    "Primary adjustment + maximum temperature",
    c("M-TEMP0", "M-TEMP1"),
    "Markum temperature sensitivity: add overall solar"
  )
)

sensitivity_registry <- list(
  list(
    site = "Otto",
    family = "Drainage-adjusted full model",
    id = "O-D1-REML",
    scale = "sqrt(NTU)",
    ml_model = otto_drainage_pair$solar
  ),
  list(
    site = "Otto",
    family = "Drainage + rainfall + RUSLE A",
    id = "O-RUSLE1-REML",
    scale = "sqrt(NTU)",
    ml_model = otto_rusle_pair$solar
  ),
  list(
    site = "Otto",
    family = "Drainage + rainfall + NDVI + season",
    id = "O-SEAS1-REML",
    scale = "sqrt(NTU)",
    ml_model = otto_season_pair$solar
  ),
  list(
    site = "Otto",
    family = "Drainage + rainfall + NDVI",
    id = "O-CORE1-REML",
    scale = "sqrt(NTU)",
    ml_model = otto_core_pair$solar
  ),
  list(
    site = "Otto",
    family = "Compact land-management adjustment",
    id = "O-LAND1-REML",
    scale = "sqrt(NTU)",
    ml_model = otto_compact_land_pair$solar
  ),
  list(
    site = "Otto",
    family = "Core + maximum temperature",
    id = "O-TEMP1-REML",
    scale = "sqrt(NTU)",
    ml_model = otto_temperature_pair$solar
  ),
  list(
    site = "Otto",
    family = "Reduced vegetation adjustment",
    id = "O-VEG1-REML",
    scale = "sqrt(NTU)",
    ml_model = otto_reduced_vegetation_pair$solar
  ),
  list(
    site = "Markum",
    family = "RUSLE A + NDVI sensitivity",
    id = "M-RN1-REML",
    scale = "log10(NTU)",
    ml_model = markum_rusle_ndvi_pair$solar
  ),
  list(
    site = "Markum",
    family = "Rainfall + RUSLE A",
    id = "M-RR1-REML",
    scale = "log10(NTU)",
    ml_model = markum_rusle_rain_pair$solar
  ),
  list(
    site = "Markum",
    family = "Drainage + rainfall + tree cover",
    id = "M-TREE1-REML",
    scale = "log10(NTU)",
    ml_model = markum_tree_pair$solar
  ),
  list(
    site = "Markum",
    family = "Primary adjustment + maximum temperature",
    id = "M-TEMP1-REML",
    scale = "log10(NTU)",
    ml_model = markum_temperature_pair$solar
  )
)

sensitivity_REML_models <- lapply(
  sensitivity_registry,
  function(entry) {
    model_data <- if (entry$site == "Otto") otto else markum
    refit_reml(entry$ml_model, model_data)
  }
)

sensitivity_solar_coefficients <- bind_rows(lapply(seq_along(sensitivity_registry), function(i) {
  entry <- sensitivity_registry[[i]]
  coefficient_table(
    sensitivity_REML_models[[i]],
    "B",
    entry$site,
    entry$family,
    entry$id,
    entry$scale
  ) %>%
    filter(model_term == "solar_panel_pct_z") %>%
    mutate(sensitivity_family = entry$family, .before = model_id)
}))

sensitivity_model_summaries <- bind_rows(lapply(seq_along(sensitivity_registry), function(i) {
  entry <- sensitivity_registry[[i]]
  model_summary(
    sensitivity_REML_models[[i]],
    "C",
    entry$site,
    entry$family,
    entry$id,
    entry$scale
  ) %>%
    mutate(sensitivity_family = entry$family, .before = model_id)
}))

sensitivity_notes <- tibble(
  note_id = 1:3,
  note = c(
    "Panel A contains matched ML comparisons.",
    "Panel B contains REML solar coefficients.",
    "Panel C contains REML model summaries and fit checks."
  )
)

sensitivity_path <- write_table_package(
  "sensitivity_analyses",
  list(
    Panel_A_additional_ML = additional_model_comparisons,
    Panel_B_solar_coefficients = sensitivity_solar_coefficients,
    Panel_C_model_summaries = sensitivity_model_summaries,
    Notes = sensitivity_notes
  )
)

# 12. Save one combined workbook --------------------------------------------
combined_table_path <- file.path(table_dir, "combined_model_results.xlsx")

writexl::write_xlsx(
  list(
    model_comparisons = comparison_all_panels,
    comparison_notes = comparison_notes,
    final_coefficients = final_coefficients,
    final_model_summary = final_model_summary,
    coefficient_notes = coefficient_notes,
    within_between_models = within_between_model_comparison,
    within_between_LRT = within_between_lrt,
    within_between_coefficients = within_between_coefficients,
    within_between_summary = within_between_model_summary,
    within_between_notes = within_between_notes,
    additional_comparisons = additional_model_comparisons,
    sensitivity_coefficients = sensitivity_solar_coefficients,
    sensitivity_summaries = sensitivity_model_summaries,
    sensitivity_notes = sensitivity_notes
  ),
  combined_table_path
)

# 13. Create Figure S5 coefficient plot -------------------------------------
si_theme <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = 13, margin = margin(b = 5)),
      axis.title = element_text(size = 11.5, color = "#1B2630"),
      axis.text = element_text(size = 10.5, color = "#1B2630"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.border = element_rect(linewidth = 0.65, color = "#53616E"),
      plot.margin = margin(10, 14, 10, 10)
    )
}

make_coefficient_panel <- function(coefficient_data, title, x_label, point_color) {
  plot_data <- coefficient_data %>%
    filter(model_term != "(Intercept)") %>%
    mutate(predictor = factor(predictor, levels = rev(unique(predictor))))

  ggplot(plot_data, aes(x = estimate, y = predictor)) +
    geom_vline(xintercept = 0, linewidth = 0.65, linetype = "dashed",
               color = "#59636C") +
    geom_segment(aes(x = CI_95_low, xend = CI_95_high, yend = predictor),
                 linewidth = 1.3, color = point_color, lineend = "round") +
    geom_point(size = 3.8, color = point_color) +
    scale_x_continuous(expand = expansion(mult = c(0.10, 0.10))) +
    labs(title = title, x = x_label, y = NULL) +
    si_theme(base_size = 12)
}

coefficient_panel_otto <- make_coefficient_panel(
  otto_primary_coefficients, "Otto primary model",
  "Fixed-effect estimate on sqrt(NTU) scale\n(per 1-SD predictor increase)", "#155A9C"
)
coefficient_panel_markum <- make_coefficient_panel(
  markum_primary_coefficients, "Markum primary model",
  "Fixed-effect estimate on log10(NTU) scale\n(per 1-SD predictor increase)", "#B53A32"
)

figure_s5 <- coefficient_panel_otto / coefficient_panel_markum +
  patchwork::plot_layout(heights = c(1.45, 1))

figure_s5_png <- file.path(figure_dir, "Figure_S5_primary_model_coefficients.png")
figure_s5_tiff <- file.path(figure_dir, "Figure_S5_primary_model_coefficients.tiff")
figure_s5_pdf <- file.path(figure_dir, "Figure_S5_primary_model_coefficients.pdf")

ggsave(figure_s5_png, figure_s5, width = 7.5, height = 8.8, units = "in",
       dpi = 600, bg = "white")
ggsave(figure_s5_tiff, figure_s5, width = 7.5, height = 8.8,
       units = "in", dpi = 600, compression = "lzw", bg = "white")
ggsave(figure_s5_pdf, figure_s5, width = 7.5, height = 8.8,
       units = "in", bg = "white")

readr::write_csv(
  bind_rows(otto_primary_coefficients, markum_primary_coefficients) %>% filter(model_term != "(Intercept)"),
  file.path(figure_dir, "Figure_S5_primary_REML_coefficients_source_data.csv")
)

# 14. Create Figures S7A and S7B diagnostics --------------------------------
diagnostic_data <- function(model, site, response_scale) {
  model_frame <- model.frame(model)
  response <- model.response(model_frame)
  fitted_values <- fitted(model)
  raw_residuals <- residuals(model)
  standardized_residuals <- raw_residuals / sigma(model)

  tibble(
    site = site,
    response_scale = response_scale,
    observed = response,
    fitted = fitted_values,
    raw_residual = raw_residuals,
    standardized_residual = standardized_residuals
  )
}

make_diagnostic_plots <- function(data, site, response_scale, labels, color) {
  fitted_label <- paste0("Conditional fitted value [", response_scale, "]")
  observed_label <- paste0("Observed value [", response_scale, "]")
  common_points <- geom_point(alpha = 0.78, size = 2.2, color = color)

  residual_plot <- ggplot(data, aes(x = fitted, y = standardized_residual)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.65,
               color = "#59636C") +
    common_points +
    geom_smooth(method = "loess", formula = y ~ x, span = 1,
                se = FALSE, linewidth = 0.85, color = "#B23A48") +
    labs(title = labels[[1]], x = fitted_label,
         y = "Residual / model residual SD") +
    si_theme()

  qq_plot <- ggplot(data, aes(sample = standardized_residual)) +
    stat_qq_line(linewidth = 0.7, linetype = "dashed", color = "#59636C") +
    stat_qq(alpha = 0.78, size = 2.2, color = color) +
    labs(title = labels[[2]], x = "Theoretical normal quantile",
         y = "Residual / model residual SD") +
    si_theme()

  observed_plot <- ggplot(data, aes(x = fitted, y = observed)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                linewidth = 0.7, color = "#59636C") +
    common_points +
    coord_cartesian() +
    labs(title = labels[[3]], x = fitted_label, y = observed_label) +
    si_theme()

  list(residual_plot, qq_plot, observed_plot)
}

otto_diagnostics_data <- diagnostic_data(otto_primary_REML, "Otto", "sqrt(NTU)")
markum_diagnostics_data <- diagnostic_data(markum_primary_REML, "Markum", "log10(NTU)")

otto_diagnostic_plots <- make_diagnostic_plots(
  otto_diagnostics_data, "Otto", "sqrt(NTU)",
  c("Otto: residuals versus fitted", "Otto: normal Q-Q",
    "Otto: observed versus fitted"), "#1E5079"
)
markum_diagnostic_plots <- make_diagnostic_plots(
  markum_diagnostics_data, "Markum", "log10(NTU)",
  c("Markum: residuals versus fitted", "Markum: normal Q-Q",
    "Markum: observed versus fitted"), "#1E5079"
)

diagnostics_otto <- patchwork::wrap_plots(otto_diagnostic_plots, ncol = 1)
diagnostics_markum <- patchwork::wrap_plots(markum_diagnostic_plots, ncol = 1)

save_diagnostic_page <- function(plot, stem) {
  paths <- setNames(file.path(figure_dir, paste0(stem, c(".png", ".tiff", ".pdf"))),
                    c("png", "tiff", "pdf"))
  ggsave(paths[["png"]], plot, width = 7.5, height = 9.3, units = "in",
         dpi = 600, bg = "white")
  ggsave(paths[["tiff"]], plot, width = 7.5, height = 9.3, units = "in",
         dpi = 600, compression = "lzw", bg = "white")
  ggsave(paths[["pdf"]], plot, width = 7.5, height = 9.3,
         units = "in", bg = "white")
  paths
}

diagnostics_otto_files <- save_diagnostic_page(
  diagnostics_otto, "Figure_S7A_Otto_primary_model_diagnostics")
diagnostics_markum_files <- save_diagnostic_page(
  diagnostics_markum, "Figure_S7B_Markum_primary_model_diagnostics")

readr::write_csv(
  bind_rows(otto_diagnostics_data, markum_diagnostics_data),
  file.path(figure_dir, "Figure_S7_primary_model_diagnostics_source_data.csv")
)

# 15. Save models and run information ---------------------------------------
model_objects <- list(
  otto_primary_pair_ML = otto_primary_pair,
  otto_drainage_pair_ML = otto_drainage_pair,
  markum_primary_pair_ML = markum_primary_pair,
  otto_solar_formulations_ML = otto_solar_models,
  markum_solar_formulations_ML = markum_solar_models,
  otto_primary_REML = otto_primary_REML,
  markum_primary_REML = markum_primary_REML,
  otto_drainage_REML = otto_drainage_REML,
  otto_within_between_REML = otto_within_between_REML,
  markum_within_between_REML = markum_within_between_REML,
  matched_sensitivity_pairs = list(
    otto_rusle = otto_rusle_pair,
    otto_season = otto_season_pair,
    markum_rusle_ndvi = markum_rusle_ndvi_pair,
    markum_rusle_rain = markum_rusle_rain_pair
  ),
  additional_sensitivity_pairs = list(
    otto_core = otto_core_pair,
    otto_compact_land = otto_compact_land_pair,
    otto_temperature = otto_temperature_pair,
    otto_reduced_vegetation = otto_reduced_vegetation_pair,
    markum_tree = markum_tree_pair,
    markum_temperature = markum_temperature_pair
  ),
  sensitivity_REML_models = sensitivity_REML_models
)

saveRDS(model_objects, file.path(model_dir, "combined_final_model_objects.rds"))

writeLines(
  capture.output(summary(otto_primary_REML)),
  file.path(log_dir, "Otto_primary_no_drainage_REML_summary.txt")
)
writeLines(
  capture.output(summary(markum_primary_REML)),
  file.path(log_dir, "Markum_primary_drainage_REML_summary.txt")
)
writeLines(
  capture.output(summary(otto_drainage_REML)),
  file.path(log_dir, "Otto_drainage_sensitivity_REML_summary.txt")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(log_dir, "R_sessionInfo.txt")
)

autocorrelation_checks <- list(
  Otto_primary = tryCatch(
    performance::check_autocorrelation(otto_primary_REML),
    error = function(e) e
  ),
  Markum_primary = tryCatch(
    performance::check_autocorrelation(markum_primary_REML),
    error = function(e) e
  )
)

saveRDS(
  autocorrelation_checks,
  file.path(model_dir, "primary_model_autocorrelation_checks.rds")
)
writeLines(
  capture.output(print(autocorrelation_checks$Otto_primary)),
  file.path(log_dir, "Otto_primary_autocorrelation_check.txt")
)
writeLines(
  capture.output(print(autocorrelation_checks$Markum_primary)),
  file.path(log_dir, "Markum_primary_autocorrelation_check.txt")
)

analysis_decisions <- c(
  "ML was used for model comparisons; REML was used for coefficient estimates.",
  "The Otto primary model did not include drainage area.",
  "The Otto drainage model is a sensitivity analysis.",
  "The Markum primary model includes drainage area.",
  "Figure S5 excludes model intercepts.",
  "Akaike weights were calculated within each model set.",
  "Singular fits are retained and flagged.",
  "The Markum model contains 27 sampling-event levels."
)
writeLines(analysis_decisions, file.path(output_dir, "ANALYSIS_DECISIONS.txt"))

# 16. Check outputs and write manifest --------------------------------------
expected_key_files <- c(
  comparison_path,
  coefficient_path,
  within_between_path,
  sensitivity_path,
  combined_table_path,
  figure_s5_png,
  figure_s5_tiff,
  figure_s5_pdf,
  unname(diagnostics_otto_files),
  unname(diagnostics_markum_files),
  file.path(model_dir, "combined_final_model_objects.rds"),
  file.path(output_dir, "ANALYSIS_DECISIONS.txt")
)

if (!all(file.exists(expected_key_files))) {
  stop(
    "One or more required outputs were not created:\n",
    paste(expected_key_files[!file.exists(expected_key_files)], collapse = "\n")
  )
}

manifest <- tibble(
  output_type = c(
    "Model comparisons", "Final coefficients", "Within-between analysis",
    "Sensitivity analyses", "Combined workbook", "Figure S5 PNG",
    "Figure S5 TIFF", "Figure S5 PDF", "Figure S7A Otto PNG",
    "Figure S7A Otto TIFF", "Figure S7A Otto PDF", "Figure S7B Markum PNG",
    "Figure S7B Markum TIFF", "Figure S7B Markum PDF", "Model archive",
    "Analysis summary"
  ),
  file = expected_key_files,
  status = "Created"
)

readr::write_csv(manifest, file.path(output_dir, "OUTPUT_MANIFEST.csv"))

message("Analysis complete.")
message("Output folder: ", output_dir)
