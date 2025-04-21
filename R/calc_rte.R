#' Calculate Rank-based Triage Effectiveness (RTE)
#'
#' @description
#' Calculates Rank-based Triage Effectiveness (RTE) metrics, which evaluate triage performance
#' based on patient queue positions rather than waiting times. RTE provides a patient-level
#' metric that complements the waiting time-based Triage Effectiveness (WTE) approach.
#'
#' @param df Data frame containing patient data (must be initialized through init())
#' @param subgroup Optional list for subgroup analysis - an alternative to var1/2
#' @param var1 Optional string for first comparison variable
#' @param var2 Optional string for second comparison variable
#' @param bootstrap Logical. Whether to perform bootstrap calculations. Default is FALSE.
#' @param bootstrap_params List of bootstrap parameters:
#'   \itemize{
#'     \item sample_percentage (default: 1) How large of the sample is used? 1 = full sample.
#'     \item n_iterations (default: 2000) How many bootstrap iterations to use
#'     \item distribution_span (default: 0.95) What is the confidence interval used?
#'   }
#' @param n_workers Number of workers for parallel processing (default: detectCores() - 1)
#' @param verbose Logical. If TRUE, returns detailed results with intermediate calculation
#'   fields from calculate_queue_metrics. Default is FALSE.
#' @param overall_only Logical. If TRUE, only overall metrics are returned (default: FALSE)
#' @param min_loset_warning Numerical, when should the function warn for low loset prevalence? (default: 5)
#' @param include_distributions Logical. If TRUE, returns patient-level RTE values for
#'   further analysis and bootstrap distrobutions. Default is FALSE.
#' @param check_convergence Logical, Default is TRUE. Will check convergence and produce convergence plots if true + bootstrap == TRUE.
#' @param seed Seed for reproducible bootstrapping (default: NULL - no seed)
#' @param quite Logical. Default is TRUE. If false outputs the recommended/used confidence interval choice.
#'
#' @return A list containing:
#'   \itemize{
#'     \item results: Tibble with RTE metrics
#'     \item metadata: List containing calculation parameters
#'     \item distributions: (If include_distributions=TRUE) Patient-level RTE values
#'   }
#'
#' @importFrom dplyr filter select mutate arrange left_join bind_rows group_by
#' @importFrom dplyr summarise ungroup across group_vars group_modify n rename first
#' @importFrom tibble as_tibble
#' @importFrom data.table as.data.table
#' @importFrom stats sd var shapiro.test qt quantile median weighted.mean
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline labs theme_minimal annotate
#' @importFrom rlang sym !!
#' @importFrom furrr future_map furrr_options
#' @importFrom progressr with_progress progressor
#'
#' @export
calc_rte <- function(df,
                     subgroup = NULL,
                     var1 = NULL,
                     var2 = NULL,
                     bootstrap = FALSE,
                     bootstrap_params = list(
                       sample_percentage = 1,
                       n_iterations = 2000,
                       distribution_span = 0.95
                     ),
                     n_workers = parallel::detectCores() - 1,
                     verbose = FALSE,
                     overall_only = FALSE,
                     min_loset_warning = 5,
                     include_distributions = TRUE,
                     check_convergence = TRUE,
                     seed = NULL,
                     quite = TRUE) {

  # Validate input data
  validate_te_data(df, subgroup, var1, var2)  # Reuse from calc_wte

  n_units <- length(unique(df$unit))
  if (n_units == 1 && overall_only) {
    # Get the unit name for the message
    unit_name <- unique(df$unit)[1]
    # Override overall_only setting
    overall_only <- FALSE
    # Inform the user
    message(sprintf("Only one unit ('%s') detected. Setting overall_only=FALSE as overall statistics are equivalent to unit statistics.", unit_name))
  }


  # Check if RTE metrics already exist
  if (!any(c("observed_RTE") %in% names(df))) {
    # Calculate RTE metrics
    df <- calculate_queue_metrics(df, verbose = verbose)
  }

  # Create grouping and aggregate results
  grouping <- create_grouping(df, subgroup, var1, var2)  # Reuse from calc_wte
  unit_results <- calculate_rte_unit_metrics(grouping)

  # Only calculate overall if more than one unit or subgroup is not used
  if (length(unique(unit_results$unit)) > 1 & is.null(subgroup)) {
    # Create overall grouping for overall calculations
    overall_grouping <- create_overall_grouping(unit_results, subgroup, var1, var2)

    # Calculate overall metrics using the grouping
    overall_results <- calculate_rte_overall_metrics(unit_results, df, overall_grouping)
    results <- bind_rows(overall_results, unit_results)
  } else {
    results <- unit_results
  }

  results <- results %>%
    ungroup()

  # Determine if non-normal distributions detected from direct calculation
  non_normal_detected <- validate_rte_results(results, min_loset_warning, quite)

  # Perform bootstrap calculations if requested
  bootstrap_results <- NULL
  if (bootstrap) {
    bootstrap_results <- calculate_rte_bootstrap(df, subgroup, var1, var2,
                                                 bootstrap_params, n_workers, seed)

    # Combine bootstrap results with direct calculation results
    results <- combine_rte_results(results, bootstrap_results$metrics)
  }

  # Create patient-level distributions if requested
  distributions_data <- NULL
  if (include_distributions) {
    # Create a base selection of columns that always exist
    select_cols <- c("unit", "id", "observed_RTE", "valid", "loset")

    # Add grouping variables if they exist
    if (!is.null(var1)) select_cols <- c(select_cols, var1)
    if (!is.null(var2)) select_cols <- c(select_cols, var2)

    # Add theoretical RTE columns only if they exist
    if ("theoretical_RTE" %in% names(df)) select_cols <- c(select_cols, "theoretical_RTE")
    if ("binary_theoretical_RTE" %in% names(df)) select_cols <- c(select_cols, "binary_theoretical_RTE")

    # Now perform the selection with only columns that exist
    distributions_data <- df %>%
      filter(loset == TRUE, valid == TRUE) %>%
      select(all_of(select_cols))
  }

  # Filter results for overall_only if requested
  if (overall_only) {
    results <- results %>%
      filter(unit == "overall")

    # ALSO filter the distributions to only include overall
    if (!is.null(distributions_data)) {
      distributions_data <- distributions_data %>%
        filter(unit == "overall")
    }

    if (!is.null(bootstrap_results) && !is.null(bootstrap_results$distributions)) {
      bootstrap_results$distributions <- bootstrap_results$distributions %>%
        filter(unit == "overall")
    }
  }

  # Check if any unit recommends nonparametric CIs
  any_nonparametric <- any(results$recommended_ci_type == "nonparametric", na.rm = TRUE)
  recommended_ci_vars <- if(bootstrap) {
    "boot"
  } else if(any_nonparametric) {
    "nonparametric"
  } else {
    "parametric"
  }

  # Create metadata
  metadata <- list(
    calculation_time = Sys.time(),
    calculation_method = "rank-based",
    group_var1 = var1,
    group_var2 = var2,
    subgr = subgroup,
    ci_method = if(bootstrap) "bootstrap" else "direct",
    recommended_ci_vars = recommended_ci_vars,
    bootstrap_params = if(bootstrap) bootstrap_params else NULL
  )

  output <- list(
    results = results,
    metadata = metadata
  )

  # Add bootstrap distributions if available
  if (bootstrap && !is.null(bootstrap_results$distributions)) {
    output$bootstrap_distributions <- bootstrap_results$distributions
  }

  # Add individual patient distribution data if requested
  if (include_distributions && !is.null(distributions_data)) {
    output$distributions <- distributions_data
  }

  # Add convergence analysis if requested and bootstrap was performed
  if (bootstrap && check_convergence) {
    output$convergence <- analyze_bootstrap_convergence(output,
                                                        subgroup, var1, var2)
  }

  # Set class for custom print method
  class(output) <- c("calc_te", class(output))

  return(output)
}


#' Calculate Unit-Level RTE Metrics
#'
#' @param grouping Grouped data frame including unit and comparison variables
#'
#' @return Tibble with unit-level RTE metrics
#'
#' @details
#' Calculates unit-level metrics including:
#' - Patient counts and valid RTE counts
#' - Mean RTE values for observed, theoretical, and binary theoretical
#' - Parametric and non-parametric confidence intervals
#' - Normality test results
#'
#' @keywords internal
calculate_rte_unit_metrics <- function(grouping) {
  # Extract the grouping variables
  group_vars <- dplyr::group_vars(grouping)

  # Calculate all metrics in one step, maintaining the grouping structure
  result <- grouping %>%
    # First calculate counts and sensitivity/specificity for all patients
    summarise(
      # Patient counts
      n_patients = n(),
      n_patients_loset = sum(loset, na.rm = TRUE),
      loset_prevalence = mean(loset, na.rm = TRUE),

      # Only calculate sensitivity/specificity if we have LOSET cases and non-LOSET cases
      sensitivity = if(sum(loset, na.rm = TRUE) > 0)
        sum(tp, na.rm = TRUE) / (sum(tp, na.rm = TRUE) + sum(fn, na.rm = TRUE))
      else NA_real_,
      specificity = if(sum(!loset, na.rm = TRUE) > 0)
        sum(tn, na.rm = TRUE) / (sum(tn, na.rm = TRUE) + sum(fp, na.rm = TRUE))
      else NA_real_,

      # Then calculate RTE metrics using only valid LOSET patients
      n_valid_rte = sum(loset == TRUE & valid == TRUE, na.rm = TRUE),

      # Temporary containers for calculations
      .valid_observed_RTE = if(sum(loset == TRUE & valid == TRUE) > 0)
        list(observed_RTE[loset == TRUE & valid == TRUE]) else list(NULL),
      .valid_theoretical_RTE = if("theoretical_RTE" %in% names(.) && sum(loset == TRUE & valid == TRUE) > 0)
        list(theoretical_RTE[loset == TRUE & valid == TRUE]) else list(NULL),
      .valid_binary_theoretical_RTE = if("binary_theoretical_RTE" %in% names(.) && sum(loset == TRUE & valid == TRUE) > 0)
        list(binary_theoretical_RTE[loset == TRUE & valid == TRUE]) else list(NULL),

      .groups = 'keep'
    ) %>%
    # Now use these containers to calculate the actual metrics
    mutate(
      # OTE calculations - only if we have valid data
      ote_te = sapply(.valid_observed_RTE, function(x) if(!is.null(x) && length(x) > 0) mean(x, na.rm = TRUE) else NA_real_),
      parametric_ote_sd = sapply(.valid_observed_RTE, function(x) if(!is.null(x) && length(x) > 0) sd(x, na.rm = TRUE) else NA_real_),
      ote_se = if(n_valid_rte > 1) parametric_ote_sd / sqrt(n_valid_rte) else NA_real_,
      t_value = if(n_valid_rte > 1) qt(0.975, df = n_valid_rte - 1) else NA_real_,
      parametric_ote_mean = ote_te,
      parametric_ote_sd_q = sapply(.valid_observed_RTE, function(x) {
        if(!is.null(x) && length(x) > 3) {
          (quantile(x, 0.75, na.rm = TRUE) - quantile(x, 0.25, na.rm = TRUE)) / 1.349
        } else NA_real_
      }),
      parametric_ote_var_lower = if(n_valid_rte > 1) ote_te - t_value * ote_se else NA_real_,
      parametric_ote_var_upper = if(n_valid_rte > 1) ote_te + t_value * ote_se else NA_real_,
      nonparametric_ote_var_lower = sapply(.valid_observed_RTE, function(x) {
        if(!is.null(x) && length(x) > 3) quantile(x, 0.025, na.rm = TRUE) else NA_real_
      }),
      nonparametric_ote_var_upper = sapply(.valid_observed_RTE, function(x) {
        if(!is.null(x) && length(x) > 3) quantile(x, 0.975, na.rm = TRUE) else NA_real_
      }),
      ote_shapiro_p = sapply(.valid_observed_RTE, function(x) {
        if(!is.null(x) && length(x) > 2 && length(x) <= 5000) {
          tryCatch(shapiro.test(x)$p.value, error = function(e) NA_real_)
        } else NA_real_
      }),

      # TTE calculations - with similar checks
      tte_te = sapply(.valid_theoretical_RTE, function(x) if(!is.null(x) && length(x) > 0) mean(x, na.rm = TRUE) else NA_real_),
      parametric_tte_sd = sapply(.valid_theoretical_RTE, function(x) if(!is.null(x) && length(x) > 0) sd(x, na.rm = TRUE) else NA_real_),
      tte_se = if(n_valid_rte > 1 && !is.na(parametric_tte_sd)) parametric_tte_sd / sqrt(n_valid_rte) else NA_real_,
      parametric_tte_mean = tte_te,
      parametric_tte_sd_q = sapply(.valid_theoretical_RTE, function(x) {
        if(!is.null(x) && length(x) > 3) {
          (quantile(x, 0.75, na.rm = TRUE) - quantile(x, 0.25, na.rm = TRUE)) / 1.349
        } else NA_real_
      }),
      parametric_tte_var_lower = if(n_valid_rte > 1 && !is.na(tte_se)) tte_te - t_value * tte_se else NA_real_,
      parametric_tte_var_upper = if(n_valid_rte > 1 && !is.na(tte_se)) tte_te + t_value * tte_se else NA_real_,
      nonparametric_tte_var_lower = sapply(.valid_theoretical_RTE, function(x) {
        if(!is.null(x) && length(x) > 3) quantile(x, 0.025, na.rm = TRUE) else NA_real_
      }),
      nonparametric_tte_var_upper = sapply(.valid_theoretical_RTE, function(x) {
        if(!is.null(x) && length(x) > 3) quantile(x, 0.975, na.rm = TRUE) else NA_real_
      }),
      tte_shapiro_p = sapply(.valid_theoretical_RTE, function(x) {
        if(!is.null(x) && length(x) > 2 && length(x) <= 5000) {
          tryCatch(shapiro.test(x)$p.value, error = function(e) NA_real_)
        } else NA_real_
      }),

      # BTTE calculations - with similar checks
      btte_te = sapply(.valid_binary_theoretical_RTE, function(x) if(!is.null(x) && length(x) > 0) mean(x, na.rm = TRUE) else NA_real_),
      parametric_btte_sd = sapply(.valid_binary_theoretical_RTE, function(x) if(!is.null(x) && length(x) > 0) sd(x, na.rm = TRUE) else NA_real_),
      btte_se = if(n_valid_rte > 1 && !is.na(parametric_btte_sd)) parametric_btte_sd / sqrt(n_valid_rte) else NA_real_,
      parametric_btte_mean = btte_te,
      parametric_btte_sd_q = sapply(.valid_binary_theoretical_RTE, function(x) {
        if(!is.null(x) && length(x) > 3) {
          (quantile(x, 0.75, na.rm = TRUE) - quantile(x, 0.25, na.rm = TRUE)) / 1.349
        } else NA_real_
      }),
      parametric_btte_var_lower = if(n_valid_rte > 1 && !is.na(btte_se)) btte_te - t_value * btte_se else NA_real_,
      parametric_btte_var_upper = if(n_valid_rte > 1 && !is.na(btte_se)) btte_te + t_value * btte_se else NA_real_,
      nonparametric_btte_var_lower = sapply(.valid_binary_theoretical_RTE, function(x) {
        if(!is.null(x) && length(x) > 3) quantile(x, 0.025, na.rm = TRUE) else NA_real_
      }),
      nonparametric_btte_var_upper = sapply(.valid_binary_theoretical_RTE, function(x) {
        if(!is.null(x) && length(x) > 3) quantile(x, 0.975, na.rm = TRUE) else NA_real_
      }),
      btte_shapiro_p = sapply(.valid_binary_theoretical_RTE, function(x) {
        if(!is.null(x) && length(x) > 2 && length(x) <= 5000) {
          tryCatch(shapiro.test(x)$p.value, error = function(e) NA_real_)
        } else NA_real_
      }),

      # Calculate OTG
      otg_te = ote_te - tte_te,

      # Determine recommended CI prefix based on normality tests
      recommended_ci_type = ifelse(!is.na(ote_shapiro_p) && ote_shapiro_p < 0.05,
                                   "nonparametric", "parametric")
    ) %>%
    # Remove temporary list columns
    select(-.valid_observed_RTE, -.valid_theoretical_RTE, -.valid_binary_theoretical_RTE)

  return(result)
}

#' Calculate RTE Overall Metrics
#'
#' @param unit_results Results from unit-level calculations
#' @param df Original patient data frame with RTE values
#' @param overall_grouping Grouped data frame for overall calculations
#'
#' @return Tibble with overall RTE metrics
#'
#' @details
#' Calculates overall metrics by combining unit results with appropriate weighting.
#' For bootstrap results, uses patient-weighted means for the most accurate estimate.
#'
#' @keywords internal
calculate_rte_overall_metrics <- function(unit_results, df, overall_grouping = NULL) {
  # If no overall_grouping is provided, calculate a single overall
  if (is.null(overall_grouping)) {
    # Standard all-in-one overall calculation
    return(calculate_single_overall(df))
  }

  # Process each group in the overall_grouping
  overall_results <- overall_grouping %>%
    group_modify(~{
      # Get the current group from the group data
      current_group <- .y

      # For each grouping variable in the current group, filter df accordingly
      filtered_df <- df
      for (var_name in names(current_group)) {
        var_value <- current_group[[var_name]]
        filtered_df <- filtered_df %>%
          filter(.data[[var_name]] == var_value)
      }

      # Calculate overall metrics for this filtered data
      result <- calculate_single_overall(filtered_df)

      # Return result with "overall" as unit name
      result
    })

  return(overall_results)
}

# Helper function to calculate a single overall result
calculate_single_overall <- function(df) {
  # Extract only valid RTE values
  valid_df <- df %>%
    filter(loset == TRUE, valid == TRUE)

  # Extract RTE values
  ote_values <- valid_df$observed_RTE[!is.na(valid_df$observed_RTE)]
  tte_values <- if("theoretical_RTE" %in% names(valid_df)) {
    valid_df$theoretical_RTE[!is.na(valid_df$theoretical_RTE)]
  } else {
    NULL
  }
  btte_values <- if("binary_theoretical_RTE" %in% names(valid_df)) {
    valid_df$binary_theoretical_RTE[!is.na(valid_df$binary_theoretical_RTE)]
  } else {
    NULL
  }

  # Run Shapiro-Wilk on pooled data - but only if sample size is appropriate
  # Shapiro-Wilk is limited to samples between 3 and 5000
  ote_shapiro <- if(length(ote_values) > 2 && length(ote_values) <= 5000) {
    # Use sample if too large for full test
    if(length(ote_values) == 5000) {
      test_sample <- sample(ote_values, 5000)
    } else {
      test_sample <- ote_values
    }
    tryCatch(shapiro.test(test_sample)$p.value, error = function(e) NA_real_)
  } else {
    NA_real_
  }

  tte_shapiro <- if(!is.null(tte_values) && length(tte_values) > 2 && length(tte_values) <= 5000) {
    # Use sample if too large for full test
    if(length(tte_values) == 5000) {
      test_sample <- sample(tte_values, 5000)
    } else {
      test_sample <- tte_values
    }
    tryCatch(shapiro.test(test_sample)$p.value, error = function(e) NA_real_)
  } else {
    NA_real_
  }

  btte_shapiro <- if(!is.null(btte_values) && length(btte_values) > 2 && length(btte_values) <= 5000) {
    # Use sample if too large for full test
    if(length(btte_values) == 5000) {
      test_sample <- sample(btte_values, 5000)
    } else {
      test_sample <- btte_values
    }
    tryCatch(shapiro.test(test_sample)$p.value, error = function(e) NA_real_)
  } else {
    NA_real_
  }

  # Determine CI type based on direct tests
  # If any test indicates non-normality or if tests are NA, default to nonparametric
  recommended_ci_type <- if((!is.na(ote_shapiro) && ote_shapiro < 0.05) ||
                            (!is.null(tte_values) && !is.na(tte_shapiro) && tte_shapiro < 0.05)) {
    "nonparametric"
  } else if(is.na(ote_shapiro) || (!is.null(tte_values) && is.na(tte_shapiro))) {
    # When shapiro tests are NA, be conservative and use nonparametric
    "nonparametric"
  } else {
    "parametric"
  }

  # Check if we have iteration column (bootstrap data)
  is_bootstrap <- "iteration" %in% names(df)

  if (is_bootstrap) {
    # For bootstrap data, use weighted mean across iterations
    overall_stats <- df %>%
      group_by(unit) %>%
      summarise(
        n_patients = sum(n_patients) / n_distinct(iteration),
        n_patients_loset = sum(n_patients_loset) / n_distinct(iteration),
        loset_prevalence = n_patients_loset / n_patients,

        # Calculate TE metrics using weighted means
        ote_te = weighted.mean(ote_te, n_patients_loset, na.rm = TRUE),
        tte_te = if("tte_te" %in% names(.))
          weighted.mean(tte_te, n_patients_loset, na.rm = TRUE) else NA_real_,
        btte_te = if("btte_te" %in% names(.))
          weighted.mean(btte_te, n_patients_loset, na.rm = TRUE) else NA_real_,

        # Calculate sensitivity and specificity
        sensitivity = weighted.mean(sensitivity, n_patients_loset, na.rm = TRUE),
        specificity = weighted.mean(specificity, n_patients, na.rm = TRUE),

        # Get total valid observations
        n_valid_rte = sum(n_valid_rte) / n_distinct(iteration),

        .groups = 'drop'
      ) %>%
      mutate(
        otg_te = if_else(!is.na(tte_te), ote_te - tte_te, NA_real_),
        unit = "overall",
        recommended_ci_type = recommended_ci_type
      )

    return(overall_stats)
  } else {
    # Standard direct calculation of overall stats
    tibble(
      unit = "overall",
      n_patients = nrow(df),
      n_patients_loset = sum(df$loset),
      n_valid_rte = length(ote_values),
      loset_prevalence = sum(df$loset) / nrow(df),

      # Direct calculation of RTE metrics
      ote_te = mean(ote_values, na.rm = TRUE),
      tte_te = if(!is.null(tte_values)) mean(tte_values, na.rm = TRUE) else NA_real_,
      btte_te = if(!is.null(btte_values)) mean(btte_values, na.rm = TRUE) else NA_real_,

      # Standard errors
      ote_se = sd(ote_values, na.rm = TRUE) / sqrt(length(ote_values)),
      tte_se = if(!is.null(tte_values)) sd(tte_values, na.rm = TRUE) / sqrt(length(tte_values)) else NA_real_,
      btte_se = if(!is.null(btte_values)) sd(btte_values, na.rm = TRUE) / sqrt(length(btte_values)) else NA_real_,

      # t-value
      t_value = qt(0.975, df = length(ote_values) - 1),

      # Mean values
      parametric_ote_mean = mean(ote_values, na.rm = TRUE),
      parametric_tte_mean = if(!is.null(tte_values)) mean(tte_values, na.rm = TRUE) else NA_real_,
      parametric_btte_mean = if(!is.null(btte_values)) mean(btte_values, na.rm = TRUE) else NA_real_,

      # Standard deviations
      parametric_ote_sd = sd(ote_values, na.rm = TRUE),
      parametric_tte_sd = if(!is.null(tte_values)) sd(tte_values, na.rm = TRUE) else NA_real_,
      parametric_btte_sd = if(!is.null(btte_values)) sd(btte_values, na.rm = TRUE) else NA_real_,

      # Robust standard deviations using IQR
      parametric_ote_sd_q = (quantile(ote_values, 0.75, na.rm = TRUE) -
                               quantile(ote_values, 0.25, na.rm = TRUE)) / 1.349,
      parametric_tte_sd_q = if(!is.null(tte_values))
        (quantile(tte_values, 0.75, na.rm = TRUE) -
           quantile(tte_values, 0.25, na.rm = TRUE)) / 1.349 else NA_real_,
      parametric_btte_sd_q = if(!is.null(btte_values))
        (quantile(btte_values, 0.75, na.rm = TRUE) -
           quantile(btte_values, 0.25, na.rm = TRUE)) / 1.349 else NA_real_,

      # Parametric confidence intervals
      parametric_ote_var_lower = mean(ote_values, na.rm = TRUE) - qt(0.975, df = length(ote_values) - 1) *
        sd(ote_values, na.rm = TRUE) / sqrt(length(ote_values)),
      parametric_ote_var_upper = mean(ote_values, na.rm = TRUE) + qt(0.975, df = length(ote_values) - 1) *
        sd(ote_values, na.rm = TRUE) / sqrt(length(ote_values)),
      parametric_tte_var_lower = if(!is.null(tte_values))
        mean(tte_values, na.rm = TRUE) - qt(0.975, df = length(tte_values) - 1) *
        sd(tte_values, na.rm = TRUE) / sqrt(length(tte_values)) else NA_real_,
      parametric_tte_var_upper = if(!is.null(tte_values))
        mean(tte_values, na.rm = TRUE) + qt(0.975, df = length(tte_values) - 1) *
        sd(tte_values, na.rm = TRUE) / sqrt(length(tte_values)) else NA_real_,
      parametric_btte_var_lower = if(!is.null(btte_values))
        mean(btte_values, na.rm = TRUE) - qt(0.975, df = length(btte_values) - 1) *
        sd(btte_values, na.rm = TRUE) / sqrt(length(btte_values)) else NA_real_,
      parametric_btte_var_upper = if(!is.null(btte_values))
        mean(btte_values, na.rm = TRUE) + qt(0.975, df = length(btte_values) - 1) *
        sd(btte_values, na.rm = TRUE) / sqrt(length(btte_values)) else NA_real_,

      # Non-parametric confidence intervals
      nonparametric_ote_var_lower = quantile(ote_values, 0.025, na.rm = TRUE),
      nonparametric_ote_var_upper = quantile(ote_values, 0.975, na.rm = TRUE),
      nonparametric_tte_var_lower = if(!is.null(tte_values))
        quantile(tte_values, 0.025, na.rm = TRUE) else NA_real_,
      nonparametric_tte_var_upper = if(!is.null(tte_values))
        quantile(tte_values, 0.975, na.rm = TRUE) else NA_real_,
      nonparametric_btte_var_lower = if(!is.null(btte_values))
        quantile(btte_values, 0.025, na.rm = TRUE) else NA_real_,
      nonparametric_btte_var_upper = if(!is.null(btte_values))
        quantile(btte_values, 0.975, na.rm = TRUE) else NA_real_,

      # Shapiro-Wilk p-values
      ote_shapiro_p = ote_shapiro,
      tte_shapiro_p = tte_shapiro,
      btte_shapiro_p = btte_shapiro,

      # OTG
      otg_te = ote_te - tte_te,

      # Classification metrics (calculate directly)
      sensitivity = sum(df$tp, na.rm = TRUE) / (sum(df$tp, na.rm = TRUE) + sum(df$fn, na.rm = TRUE)),
      specificity = sum(df$tn, na.rm = TRUE) / (sum(df$tn, na.rm = TRUE) + sum(df$fp, na.rm = TRUE)),

      # Recommended CI type
      recommended_ci_type = recommended_ci_type
    )
  }
}

#' Calculate RTE Bootstrap Metrics
#'
#' @param df Data frame containing patient data
#' @param subgroup Optional list for subgroup analysis
#' @param var1 Optional string for first comparison variable
#' @param var2 Optional string for second comparison variable
#' @param bootstrap_params List of bootstrap parameters
#' @param n_workers Number of workers for parallel processing
#' @param seed Seed for reproducible bootstrapping
#'
#' @return List containing bootstrap metrics and distributions
#'
#' @details
#' Performs bootstrap resampling of patients and recalculates RTE metrics for each
#' bootstrap sample. This approach mirrors the WTE bootstrap method by capturing
#' uncertainty in the overall metric rather than just patient-level values.
#'
#' @keywords internal
calculate_rte_bootstrap <- function(df, subgroup = NULL, var1 = NULL, var2 = NULL,
                                    bootstrap_params, n_workers, seed = NULL) {
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Calculate sample size based on percentage
  sample_size <- floor(nrow(df) * bootstrap_params$sample_percentage)

  # Ensure sample size is at least 1
  sample_size <- max(1, sample_size)

  # Set up parallel processing
  setup_parallel(n_workers = n_workers)
  on.exit(cleanup_parallel())

  # Initialize progress reporting
  init_progressr()

  print(paste("Starting RTE bootstrap with", bootstrap_params$n_iterations,
              "iterations, sample size", sample_size, "at", Sys.time()))

  # Perform bootstrap iterations with progress bar
  bootstrap_results <- with_progress({
    p <- progressor(steps = bootstrap_params$n_iterations)

    future_map(1:bootstrap_params$n_iterations, function(i) {
      # Sample with replacement
      boot_indices <- sample(nrow(df), size = sample_size, replace = TRUE)
      boot_sample <- df[boot_indices, ]

      # Recalculate queue metrics for this bootstrap sample if they don't exist
      if (!any(c("observed_RTE") %in% names(boot_sample))) {
        boot_sample <- calculate_queue_metrics(boot_sample, verbose = FALSE)
      }

      # Create grouping and calculate RTE metrics for this bootstrap iteration
      boot_grouping <- create_grouping(boot_sample, subgroup, var1, var2)
      result <- calculate_rte_unit_metrics(boot_grouping)

      # Add overall metrics if more than one unit
      if (length(unique(result$unit)) > 1 & is.null(subgroup)) {
        overall_grouping <- create_overall_grouping(result, subgroup, var1, var2)
        overall_result <- calculate_rte_overall_metrics(result, boot_sample, overall_grouping)
        result <- bind_rows(overall_result, result)
      }

      # Add iteration number
      result$iteration <- i

      p()
      result
    }, .options = furrr_options(seed = TRUE))
  })

  print(paste("RTE bootstrap completed at", Sys.time()))

  # Combine results into a single data frame
  bootstrap_df <- bind_rows(bootstrap_results)

  # Get all grouping variables - THIS IS THE FIX
  group_vars <- c("unit")
  if (!is.null(var1)) group_vars <- c(group_vars, var1)
  if (!is.null(var2)) group_vars <- c(group_vars, var2)

  # Calculate bootstrap metrics with proper grouping
  bootstrap_metrics <- bootstrap_df %>%
    group_by(across(all_of(group_vars))) %>%  # Use all grouping variables
    summarise(
      # OTE metrics
      boot_ote_mean = mean(ote_te, na.rm = TRUE),
      boot_ote_sd = sd(ote_te, na.rm = TRUE),
      boot_ote_sd_q = (quantile(ote_te, 0.75, na.rm = TRUE) -
                         quantile(ote_te, 0.25, na.rm = TRUE)) / 1.349,
      boot_ote_var_lower = quantile(ote_te, probs = (1 - bootstrap_params$distribution_span) / 2, na.rm = TRUE),
      boot_ote_var_upper = quantile(ote_te, probs = 1 - (1 - bootstrap_params$distribution_span) / 2, na.rm = TRUE),

      # TTE metrics if available
      boot_tte_mean = if("tte_te" %in% names(.)) mean(tte_te, na.rm = TRUE) else NA_real_,
      boot_tte_sd = if("tte_te" %in% names(.)) sd(tte_te, na.rm = TRUE) else NA_real_,
      boot_tte_sd_q = if("tte_te" %in% names(.))
        (quantile(tte_te, 0.75, na.rm = TRUE) -
           quantile(tte_te, 0.25, na.rm = TRUE)) / 1.349 else NA_real_,
      boot_tte_var_lower = if("tte_te" %in% names(.))
        quantile(tte_te, probs = (1 - bootstrap_params$distribution_span) / 2, na.rm = TRUE) else NA_real_,
      boot_tte_var_upper = if("tte_te" %in% names(.))
        quantile(tte_te, probs = 1 - (1 - bootstrap_params$distribution_span) / 2, na.rm = TRUE) else NA_real_,

      # BTTE metrics if available
      boot_btte_mean = if("btte_te" %in% names(.)) mean(btte_te, na.rm = TRUE) else NA_real_,
      boot_btte_sd = if("btte_te" %in% names(.)) sd(btte_te, na.rm = TRUE) else NA_real_,
      boot_btte_sd_q = if("btte_te" %in% names(.))
        (quantile(btte_te, 0.75, na.rm = TRUE) -
           quantile(btte_te, 0.25, na.rm = TRUE)) / 1.349 else NA_real_,
      boot_btte_var_lower = if("btte_te" %in% names(.))
        quantile(btte_te, probs = (1 - bootstrap_params$distribution_span) / 2, na.rm = TRUE) else NA_real_,
      boot_btte_var_upper = if("btte_te" %in% names(.))
        quantile(btte_te, probs = 1 - (1 - bootstrap_params$distribution_span) / 2, na.rm = TRUE) else NA_real_,

      # OTG metrics if we have both OTE and TTE
      boot_otg_mean = if("tte_te" %in% names(.)) mean(ote_te - tte_te, na.rm = TRUE) else NA_real_,
      boot_otg_sd = if("tte_te" %in% names(.)) sd(ote_te - tte_te, na.rm = TRUE) else NA_real_,
      boot_otg_sd_q = if("tte_te" %in% names(.))
        (quantile(ote_te - tte_te, 0.75, na.rm = TRUE) -
           quantile(ote_te - tte_te, 0.25, na.rm = TRUE)) / 1.349 else NA_real_,
      boot_otg_var_lower = if("tte_te" %in% names(.))
        quantile(ote_te - tte_te, probs = (1 - bootstrap_params$distribution_span) / 2, na.rm = TRUE) else NA_real_,
      boot_otg_var_upper = if("tte_te" %in% names(.))
        quantile(ote_te - tte_te, probs = 1 - (1 - bootstrap_params$distribution_span) / 2, na.rm = TRUE) else NA_real_,

      # Sample information
      boot_mean_n_patients = mean(n_patients, na.rm = TRUE),
      boot_mean_n_loset = mean(n_patients_loset, na.rm = TRUE),

      .groups = 'drop'
    )

  # Return both the metrics and full distribution
  return(list(
    metrics = bootstrap_metrics,
    distributions = bootstrap_df,
    params = bootstrap_params
  ))
}

#' Calculate RTE with Segment Bootstrap
#'
#' @param df Data frame containing patient data
#' @param subgroup Optional list for subgroup analysis
#' @param var1 Optional string for first comparison variable
#' @param var2 Optional string for second comparison variable
#' @param bootstrap_params List of bootstrap parameters
#' @param n_workers Number of workers for parallel processing
#' @param seed Seed for reproducible bootstrapping
#'
#' @return List containing bootstrap metrics and distributions
#'
#' @details
#' Performs segment-based bootstrap for RTE calculations, where entire queue
#' segments are resampled rather than individual patients. This preserves the
#' queue dynamics while estimating the uncertainty in the RTE metrics.
#'
#' @keywords internal
calculate_rte_segment_bootstrap <- function(df, subgroup = NULL, var1 = NULL, var2 = NULL,
                                            bootstrap_params, n_workers, seed = NULL) {
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Check if segment column exists; if not, create segments
  if (!"segment" %in% names(df)) {
    print(paste("Creating segments for segment bootstrap", Sys.time()))
    df <- create_segments(df, n_workers)
  } else {
    print(paste("Using existing segments for bootstrap", Sys.time()))
  }

  # Get unit-segment combinations to preserve unit separation during sampling
  # This ensures we maintain queue independence between units
  unit_segment_pairs <- df %>%
    select(unit, segment) %>%
    distinct()

  n_segments <- nrow(unit_segment_pairs)
  n_samples <- floor(n_segments * bootstrap_params$sample_percentage)

  print(paste("Starting RTE segment bootstrap with", n_segments, "unit-segment pairs,",
              n_samples, "samples per iteration,",
              bootstrap_params$n_iterations, "iterations", Sys.time()))

  # Set up parallel processing
  setup_parallel(n_workers = n_workers)
  on.exit(cleanup_parallel())

  # Initialize progress reporting
  init_progressr()

  # Perform bootstrap iterations with progress bar
  bootstrap_results <- with_progress({
    p <- progressor(steps = bootstrap_params$n_iterations)
    future_map(1:bootstrap_params$n_iterations, function(i) {
      # Sample unit-segment pairs with replacement
      sample_indices <- sample(n_segments, size = n_samples, replace = TRUE)
      sampled_pairs <- unit_segment_pairs[sample_indices, ]

      # Create bootstrap sample by filtering rows matching sampled unit-segment combinations
      boot_sample <- df %>%
        dplyr::inner_join(sampled_pairs, by = c("unit", "segment"), relationship = "many-to-many")

      # Recalculate queue metrics for this bootstrap sample if they don't exist
      if (!any(c("observed_RTE") %in% names(boot_sample))) {
        boot_sample <- calculate_queue_metrics(boot_sample, verbose = FALSE)
      }

      # Create grouping and calculate RTE metrics for this bootstrap iteration
      boot_grouping <- create_grouping(boot_sample, subgroup, var1, var2)
      result <- calculate_rte_unit_metrics(boot_grouping)

      # Add overall metrics if more than one unit
      if (length(unique(result$unit)) > 1 & is.null(subgroup)) {
        overall_grouping <- create_overall_grouping(result, subgroup, var1, var2)
        overall_result <- calculate_rte_overall_metrics(result, boot_sample, overall_grouping)
        result <- bind_rows(overall_result, result)
      }

      # Add iteration number
      result$iteration <- i

      p()
      result
    }, .options = furrr_options(seed = TRUE))
  })

  print(paste("RTE segment bootstrap completed at", Sys.time()))

  # Combine results into a single data frame
  bootstrap_df <- bind_rows(bootstrap_results)

  # Get all grouping variables - THIS IS THE FIX
  group_vars <- c("unit")
  if (!is.null(var1)) group_vars <- c(group_vars, var1)
  if (!is.null(var2)) group_vars <- c(group_vars, var2)

  # Calculate bootstrap metrics using weighted means for overall calculations with proper grouping
  bootstrap_metrics <- bootstrap_df %>%
    group_by(across(all_of(group_vars))) %>%  # Use all grouping variables
    summarise(
      # OTE metrics
      boot_ote_mean = mean(ote_te, na.rm = TRUE),
      boot_ote_sd = sd(ote_te, na.rm = TRUE),
      boot_ote_sd_q = (quantile(ote_te, 0.75, na.rm = TRUE) -
                         quantile(ote_te, 0.25, na.rm = TRUE)) / 1.349,
      boot_ote_var_lower = quantile(ote_te, probs = (1 - bootstrap_params$distribution_span) / 2, na.rm = TRUE),
      boot_ote_var_upper = quantile(ote_te, probs = 1 - (1 - bootstrap_params$distribution_span) / 2, na.rm = TRUE),

      # TTE metrics if available
      boot_tte_mean = if("tte_te" %in% names(.)) mean(tte_te, na.rm = TRUE) else NA_real_,
      boot_tte_sd = if("tte_te" %in% names(.)) sd(tte_te, na.rm = TRUE) else NA_real_,
      boot_tte_sd_q = if("tte_te" %in% names(.))
        (quantile(tte_te, 0.75, na.rm = TRUE) -
           quantile(tte_te, 0.25, na.rm = TRUE)) / 1.349 else NA_real_,
      boot_tte_var_lower = if("tte_te" %in% names(.))
        quantile(tte_te, probs = (1 - bootstrap_params$distribution_span) / 2, na.rm = TRUE) else NA_real_,
      boot_tte_var_upper = if("tte_te" %in% names(.))
        quantile(tte_te, probs = 1 - (1 - bootstrap_params$distribution_span) / 2, na.rm = TRUE) else NA_real_,

      # BTTE metrics if available
      boot_btte_mean = if("btte_te" %in% names(.)) mean(btte_te, na.rm = TRUE) else NA_real_,
      boot_btte_sd = if("btte_te" %in% names(.)) sd(btte_te, na.rm = TRUE) else NA_real_,
      boot_btte_sd_q = if("btte_te" %in% names(.))
        (quantile(btte_te, 0.75, na.rm = TRUE) -
           quantile(btte_te, 0.25, na.rm = TRUE)) / 1.349 else NA_real_,
      boot_btte_var_lower = if("btte_te" %in% names(.))
        quantile(btte_te, probs = (1 - bootstrap_params$distribution_span) / 2, na.rm = TRUE) else NA_real_,
      boot_btte_var_upper = if("btte_te" %in% names(.))
        quantile(btte_te, probs = 1 - (1 - bootstrap_params$distribution_span) / 2, na.rm = TRUE) else NA_real_,

      # OTG metrics if we have both OTE and TTE
      boot_otg_mean = if("tte_te" %in% names(.)) mean(ote_te - tte_te, na.rm = TRUE) else NA_real_,
      boot_otg_sd = if("tte_te" %in% names(.)) sd(ote_te - tte_te, na.rm = TRUE) else NA_real_,
      boot_otg_sd_q = if("tte_te" %in% names(.))
        (quantile(ote_te - tte_te, 0.75, na.rm = TRUE) -
           quantile(ote_te - tte_te, 0.25, na.rm = TRUE)) / 1.349 else NA_real_,
      boot_otg_var_lower = if("tte_te" %in% names(.))
        quantile(ote_te - tte_te, probs = (1 - bootstrap_params$distribution_span) / 2, na.rm = TRUE) else NA_real_,
      boot_otg_var_upper = if("tte_te" %in% names(.))
        quantile(ote_te - tte_te, probs = 1 - (1 - bootstrap_params$distribution_span) / 2, na.rm = TRUE) else NA_real_,

      # Sample information
      boot_mean_n_patients = mean(n_patients, na.rm = TRUE),
      boot_mean_n_loset = mean(n_patients_loset, na.rm = TRUE),

      .groups = 'drop'
    )

  # Return both the metrics and full distribution
  return(list(
    metrics = bootstrap_metrics,
    distributions = bootstrap_df,
    params = bootstrap_params,
    method = "segment"
  ))
}


#' Validate RTE Results and Generate Warnings
#'
#' @param results Combined unit and overall results
#' @param min_loset_warning Least number of LOSET cases that should exist or the function will warn.
#'
#' @return Logical indicating if non-normality was detected in any unit
#'
#' @details
#' Performs validation checks on calculated results:
#' - Low LOSET counts
#' - Units with too few valid RTE measurements
#' - Potential non-normal distributions
#'
#' @keywords internal
validate_rte_results <- function(results, min_loset_warning, quite) {
  # Check for low loset counts
  low_loset <- results %>%
    filter(unit != "overall", n_patients_loset < min_loset_warning)

  if (nrow(low_loset) > 0) {
    warning("Units with fewer than specified minimum LOSET cases found: ",
            paste(low_loset$unit, collapse = ", "))
  }

  # Check for units with insufficient valid RTE measurements
  low_valid_rte <- results %>%
    filter(unit != "overall", n_valid_rte < 3)

  if (nrow(low_valid_rte) > 0) {
    warning("Units with fewer than 3 valid RTE measurements found: ",
            paste(low_valid_rte$unit, collapse = ", "),
            ". Statistical inference may be unreliable.")
  }

  # Check normality for both OTE and TTE
  normality_issues_ote <- results %>%
    filter(unit != "overall", !is.na(ote_shapiro_p), ote_shapiro_p < 0.05)

  normality_issues_tte <- results %>%
    filter(unit != "overall", !is.na(tte_shapiro_p), tte_shapiro_p < 0.05)

  # Combine units with normality issues
  all_normality_issues <- bind_rows(
    normality_issues_ote %>% select(unit),
    normality_issues_tte %>% select(unit)
  ) %>% distinct()

  if (nrow(all_normality_issues) > 0 & !quite) {
    warning("Non-normal RTE distributions detected in: ",
            paste(all_normality_issues$unit, collapse = ", "),
            ". Using non-parametric confidence intervals.")

    # Return TRUE if non-normality was detected
    return(TRUE)
  }

  # Return FALSE if all distributions appear normal
  return(FALSE)
}


#' Combine RTE Base Results with Bootstrap Results
#'
#' @param base_results Tibble with base RTE calculations
#' @param bootstrap_metrics Tibble with bootstrap metrics
#'
#' @return Combined results tibble
#'
#' @details
#' Merges direct calculation results with bootstrap metrics.
#'
#' @keywords internal
combine_rte_results <- function(base_results, bootstrap_metrics) {
  # Get the key columns for joining
  key_cols <- intersect(
    names(base_results),
    names(bootstrap_metrics)
  )

  # If key_cols only contains "unit", use that
  # Otherwise filter to include only unit and any grouping vars
  join_cols <- if(length(key_cols) == 1) "unit" else key_cols

  # Join base results with bootstrap metrics
  combined_metrics <- base_results %>%
    left_join(
      bootstrap_metrics %>%
        select(starts_with("boot_"), all_of(join_cols)),
      by = join_cols
    )

  return(combined_metrics)
}

#' Plot RTE Distribution
#'
#' @description
#' Creates a histogram visualization of RTE values with statistical annotations.
#' The function filters to valid time-critical patients, removes extreme outliers,
#' and can display normality test results and basic statistics.
#'
#' @param data Data frame containing RTE values
#' @param metric Character. Which RTE metric to plot: "observed_RTE", "theoretical_RTE",
#'   or "binary_theoretical_RTE". Default is "observed_RTE".
#' @param cutoff_percentile Numeric. Percentile cutoff for removing extreme outliers.
#'   Default is 0.01 (removes the lowest 1%).
#' @param bin_width Numeric. Width of histogram bins. Default is 0.1.
#' @param include_stats Logical. Whether to display statistical annotations. Default is TRUE.
#' @param color Character. Color for histogram bars. Default is "#D4AF37" (gold).
#' @param facet_by_unit Logical. Default is FALSE. If TRUE, creates separate facet panels for each unit, allowing comparison of RTE distributions across different units.
#'
#' @return A ggplot object
#'
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline labs theme_minimal annotate
#' @importFrom stats sd median shapiro.test quantile
#' @importFrom rlang sym !!
#'
#' @export
plot_rte_distribution <- function(data, metric = "observed_RTE",
                                  cutoff_percentile = 0.01, bin_width = 0.1,
                                  include_stats = TRUE, color = "#D4AF37",
                                  facet_by_unit = FALSE) {
  # Handle both data frame and calc_te inputs
  if (inherits(data, "calc_te")) {
    if (is.null(data$distributions)) {
      stop("No distribution data available. Re-run calc_rte with include_distributions = TRUE")
    }
    filtered_data <- data$distributions
  } else {
    filtered_data <- data
  }

  # Filter valid RTE values
  filtered_data <- filtered_data %>%
    filter(valid == TRUE, !is.na(!!sym(metric)))

  # Apply outlier trimming
  if (nrow(filtered_data) > 0) {
    min_threshold <- quantile(filtered_data[[metric]], cutoff_percentile, na.rm = TRUE)
    max_threshold <- quantile(filtered_data[[metric]], 1 - cutoff_percentile, na.rm = TRUE)
    filtered_data <- filtered_data %>%
      filter(!!sym(metric) >= min_threshold & !!sym(metric) <= max_threshold)
  }

  # Create base plot
  p <- ggplot(filtered_data, aes(x = !!sym(metric))) +
    geom_histogram(binwidth = bin_width, fill = color, alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
    geom_vline(xintercept = 1, linetype = "dotted", color = "darkgray") +
    labs(title = paste("Distribution of", gsub("_", " ", metric)),
         x = "Value", y = "Count") +
    theme_minimal()

  # Add faceting by unit if requested
  if (facet_by_unit && "unit" %in% names(filtered_data)) {
    p <- p + facet_wrap(~unit, scales = "free_y")
  }

  # Add statistics
  if (include_stats && nrow(filtered_data) > 2) {
    # Calculate statistics by unit or overall
    if (facet_by_unit && "unit" %in% names(filtered_data)) {
      stats_data <- filtered_data %>%
        group_by(unit) %>%
        summarize(
          n = n(),
          mean_val = mean(!!sym(metric), na.rm = TRUE),
          median_val = median(!!sym(metric), na.rm = TRUE),
          sd_val = sd(!!sym(metric), na.rm = TRUE),
          shapiro_p = if(n() > 2 && n() <= 5000) {
            tryCatch(shapiro.test(!!sym(metric))$p.value,
                     error = function(e) NA_real_)
          } else NA_real_
        )

      # Add text annotations to each facet
      for (u in unique(filtered_data$unit)) {
        unit_stats <- stats_data %>% filter(unit == u)
        stats_text <- sprintf("n = %d\nMean: %.2f\nMedian: %.2f\nSD: %.2f",
                              unit_stats$n, unit_stats$mean_val,
                              unit_stats$median_val, unit_stats$sd_val)

        if (!is.na(unit_stats$shapiro_p)) {
          stats_text <- paste0(stats_text,
                               sprintf("\nShapiro-Wilk p: %.4g %s",
                                       unit_stats$shapiro_p,
                                       ifelse(unit_stats$shapiro_p < 0.05,
                                              "(non-normal)", "(normal)")))
        }

        p <- p + annotate("text", x = Inf, y = Inf, label = stats_text,
                          hjust = 1.1, vjust = 1.1, size = 3.5,
                          data = data.frame(unit = u))
      }
    } else {
      # Overall statistics
      mean_val <- mean(filtered_data[[metric]], na.rm = TRUE)
      median_val <- median(filtered_data[[metric]], na.rm = TRUE)
      sd_val <- sd(filtered_data[[metric]], na.rm = TRUE)
      n_val <- nrow(filtered_data)

      shapiro_result <- if(n_val > 2 && n_val <= 5000) {
        tryCatch(shapiro.test(filtered_data[[metric]])$p.value,
                 error = function(e) NA_real_)
      } else {
        NA_real_
      }

      stats_text <- sprintf("n = %d\nMean: %.2f\nMedian: %.2f\nSD: %.2f",
                            n_val, mean_val, median_val, sd_val)

      if (!is.na(shapiro_result)) {
        stats_text <- paste0(stats_text,
                             sprintf("\nShapiro-Wilk p: %.4g %s",
                                     shapiro_result,
                                     ifelse(shapiro_result < 0.05,
                                            "(non-normal)", "(normal)")))
      }

      p <- p + annotate("text", x = Inf, y = Inf, label = stats_text,
                        hjust = 1.1, vjust = 1.1, size = 3.5)
    }
  }

  return(p)
}



#' Calculate Queue Metrics for Rank-based Triage Effectiveness (RTE) Analysis
#'
#' @description
#' This function implements the Rank-based Triage Effectiveness (RTE) methodology to evaluate
#' triage system performance based on patient queue positioning rather than waiting times.
#' RTE measures how effectively a triage system prioritizes time-critical patients relative
#' to their positions in a theoretical first-come-first-served system, providing a patient-level
#' metric that complements the waiting time-based Triage Effectiveness (WTE) approach.
#'
#' @param data A data frame that has been initiated with init() with the following required columns:
#'   \itemize{
#'     \item segment: Created in init(), a grouping variable identifying distinct queue segments.
#'     \item arrival_minute: Time of patient arrival (numeric).
#'     \item observed_resolve_time: Time when patient was actually resolved (numeric), created by init().
#'     \item theoretica_resolve_time: Theoretical resolution time (numeric) from sim_te(). Optional.
#'     \item binary_theoretical_resolve_time: Theoretical resolution time when using binary priorities (numeric). Optional.
#'     \item loset: Logical flag indicating if patient is time-critical (TRUE) or not (FALSE).
#'   }
#' @param verbose logical, determines if intermediate calculation fields are returned (TRUE) or only the final RTE results (FALSE).
#' @param ambiguity_invalid logical, default is TRUE, sets if cases where n_tc => L are invalid since they are ambigous and thus cannot be used to draw conclusions. I.e. p = 1 theoretically gives both TE = 100% and 0%.
#' @param ambiguity_perfect logical, default is TRUE, used if ambiguity_invalit == FALSE, if TRUE, cases where n_tc => L gives 100% TE, if false, gives 0% TE.
#'
#' @return A data frame with the original data plus additional columns:
#'   \itemize{
#'     \item arrivals_up_to_now: (if verbose == TRUE) Cumulative count of arrivals at each patient's arrival time.
#'     \item observed_resolves_before: (if verbose == TRUE) Count of patients resolved before each patient's arrival.
#'     \item L: (if verbose == TRUE) Queue length at patient arrival (including the patient).
#'     \item observed_p: (if verbose == TRUE) Patient's position in the resolution sequence, counting from arrival to observed resolution.
#'     \item theoretical_p: (if verbose == TRUE) Patient's position in the theoretical resolution sequence.
#'     \item n_tc: (if verbose == TRUE) Number of time-critical patients in the concurrent window (from last resolve before arrival to first resolve after).
#'     \item is_cluster: Logical flag indicating if the patient is part of a time-critical cluster.
#'     \item cluster_id: Identifier for patients belonging to the same time-critical cluster.
#'     \item valid: Logical flag indicating if the time-critical patient's RTE calculation is valid (see details).
#'     \item observed_RTE: Rank-based Triage Effectiveness based on observed resolution times (NA for non-time-critical patients).
#'     \item theoretical_RTE: Rank-based Triage Effectiveness based on theoretical resolution times (NA for non-time-critical patients).
#'   }
#'
#' @note
#' This function makes specific assumptions when handling multiple patients resolved in
#' the same minute:
#'
#' For observed_p calculation: When multiple patients are resolved in the same minute,
#' the function assumes that higher priority patients were resolved first, followed by
#' patients with equal priority in arrival order. This may not reflect the actual
#' sequence that occurred in practice, as this level of detail is often not captured
#' in typical ED datasets with minute-level resolution. The calculation therefore
#' interprets the observed data in the "best light" - assuming ideal priority-based
#' ordering within each minute.
#'
#' For theoretical_p and binary_theoretical_p: The priority-based ordering is
#' consistent with the simulation used in sim_te() and sim_heat(), ensuring that
#' within the same minute, higher priority patients are always processed before
#' lower priority patients, and equal priority patients are processed in arrival order.
#'
#' Cluster handling: The function identifies time-critical clusters when there are
#' multiple time-critical cases (n_tc > 1) and at least one has n_tc > L. Patients
#' in these clusters are evaluated using a modified RTE formula that better accounts
#' for the unique dynamics of clustered time-critical arrivals.
#'
#' This handling is necessary to achieve consistent 100% TE scores in scenarios
#' with perfect triage (100% sensitivity and specificity), particularly when
#' multiple patients are resolved within the same minute or when time-critical
#' patients arrive in clusters.
#'
#' @details
#' ## Conceptual Framework
#'
#' RTE quantifies triage performance through queue positioning rather than waiting times. It evaluates
#' how effectively a triage system moves time-critical patients forward in the resolution queue compared
#' to their position in a first-come-first-served (FCFS) system, where:
#'
#' * RTE = 0 means the patient's position was equivalent to FCFS (no triage benefit)
#' * RTE = 1 means optimal positioning (perfect triage)
#' * RTE < 0 means worse positioning than FCFS (triage was counterproductive)
#'
#' ## Key Metrics Calculation
#'
#' \strong{L (Queue Length)}: Represents the number of patients in queue when a patient arrives,
#' including the patient themselves. This value indicates what position the patient would have in a
#' FCFS system.
#'
#' \strong{p (Position)}: The patient's position in the resolution sequence, counting from their arrival
#' time to their resolution time. A p value of 1 means the patient was the first to be resolved among
#' all patients in queue at their arrival time.
#'
#' \strong{n_tc}: Number of time-critical patients (loset=TRUE) in the concurrent window between the last resolve
#' before the patient's arrival and the first resolve after their arrival, including the patient themselves.
#' This represents the "cluster" of time-critical cases that arrived in close temporal proximity.
#'
#' \strong{valid}: A time-critical patient's RTE calculation is considered valid only when triage can
#' make a meaningful difference. Specifically, the calculation is valid if any of these conditions are met:
#'   \itemize{
#'     \item There's at least one non-time-critical patient already in the queue at arrival time, OR
#'     \item A non-time-critical patient will arrive before L resolves happen (if L = n_tc), OR
#'     \item L ≠ n_tc (queue composition allows for meaningful prioritization)
#'   }
#'   This validation prevents analysis of scenarios where the order doesn't matter (e.g., when all patients
#'   in queue are time-critical and no new non-time-critical patients arrive before they're all resolved).
#'
#' @details
#' ## RTE Calculation Logic
#'
#' The RTE formula varies depending on the queue dynamics to ensure consistent interpretation across different scenarios:
#'
#' 1. \strong{Perfect Triage} (when p ≤ n_tc):
#'    * The patient is resolved within the first n_tc positions
#'    * If n_tc ≥ L and ambiguity_perfect=TRUE: RTE = 1 (perfect score)
#'    * If n_tc ≥ L and ambiguity_perfect=FALSE: RTE = 0
#'    * Otherwise (when n_tc < L): RTE = 1 (perfect score)
#'
#' 2. \strong{Positive Triage} (when n_tc < p ≤ L):
#'    * The patient is resolved after the first n_tc positions but within or at L
#'    * RTE = (L - p) / (L - n_tc)
#'
#' 3. \strong{Negative Triage} (when p > L):
#'    * The patient is resolved after L positions (worse than first-come, first-served)
#'    * RTE = (L - p) / L (this will be negative)
#'
#' ## Ambiguity Handling
#'
#' An ambiguous case occurs when n_tc ≥ L, meaning there are as many or more time-critical patients
#' than the queue length. In this situation, the position would be the same with perfect prioritization
#' as in a first-come, first-serve system.
#'
#' By default (ambiguity_invalid=TRUE), these cases are considered invalid for RTE calculation.
#' When ambiguity_invalid=FALSE, the result depends on ambiguity_perfect:
#' * When ambiguity_perfect=TRUE (default): These cases are scored as "perfect" (RTE = 1)
#' * When ambiguity_perfect=FALSE: These cases are scored as "zero" (RTE = 0)
#' ## Examples
#'
#' \strong{Example 1:} A time-critical patient arrives and finds 5 patients already in queue (L = 6).
#' They are the only time-critical patient in this window (n_tc = 1). They get resolved as the 2nd patient (p = 2).
#'
#' Using case 3 formula: RTE = (6 - 2)/(6 - (1+1)/2) = 4/5 = 0.8
#'
#' This indicates very good triage performance, moving the patient nearly to optimal position.
#'
#' \strong{Example 2:} A time-critical patient arrives at an empty queue (L = 1) and is the only time-critical
#' patient (n_tc = 1). Before this patient is resolved, 3 non-time-critical patients arrive. The time-critical
#' patient gets resolved as the 3rd patient (p = 3).
#'
#' Using case 4 formula: RTE = (1 + ((1+1)/2-1) - 3)/1 = -1.5
#'
#' This indicates poor triage, as the time-critical patient was significantly delayed despite arriving first.
#'
#' \strong{Example 3:} During a major incident, 4 time-critical patients arrive in rapid succession to an ED
#' with 2 patients already waiting (L = 3 for the first arriving time-critical patient, n_tc = 4). The first
#' time-critical patient gets resolved as the 2nd patient (p = 2).
#'
#' Using case 2 formula: RTE = ((4+1)/2 - 2)/((4+1)/2 - 1) = 0.5/1.5 = 0.33
#'
#' This indicates moderate triage performance within the context of a time-critical patient cluster.
#'
#' ## Differences between Rank-based TE (RTE) and Waiting time-based TE (WTE)
#'
#' 1. \strong{Conceptual Focus:}
#'    * RTE measures position improvement in the queue (process metric)
#'    * WTE measures waiting time reduction (outcome metric)
#'
#' 2. \strong{Analytical Granularity:}
#'    * RTE provides patient-level measurements, enabling standard statistical analyses without bootstrapping
#'    * WTE produces aggregate measurements requiring bootstrapping for confidence intervals
#'
#' 3. \strong{Queue Dynamics Handling:}
#'    * RTE handles non-coincident arrivals and resolves better than WTE
#'    * WTE better reflects the actual time impact experienced by patients
#'
#' 4. \strong{Edge Case Behavior:}
#'    * RTE can punish placement severely in short queues where the time effect is minimal
#'    * WTE can punish EDs with lower resolve rates and fewer priority 1 cases
#'
#' 5. \strong{System Evaluation:}
#'    * RTE focuses specifically on triage decision quality regardless of system throughput
#'    * WTE combines effects of both triage decisions and system throughput capacity
#'
#' 6. \strong{Interpretation Differences:}
#'    * RTE measures position advantage relative to FCFS position and cluster composition
#'    * WTE measures time advantage relative to mean waiting time for all patients
#'
#' 7. \strong{Cross-System Comparison:}
#'    * RTE may be more stable across different ED environments with varying throughput capacities
#'    * WTE more directly reflects patient experience in terms of actual waiting time
#'
#' @note
#' This function makes specific assumptions when handling multiple patients resolved in
#' the same minute:
#'
#' For observed_p calculation: When multiple patients are resolved in the same minute,
#' the function assumes that higher priority patients were resolved first, followed by
#' patients with equal priority in arrival order. This may not reflect the actual
#' sequence that occurred in practice, as this level of detail is often not captured
#' in typical ED datasets with minute-level resolution. The calculation therefore
#' interprets the observed data in the "best light" - assuming ideal priority-based
#' ordering within each minute.
#'
#' For theoretical_p and binary_theoretical_p: The priority-based ordering is
#' consistent with the simulation used in sim_te() and sim_heat(), ensuring that
#' within the same minute, higher priority patients are always processed before
#' lower priority patients, and equal priority patients are processed in arrival order.
#'
#' This handling is necessary to achieve consistent 100% TE scores in scenarios
#' with perfect triage (100% sensitivity and specificity), particularly when
#' multiple patients are resolved within the same minute.
#'
#' @seealso
#' \code{\link{init}} for initializing data for triage analysis
#' \code{\link{sim_te}} for simulating theoretical resolve times
#' Calculate Queue Metrics for Rank-based Triage Effectiveness (RTE) Analysis using data.table
#'
#' @importFrom data.table data.table as.data.table setDT setkey setorder :=
#' @importFrom tibble as_tibble
#' @export
calculate_queue_metrics <- function(data, verbose = TRUE, ambiguity_invalid = TRUE, ambiguity_perfect = TRUE) {
  RTE_method <- function(L, p, n_tc) {
    # Perfect triage: patient is seen within n_tc positions
    if (p <= n_tc) {
      if (n_tc >= L & ambiguity_perfect) {
        return("perfect")
      } else if(n_tc >= L & !ambiguity_perfect) {
        return("zero")
      } else {
        return("perfect")
      }
    }
    # Positive triage: patient is seen after n_tc but within or at L
    else if (p <= L) {
      return("positive")
    }
    # Negative triage: patient is seen after L
    else {
      return("negative")
    }
  }

  RTE_calculate <- function(L, p, n_tc, method) {
    switch(method,
           "perfect" = 1,
           "positive" = return((L - p) / (L - n_tc)),
           "negative" = (L - p) / L,
           "zero" = 0,
           NA_real_  # Default case
    )
  }

  # Check data columns and prepare the data
  has_observed <- "priority" %in% names(data)
  has_theoretical <- "theoretical_wait_time" %in% names(data)
  has_binary_theoretical <- "binary_theoretical_wait_time" %in% names(data)

  # Convert to data.table - more efficient for large operations
  dt <- as.data.table(data)

  # PRE-CALCULATE THEORETICAL RESOLVE TIMES if needed
  if (has_theoretical) {
    dt[, theoretical_resolve_time := arrival_minute + theoretical_wait_time]
  }

  if (has_binary_theoretical) {
    dt[, binary_theoretical_resolve_time := arrival_minute + binary_theoretical_wait_time]
  }

  # Sort by segment and arrival_minute for more efficient processing
  setkey(dt, segment, arrival_minute)

  # Calculate basic metrics in one pass
  dt[, `:=`(
    # Count arrivals within segment
    arrivals_up_to_now = seq_len(.N),
    # Store all resolve times for each segment
    segment_resolves = list(resolve_minute)
  ), by = segment]

  # Calculate resolves before each arrival
  dt[, observed_resolves_before := {
    resolves <- segment_resolves[[1]]
    vapply(arrival_minute, function(t) sum(resolves < t), integer(1))
  }, by = segment]

  # Calculate queue length
  dt[, L := arrivals_up_to_now - observed_resolves_before]

  # Remove temporary column
  dt[, segment_resolves := NULL]

  # Calculate window boundaries
  dt[, `:=`(
    # Find last resolve before each arrival
    last_resolve_before = {
      vapply(arrival_minute, function(t) {
        resolves_before <- resolve_minute[resolve_minute < t]
        ifelse(length(resolves_before) > 0, max(resolves_before), -Inf)
      }, numeric(1))
    },
    # Find first resolve after or at arrival
    first_resolve_after = {
      vapply(arrival_minute, function(t) {
        resolves_after <- resolve_minute[resolve_minute >= t]
        ifelse(length(resolves_after) > 0, min(resolves_after), Inf)
      }, numeric(1))
    }
  ), by = segment]

  # Create window identifier
  dt[, window_id := paste0(last_resolve_before, "_", first_resolve_after)]

  # Create window statistics table
  window_stats <- dt[, .(
    # Count loset patients in window
    loset_arrivals = sum(loset),
    # For ordering
    min_arrival = min(arrival_minute),
    window_end = first(first_resolve_after)
  ), by = .(segment, window_id)]

  # Sort window stats by segment and arrival time
  setorder(window_stats, segment, min_arrival)

  # Get resolve times for each segment
  segment_resolves <- dt[, .(resolve_times = list(resolve_minute)), by = segment]

  # Join window stats with segment resolves
  window_stats <- merge(window_stats, segment_resolves, by = "segment", sort = FALSE)

  # Calculate number of resolves at window end
  window_stats[, n_resolves := sapply(seq_len(.N), function(i) {
    sum(resolve_times[[i]] == window_end[i])
  })]

  # Clean up
  window_stats[, resolve_times := NULL]

  # Calculate carryover and adjusted n_tc by segment
  window_stats[, `:=`(
    carryover = {
      n_rows <- .N
      carryover <- integer(n_rows)
      unresolved_tc <- 0

      for (i in 1:n_rows) {
        carryover[i] <- unresolved_tc
        unresolved_tc <- unresolved_tc + loset_arrivals[i] - n_resolves[i]
        unresolved_tc <- max(0, unresolved_tc)
      }
      carryover
    },
    adjusted_n_tc = {
      n_rows <- .N
      adjusted <- integer(n_rows)
      unresolved_tc <- 0

      for (i in 1:n_rows) {
        adjusted[i] <- loset_arrivals[i] + unresolved_tc
        unresolved_tc <- unresolved_tc + loset_arrivals[i] - n_resolves[i]
        unresolved_tc <- max(0, unresolved_tc)
      }
      adjusted
    }
  ), by = segment]

  # Join adjusted n_tc back to main data
  dt <- merge(
    dt,
    window_stats[, .(segment, window_id, adjusted_n_tc, carryover)],
    by = c("segment", "window_id"),
    all.x = TRUE
  )

  # Update n_tc with adjusted value and clean up
  dt[, `:=`(
    n_tc = adjusted_n_tc,
    adjusted_n_tc = NULL,
    carryover = NULL,
    window_id = NULL
  )]


  # Calculate observed_p
  if (has_observed) {
    # Process each segment separately
    dt[, observed_p := {
      result <- integer(.N)

      # For each patient in this segment
      for (i in 1:.N) {
        # Skip single-patient queues
        if (L[i] <= 1) {
          result[i] <- 1
          next
        }

        arrival_t <- arrival_minute[i]
        resolve_t <- resolve_minute[i]

        # Find candidates resolved during this wait
        candidates <- which(
          resolve_minute >= arrival_t &
            resolve_minute <= resolve_t
        )

        # If no eligible patients or only me, position is 1
        if (length(candidates) <= 1) {
          result[i] <- 1
          next
        }

        # Get candidate values for sorting
        candidate_times <- resolve_minute[candidates]
        candidate_priorities <- priority[candidates]
        candidate_arrivals <- arrival_minute[candidates]

        # Create ordering based on resolve time, priority, arrival
        order_index <- order(
          candidate_times,
          candidate_priorities,
          candidate_arrivals
        )

        # Find position in ordered list
        my_pos <- which(candidates[order_index] == i)

        # Set result
        result[i] <- if (length(my_pos) == 0) 1 else my_pos
      }

      result
    }, by = segment]
  }


  # Calculate theoretical_p
  if (has_theoretical) {
    # Process each segment separately
    dt[, theoretical_p := {
      result <- integer(.N)

      # For each patient in this segment
      for (i in 1:.N) {
        # Skip single-patient queues
        if (L[i] <= 1) {
          result[i] <- 1
          next
        }

        arrival_t <- arrival_minute[i]
        resolve_t <- theoretical_resolve_time[i]

        # Find candidates resolved during this wait
        candidates <- which(
          theoretical_resolve_time >= arrival_t &
            theoretical_resolve_time <= resolve_t
        )

        # If no eligible patients or only me, position is 1
        if (length(candidates) <= 1) {
          result[i] <- 1
          next
        }

        # Get candidate values for sorting
        candidate_times <- theoretical_resolve_time[candidates]
        candidate_priorities <- priority[candidates]
        candidate_arrivals <- arrival_minute[candidates]

        # Create ordering based on resolve time, priority, arrival
        order_index <- order(
          candidate_times,
          candidate_priorities,
          candidate_arrivals
        )

        # Find position in ordered list
        my_pos <- which(candidates[order_index] == i)

        # Set result
        result[i] <- if (length(my_pos) == 0) 1 else my_pos
      }

      result
    }, by = segment]
  }

  # Calculate binary_theoretical_p
  if (has_binary_theoretical) {
    # Process each segment separately
    dt[, binary_theoretical_p := {
      result <- integer(.N)

      # For each patient in this segment
      for (i in 1:.N) {
        # Skip single-patient queues
        if (L[i] <= 1) {
          result[i] <- 1
          next
        }

        arrival_t <- arrival_minute[i]
        resolve_t <- binary_theoretical_resolve_time[i]

        # Find candidates resolved during this wait
        candidates <- which(
          binary_theoretical_resolve_time >= arrival_t &
            binary_theoretical_resolve_time <= resolve_t
        )

        # If no eligible patients or only me, position is 1
        if (length(candidates) <= 1) {
          result[i] <- 1
          next
        }

        # Get candidate values for sorting
        candidate_times <- binary_theoretical_resolve_time[candidates]
        candidate_priorities <- priority_binary[candidates]
        candidate_arrivals <- arrival_minute[candidates]

        # Create ordering based on resolve time, priority, arrival
        order_index <- order(
          candidate_times,
          candidate_priorities,
          candidate_arrivals
        )

        # Find position in ordered list
        my_pos <- which(candidates[order_index] == i)

        # Set result
        result[i] <- if (length(my_pos) == 0) 1 else my_pos
      }

      result
    }, by = segment]
  }


  # Calculate valid flag
  dt[, `:=`(
    segment_has_non_loset = {
      any(!loset)
    },
    # Only time-critical patients need the valid flag to be meaningful
    valid = TRUE
  ), by = segment]

  # Update valid flag with segment_has_non_loset info
  dt[, valid := ifelse(!loset, TRUE, segment_has_non_loset)]

  # Apply ambiguity invalidation
  dt[n_tc >= L, valid := !ambiguity_invalid]

  # Clean up
  dt[, segment_has_non_loset := NULL]


  if (has_observed) {
    # Calculate RTE metrics
    dt[loset == TRUE & valid == TRUE, `:=`(
      # Determine RTE calculation method
      observed_RTE_case = mapply(RTE_method, L, observed_p, n_tc),
      # Calculate RTE
      observed_RTE = mapply(RTE_calculate, L, observed_p, n_tc,
                            mapply(RTE_method, L, observed_p, n_tc))
    )]
  }

  # Add theoretical RTE if available
  if (has_theoretical) {
    dt[loset == TRUE & valid == TRUE, `:=`(
      # Determine theoretical RTE calculation method
      theoretical_RTE_case = mapply(RTE_method, L, theoretical_p, n_tc),
      # Calculate theoretical RTE
      theoretical_RTE = mapply(RTE_calculate, L, theoretical_p, n_tc,
                               mapply(RTE_method, L, theoretical_p, n_tc))
    )]
  }

  # Add binary theoretical RTE if available
  if (has_binary_theoretical) {
    dt[loset == TRUE & valid == TRUE, `:=`(
      # Determine binary theoretical RTE calculation method
      binary_theoretical_RTE_case = mapply(RTE_method, L, binary_theoretical_p, n_tc),
      # Calculate binary theoretical RTE
      binary_theoretical_RTE = mapply(RTE_calculate, L, binary_theoretical_p, n_tc,
                                      mapply(RTE_method, L, binary_theoretical_p, n_tc))
    )]
  }

  # If verbose is FALSE, remove intermediate calculation columns
  if (!verbose) {
    # Identify intermediate calculation columns to exclude
    intermediate_cols <- c("arrivals_up_to_now", "observed_resolves_before", "L",
                           "n_tc", "last_resolve_before", "first_resolve_after")

    # Add observed columns only if they exist
    if (has_observed) {
      intermediate_cols <- c(intermediate_cols, "observed_p", "observed_RTE_case")
    }

    # Add theoretical columns if they exist
    if (has_theoretical) {
      intermediate_cols <- c(intermediate_cols, "theoretical_p",
                             "theoretical_resolve_time", "theoretical_RTE_case")
    }

    # Add binary theoretical columns if they exist
    if (has_binary_theoretical) {
      intermediate_cols <- c(intermediate_cols, "binary_theoretical_p",
                             "binary_theoretical_resolve_time", "binary_theoretical_RTE_case")
    }

    # Only try to remove columns that actually exist in the data
    cols_to_remove <- intersect(intermediate_cols, names(dt))
    if (length(cols_to_remove) > 0) {
      dt[, (cols_to_remove) := NULL]
    }
  }

  # Convert back to tibble for consistent output format
  return(as_tibble(dt))
}
