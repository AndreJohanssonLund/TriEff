#' Calculate Triage Effectiveness
#'
#' @description
#' Calculates triage effectiveness metrics including Observed Triage Effectiveness (OTE),
#' Theoretical Triage Effectiveness (TTE), and their bootstrap confidence intervals if requested.
#' Note the separate help section on convergence documentation.
#'
#' @param df Data frame containing patient data (must be initialized through init())
#' @param subgroup Optional list for subgroup analysis - an alternative to var1/2
#' @param var1 Optional string for first comparison variable
#' @param var2 Optional string for second comparison variable
#' @param bootstrap Logical or character. Whether to perform bootstrap calculations and what method to use:
#'   \itemize{
#'     \item FALSE: No bootstrapping (default)
#'     \item TRUE or "standard": Standard patient-level bootstrapping
#'     \item "segment": A custom type of Block bootstrapping using queue segments
#'   }
#' @param bootstrap_params List of bootstrap parameters:
#'   \itemize{
#'     \item sample_percentage (default: 1) How large of the sample is used? 1 = full sample. In segment bootstrapping this reffers to the percentage of samples that are used
#'     \item n_iterations (default: 2000) How many bootstrap iterations to use
#'     \item distribution_span (default: 0.95) What is the confidence interval used?
#'   }
#' @param n_workers Number of workers for parallel processing (default: detectCores() - 1)
#' @param overall_only Logical. If TRUE, only overall metrics are returned (default: FALSE)
#' @param check_convergence Logical. If true will return a list with convergence plots, if bootstrap is done (default: TRUE)
#' @param min_loset_warning Numerical, when should the function warn for low loset prevalence? (default: 5)
#' @param seed Seed for reproducible bootstrapping (default: NULL - no seed)
#'
#' @return A list containing:
#'   \itemize{
#'     \item results: Tibble with TE metrics:
#'       \itemize{
#'         \item unit: Medical specialty unit or "overall" for combined results
#'         \item n_patients: Total number of patients in the group
#'         \item n_patients_loset: Number of LOSET-positive (time-critical) patients
#'         \item loset_prevalence: Proportion of LOSET-positive patients
#'         \item mean_all: Mean waiting time for all patients
#'         \item ote_mean_loset: Mean observed waiting time for LOSET patients
#'         \item ote_te: Observed Triage Effectiveness
#'         \item tte_mean_loset: Mean theoretical waiting time for LOSET patients
#'         \item tte_te: Theoretical Triage Effectiveness
#'         \item btte_mean_loset: Mean binary theoretical waiting time for LOSET patients
#'         \item btte_te: Binary Theoretical Triage Effectiveness
#'         \item otg_te: Observed-Theoretical Gap (ote_te - tte_te)
#'         \item sensitivity: True positive rate for identifying time-critical patients
#'         \item specificity: True negative rate for identifying non-time-critical patients
#'         \item boot_ote_mean: Bootstrap mean for OTE (if bootstrap=TRUE)
#'         \item boot_ote_sd: Bootstrap standard deviation for OTE
#'         \item boot_ote_sd_q: Robust standard deviation for OTE based on IQR, Calculated as IQR/1.349
#'         \item boot_ote_var_lower: Lower bound of OTE confidence interval
#'         \item boot_ote_var_upper: Upper bound of OTE confidence interval
#'         \item boot_tte_mean: Bootstrap mean for TTE
#'         \item boot_tte_sd: Bootstrap standard deviation for TTE
#'         \item boot_tte_sd_q: Robust standard deviation for TTE based on IQR, Calculated as IQR/1.349
#'         \item boot_tte_var_lower: Lower bound of TTE confidence interval
#'         \item boot_tte_var_upper: Upper bound of TTE confidence interval
#'         \item boot_btte_mean: Bootstrap mean for BTTE
#'         \item boot_btte_sd: Bootstrap standard deviation for BTTE
#'         \item boot_btte_sd_q: Robust standard deviation for BTTE based on IQR, Calculated as IQR/1.349
#'         \item boot_btte_var_lower: Lower bound of BTTE confidence interval
#'         \item boot_btte_var_upper: Upper bound of BTTE confidence interval
#'         \item boot_mean_n_patients: Mean number of patients across bootstrap samples
#'         \item boot_mean_n_loset: Mean number of LOSET patients across bootstrap samples
#'       }
#'     \item metadata: List containing calculation parameters:
#'       \itemize{
#'         \item calculation_time: Timestamp of when calculation was performed
#'         \item bootstrap_params: List of parameters used for bootstrap (if applicable):
#'           \itemize{
#'             \item sample_percentage: Fraction of data used in each bootstrap iteration
#'             \item n_iterations: Number of bootstrap iterations performed
#'             \item distribution_span: Width of confidence intervals (e.g., 0.95 for 95% CI)
#'           }
#'         \item bootstrap_method: Type of bootstrap method use - standard or segment (if applicable)
#'         \item group_var1: First grouping variable used (if any)
#'         \item group_var2: Second grouping variable used (if any)
#'         \item subgr: Subgroup criteria used (if any)
#'       }
#'     \item bootstrap_distributions: Full bootstrap iteration data (if bootstrap=TRUE)
#'     \item convergence: Convergence analysis results (if check_convergence=TRUE)
#'   }
#'
#' @importFrom dplyr filter select mutate arrange left_join bind_rows group_by group_split
#' @importFrom dplyr summarise ungroup across group_vars group_modify n rename distinct inner_join
#' @importFrom magrittr %>%
#' @importFrom purrr map map2 map_dfr
#' @importFrom stats var sd quantile median weighted.mean
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom scales percent_format
#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot aes geom_line labs theme_minimal geom_vline geom_segment
#' @importFrom ggplot2 geom_point geom_text scale_color_manual coord_cartesian
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous ggtitle theme element_blank
#' @importFrom ggplot2 element_text element_line margin
#' @importFrom rlang sym !!
#'
#' @examples
#' \dontrun{
#' # Basic calculation
#' te_results <- calc_te(patient_data)
#'
#' # With bootstrap
#' te_results <- calc_te(patient_data,
#'                      bootstrap = TRUE)
#'
#' # With subgroup analysis
#' te_results <- calc_te(patient_data,
#'                      subgroup = list(age_group = "elderly"))
#' }
#'
#' @export
# Create metadata about the calculation
calc_wte <- function(df,
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
                    overall_only = FALSE,
                    check_convergence = TRUE,
                    min_loset_warning = 5,
                    seed = NULL
) {


  # Validate input data
  validate_te_data(df, subgroup, var1, var2)

  # Check if there's only one unit and overall_only is TRUE
  n_units <- length(unique(df$unit))
  if (n_units == 1 && overall_only) {
    # Get the unit name for the message
    unit_name <- unique(df$unit)[1]
    # Override overall_only setting
    overall_only <- FALSE
    # Inform the user
    message(sprintf("Only one unit ('%s') detected. Setting overall_only=FALSE as overall statistics are equivalent to unit statistics.", unit_name))
  }

  # Validate bootstrap parameter
  if (!is.null(bootstrap)) {
    valid_bootstrap_values <- c(TRUE, FALSE, "standard", "segment")

    if (!(is.logical(bootstrap) || is.character(bootstrap))) {
      stop("Invalid bootstrap parameter type. Must be logical (TRUE/FALSE) or character ('standard'/'segment').")
    }

    if (is.character(bootstrap) && !(bootstrap %in% valid_bootstrap_values)) {
      stop(paste0("Invalid bootstrap value: '", bootstrap, "'. Valid values are: TRUE, FALSE, 'standard', or 'segment'."))
    }

    # Handle backward compatibility - normalize values
    if (bootstrap == TRUE) {
      bootstrap <- "standard"  # Normalize to string version for consistency
    }
  }

  # Base calculations
  base_results <- calculate_te(df, subgroup, var1, var2, min_loset_warning)

  bootstrap_results = NULL
  if (bootstrap == "standard") {
    bootstrap_results <- calculate_te_bootstrap(df, bootstrap_params, n_workers,
                                                subgroup, var1, var2,
                                                min_loset_warning, seed)
  } else if (bootstrap == "segment") {
    bootstrap_results <- calculate_te_segment_bootstrap(df, bootstrap_params, n_workers,
                                                        subgroup, var1, var2,
                                                        min_loset_warning, seed)
  }


  final_results <- combine_results(base_results, bootstrap_results, var1, var2, subgroup)


  if (overall_only) {
    final_results$results <- final_results$results %>%
      filter(unit == "overall")

    if (!is.null(final_results$bootstrap_distributions)) {
      final_results$bootstrap_distributions <- final_results$bootstrap_distributions %>%
        filter(unit == "overall")
    }
  }

  if ((bootstrap == "standard" || bootstrap == "segment") && check_convergence) {
    final_results$convergence <- analyze_bootstrap_convergence(final_results,
                                                               subgroup, var1, var2)
  }

  class(final_results) <- c("calc_te", class(final_results))

  return(final_results)
}



#' Validate Triage Effectiveness Data
#'
#' @param df Data frame containing patient data
#' @param subgroup Optional list for subgroup analysis
#' @param var1 Optional string for first comparison variable
#' @param var2 Optional string for second comparison variable
#'
#' @details
#' Performs validation checks on input data including:
#' - Required columns presence
#' - Data type validation
#' - Value range validation
#' - Comparison variable validity
#'
#' @keywords internal
validate_te_data <- function(df, subgroup = NULL, var1 = NULL, var2 = NULL) {
  # First check for init
  if (!all(c("arrival_minute", "resolve_minute") %in% names(df))) {
    stop("Data frame appears uninitialized. Please run init() first.")
  }

  # Required columns
  required_cols <- c("unit", "loset", "observed_wait_time")
  if (!all(required_cols %in% names(df))) {
    stop("Missing required columns: ",
         paste(setdiff(required_cols, names(df)), collapse = ", "))
  }

  # Data type validation
  if (!is.numeric(df$observed_wait_time)) {
    stop("observed_wait_time must be numeric")
  }
  if (!is.logical(df$loset)) {
    stop("loset must be logical")
  }

  # Value validation
  if (any(df$observed_wait_time < 0, na.rm = TRUE)) {
    stop("Negative waiting times detected")
  }

  # Validate so no unit totally lacks LOSET
  n_loset <- df %>%
    group_by(unit) %>%
    summarise(
      n_loset = sum(loset)
    )

  units_no_loset <- n_loset %>%
    filter(n_loset == 0) %>%
    pull(unit)

  if (length(units_no_loset) > 0) {
    warning(sprintf("The following units have no LOSET cases: %s\nThis will generate incomplete results that may interfere with plots.\nIt is recommended to filter out these units.",
                    paste(units_no_loset, collapse = ", ")))
  }

  # Validate comparison variables
  if (!is.null(var1)) {
    if (!is.character(var1) || length(var1) != 1) {
      stop("var1 must be a single character string")
    }
    if (!var1 %in% names(df)) {
      stop(sprintf("Column '%s' specified in var1 not found in data frame", var1))
    }
  }

  if (!is.null(var2)) {
    if (is.null(var1)) {
      stop("var2 cannot be specified without var1")
    }
    if (!is.character(var2) || length(var2) != 1) {
      stop("var2 must be a single character string")
    }
    if (!var2 %in% names(df)) {
      stop(sprintf("Column '%s' specified in var2 not found in data frame", var2))
    }
  }

  # Validate subgroup
  if (!is.null(subgroup)) {
    if (!is.list(subgroup)) {
      stop("subgroup must be a list")
    }
    for (var_name in names(subgroup)) {
      if (!var_name %in% names(df)) {
        stop(sprintf("Subgroup variable '%s' not found in data frame", var_name))
      }
      if (!is.vector(subgroup[[var_name]])) {
        stop(sprintf("Subgroup values for '%s' must be a vector", var_name))
      }
    }
  }
}





#' Create Grouped Data Structure for TE Calculations
#'
#' @param df Data frame containing patient data
#' @param subgroup Optional list for subgroup analysis
#' @param var1 Optional string for first comparison variable
#' @param var2 Optional string for second comparison variable
#'
#' @return A grouped dataframe with appropriate grouping structure for TE calculations
#'
#' @details
#' Creates a grouped structure based on input parameters:
#' - Basic: Groups only by unit
#' - Subgroup: Filters by subgroup criteria, then groups by unit
#' - Single comparison: Groups by unit and one comparison variable
#' - Double comparison: Groups by unit and two comparison variables
#'
#' @keywords internal
create_grouping <- function(df, subgroup = NULL, var1 = NULL, var2 = NULL) {
  # Input validation - check for incompatible combinations
  if (!is.null(subgroup) && (!is.null(var1) || !is.null(var2))) {
    stop("Cannot combine subgroup analysis with comparison variables")
  }
  # Start with base unit grouping
  grouped_df <- df

  # Handle subgroup if specified
  if (!is.null(subgroup)) {
    # Apply each subgroup filter
    for (var in names(subgroup)) {
      grouped_df <- grouped_df %>%
        filter(.data[[var]] %in% subgroup[[var]])
    }
    # Group by unit after filtering
    return(group_by(grouped_df, unit))
  }
  # Handle comparison variable(s)
  if (!is.null(var1) && !is.null(var2)) {
    # Check variables exist in dataframe
    if (!all(c(var1, var2) %in% names(df))) {
      stop("Comparison variables not found in data frame")
    }
    return(group_by(grouped_df, unit, !!sym(var1), !!sym(var2)))
  }
  if (!is.null(var1)) {
    # Check variable exists in dataframe
    if (!var1 %in% names(df)) {
      stop("Comparison variable not found in data frame")

    }
    return(group_by(grouped_df, unit, !!sym(var1)))
  }

  # Default case: just group by unit
  return(group_by(grouped_df, unit))
}



#' Create Overall Grouping Structure
#'
#' @param results Unit-level results
#' @param subgroup Optional list for subgroup analysis
#' @param var1 Optional string for first comparison variable
#' @param var2 Optional string for second comparison variable
#'
#' @return A grouped data frame for overall calculations
#'
#' @details
#' Creates appropriate grouping structure for calculating overall metrics
#' based on the analysis type (basic, subgroup, or comparison variables).
#'
#' @keywords internal
create_overall_grouping <- function(results, subgroup = NULL, var1 = NULL, var2 = NULL) {
  # Start with ungrouped data
  grouped_data <- ungroup(results)

  # Create grouping based on comparison variables
  if (!is.null(var1) && !is.null(var2)) {
    return(group_by(grouped_data, !!sym(var1), !!sym(var2)))
  }
  if (!is.null(var1)) {
    return(group_by(grouped_data, !!sym(var1)))
  }
  if (!is.null(subgroup)) {
    return(grouped_data)  # Single group for subgroup analysis
  }

  return(grouped_data)  # No grouping for basic calculation
}





#' Calculate Reference Values for TE Calculations
#'
#' @param df Data frame containing patient data
#' @param subgroup Optional list for subgroup analysis
#' @param var1 Optional string for first comparison variable
#' @param var2 Optional string for second comparison variable
#'
#' @return A tibble containing:
#' - Grouping variables (unit and any comparison variables)
#' - mean_all: Mean waiting time for the group
#' - n_loset: Number of loset cases in the group
#'
#' @details
#' Calculates essential reference values needed for TE calculations using the
#' same grouping structure as the main calculations.
#'
#' @keywords internal
calculate_reference_values <- function(df, subgroup = NULL, var1 = NULL, var2 = NULL) {
  # Use the same grouping mechanism we created
  grouped_df <- create_grouping(df, subgroup = subgroup, var1 = var1, var2 = var2)
  # Calculate reference values for each group
  reference_values <- grouped_df %>%
    summarise(
      mean_all = mean(observed_wait_time),
      n_loset = sum(loset),
      n_patients = n(),
      loset_prev = mean(loset),
      .groups = 'drop' # Drop grouping to return a clean tibble
    )
  return(reference_values)
}



### _------------------------------------------------------------------------BASE CALC------------------------------------

#' Calculate Triage Effectiveness Metrics
#'
#' @param df A data frame containing patient data (must be initialized through init())
#' @param subgroup Optional list defining subgroup criteria
#' @param var1 Optional string naming first comparison variable
#' @param var2 Optional string naming second comparison variable
#' @param min_loset_warning Will be sent on to check if a valid amount of loset cases exists.
#'
#' @return A tibble containing:
#'   \itemize{
#'     \item Unit-level metrics (priority counts, patient counts, LOSET prevalence)
#'     \item OTE (Observed Triage Effectiveness) metrics
#'     \item TTE (Theoretical Triage Effectiveness) metrics if available
#'     \item BTTE (Binary Theoretical Triage Effectiveness) metrics if available
#'     \item Sensitivity and specificity values
#'     \item Overall metrics across all units
#'   }
#'
#' @details
#' Internal implementation function that performs the core TE calculations.
#' Creates unit-level metrics and overall metrics using weighted means
#' based on number of LOSET patients. This function is called by the
#' main exported `calc_te` function and should not be used directly.
#'
#' @keywords internal
calculate_te <- function(df, subgroup = NULL, var1 = NULL, var2 = NULL, min_loset_warning) {

  # Create groupings
  grouping <- create_grouping(df, subgroup, var1, var2)

  # Calculate unit-level metrics
  unit_results <- calculate_unit_metrics(grouping)

  # Only calculate overall if more than one unit or subgroup is not used
  if (length(unique(unit_results$unit)) > 1 & is.null(subgroup)) {
    # Create overall grouping and calculate overall metrics
    overall_grouping <- create_overall_grouping(unit_results, subgroup, var1, var2)
    overall_results <- calculate_overall_metrics(unit_results, overall_grouping)
    results <- bind_rows(overall_results, unit_results)
  } else {
    results <- unit_results
  }

  results <- results %>%
    ungroup()

  # Validate results and generate warnings
  validate_results(results, min_loset_warning)

  return(results)
}

#' Calculate Unit-Level TE Metrics
#'
#' @param grouping Grouped data frame including unit and comparison variables
#'
#' @return Tibble with unit-level TE metrics
#'
#' @details
#' Calculates unit-level metrics including:
#' - Patient counts and LOSET prevalence
#' - Mean waiting times
#' - OTE, TTE, and BTTE metrics
#' - Sensitivity and specificity
#'
#' @keywords internal
calculate_unit_metrics <- function(grouping) {
  grouping %>%
    summarise(
      # Patient counts
      n_patients = n(),
      n_patients_loset = sum(loset, na.rm = TRUE),
      loset_prevalence = mean(loset, na.rm = TRUE),

      # Calculate mean waiting times
      mean_all = mean(observed_wait_time, na.rm = TRUE),

      # Calculate OTE
      ote_mean_loset = mean(observed_wait_time[loset], na.rm = TRUE),
      ote_te = 1 - (ote_mean_loset / mean_all),

      # Calculate TTE if available
      tte_mean_loset = if("theoretical_wait_time" %in% names(.) &&
                          all(!is.na(theoretical_wait_time)))
        mean(theoretical_wait_time[loset], na.rm = TRUE) else NA_real_,
      tte_te = if(!is.na(tte_mean_loset))
        1 - (tte_mean_loset / mean_all) else NA_real_,

      # Calculate BTTE if available
      btte_mean_loset = if("binary_theoretical_wait_time" %in% names(.) &&
                           all(!is.na(binary_theoretical_wait_time)))
        mean(binary_theoretical_wait_time[loset], na.rm = TRUE) else NA_real_,
      btte_te = if(!is.na(btte_mean_loset))
        1 - (btte_mean_loset / mean_all) else NA_real_,

      # Calculate sensitivity and specificity
      sensitivity = sum(tp, na.rm = TRUE) / (sum(tp, na.rm = TRUE) + sum(fn, na.rm = TRUE)),
      specificity = sum(tn, na.rm = TRUE) / (sum(tn, na.rm = TRUE) + sum(fp, na.rm = TRUE)),

      .groups = 'keep'
    ) %>%
    mutate(otg_te = if_else(!is.na(tte_te), ote_te - tte_te, NA_real_))
}


#' Calculate Overall TE Metrics
#'
#' @param unit_results Results from unit-level calculations
#' @param overall_grouping Grouped data frame for overall calculations
#'
#' @return Tibble with overall TE metrics
#'
#' @details
#' Calculates overall metrics using weighted means based on
#' number of LOSET patients in each unit.
#'
#' @keywords internal
calculate_overall_metrics <- function(unit_results, overall_grouping) {
  overall_grouping %>%
    summarise(
      unit = "overall",
      n_patients_overall = sum(n_patients),
      n_patients_loset_overall = sum(n_patients_loset),
      loset_prevalence_overall = sum(n_patients_loset) / sum(n_patients),

      # Weighted means for TE metrics using n_patients_loset as weights
      mean_all_overall = weighted.mean(mean_all, n_patients, na.rm = TRUE),
      ote_mean_loset_overall = weighted.mean(ote_mean_loset, n_patients_loset, na.rm = TRUE),
      ote_te_overall = weighted.mean(ote_te, n_patients_loset, na.rm = TRUE),

      # Weighted means for TE metrics
      tte_mean_loset_overall = if("tte_mean_loset" %in% names(.) &&
                                  all(!is.na(tte_mean_loset)))
        weighted.mean(tte_mean_loset, n_patients_loset, na.rm = TRUE) else NA_real_,
      tte_te_overall = if(!is.na(tte_mean_loset_overall))
        weighted.mean(tte_te, n_patients_loset, na.rm = TRUE) else NA_real_,

      btte_mean_loset_overall = if("btte_mean_loset" %in% names(.) &&
                                   all(!is.na(btte_mean_loset)))
        weighted.mean(btte_mean_loset, n_patients_loset, na.rm = TRUE) else NA_real_,
      btte_te_overall = if(!is.na(btte_mean_loset_overall))
        weighted.mean(btte_te, n_patients_loset, na.rm = TRUE) else NA_real_,

      sensitivity_overall = weighted.mean(sensitivity, n_patients_loset, na.rm = TRUE),
      specificity_overall = weighted.mean(specificity, n_patients_loset, na.rm = TRUE),

      .groups = 'drop'
    ) %>%
    mutate(
      otg_te_overall = if_else(!is.na(tte_te_overall),
                               ote_te_overall - tte_te_overall,
                               NA_real_)
    ) %>%
    rename(
      n_patients = n_patients_overall,
      n_patients_loset = n_patients_loset_overall,
      loset_prevalence = loset_prevalence_overall,
      mean_all = mean_all_overall,
      ote_mean_loset = ote_mean_loset_overall,
      ote_te = ote_te_overall,
      tte_mean_loset = tte_mean_loset_overall,
      tte_te = tte_te_overall,
      btte_mean_loset = btte_mean_loset_overall,
      btte_te = btte_te_overall,
      otg_te = otg_te_overall,
      sensitivity = sensitivity_overall,
      specificity = specificity_overall
    )
}



#' Validate Results and Generate Warnings
#'
#' @param results Combined unit and overall results
#' @param min_loset_warning Least number of LOSET cases that should exists or the function will warn.
#'
#' @return Invisible NULL, generates warnings if issues found
#'
#' @details
#' Performs validation checks on calculated results:
#' - Low LOSET counts
#' - Extreme TE values
#'
#' @keywords internal
validate_results <- function(results, min_loset_warning) {
  # Check for low loset counts
  low_loset <- results %>%
    filter(unit != "overall", n_patients_loset < min_loset_warning)

  if (nrow(low_loset) > 0) {
    warning("Units with fewer than specifiedd low_loset found: ",
            paste(low_loset$unit, collapse = ", "))
  }

  # Check for extreme TE values
  if (any(results$ote_te > 1, na.rm = TRUE)) {
    warning("OTE values greater than 1 detected")
  }
  if (any(results$tte_te > 1, na.rm = TRUE)) {
    warning("TTE values greater than 1 detected")
  }
  if (any(results$btte_te > 1, na.rm = TRUE)) {
    warning("BTTE values greater than 1 detected")
  }

  invisible(NULL)
}


#--------------------------------------------------------- boootstrapping-----------------------------

#' Calculate Triage Effectiveness with Bootstrap
#'
#' @param df Data frame containing patient data
#' @param bootstrap_params List of bootstrap parameters
#' @param n_workers Number of workers for parallel processing
#' @param subgroup Optional list for subgroup analysis
#' @param var1 Optional string for first comparison variable
#' @param var2 Optional string for second comparison variable
#' @param min_loset_warning will be sent onward to check if a valid amount of loset cases exists
#' @param seed Seed for reproducible results
#'
#' @return List containing metrics, distributions, and convergence information
#'
#' @details
#' Performs bootstrap analysis of TE metrics with parallel processing
#' and progress tracking.
#'
#' @keywords internal
calculate_te_bootstrap <- function(df, bootstrap_params, n_workers,
                                   subgroup, var1, var2, min_loset_warning,
                                   seed) {

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Calculate sample size based on percentage
  n_samples <- floor(nrow(df) * bootstrap_params$sample_percentage)

  # Set up parallel processing if needed
  setup_parallel(n_workers = n_workers)
  on.exit(cleanup_parallel())

  # Initialize progress reporting
  init_progressr()


  # Perform bootstrap iterations with progress bar
  print(paste("Starting bootstrap iterations", Sys.time()))
  bootstrap_results <- with_progress({
    p <- progressor(steps = bootstrap_params$n_iterations)
    future_map(1:bootstrap_params$n_iterations, function(i) {
      # Sample with replacement
      boot_sample <- df[sample(nrow(df), size = n_samples, replace = TRUE), ]
      # Calculate TE using existing functions
      result <- calculate_te(boot_sample, subgroup, var1, var2, min_loset_warning)
      p()
      result
    }, .options = furrr_options(seed = TRUE))
  })
  print(paste("Bootstrap iterations done", Sys.time()))

  # Convert list of results to a more manageable format
  results_df <- bind_rows(bootstrap_results, .id = "iteration")

  # Create grouped structure based on input parameters
  results <- create_grouping(results_df, subgroup, var1, var2)
  return(calculate_bootstrap_metrics(results, results_df, bootstrap_params, "standard"))

}


#' Calculate Triage Effectiveness with Segment Bootstrap
#'
#' @param df Data frame containing patient data
#' @param bootstrap_params List of bootstrap parameters
#' @param n_workers Number of workers for parallel processing
#' @param subgroup Optional list for subgroup analysis
#' @param var1 Optional string for first comparison variable
#' @param var2 Optional string for second comparison variable
#' @param min_loset_warning Minimum number of LOSET cases before warning is issued
#' @param seed Seed for reproducible results
#'
#' @return List containing metrics, distributions, and method information
#'
#' @details
#' Performs segment-based bootstrap analysis of TE metrics, where entire queue
#' segments are resampled rather than individual patients.
#'
#' @keywords internal
calculate_te_segment_bootstrap <- function(df, bootstrap_params, n_workers,
                                           subgroup, var1, var2, min_loset_warning,
                                           seed) {

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

  print(paste("Starting segment bootstrap with", n_segments, "unit-segment pairs,",
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

      # Calculate TE using existing functions
      result <- calculate_te(boot_sample, subgroup, var1, var2, min_loset_warning)
      p()
      result
    }, .options = furrr_options(seed = TRUE))
  })

  print(paste("Segment bootstrap iterations done", Sys.time()))

  # Convert list of results to a more manageable format
  results_df <- bind_rows(bootstrap_results, .id = "iteration")

  # Create grouped structure based on input parameters
  results <- create_grouping(results_df, subgroup, var1, var2)

  # Calculate variation metrics using existing grouping
  return(calculate_bootstrap_metrics(results, results_df, bootstrap_params, "segment"))
}



#' Calculate Bootstrap Metrics for Triage Effectiveness
#'
#' @description
#' This function performs bootstrap analysis on Triage Effectiveness data, calculating
#' confidence intervals and statistical properties of TE metrics.
#'
#' @param results A list containing the complete calculation results
#' @param results_df A data frame containing the summarized TE results
#' @param bootstrap_params A list containing bootstrap parameters:
#'   \itemize{
#'     \item sample_percentage: Percentage of data to sample in each iteration
#'     \item n_iterations: Number of bootstrap iterations
#'     \item distribution_span: Width of confidence interval
#'   }
#' @param method Character string specifying the bootstrap method ("standard" or "segment")
#'
#' @return A modified results data frame with added bootstrap metrics
#'
#' @keywords internal
calculate_bootstrap_metrics <- function(results, results_df, bootstrap_params, method) {
  # Calculate confidence interval boundaries
  alpha <- (1 - bootstrap_params$distribution_span) / 2

  # Calculate variation metrics
  variation_metrics <- results %>%
    summarise(
      # OTE metrics
      boot_ote_mean = mean(ote_te, na.rm = TRUE),
      boot_ote_sd = sd(ote_te, na.rm = TRUE),
      boot_ote_sd_q = (quantile(ote_te, 0.75, na.rm = TRUE) -
                         quantile(ote_te, 0.25, na.rm = TRUE)) / 1.349,
      boot_ote_var_lower = quantile(ote_te, probs = alpha, na.rm = TRUE),
      boot_ote_var_upper = quantile(ote_te, probs = 1 - alpha, na.rm = TRUE),

      # TTE metrics (if available)
      boot_tte_mean = mean(tte_te, na.rm = TRUE),
      boot_tte_sd = sd(tte_te, na.rm = TRUE),
      boot_tte_sd_q = (quantile(tte_te, 0.75, na.rm = TRUE) -
                         quantile(tte_te, 0.25, na.rm = TRUE)) / 1.349,
      boot_tte_var_lower = quantile(tte_te, probs = alpha, na.rm = TRUE),
      boot_tte_var_upper = quantile(tte_te, probs = 1 - alpha, na.rm = TRUE),

      # BTTE metrics (if available)
      boot_btte_mean = mean(btte_te, na.rm = TRUE),
      boot_btte_sd = sd(btte_te, na.rm = TRUE),
      boot_btte_sd_q = (quantile(btte_te, 0.75, na.rm = TRUE) -
                          quantile(btte_te, 0.25, na.rm = TRUE)) / 1.349,
      boot_btte_var_lower = quantile(btte_te, probs = alpha, na.rm = TRUE),
      boot_btte_var_upper = quantile(btte_te, probs = 1 - alpha, na.rm = TRUE),

      # OTE metrics (if available)
      boot_otg_mean = mean(ote_te - tte_te, na.rm = TRUE),
      boot_otg_sd = sd(ote_te - tte_te, na.rm = TRUE),
      boot_otg_sd_q = (quantile(ote_te - tte_te, 0.75, na.rm = TRUE) -
                         quantile(ote_te - tte_te, 0.25, na.rm = TRUE)) / 1.349,
      boot_otg_var_lower = quantile(ote_te - tte_te, probs = alpha, na.rm = TRUE),
      boot_otg_var_upper = quantile(ote_te - tte_te, probs = 1 - alpha, na.rm = TRUE),

      # Basic statistics
      boot_mean_n_patients = mean(n_patients, na.rm = TRUE),
      boot_mean_n_loset = mean(n_patients_loset, na.rm = TRUE),
      .groups = 'drop'
    )

  variation_metrics <- variation_metrics %>%
    ungroup()

  # Return structured results
  list(
    metrics = variation_metrics,
    distributions = results_df,
    params = bootstrap_params,
    method = method
  )
}



#' Combine Base Results with Bootstrap Results
#'
#' @param base_results Tibble with base TE calculations
#' @param bootstrap_results List containing bootstrap metrics and distributions
#' @param var1 If applied, var1 is added.
#' @param var2 If applied, var2 is added.
#' @param subgroup If applied, the subgroup chain is added.
#'
#' @return List with combined results and metadata
#'
#' @details
#' Merges base calculations with bootstrap results and creates
#' metadata about the calculation process.
#'
#' @keywords internal
combine_results <- function(base_results, bootstrap_results, var1, var2, subgroup) {

  if(inherits(bootstrap_results$metrics, "tbl_df")) {
    # Get the key columns for joining
    key_cols <- intersect(
      names(base_results),
      names(bootstrap_results$metrics)
    )

    # Join base results with bootstrap metrics
    combined_metrics <- base_results %>%
      left_join(
        bootstrap_results$metrics %>%
          select(starts_with("boot_"), all_of(key_cols)),
        by = key_cols
      )
  } else {
    combined_metrics <- base_results
  }

  # Create metadata about the calculation
  metadata <- list(
    calculation_time = Sys.time(),
    calculation_method = "waiting-time-based",
    bootstrap_params = bootstrap_results$params,
    bootstrap_method = bootstrap_results$method,
    group_var1 = var1,
    group_var2 = var2,
    subgr = subgroup,
    ci_method = if(inherits(bootstrap_results$metrics, "tbl_df")) "bootstrap" else "none",
    recommended_ci_vars = if(inherits(bootstrap_results$metrics, "tbl_df")) "boot" else "none" # Change here
  )

  # Return structured output
  if(inherits(bootstrap_results$metrics, "tbl_df")) {
    combined_results <- list(
      results = combined_metrics,
      bootstrap_distributions =  bootstrap_results$distributions,
      metadata = metadata
    )
  } else {
    combined_results <- list(
      results = combined_metrics,
      metadata = metadata
    )
  }

  return(combined_results)
}


#### ------------------------------------------------------- convergence check!-------------------------------------


#' Analyze Bootstrap Convergence
#'
#' @param results List containing bootstrap results from calc_te
#' @param subgroup Optional list defining subgroup criteria
#' @param var1 Optional string for first comparison variable
#' @param var2 Optional string for second comparison variable
#'
#' @return List containing:
#'   \itemize{
#'     \item groups: Tibble with group-level convergence results
#'     \item summary: Overall convergence statistics
#'     \item metrics_analyzed: Vector of analyzed TE metrics
#'     \item calculation_time: Timestamp of analysis
#'   }
#'
#' @details
#' Analyzes convergence stability of bootstrap iterations by calculating
#' running means, standard errors, and confidence intervals for each group.
#'
#' @keywords internal
analyze_bootstrap_convergence <- function(results,
                                          subgroup = NULL, var1 = NULL, var2 = NULL) {

  # Extract bootstrap distributions
  dist_df <- results$bootstrap_distributions

  # Use create_grouping to get the correct grouping structure
  grouped_df <- create_grouping(dist_df, subgroup, var1, var2)

  # Get available TE metrics
  metric_cols <- c("ote_te", "tte_te", "btte_te")
  available_metrics <- intersect(metric_cols, names(dist_df))

  # Calculate convergence for each group
  convergence_by_group <- grouped_df %>%
    group_by(across(group_vars(.))) %>%
    group_modify(~{
      iter_data <- .x %>% arrange(iteration)

      # Results structure for this group
      group_results <- list()

      # Get group info from .y
      group_info <- paste(
        map2(names(.y), .y, ~paste(.x, ":", .y)),
        collapse = ", "
      )

      for (metric in available_metrics) {
        if (all(is.na(iter_data[[metric]]))) next

        values <- iter_data[[metric]]
        iterations <- seq_len(length(values))

        running_stats <- tibble(
          iteration = iterations,
          mean = cumsum(values) / iterations,
          var = map_dbl(iterations, ~if(.x == 1) 0 else var(values[1:.x], na.rm = TRUE))
        ) %>%
          mutate(
            se = sqrt(var / iteration),
            ci_width = se * 1.96
          )

        # Then when creating plots:
        means_plot <- ggplot2::ggplot(running_stats, ggplot2::aes(x = iteration, y = mean)) +
          ggplot2::geom_line(color = "blue") +
          ggplot2::labs(x = "Iteration", y = "Running Mean",
                        title = paste(metric, "- Running Mean, for"),
                        subtitle = group_info)

        ci_plot <- ggplot2::ggplot(running_stats, ggplot2::aes(x = iteration, y = ci_width)) +
          ggplot2::geom_line(color = "red") +
          ggplot2::labs(x = "Iteration", y = "CI Width",
                        title = paste(metric, "- CI Width, for"),
                        subtitle = group_info)

        se_plot <- ggplot2::ggplot(running_stats, ggplot2::aes(x = iteration, y = se)) +
          ggplot2::geom_line(color = "purple") +
          ggplot2::labs(x = "Iteration", y = "Standard Error",
                        title = paste(metric, "- Standard Error for"),
                        subtitle = group_info)

        # Calculate convergence metrics
        final_values <- tail(running_stats, 100)  # Look at last 100 iterations
        mean_change <- abs(diff(range(final_values$mean))) / mean(final_values$mean)
        se_change <- abs(diff(range(final_values$se))) / mean(final_values$se)


        group_results[[metric]] <- list(
          running_stats = running_stats,
          final_mean = tail(running_stats$mean, 1),
          final_se = tail(running_stats$se, 1),
          final_ci_width = tail(running_stats$ci_width, 1),
          n_iterations = length(iterations),
          mean_change = mean_change,
          se_change = se_change,
          plots = list(
            means = means_plot,
            ci_width = ci_plot,
            se = se_plot
          )
        )
      }

      # Return results for this group
      tibble(
        metric_results = list(group_results),
        n_iterations = length(unique(iter_data$iteration)),
        n_patients = nrow(iter_data)
      )
    })

  # Calculate summary statistics
  summary_stats <- convergence_by_group %>%
    ungroup() %>%
    summarise(
      total_groups = n(),
      mean_iterations = mean(n_iterations),
      total_patients = sum(n_patients)
    )

  # Return structured results
  conv_results <- list(
    groups = convergence_by_group,
    summary = summary_stats,
    metrics_analyzed = available_metrics,
    calculation_time = Sys.time()
  )

  # Add class for convergence results
  structure(conv_results, class = c("te_convergence", "list"))

}


#' Print Method for TE Convergence Results
#'
#' @description
#' Displays convergence analysis results, showing key metadata and all diagnostic plots
#' for assessing bootstrap convergence stability.
#'
#' @param x A te_convergence object containing convergence analysis results
#' @param ... Additional arguments passed to print (not currently used)
#'
#' @return Invisibly returns the input object
#'
#' @importFrom ggplot2 ggtitle theme_minimal
#' @importFrom patchwork wrap_plots plot_annotation
#'
#' @method print te_convergence
#' @export
print.te_convergence <- function(x, ...) {
  # Minimal useful metadata or we could remove this entirely
  cat("\nConvergence Analysis Summary\n")
  cat("=========================\n")
  cat(sprintf("Metrics analyzed: %s\n\n", paste(x$metrics_analyzed, collapse = ", ")))

  # Display plots with better labeling
  for (i in seq_len(nrow(x$groups))) {
    group_row <- x$groups[i, ]
    metric_results <- group_row$metric_results[[1]]

    for (metric in names(metric_results)) {
      # Add metric type to plot titles
      means_plot <- metric_results[[metric]]$plots$means
      ci_plot <- metric_results[[metric]]$plots$ci_width
      se_plot <- metric_results[[metric]]$plots$se

      # Combine and print plots without showing list structure
      combined_plots <- wrap_plots(means_plot, ci_plot, se_plot)
      print(combined_plots)
    }
  }

  invisible(x)
}

#' Plot Method for TE Convergence Results
#'
#' @description
#' Alternative method to display convergence analysis results. Functions identically
#' to print.te_convergence for consistency in output display.
#'
#' @param x A te_convergence object containing convergence analysis results
#' @param ... Additional arguments passed to plot (not currently used)
#'
#' @return Invisibly returns the input object
#'
#' @method plot te_convergence
#' @export
plot.te_convergence <- function(x, ...) {
  # Explicitly call our print method
  print.te_convergence(x)
  invisible(x)
}



####----------------------------------------- support function to print the result
#' Print Method for Triage Effectiveness Results
#'
#' @description
#' Formats and displays triage effectiveness analysis results in a readable format.
#' Shows metrics for each analysis group including sample sizes, classification metrics,
#' and triage effectiveness values. If bootstrap analysis was performed, also displays
#' confidence intervals.
#'
#' @param x A calc_te object containing:
#'   \itemize{
#'     \item results: Tibble with TE metrics
#'     \item metadata: List with calculation parameters
#'     \item bootstrap_distributions: Optional bootstrap results
#'   }
#' @param ... Additional arguments passed to print (not currently used)
#'
#' @details
#' The output is organized in sections:
#' 1. Results for each group showing:
#'    - Total patients and LOSET prevalence
#'    - Classification metrics (sensitivity, specificity)
#'    - Triage Effectiveness metrics (OTE, TTE, OTG)
#'    - Bootstrap variation intervals (if available)
#' 2. Computation information including:
#'    - Calculation method (bootstrap or direct)
#'    - Number of iterations (if bootstrap)
#'    - Sample percentage (if bootstrap)
#'    - Calculation timestamp
#'
#' @method print calc_te
#' @export
print.calc_te <- function(x, ...) {
  # Helper function to format percentages with CI
  format_metric_ci <- function(mean, lower, upper, lower_np = NA, upper_np = NA, shapiro_p = NA) {
    if (is.na(mean)) return("NA")

    # Determine whether to use non-parametric CI based on normality test
    use_nonparametric <- FALSE
    if (!is.na(shapiro_p) && !is.na(lower_np) && !is.na(upper_np)) {
      use_nonparametric <- shapiro_p < 0.05
    }

    formatted_mean <- sprintf("%.1f%%", mean * 100)

    if (use_nonparametric) {
      # Use non-parametric CI
      return(sprintf("%s (%.1f%% to %.1f%%) (non-parametric)",
                     formatted_mean, lower_np * 100, upper_np * 100))
    } else if (!is.na(lower) && !is.na(upper)) {
      # Use parametric CI
      return(sprintf("%s (%.1f%% to %.1f%%)",
                     formatted_mean, lower * 100, upper * 100))
    }

    # Fallback to just the mean if CI not available
    return(formatted_mean)
  }

  # Extract results and metadata
  results_df <- x$results
  metadata <- x$metadata

  cat("\nTriage Effectiveness Analysis Results\n")
  cat("===================================\n\n")

  # Determine if we're displaying RTE or WTE
  is_rte <- !is.null(metadata$calculation_method) && metadata$calculation_method == "rank-based"

  # Process rows in their original order
  for(i in 1:nrow(results_df)) {
    row <- results_df[i,]

    # Get grouping variables for this row
    group_vars <- c("unit")
    if (!is.null(metadata$group_var1)) group_vars <- c(group_vars, metadata$group_var1)
    if (!is.null(metadata$group_var2)) group_vars <- c(group_vars, metadata$group_var2)

    group_values <- row[, intersect(group_vars, names(row)), drop = FALSE]

    # Create header label
    if (group_values$unit == "overall") {
      header <- "Results for Overall"

      # Add grouping variables if present
      if (!is.null(metadata$group_var1) && metadata$group_var1 %in% names(group_values)) {
        header <- paste(header, group_values[[metadata$group_var1]], sep = ", ")
      }
      if (!is.null(metadata$group_var2) && metadata$group_var2 %in% names(group_values)) {
        header <- paste(header, group_values[[metadata$group_var2]], sep = ", ")
      }
    } else {
      header <- paste("Results for unit:", group_values$unit)

      # Add grouping variables if present
      if (!is.null(metadata$group_var1) && metadata$group_var1 %in% names(group_values)) {
        header <- paste(header, group_values[[metadata$group_var1]], sep = ", ")
      }
      if (!is.null(metadata$group_var2) && metadata$group_var2 %in% names(group_values)) {
        header <- paste(header, group_values[[metadata$group_var2]], sep = ", ")
      }
    }

    cat(header, "\n")
    cat(paste(rep("-", nchar(header)), collapse = ""), "\n")

    # Print sample sizes with correct LOSET percentage
    if (is_rte && "n_valid_rte" %in% names(row)) {
      cat(sprintf("Total patients: %d (%.1f%% LOSET positive, %d valid for RTE)\n",
                  row$n_patients,
                  row$loset_prevalence * 100,
                  row$n_valid_rte))
    } else {
      cat(sprintf("Total patients: %d (%.1f%% LOSET positive)\n",
                  row$n_patients,
                  row$loset_prevalence * 100))
    }

    # Print classification metrics
    cat("\nClassification Metrics:\n")
    cat(sprintf("  Sensitivity: %.1f%%\n", row$sensitivity * 100))
    cat(sprintf("  Specificity: %.1f%%\n", row$specificity * 100))

    # Print TE metrics with appropriate header
    if (is_rte) {
      cat("\nRank-based Triage Effectiveness Metrics:\n")
    } else {
      cat("\nTriage Effectiveness Metrics:\n")
    }

    cat(sprintf("  OTE: %.1f%%\n", row$ote_te * 100))

    if (!is.na(row$tte_te)) {
      cat(sprintf("  TTE: %.1f%%\n", row$tte_te * 100))
      cat(sprintf("  OTG: %.1f%%\n", row$otg_te * 100))
    }

    if ("btte_te" %in% names(row) && !is.na(row$btte_te)) {
      cat(sprintf("  BTTE: %.1f%%\n", row$btte_te * 100))
    }

    # Determine if we should show confidence intervals
    show_ci <- !(is.null(metadata$recommended_ci_vars) || metadata$recommended_ci_vars == "none")

    if (show_ci) {
      # Print confidence intervals
      cat("\nConfidence Intervals (95%):\n")

      # Determine which CI variables to use based on metadata
      ci_vars_prefix <- if (!is.null(metadata$recommended_ci_vars)) {
        metadata$recommended_ci_vars
      } else {
        "parametric"  # Default fallback
      }

      # Determine parametric CI variables
      par_lower_ote <- paste0(ci_vars_prefix, "_ote_var_lower")
      par_upper_ote <- paste0(ci_vars_prefix, "_ote_var_upper")
      par_lower_tte <- paste0(ci_vars_prefix, "_tte_var_lower")
      par_upper_tte <- paste0(ci_vars_prefix, "_tte_var_upper")
      par_lower_btte <- paste0(ci_vars_prefix, "_btte_var_lower")
      par_upper_btte <- paste0(ci_vars_prefix, "_btte_var_upper")

      # Default for non-parametric CI variables
      np_lower_ote <- "nonparametric_ote_var_lower"
      np_upper_ote <- "nonparametric_ote_var_upper"
      np_lower_tte <- "nonparametric_tte_var_lower"
      np_upper_tte <- "nonparametric_tte_var_upper"
      np_lower_btte <- "nonparametric_btte_var_lower"
      np_upper_btte <- "nonparametric_btte_var_upper"

      # Get values or NA if not present
      get_value <- function(df, col) {
        if (col %in% names(df)) df[[col]] else NA
      }

      # Format OTE CI
      cat("  OTE: ", format_metric_ci(
        row$ote_te,
        get_value(row, par_lower_ote),
        get_value(row, par_upper_ote),
        get_value(row, np_lower_ote),
        get_value(row, np_upper_ote),
        get_value(row, "ote_shapiro_p")
      ), "\n")

      # Format TTE CI if available
      if (!is.na(row$tte_te)) {
        cat("  TTE: ", format_metric_ci(
          row$tte_te,
          get_value(row, par_lower_tte),
          get_value(row, par_upper_tte),
          get_value(row, np_lower_tte),
          get_value(row, np_upper_tte),
          get_value(row, "tte_shapiro_p")
        ), "\n")

        # OTG has no direct non-parametric CI
        cat("  OTG: ", format_metric_ci(
          row$otg_te,
          NA, NA, NA, NA, NA  # No CIs for OTG
        ), "\n")
      }

      # Format BTTE CI if available
      if ("btte_te" %in% names(row) && !is.na(row$btte_te)) {
        cat("  BTTE: ", format_metric_ci(
          row$btte_te,
          get_value(row, par_lower_btte),
          get_value(row, par_upper_btte),
          get_value(row, np_lower_btte),
          get_value(row, np_upper_btte),
          get_value(row, "btte_shapiro_p")
        ), "\n")
      }
    }

    cat("\n")  # Add spacing between rows
  }

  # Print computation information footer
  cat("Computation Information\n")
  cat("=====================\n")

  # Check if this is rank-based or waiting-time-based
  is_rank_based <- !is.null(metadata$calculation_method) &&
    metadata$calculation_method == "rank-based"

  if (is_rank_based) {
    cat("Method: Rank-based triage effectiveness calculation\n")
  } else if (!is.null(metadata$bootstrap_params)) {
    cat("Method: Bootstrap calculation (waiting-time-based)\n")
    cat(sprintf("Number of iterations: %d\n", metadata$bootstrap_params$n_iterations))
    cat(sprintf("Sample percentage: %.0f%%\n",
                metadata$bootstrap_params$sample_percentage * 100))
    cat(sprintf("Distribution span: %.0f%%\n",
                metadata$bootstrap_params$distribution_span * 100))
  } else {
    cat("Method: Direct calculation (waiting-time-based)\n")
  }

  cat(sprintf("Calculation time: %s\n", metadata$calculation_time))

  # Return invisibly
  invisible(x)
}
