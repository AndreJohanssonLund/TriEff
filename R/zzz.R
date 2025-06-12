#' Internal package utilities and initialization
#'
#' @keywords internal
#' @import progressr
#' @importFrom progressr handlers progressor
"_PACKAGE"

#' Internal package global variables
#'
#' @title Internal package global variables
#' @name trieff-globals
#' @keywords internal
NULL
utils::globalVariables(c(
  # Basic identifiers and metadata
  "id", "unit", ".", ".data", "all_of", "any_of", "everything", "lag", "matches",
  "runif", "starts_with", "segment_id",

  # Time-related variables
  "arrival", "resolve", "arrival_minute", "resolve_minute", "event_minute",

  # Priority and classification
  "priority", "priority_binary", "loset", "tp", "fp", "fn", "tn", "has_loset",

  # Wait time calculations
  "observed_wait_time", "theoretical_wait_time", "binary_theoretical_wait_time",
  "mean_all", "mean_loset", "ote_mean_loset", "tte_mean_loset", "btte_mean_loset",

  # Processing order variables
  "resolve_time", "resolve_event_order", "resolve_order", "binary_resolve_order",
  "theoretical_resolve_time", "binary_theoretical_resolve_time",

  # Statistics and calculations
  "sensitivity", "specificity", "n_patients", "n_patients_loset",
  "loset_prevalence", "ote_te", "tte_te", "btte_te", "otg_te",
  "btte_mean_loset_overall", "btte_te_overall", "loset_prevalence_overall",
  "mean_all_overall", "n_patients_loset_overall", "n_patients_overall",
  "ote_mean_loset_overall", "ote_te_overall", "otg_te_overall",
  "sensitivity_overall", "specificity_overall", "tte_mean_loset_overall",
  "tte_te_overall",

  # Bootstrap related
  "boot_ote_mean", "boot_ote_sd", "boot_ote_sd_q", "boot_ote_var_lower", "boot_ote_var_upper",
  "boot_tte_mean", "boot_tte_sd", "boot_tte_sd_q", "boot_tte_var_lower", "boot_tte_var_upper",
  "boot_btte_mean", "boot_btte_sd", "boot_btte_sd_q", "boot_btte_var_lower", "boot_btte_var_upper",
  "boot_mean_n_patients", "boot_mean_n_loset", "iteration", "n_iterations", "include_distributions",
  "boot_sample_size",

  # RTE specific variables
  "observed_RTE", "theoretical_RTE", "binary_theoretical_RTE", "valid", "n_distinct",
  "n_valid_rte", "observed_p", "theoretical_p", "binary_theoretical_p", "n_tc", "original_df",

  # RTE calculation components
  "arrivals_up_to_now", "observed_resolves_before", "L", "is_cluster", "cluster_id",
  "last_resolve_before", "first_resolve_after", "window_id", "carryover",
  "adjusted_n_tc", "window_end", "loset_arrivals", "min_arrival", "resolve_times",
  "n_resolves", "segment_has_non_loset", "facet_wrap",

  # RTE statistical metrics
  "ote_se", "tte_se", "btte_se", "t_value",
  "observed_RTE_mean", "observed_RTE_sd", "observed_RTE_q025", "observed_RTE_q975", "observed_RTE_sd_q",
  "theoretical_RTE_mean", "theoretical_RTE_sd", "theoretical_RTE_q025", "theoretical_RTE_q975", "theoretical_RTE_sd_q",
  "binary_theoretical_RTE_mean", "binary_theoretical_RTE_sd", "binary_theoretical_RTE_q025", "binary_theoretical_RTE_q975", "binary_theoretical_RTE_sd_q",
  "n_samples", "observed_RTE_case", "theoretical_RTE_case", "binary_theoretical_RTE_case", ".valid_observed_RTE", ".valid_theoretical_RTE", ".valid_binary_theoretical_RTE",

  # Parametric and nonparametric stats
  "parametric_ote_mean", "parametric_ote_sd", "parametric_ote_sd_q", "parametric_ote_var_lower", "parametric_ote_var_upper",
  "parametric_tte_mean", "parametric_tte_sd", "parametric_tte_sd_q", "parametric_tte_var_lower", "parametric_tte_var_upper",
  "parametric_btte_mean", "parametric_btte_sd", "parametric_btte_sd_q", "parametric_btte_var_lower", "parametric_btte_var_upper",
  "nonparametric_ote_var_lower", "nonparametric_ote_var_upper",
  "nonparametric_tte_var_lower", "nonparametric_tte_var_upper",
  "nonparametric_btte_var_lower", "nonparametric_btte_var_upper",
  "ote_shapiro_p", "tte_shapiro_p", "btte_shapiro_p", "recommended_ci_type",

  # All window related
  "all_resolves", "all_patients", "candidates", "arrival_t", "resolve_t", "current_priority",
  "my_segment", "candidate_times", "candidate_priorities", "candidate_arrivals",
  "order_index", "my_pos", "window_stats", "seg_windows", "unresolved_tc", "L_min_p",
  ".N",

  # Event processing
  "event_type", "cum_arrivals", "cum_resolves", "queue_size", "segment_start", "segment",

  # Plotting related
  "x", "y", "y_pos", "te", "label", "color", "group_label", "random_value",

  # Transferability-related variables
  "applied_unit", "original_unit", "combination", "diff_from_original",
  "loset_count", "non_loset_count", "mean_rte", "mean_diff", "lower_ci", "upper_ci",
  "running_mean", "running_var", "running_se", "ci_width",

  # Heat map related
  "te_label", "normalized_te", "un_normalized_te", "target_100_100", "original_100_100"
))



# Package environment for storing internal state
.trieff_env <- new.env(parent = emptyenv())

#' Determine Optimal Parallel Backend
#'
#' @description
#' Determines the most suitable parallel processing backend based on the system environment.
#' Returns "multisession" for RStudio environments and Windows systems, and "multicore" for
#' other systems where available.
#'
#' @return Character string specifying the optimal backend ("multisession" or "multicore")
#' @keywords internal
get_optimal_backend <- function() {
  # Always use multisession in RStudio
  if (Sys.getenv("RSTUDIO") == "1") {
    return("multisession")
  }

  # For non-RStudio environments
  if (.Platform$OS.type != "windows" && requireNamespace("parallel", quietly = TRUE)) {
    "multicore"
  } else {
    "multisession"
  }
}

#' Setup Parallel Processing Environment
#'
#' @description
#' Configures the parallel processing environment using the future package.
#' Automatically selects and sets up the appropriate backend based on the system
#' environment and user settings. When using multisession, uses the callr backend
#' for better cleanup behavior on forced termination.
#'
#' @param n_workers Integer. Number of cores to use for parallel processing.
#'   Defaults to the value set in trieff.cores option.
#'
#' @return Character string indicating which backend was chosen
#' @importFrom parallel detectCores
#' @importFrom future plan multicore multisession
#' @importFrom callr r_bg
#' @keywords internal
setup_parallel <- function(n_workers = getOption("trieff.cores")) {
  if (!requireNamespace("future", quietly = TRUE)) {
    stop("Package 'future' is required for parallel processing")
  }

  if (is.null(n_workers)) {
    n_workers <- max(1, parallel::detectCores() - 1)
  }

  # Get the backend type
  backend <- getOption("trieff.backend", get_optimal_backend())

  # Force multisession in RStudio
  if (Sys.getenv("RSTUDIO") == "1") {
    backend <- "multisession"
  }

  # Set up the parallel backend
  if (backend == "multicore") {
    future::plan(future::multicore, workers = n_workers)
  } else {
    # Use callr backend for multisession for better cleanup on force stop
    options(future.fork.enable = FALSE)  # Ensure fork is disabled
    future::plan(future::multisession, workers = n_workers)
  }

  return(backend)
}



#' Clean Up Parallel Processing Environment
#'
#' @description
#' Resets the parallel processing environment to sequential processing.
#' Should be called when parallel processing is complete or on package unload.
#'
#' @importFrom future plan sequential
#' @keywords internal
cleanup_parallel <- function() {
  if (requireNamespace("future", quietly = TRUE)) {
    future::plan(future::sequential)
  }
}

#' Initialize Progress Reporting Handlers
#'
#' @description
#' Sets up progress reporting handlers if progress reporting is enabled in package options.
#'
#' @importFrom progressr handlers
#' @keywords internal
init_progressr <- function() {
  if (getOption("trieff.progress", default = TRUE)) {
    progressr::handlers("progress")
  }
}

#' Package Loading Hook
#'
#' @description
#' Initializes package options and environment when the package is loaded.
#' Sets default values for cores, backend, progress reporting, and logging.
#'
#' @param libname Character string giving the library where the package is installed
#' @param pkgname Character string giving the name of the package
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  # Set default options if they don't exist
  op <- options()
  op.trieff <- list(
    trieff.cores = max(1, parallel::detectCores() - 1),
    trieff.backend = get_optimal_backend(),
    trieff.progress = TRUE,
    trieff.log_level = "INFO"
  )
  toset <- !(names(op.trieff) %in% names(op))
  if(any(toset)) options(op.trieff[toset])

  # Initialize package environment
  .trieff_env$initialized <- TRUE
  .trieff_env$last_calc <- NULL
}

#' Package Attachment Hook
#'
#' @description
#' Displays startup messages when the package is attached.
#' Shows version information and parallel processing configuration.
#'
#' @param libname Character string giving the library where the package is installed
#' @param pkgname Character string giving the name of the package
#' @importFrom utils packageVersion
#' @keywords internal
.onAttach <- function(libname, pkgname) {
  # Version message
  ver <- utils::packageVersion("trieff")
  packageStartupMessage(sprintf("Loaded trieff %s", ver))

  # System capability info
  backend <- get_optimal_backend()
  max_workers <- parallel::detectCores()

}

#' Package Unload Hook
#'
#' @description
#' Performs cleanup operations when the package is unloaded.
#' Ensures proper cleanup of parallel processing and package environment.
#'
#' @param libpath Character string giving the library where the package is installed
#' @keywords internal
.onUnload <- function(libpath) {
  # Clean up any parallel processing
  try(cleanup_parallel(), silent = TRUE)

  # Clear package environment
  rm(list = ls(.trieff_env), envir = .trieff_env)
}





#' Synthetic Emergency Department Data from Skåne Emergency Medicine (SEM) Cohort
#'
#' A synthetic dataset generated based on statistical properties of ED visits across all eight
#' hospitals in the Skåne Emergency Medicine (SEM) database from 2017-2018. Generated using
#' synthpop methodology with complaint-based grouping, this dataset mimics real patterns while
#' containing no actual patient data.
#'
#' @format A data frame with 460,986 observations and 15 variables:
#' \describe{
#'   \item{priority}{Factor with 5 levels (1-5) representing triage priority (RETTS system), where 1 is highest priority}
#'   \item{ambulance}{Logical indicating if patient arrived by ambulance}
#'   \item{unit}{Factor with 8 levels representing the treating hospital: "Helsingborg", "Hässleholm", "Kristianstad", "Landskrona", "Lund", "Malmö", "Trelleborg", "Ystad"}
#'   \item{chief_complaint}{Factor representing the main reason for ED visit (translated to English and grouped into clinical categories)}
#'   \item{age_at_arrival}{Numeric age of patient at time of arrival}
#'   \item{gender}{Factor with 2 levels: "F" (female), "M" (male)}
#'   \item{discharged_to}{Factor with 5 levels: "Admitted", "Died", "Home", "Moved to other hospital", indicating ED visit outcome}
#'   \item{admitted_to_hospital}{Logical indicating if patient was admitted to hospital}
#'   \item{ward}{Factor indicating type of ward for admitted patients (e.g., "ICU", "Cardiology", "MedSurg", etc.)}
#'   \item{loset}{Logical indicating if case was time-critical according to LOSET criteria}
#'   \item{arrival_at_hospital}{POSIXct timestamp of arrival at hospital}
#'   \item{exit_ed}{POSIXct timestamp of exit from ED}
#'   \item{resolve}{POSIXct timestamp of first physician contact}
#'   \item{arrival}{POSIXct timestamp of triage assessment}
#'   \item{chief_complaint_category}{Factor grouping chief complaints into broader clinical categories}
#' }
#'
#' @details
#' This synthetic dataset is based on ED visits from 2017 across two hospitals in
#' the Skåne region of Sweden. It preserves the statistical relationships and patterns of the
#' original SEM cohort while containing no real patient information. The dataset includes adult
#' patients (age >= 18) with medical, surgical, or orthopedic complaints.
#'
#' While the synthetic data maintains similar priority distributions and demographic patterns
#' as the original, triage effectiveness values are consistently lower (e.g., overall Observed
#' TE: 39.1% vs 45.8% in original). However, relative performance patterns across hospitals
#' are preserved, making it suitable for testing the TriEff statistical package.
#'
#' Time variables follow this sequence:
#' arrival_at_hospital <= arrival (triage) <= resolve (physician) <= exit_ed
#'
#' @source Based on statistical properties of the Skåne Emergency Medicine (SEM) database.
#' Reference: Ekelund, U., et al. (2024). The Skåne Emergency Medicine (SEM) cohort.
#' Scandinavian Journal of Trauma Resuscitation and Emergency Medicine, 32(1).
#' doi: https://doi.org/10.1186/s13049-024-01206-0
#'
#' @export
load_sem_synth <- function() {
  path <- system.file("extdata", "sem_synth.rds", package = "trieff")
  readRDS(path)
}

