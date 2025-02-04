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
  "runif", "starts_with",

  # Time-related variables
  "arrival", "resolve", "arrival_minute", "resolve_minute", "event_minute",

  # Priority and classification
  "priority", "priority_binary", "loset", "tp", "fp", "fn", "tn",

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

  # Event processing
  "event_type", "cum_arrivals", "cum_resolves", "queue_size", "segment_start", "segment",

  # Plotting related
  "x", "y", "y_pos", "te", "label", "color", "group_label", "random_value",

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
  packageStartupMessage(sprintf("System capable of parallel processing with %s backend", backend))
  packageStartupMessage(sprintf("Maximum available workers: %d", max_workers))
  packageStartupMessage("Note: Individual functions will determine optimal core usage at runtime")
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
