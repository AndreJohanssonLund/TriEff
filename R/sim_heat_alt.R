#' Simulate Triage Effectiveness including potential normalization techniques - mainly for validation and methodological research
#'
#' @description
#' Performs simulations to analyze how Triage Effectiveness (TE) varies with different
#' combinations of sensitivity and specificity. The function provides options for
#' alternative TE calculations and can store complete simulation data for detailed analysis.
#'
#' @param df Data frame containing patient data (must be initialized through init())
#' @param step_size Numeric. The step size for sensitivity and specificity values.
#'   Must be one of: 25, 10, 5, 2.5, or 1.
#' @param n_workers Integer. Number of cores to use for parallel processing.
#'   Default is (detectCores() - 1).
#' @param include_dataframes Logical. If TRUE, includes complete simulation data frames
#'   in output for detailed analysis. Default is FALSE to minimize memory usage.
#' @param alt_calc Logical. If TRUE, calculates alternative TE metrics using cube root,
#'   log, and median-based normalization. Default is FALSE.
#' @param pos_values_only Logical. If TRUE, only includes combinations where
#'   sensitivity + specificity >= 1. Default is FALSE.
#' @param n_loset Integer, least number of LOSET cases in the datafram. Will duplicate
#' the data to reach this.
#' @param seed Seed for reprorucible heatmaps (default: null - no seed)
#'
#'
#' @return A tibble containing:
#'   \itemize{
#'     \item sensitivity: Sensitivity value for the simulation
#'     \item specificity: Specificity value for the simulation
#'     \item te: Standard Triage Effectiveness value
#'     \item te_cube: (If alt_calc=TRUE) TE with cube root normalization
#'     \item te_log: (If alt_calc=TRUE) TE with log normalization
#'     \item te_median_local: (If alt_calc=TRUE) TE using median-based calculation using the simulations median
#'     \item te_median_global: (If alt_calc=TRUE) TE using median-based calculations using the median_all from 100% sensitivity/0%specificity
#'     \item sim_data: (If include_dataframes=TRUE) Complete simulation data frames
#'   }
#'
#' @details
#' The function performs these key steps:
#' 0. Optional, (if alt_calc=TRUE) runs a simulation to evaluate median TE at
#'  sensitivity = 100%, specificity = 0% - i.e a first come, first serve ED. This
#'  acts like the denominator for the global median calculation.
#' 1. Creates combinations of sensitivity/specificity values based on step_size
#' 2. Identifies segments containing time-critical (LOSET) cases
#' 3. Processes each combination in parallel:
#'    - Assigns priorities based on sensitivity/specificity
#'    - Simulates queue processing
#'    - Calculates TE metrics
#' 4. Returns combined results in a tibble format
#'
#' Alternative TE calculations (when alt_calc=TRUE):
#' - Cube root: Normalizes wait times using cube root transformation
#' - Log: Normalizes wait times using log transformation
#' - Median: Uses median instead of mean for central tendency
#'
#' Performance considerations:
#' - Large step_sizes (e.g., 25) are suitable for initial exploration
#' - Smaller step_sizes provide more detailed analysis but increase computation time
#' - Setting include_dataframes=TRUE significantly increases memory usage
#' - Processing time increases with dataset size and number of combinations
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' results <- sim_heat(patient_data, step_size = 10)
#'
#' # With alternative calculations
#' results <- sim_heat(patient_data,
#'                     step_size = 10,
#'                     alt_calc = TRUE)
#'
#' # Include complete simulation data
#' results <- sim_heat(patient_data,
#'                     step_size = 10,
#'                     include_dataframes = TRUE)
#'
#' # Only positive TE combinations
#' results <- sim_heat(patient_data,
#'                     step_size = 10,
#'                     pos_values_only = TRUE)
#' }
#'
#' @importFrom parallel detectCores
#' @importFrom dplyr filter select mutate arrange bind_rows
#' @importFrom tidyr expand_grid
#' @importFrom purrr map map_dbl
#' @importFrom data.table as.data.table fifelse
#' @importFrom furrr future_map furrr_options
#' @importFrom future plan sequential multisession
#' @importFrom progressr progressor handlers
#'
#' @export
sim_heat_alt <- function(df, step_size, n_workers = detectCores() - 1,
                     include_dataframes = FALSE, alt_calc = FALSE,
                     pos_values_only = FALSE, n_loset = NULL, seed = NULL) {
  # Input validation
  valid_step_sizes <- c(25, 10, 5, 2.5, 1)
  if (!step_size %in% valid_step_sizes) {
    stop("Invalid step_size. Must be one of: 25, 10, 5, 2.5, or 1.")
  }

  validate_sim_heat_data(df)

  # Create sensitivity and specificity combinations
  steps <- seq(0, 1, by = step_size / 100)
  combinations <- expand_grid(sensitivity = steps, specificity = steps)

  if (pos_values_only) {
    combinations <- combinations %>%
      filter(sensitivity + specificity >= 1)
  }

  # Create segments
  print(paste("Creating segments out of dataframe", Sys.time()))
  segments <- create_segments(df, n_workers)
  print(paste("Processing", length(segments), "segments."))

  # Duplicate segments if n_loset is specified
  if (!is.null(n_loset)) {
    # Count current LOSET cases
    total_loset <- sum(sapply(segments, function(segment) sum(segment$loset)))

    # Calculate needed duplications (rounding up)
    n_duplications <- ceiling(n_loset / total_loset)

    if (n_duplications > 1) {
      print(paste("Duplicating segments", n_duplications, "times to reach minimum", n_loset, "LOSET cases"))
      segments <- rep(segments, n_duplications)
      new_total_loset <- sum(sapply(segments, function(segment) sum(segment$loset)))
      print(paste("New total LOSET cases:", new_total_loset))
    }
  }


  # If alt_calc is TRUE, calculate global median from reference simulation
  median_all_global <- if(alt_calc) {
    print(paste("Calculating global median from reference simulation", Sys.time()))
    all_wait_times <- vector("numeric")

    # Process each segment with all priorities set to 1 (equivalent to 100% sens/0% spec)
    for(segment in segments) {
      dt <- data.table::as.data.table(segment)
      # Set all priorities to 1 since this represents 100% sensitivity/0% specificity
      dt[, priority_binary := 1]

      sim_result <- run_simulation(dt, tte = FALSE, btte = TRUE)
      all_wait_times <- c(all_wait_times, sim_result$binary_theoretical_wait_time)
      rm(sim_result)
    }

    median(all_wait_times)
  } else {
    NULL
  }

  combinations <- split(combinations, seq(nrow(combinations)))

  # Set up parallel processing
  setup_parallel(n_workers = n_workers)
  on.exit(cleanup_parallel())

  # Initialize progress reporting
  init_progressr()

  print(paste("Starting simulations with", nrow(combinations), "combinations", Sys.time()))

  # Run parallel simulations by sesitivity/specificity combination
  results <- with_progress({
    p <- progressor(steps = length(combinations))

    # If seed provided, set it
    if (!is.null(seed)) {
      set.seed(seed)
    }

    future_map(combinations, function(combinations) {
      # Process each combination
      batch_results <- map(1:nrow(combinations), function(i) {
        sens <- combinations$sensitivity[i]
        spec <- combinations$specificity[i]

        # Storage for simulation results
        all_wait_times <- vector("numeric")
        loset_indices <- vector("logical")

        # Process each segment
        for(segment in segments) {
          dt <- data.table::as.data.table(segment)
          random_values <- runif(nrow(dt))
          dt[, priority_binary := data.table::fifelse(
            loset == TRUE,
            data.table::fifelse(random_values <= sens, 1, 2),
            data.table::fifelse(random_values <= spec, 2, 1)
          )]

          sim_result <- run_simulation(dt, tte = FALSE, btte = TRUE)
          all_wait_times <- c(all_wait_times, sim_result$binary_theoretical_wait_time)
          loset_indices <- c(loset_indices, sim_result$loset)

          # Store complete simulation results if requested
          if (include_dataframes) {
            dt$sim_wait_times <- sim_result$binary_theoretical_wait_time
          }

          rm(sim_result)
        }

        # Calculate mean values for this simulation
        mean_all <- mean(all_wait_times)
        mean_loset <- mean(all_wait_times[loset_indices])

        # Calculate TE values and store reference values
        result <- list(
          sensitivity = sens,
          specificity = spec,
          te = 1 - (mean_loset / mean_all),
          mean_all = mean_all
        )

        # Calculate alternative TE metrics if requested
        if (alt_calc) {
          # Cube root normalization
          cube_all <- mean(all_wait_times^(1/3))
          cube_loset <- mean((all_wait_times[loset_indices])^(1/3))
          result$te_cube <- 1 - (cube_loset / cube_all)
          result$cube_all <- cube_all

          # Log normalization
          log_all <- mean(log(all_wait_times + 1))
          log_loset <- mean(log(all_wait_times[loset_indices] + 1))
          result$te_log <- 1 - (log_loset / log_all)
          result$log_all <- log_all

          # Median based (local)
          median_all_local <- median(all_wait_times)
          median_loset <- median(all_wait_times[loset_indices])
          result$te_median_local <- 1 - (median_loset / median_all_local)
          result$median_all_local <- median_all_local

          # Median based (global reference)
          result$te_median_global <- 1 - (median_loset / median_all_global)
          result$median_all_global <- median_all_global
        }

        # Store complete data if requested
        if (include_dataframes) {
          result$sim_data <- dt
        }

        result
      })

      p()
      batch_results
    }, .options = furrr::furrr_options(seed = TRUE))
  })

  # Convert results to tibble
  results_df <- tibble(
    sensitivity = map_dbl(unlist(results, recursive = FALSE), ~.x$sensitivity),
    specificity = map_dbl(unlist(results, recursive = FALSE), ~.x$specificity),
    te = map_dbl(unlist(results, recursive = FALSE), ~.x$te),
    mean_all = map_dbl(unlist(results, recursive = FALSE), ~.x$mean_all)
  )

  # Add alternative metrics if calculated
  if (alt_calc) {
    results_df <- results_df %>%
      mutate(
        te_cube = map_dbl(unlist(results, recursive = FALSE), ~.x$te_cube),
        cube_all = map_dbl(unlist(results, recursive = FALSE), ~.x$cube_all),
        te_log = map_dbl(unlist(results, recursive = FALSE), ~.x$te_log),
        log_all = map_dbl(unlist(results, recursive = FALSE), ~.x$log_all),
        te_median_local = map_dbl(unlist(results, recursive = FALSE), ~.x$te_median_local),
        median_all_local = map_dbl(unlist(results, recursive = FALSE), ~.x$median_all_local),
        te_median_global = map_dbl(unlist(results, recursive = FALSE), ~.x$te_median_global),
        median_all_global = map_dbl(unlist(results, recursive = FALSE), ~.x$median_all_global)
      )
  }

  # Add simulation data if requested
  if (include_dataframes) {
    results_df$sim_data <- map(unlist(results, recursive = FALSE), ~.x$sim_data)
  }

  print(paste("Simulations complete", Sys.time()))

  return(results_df)
}
