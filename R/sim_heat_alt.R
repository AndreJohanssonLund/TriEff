#' Simulate Triage Effectiveness including alternative calculations for methodology research
#'
#' @description
#' Performs simulations to analyze how Triage Effectiveness (TE) varies with different
#' combinations of sensitivity and specificity. This function supports both waiting time-based TE (WTE)
#' and rank-based TE (RTE) calculations, with options for alternative transformations to explore
#' methodological aspects of triage performance measurement.
#'
#' @param df Data frame containing patient data (must be initialized through init())
#' @param step_size Numeric. The step size for sensitivity and specificity values.
#'   Must be one of: 25, 10, 5, 2.5, or 1.
#' @param n_workers Integer. Number of cores to use for parallel processing.
#'   Default is (detectCores() - 1).
#' @param include_dataframes Logical. If TRUE, includes complete simulation data frames
#'   in output for detailed analysis. Default is FALSE to minimize memory usage.
#' @param alt_calc Logical. If TRUE, calculates alternative TE metrics using transformations
#'   appropriate to the calculation method. Default is FALSE.
#' @param pos_values_only Logical. If TRUE, only includes combinations where
#'   sensitivity + specificity >= 1. Default is FALSE.
#' @param n_loset Integer. Minimum number of LOSET cases desired. Will duplicate
#'   the data to reach this threshold. Default is NULL (no duplication).
#' @param seed Integer. Seed for reproducible heatmaps. Default is NULL (no seed).
#' @param calc_method String. Method to use for calculating TE: "wte" for waiting time-based TE,
#'   or "rte" for rank-based TE (default).
#'
#' @return A tibble containing:
#'   \itemize{
#'     \item sensitivity: Sensitivity value for the simulation
#'     \item specificity: Specificity value for the simulation
#'     \item te: Standard TE value (WTE or RTE based on calc_method)
#'     \item When alt_calc=TRUE and calc_method="wte":
#'       \itemize{
#'         \item te_cube: TE with cube root normalization of waiting times
#'         \item te_log: TE with log normalization of waiting times
#'         \item te_median_local: TE using local median-based calculation
#'         \item te_median_global: TE using global median-based calculation
#'         \item Additional columns with raw metrics (mean_all, cube_all, etc.)
#'       }
#'     \item When alt_calc=TRUE and calc_method="rte":
#'       \itemize{
#'         \item te_cube: RTE with cube root transformation
#'         \item te_log: RTE with signed logarithmic transformation
#'         \item Additional raw metrics as appropriate
#'       }
#'     \item sim_data: (If include_dataframes=TRUE) Complete simulation data frames
#'   }
#'
#' @details
#' The function performs these key steps:
#' 1. Creates combinations of sensitivity/specificity values based on step_size
#' 2. Identifies segments containing time-critical (LOSET) cases
#' 3. Processes each combination in parallel:
#'    - Assigns priorities based on sensitivity/specificity
#'    - Simulates queue processing
#'    - Calculates TE metrics (WTE or RTE based on calc_method)
#'    - Applies alternative calculations if requested
#' 4. Returns combined results in a tibble format
#'
#' WTE alternative calculations (when alt_calc=TRUE and calc_method="wte"):
#' - Cube root: Normalizes wait times using cube root transformation
#' - Log: Normalizes wait times using log transformation
#' - Median: Uses median instead of mean for central tendency
#'
#' RTE alternative calculations (when alt_calc=TRUE and calc_method="rte"):
#' - Cube root: Applies cube root transformation to RTE values (preserves sign)
#' - Signed log: Applies sign-preserving log transformation: sign(x) * log(1 + abs(x))
#'
#' Performance considerations:
#' - Large step_sizes (e.g., 25) are suitable for initial exploration
#' - Smaller step_sizes provide more detailed analysis but increase computation time
#' - Setting include_dataframes=TRUE significantly increases memory usage
#' - Processing time increases with dataset size and number of combinations
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
                         pos_values_only = FALSE, n_loset = NULL, seed = NULL,
                         calc_method = "rte") {
  # Input validation
  valid_step_sizes <- c(25, 10, 5, 2.5, 1)
  if (!step_size %in% valid_step_sizes) {
    stop("Invalid step_size. Must be one of: 25, 10, 5, 2.5, or 1.")
  }

  # Validate calculation method
  valid_methods <- c("wte", "rte")
  if (!calc_method %in% valid_methods) {
    stop("Invalid calc_method. Must be one of: 'wte' or 'rte'.")
  }

  validate_sim_heat_data(df)

  df <- df %>%
    select("id", "arrival_minute", "resolve_minute", "loset", "observed_wait_time", "unit", "segment")

  # Create sensitivity and specificity combinations
  steps <- seq(0, 1, by = step_size / 100)
  combinations <- expand_grid(sensitivity = steps, specificity = steps)

  if (pos_values_only) {
    combinations <- combinations %>%
      filter(sensitivity + specificity >= 1)
  }

  # Duplicate segments if n_loset is specified
  if (!is.null(n_loset)) {
    # Count current LOSET cases
    total_loset <- sum(df$loset)

    # Calculate needed duplications (rounding up)
    n_duplications <- ceiling(n_loset / total_loset)

    if (n_duplications > 1) {
      print(paste("Duplicating segments", n_duplications, "times to reach minimum", n_loset, "LOSET cases"))

      # Create the original plus n_duplications-1 copies
      df_copies <- c(
        list(df),  # Original data
        lapply(2:n_duplications, function(i) {
          # Create copy with modified segment IDs
          copy <- df
          copy$segment <- paste0(copy$segment, "_dup", i)
          return(copy)
        })
      )

      df <- bind_rows(df_copies)
      new_total_loset <- sum(df$loset)
      new_total_rows <- nrow(df)
      print(paste("New total LOSET cases:", new_total_loset))
      print(paste("Total patient count:", new_total_rows))
    }
  }

  # If alt_calc is TRUE for WTE, handle global median calculation
  median_all_global <- if(alt_calc && calc_method == "wte") {
    print(paste("Calculating global median from reference simulation", Sys.time()))

    # Create copy of segmented data
    df_median <- df

    # Set all priorities to 1 (equivalent to 100% sens/0% spec)
    df_median$priority_binary <- 1

    # Process by segment groups
    wait_times <- df_median %>%
      group_by(segment) %>%
      group_split() %>%
      map(~ run_simulation(.x, tte = FALSE, btte = TRUE)) %>%
      bind_rows() %>%
      pull(binary_theoretical_wait_time)

    median(wait_times)
  } else {
    NULL
  }

  combinations <- split(combinations, seq(nrow(combinations)))

  # Set up parallel processing
  setup_parallel(n_workers = n_workers)
  on.exit(cleanup_parallel())

  # Initialize progress reporting
  init_progressr()

  print(paste("Starting simulations with", length(combinations), "combinations", Sys.time()))

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

        # Create a copy of the segmented data for this combination
        df_copy <- df
        random_values <- runif(nrow(df_copy))

        # Assign priorities based on sensitivity/specificity
        df_copy$priority_binary <- data.table::fifelse(
          df_copy$loset == TRUE,
          data.table::fifelse(random_values <= sens, 1, 2),
          data.table::fifelse(random_values <= spec, 2, 1)
        )

        # Process by segment groups
        sim_results <- df_copy %>%
          group_by(segment) %>%
          group_split() %>%
          map(~ run_simulation(.x, tte = FALSE, btte = TRUE)) %>%
          bind_rows()

        # Different calculation methods for TE
        if (calc_method == "wte") {
          # Extract wait times and LOSET status
          all_wait_times <- sim_results$binary_theoretical_wait_time
          loset_indices <- sim_results$loset

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
        } else if (calc_method == "rte") {
          # Apply calculate_queue_metrics to get RTE values
          rte_results <- calculate_queue_metrics(sim_results, verbose = include_dataframes)

          # Extract valid RTE values for LOSET patients
          valid_rte <- rte_results %>%
            filter(loset == TRUE, valid == TRUE) %>%
            pull(binary_theoretical_RTE)

          # Calculate mean RTE if there are valid values, otherwise NA
          te <- if (length(valid_rte) > 0) mean(valid_rte, na.rm = TRUE) else NA_real_

          # Store base results
          result <- list(
            sensitivity = sens,
            specificity = spec,
            te = te
          )

          # Calculate alternative RTE metrics if requested
          if (alt_calc && length(valid_rte) > 0) {
            # Cube root transformation (preserves sign)
            cube_transform <- function(x) sign(x) * abs(x)^(1/3)
            rte_cube <- mean(cube_transform(valid_rte), na.rm = TRUE)
            result$te_cube <- rte_cube

            # Signed log transformation
            signed_log_transform <- function(x) sign(x) * log(1 + abs(x))
            rte_log <- mean(signed_log_transform(valid_rte), na.rm = TRUE)
            result$te_log <- rte_log

            # Median-based RTE (without transformation)
            rte_median <- median(valid_rte, na.rm = TRUE)
            result$te_median <- rte_median

            # Store original RTE values for reference
            result$rte_values <- list(valid_rte)
          }
        }

        # Store complete data if requested
        if (include_dataframes) {
          if (calc_method == "rte") {
            result$sim_data <- rte_results
          } else {
            result$sim_data <- sim_results
          }
        }

        result
      })

      p()
      batch_results
    }, .options = furrr_options(seed = TRUE))
  })

  # Convert results to tibble
  results_df <- tibble(
    sensitivity = map_dbl(unlist(results, recursive = FALSE), ~.x$sensitivity),
    specificity = map_dbl(unlist(results, recursive = FALSE), ~.x$specificity),
    te = map_dbl(unlist(results, recursive = FALSE), ~.x$te)
  )

  # Add calculation method to the results
  results_df$calc_method <- calc_method

  # Add method-specific fields
  if (calc_method == "wte") {
    results_df$mean_all <- map_dbl(unlist(results, recursive = FALSE), ~.x$mean_all)

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
  } else if (calc_method == "rte") {
    # Handle RTE base case

    # Add alternative RTE calculations if requested
    if (alt_calc) {
      # Use safer extraction approach
      te_cube_values <- sapply(unlist(results, recursive = FALSE), function(x) {
        if(is.null(x$te_cube)) NA_real_ else x$te_cube
      })

      te_log_values <- sapply(unlist(results, recursive = FALSE), function(x) {
        if(is.null(x$te_log)) NA_real_ else x$te_log
      })

      te_median_values <- sapply(unlist(results, recursive = FALSE), function(x) {
        if(is.null(x$te_median)) NA_real_ else x$te_median
      })

      results_df$te_cube <- te_cube_values
      results_df$te_log <- te_log_values
      results_df$te_median <- te_median_values

      # Optionally store the original RTE values if include_dataframes is TRUE
      if (include_dataframes) {
        results_df$rte_values <- lapply(unlist(results, recursive = FALSE), function(x) {
          if(is.null(x$rte_values)) list(NULL) else x$rte_values
        })
      }
    }
  }

  # Add simulation data if requested
  if (include_dataframes) {
    results_df$sim_data <- lapply(unlist(results, recursive = FALSE), function(x) {
      if(is.null(x$sim_data)) NULL else x$sim_data
    })
  }

  print(paste("Simulations complete", Sys.time()))
  print(paste("Calculation method used:", calc_method))

  return(results_df)
}
