#' Simulate Triage Effectiveness over different levels of sensitivity/specificity.
#'
#' @description
#' Efficiently generates a heatmap analyzing the relationship between sensitivity, specificity,
#' and Triage Effectiveness (TE). This function optimizes performance by:
#' 1. Processing only segments containing time-critical (LOSET) cases and
#' 2. Parallelizing across sensitivity/specificity combinations
#'
#' @details
#' The function processes data through these key steps:
#' 1. Pre-calculation:
#'    - Computes overall mean waiting time for TE calculations
#'    - Identifies segments containing LOSET cases
#'    - Creates sensitivity/specificity combinations based on step_size
#'
#' 2. Parallel Processing:
#'    - Processes each combination independently across workers
#'    - For each combination:
#'      * Scholastically applies priority assignment based on sensitivity/specificity
#'      * Simulates queue processing for each step
#'      * Calculates resulting (binary) TE value
#'
#' 3. Progress Tracking:
#'    - Updates progress after each step
#'
#' Performance considerations:
#' - Large step_sizes (e.g., 25) are suitable for initial exploration
#' - Smaller step_sizes provide more detailed analysis but increase computation time
#' - Processing time increases with dataset size and number of combinations
#'
#' @param df A data frame containing patient data (must be initialized through init())
#' @param step_size Numeric. The step size for sensitivity and specificity values.
#'   Must be one of: 25, 10, 5, 2.5, or 1. representing the percentual step in
#'   sensitivity/specificity that is taken.
#' @param n_workers Integer. Number of cores to use for parallel processing.
#'   Default is (detectCores() - 1).
#' @param pos_values_only Logical. If TRUE, only includes combinations where
#'   sensitivity + specificity >= 1. Default is FALSE.
#' @param n_loset Integer, least number of LOSET cases in the data frame. Will duplicate
#' the data to reach this. Can by doing this get more LOSET cases.
#' @param seed Seed for reproducible heat maps (default: null - no seed)
#' @param calc_method String. Method to use for calculating TE: "wte" for waiting time-based TE,
#'   or "rte" for rank-based TE (default).
#' @param include_dataframes Logical. If TRUE, includes complete simulation data frames
#'   in output for detailed analysis. Default is FALSE to minimize memory usage.
#'
#' @return A Tibble with sensitivity, specificity and TE values
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
#' @examples
#' \dontrun{
#' # Basic usage
#' results <- sim_heat(patient_data, step_size = 10)
#'
#' # Only positive values
#' results <- sim_heat(patient_data,
#'                     step_size = 10,
#'                     pos_values_only = TRUE)
#'
#' # Adjust worker count
#' results <- sim_heat(patient_data,
#'                     step_size = 10,
#'                     n_workers = 4)
#'
#' # Include simulation data frames
#' results <- sim_heat(patient_data,
#'                     step_size = 10,
#'                     include_dataframes = TRUE)
#' }
#' @export
sim_heat <- function(df, step_size, n_workers = detectCores() - 1,
                     pos_values_only = FALSE, n_loset = NULL, seed = NULL,
                     calc_method = "rte", include_dataframes = FALSE) {
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

  # Calculate mean_all once at the start (only needed for WTE)
  mean_all <- mean(df$observed_wait_time)

  # Create sensitivity and specificity combinations
  steps <- seq(0, 1, by = step_size / 100)
  combinations <- expand_grid(sensitivity = steps, specificity = steps)

  if (pos_values_only) {
    combinations <- combinations %>%
      filter(sensitivity + specificity >= 1)
  }

  # Filter to segments containing LOSET cases
  segments_with_loset <- df %>%
    group_by(segment) %>%
    summarize(has_loset = any(loset), .groups = "drop") %>%
    filter(has_loset) %>%
    pull(segment)

  df <- df %>%
    filter(segment %in% segments_with_loset)

  print(paste("Processing", length(unique(df$segment)), "segments with LOSET cases"))

  # If n_loset specified, handle replication differently with the segment column approach
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
      print(paste("Total patient count in original df:", nrow(df)))
      print(paste("Total patient count based on duplication", nrow(df) * n_duplications))
    }
  }

  combinations <- split(combinations, seq(nrow(combinations)))

  # Set up parallel processing
  setup_parallel(n_workers = n_workers)
  on.exit(cleanup_parallel())

  # Initialize progress reporting
  init_progressr()

  print(paste("Starting simulations with", length(combinations), "combinations", Sys.time()))

  # Run parallel simulations by sensitivity/specificity combination
  results <- with_progress({
    p <- progressor(steps = length(combinations))

    # If seed provided, set it
    if (!is.null(seed)) {
      set.seed(seed)
    }

    # Then in the batch processing loop, modify the segment processing:
    future_map(combinations, function(combinations) {
      # Process each combination
      batch_results <- map(1:nrow(combinations), function(i) {
        sens <- combinations$sensitivity[i]
        spec <- combinations$specificity[i]

        # Instead of processing each segment separately, process all segments together
        # with proper grouping
        df_copy <- df
        random_values <- runif(nrow(df_copy))

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
          # Original waiting time-based TE calculation
          # Extract wait times and LOSET status
          all_wait_times <- sim_results$binary_theoretical_wait_time
          loset_indices <- sim_results$loset

          # Calculate TE for this combination
          mean_loset <- mean(all_wait_times[loset_indices])
          te <- 1 - (mean_loset / mean_all)
        } else if (calc_method == "rte") {
          # Rank-based TE calculation
          # Apply calculate_queue_metrics to get RTE values
          rte_results <- calculate_queue_metrics(sim_results,
                                                 verbose = include_dataframes)

          # Extract valid RTE values for LOSET patients
          valid_rte <- rte_results %>%
            filter(loset == TRUE, valid == TRUE) %>%
            pull(binary_theoretical_RTE)

          # Calculate mean RTE if there are valid values, otherwise 0
          te <- if (length(valid_rte) > 0) mean(valid_rte, na.rm = TRUE) else NA_integer_
        }

        # Create result structure
        result <- list(
          sensitivity = sens,
          specificity = spec,
          te = te
        )

        # Store complete data if requested
        if (include_dataframes) {
          result$sim_data <- sim_results
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
    te = map_dbl(unlist(results, recursive = FALSE), ~.x$te)
  )

  # Add simulation data if requested
  if (include_dataframes) {
    results_df$sim_data <- map(unlist(results, recursive = FALSE), ~.x$sim_data)
  }

  print(paste("Simulations complete", Sys.time()))
  print(paste("Calculation method used:", calc_method))

  return(results_df)
}
