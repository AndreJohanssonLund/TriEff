#' Full Triage Effectiveness Simulation for Research
#'
#' @description
#' Performs a complete simulation of emergency department patient flow, calculating
#' theoretical wait times for all patients regardless of their time-critical status.
#' This function is designed for research purposes when complete queue behavior analysis
#' is needed, and is more computationally intensive than the standard sim_te function.
#'
#' @param df A data frame containing patient data, including:
#'   \itemize{
#'     \item id: Unique patient identifier
#'     \item arrival_minute: Minutes since start of study period until patient arrival
#'     \item resolve_minute: Minutes since start until patient sees physician
#'     \item priority: Original triage priority
#'     \item priority_binary: Binary priority (typically 1 for high, 2 for low)
#'     \item loset: Logical indicating if patient is time-critical per LOSET criteria
#'     \item observed_wait_time: Actual waiting time in minutes
#'   }
#' @param tte Logical. If TRUE, calculates Theoretical Triage Efficacy (TTE). Default is TRUE.
#' @param btte Logical. If TRUE, calculates Binary Theoretical Triage Efficacy (BTTE). Default is FALSE.
#' @param n_workers Integer. Number of cores to use for parallel processing. Default is 1.
#'
#' @return A data frame with the original patient data and additional columns:
#'   \itemize{
#'     \item theoretical_wait_time: Simulated wait time based on original priorities (if tte=TRUE)
#'     \item binary_theoretical_wait_time: Simulated wait time based on binary priorities (if btte=TRUE)
#'   }
#'
#' @details
#' Unlike sim_te which optimizes performance by focusing only on segments containing
#' time-critical cases, sim_te_alt simulates waiting times for all patients. This makes
#' it more suitable for research purposes such as:
#'
#' - Analyzing complete queue behavior patterns
#' - Studying wait time distributions for non-time-critical patients
#' - Validating simulation methodology
#' - Investigating effects of different priority schemes
#'
#' The simulation operates in a non-preemptive manner, treating physician contacts as
#' available time slots and processing patients strictly according to priority and
#' arrival order. While this provides complete information about theoretical queue
#' behavior, it is computationally more intensive than sim_te and may not be necessary
#' for standard TE calculations.
#'
#' Note that while sim_te_alt calculates wait times for all patients, it produces
#' identical TE values to sim_te since TE metrics depend only on the wait times of
#' time-critical patients relative to the overall mean.
#'
#' @examples
#' \dontrun{
#' # Basic usage with default settings
#' results <- sim_te_alt(patient_data)
#'
#' # Calculate both TTE and BTTE with parallel processing
#' results <- sim_te_alt(patient_data,
#'                       tte = TRUE,
#'                       btte = TRUE,
#'                       n_workers = 4)
#' }
#'
#' @seealso
#' \code{\link{sim_te}} for the optimized simulation approach
#' recommended for standard TE calculations
#'
#' @importFrom dplyr select mutate arrange left_join bind_rows filter summarize
#' @importFrom purrr map_dfr map
#' @importFrom parallel detectCores
#' @importFrom furrr future_map furrr_options
#' @importFrom future plan multisession sequential
#'
#' @export
sim_te_alt <- function(df, tte = TRUE, btte = FALSE, n_workers = 1) {
  # Check data columns and verify this data hasn't been simulated before
  if ((("theoretical_wait_time" %in% names(df)) && any(!is.na(df$theoretical_wait_time))) ||
      (("binary_theoretical_wait_time" %in% names(df)) && any(!is.na(df$binary_theoretical_wait_time)))) {
    stop("Data appears to have already been simulated. Each dataset should only be simulated once. ",
         "If you need to re-simulate, please start with the original initialized data.")
  }
  # Extracting only relevant fields for simulation for RAM conservation
  df_sim <- df %>%
    select(id, unit, arrival_minute, resolve_minute,
           priority, priority_binary, loset, segment)
  df <- df %>%
    select(-c(unit, arrival_minute, resolve_minute,
              priority, priority_binary, loset))

  print(paste("Starting multisession with", n_workers, "cores.", Sys.time()))
  setup_parallel(n_workers = n_workers)
  on.exit(cleanup_parallel())

  # Initialize progress reporting
  init_progressr()

  # Process segments in parallel (no batches)
  segment_groups <- df_sim %>% group_by(segment) %>% group_split()

  results <- with_progress({
    p <- progressor(steps = length(segment_groups))
    future_map(segment_groups, function(segment_df) {
      result <- run_simulation(segment_df, tte, btte)
      p()
      result
    }, .options = furrr::furrr_options(seed = TRUE))
  })

  print(paste("Simulations are done!", Sys.time()))
  plan(sequential)

  # Combine results from all segments
  merged_results <- bind_rows(results)

  # Merge results back with original data
  final_results <- df %>%
    left_join(merged_results, by = "id")

  # Sort the results by id
  final_results <- final_results %>%
    arrange(id)

  return(final_results)
}


#' Triage Efficacy Simulation
#'
#' This function performs an optimized simulation to calculate Triage Efficacy (TE) metrics
#' by only simulating segments containing time-critical (LOSET-positive) cases. For segments
#' without LOSET cases, it uses observed wait times directly, since these segments don't
#' affect TE calculations. This optimization significantly reduces computation time while
#' maintaining identical TE results compared to the full simulation.
#'
#' @param df A data frame containing patient data, including (note that most of theese variables are created by init()):
#'   \itemize{
#'     \item id: Unique patient identifier
#'     \item arrival_minute: Minutes since start of study period until patient arrival
#'     \item resolve_minute: Minutes since start until patient sees physician
#'     \item priority: Original triage priority
#'     \item priority_binary: Binary priority (typically 1 for high, 2 for low)
#'     \item loset: Logical indicating if patient is time-critical per LOSET criteria
#'     \item observed_wait_time: Actual waiting time in minutes
#'   }
#' @param tte Logical. If TRUE, calculates Theoretical Triage Efficacy (TTE). Default is TRUE.
#' @param btte Logical. If TRUE, calculates Binary Theoretical Triage Efficacy (BTTE). Default is FALSE.
#' @param n_workers Integer. Number of cores to use for parallel processing. Default is 1.
#'
#' @return A data frame with the original patient data and additional columns:
#'   \itemize{
#'     \item theoretical_wait_time: Simulated wait time based on original priorities (if tte=TRUE)
#'     \item binary_theoretical_wait_time: Simulated wait time based on binary priorities (if btte=TRUE)
#'   }
#'
#' @details
#' This function simulates how patients would be processed if strictly following their
#' assigned priority levels. The simulation operates by:
#'
#' 1. Identifying continuous queue segments (periods where queue size goes from 0 to 0)
#' 2. Determining which segments contain time-critical (LOSET) cases
#' 3. Simulating only segments containing LOSET cases using a non-preemptive approach:
#'    - Each physician contact becomes an available time slot
#'    - Patients are queued by priority and arrival time
#'    - Wait times are calculated from slot assignment
#' 4. For segments without LOSET cases, the observed wait times are retained as theoretical
#'    wait times. This is valid because:
#'    - TE calculations depend only on LOSET patient wait times relative to mean_all
#'    - Using observed times for non-LOSET segments maintains the correct mean_all
#'    - The waiting times of non-LOSET patients do not affect TE calculations
#'
#' The simulation ensures that:
#' - Within each simulated segment, patients are seen strictly in order of priority
#' - Within the same priority level, patients are seen in order of arrival
#' - No post-triage priority changes occur
#'
#' @examples
#' \dontrun{
#' # Basic usage with default settings
#' results <- sim_te_fast(patient_data)
#'
#' # Calculate both TTE and BTTE
#' results <- sim_te_fast(patient_data, tte = TRUE, btte = TRUE)

#' }
#'
#' @seealso
#' \code{\link{sim_te}} for the full simulation approach
#'
#' @importFrom dplyr select mutate arrange left_join bind_rows filter summarize
#' @importFrom purrr map_dfr map
#' @importFrom parallel detectCores
#' @importFrom furrr future_map furrr_options
#' @importFrom future plan multisession sequential
#'
#' @export
sim_te <- function(df, tte = TRUE, btte = FALSE, n_workers = 1) {
  # Check data columns and verify this data hasn't been simulated before
  if ((("theoretical_wait_time" %in% names(df)) && any(!is.na(df$theoretical_wait_time))) ||
      (("binary_theoretical_wait_time" %in% names(df)) && any(!is.na(df$binary_theoretical_wait_time)))) {
    stop("Data appears to have already been simulated. Each dataset should only be simulated once. ",
         "If you need to re-simulate, please start with the original initialized data.")
  }
  # Extract minimal data for simulation while preserving original
  df_sim <- df %>%
    dplyr::select(id, unit, arrival_minute, resolve_minute,
           priority, priority_binary, loset, segment)

  # Filter to keep only rows in segments containing LOSET cases
  loset_segments <- df_sim %>%
    dplyr::group_by(segment) %>%
    dplyr::summarize(has_loset = any(loset)) %>%
    dplyr::filter(has_loset) %>%
    dplyr::pull(segment)

  total_patients_before <- nrow(df_sim)
  df_sim <- df_sim %>% filter(segment %in% loset_segments)
  total_loset_cases <- sum(df_sim$loset)

  print(paste("Filtered to", length(loset_segments), "segments with LOSET cases"))
  print(paste("Total LOSET cases:", total_loset_cases))

  # Process segments in parallel
  print(paste("Starting multisession with", n_workers, "cores.", Sys.time()))
  setup_parallel(n_workers = n_workers)
  on.exit(cleanup_parallel())

  # Initialize progress reporting
  init_progressr()

  # Group data by segment and process each group
  segment_groups <- df_sim %>% group_by(segment) %>% group_split()

  results <- future_map(segment_groups, function(segment_df) {
    run_simulation(segment_df, tte, btte)
  }, .options = furrr::furrr_options(seed = TRUE))

  print(paste("Simulations are done!", Sys.time()))
  plan(sequential)

  # Combine simulation results (just id and wait times)
  sim_results <- bind_rows(results) %>%
    select(id, matches("theoretical_wait_time|binary_theoretical_wait_time"))

  # Join back to original data, maintaining all columns
  final_results <- df %>%
    left_join(sim_results, by = "id") %>%
    arrange(id)

  # Fill in non-simulated cases with observed wait times
  if (tte) {
    final_results <- final_results %>%
      mutate(
        theoretical_wait_time = case_when(
          !is.na(theoretical_wait_time) ~ theoretical_wait_time,
          TRUE ~ observed_wait_time
        )
      )
  }
  if (btte) {
    final_results <- final_results %>%
      mutate(
        binary_theoretical_wait_time = case_when(
          !is.na(binary_theoretical_wait_time) ~ binary_theoretical_wait_time,
          TRUE ~ observed_wait_time
        )
      )
  }

  return(final_results)
}

