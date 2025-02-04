#' Optimized Triage Efficacy Simulation
#'
#' This function performs an optimized simulation to calculate various Triage Efficacy (TE) metrics
#' in a deterministic way. It processes the input data frame, creates segments and batches,
#' runs parallel simulations, and applies corrections if specified.
#'
#' @param df A data frame containing patient data, including arrival times, resolve times, and priorities.
#' @param tte Logical. If TRUE, calculates Theoretical Triage Efficacy (TTE). Default is TRUE.
#' @param btte Logical. If TRUE, calculates Binary Theoretical Triage Efficacy (BTTE). Default is FALSE.
#' @param n_workers Integer. Number of cores to use for parallel processing. Default is 1.
#'
#' @return A data frame with the original patient data and additional columns for calculated TE metrics.
#'
#' @details
#' This function is the optimized version of the verbose simulation. It uses segmentation and batching
#' techniques to enable efficient parallel processing. The simulation calculates wait times in a
#' deterministic world where patients are strictly processed according to their priority and wait time
#' within the same priority level - this is the basis for the TTE value.
#'
#' Key steps in the process:
#' 1. Data preparation and extraction of relevant fields.
#' 2. Creation of segments (periods where queue goes from 0 to 0) for independent processing.
#' 3. Batching of segments for balanced parallel processing.
#' 4. Parallel execution of simulations on batches.
#' 5. Aggregation of results
#'
#' @importFrom dplyr select mutate arrange left_join bind_rows filter
#' @importFrom purrr map_dfr map
#' @importFrom parallel detectCores
#' @importFrom furrr future_map furrr_options
#' @importFrom future plan multisession sequential
#'
#' @export
sim_te_alt <- function(df, tte = TRUE, btte = FALSE, n_workers = 1) {
  # Extracting only relevant fields for simulation for RAM conservation
  df_sim <- df %>%
    select(id, unit, arrival_minute, resolve_minute,
           priority, priority_binary, loset)
  df <- df %>%
    select(-c(unit, arrival_minute, resolve_minute,
              priority, priority_binary, loset))

  print(paste("Creating segments out of dataframe", Sys.time()))
  segments <- create_segments(df_sim, n_workers)

  print(paste("Creating batches", Sys.time()))
  n_batches <- n_workers * 3  # Creating 3 times as many batches as cores
  batched_segments <- create_batches(segments, n_batches)

  print(paste("Starting multisession with", n_workers, "cores.", Sys.time()))
  setup_parallel(n_workers = n_workers)
  on.exit(cleanup_parallel())

  # Initialize progress reporting
  init_progressr()

  results <- with_progress({
    p <- progressor(steps = length(batched_segments))
    future_map(batched_segments, function(batch) {
      result <- map_dfr(batch, ~run_simulation(.x, tte, btte))
      p()
      result
    }, .options = furrr::furrr_options(seed = TRUE))
  })
  print(paste("Simulations are done!", Sys.time()))
  plan(sequential)

  # Combine results from all batches
  merged_results <- bind_rows(results)

  # Merge results back with original data
  final_results <- df %>%
    left_join(merged_results, by = "id")

  # Sort the results by id
  final_results <- final_results %>%
    arrange(id)

  return(final_results)
}





#' Fast Triage Efficacy Simulation
#'
#' This function performs an optimized simulation to calculate Triage Efficacy (TE) metrics
#' by only simulating segments containing time-critical (LOSET-positive) cases. For segments
#' without LOSET cases, it uses observed wait times directly, since these segments don't
#' affect TE calculations. This optimization significantly reduces computation time while
#' maintaining identical TE results compared to the full simulation.
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
#' This function optimizes the standard sim_te simulation by:
#' 1. Only simulating segments containing LOSET-positive cases
#' 2. Using observed wait times for segments without LOSET cases
#' 3. Maintaining minimal data during simulation to reduce memory usage
#'
#' This optimization is valid because TE calculation only depends on the waiting times
#' of time-critical patients relative to the overall mean. Since segments without
#' LOSET cases contribute identically to the mean whether simulated or observed,
#' simulating these segments is unnecessary.
#'
#' Key steps in the process:
#' 1. Extract minimal required data for simulation
#' 2. Create and identify segments containing LOSET cases
#' 3. Simulate only LOSET-containing segments
#' 4. Copy observed wait times for non-simulated segments
#' 5. Combine results maintaining all original data
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
#' @importFrom dplyr select mutate arrange left_join bind_rows filter
#' @importFrom purrr map_dfr map
#' @importFrom parallel detectCores
#' @importFrom furrr future_map furrr_options
#' @importFrom future plan multisession sequential
#'
#' @export
sim_te <- function(df, tte = TRUE, btte = FALSE, n_workers = 1) {
  # Extract minimal data for simulation while preserving original
  df_sim <- df %>%
    select(id, unit, arrival_minute, resolve_minute,
           priority, priority_binary, loset)

  # Create segments with minimal data
  print(paste("Creating segments out of dataframe", Sys.time()))
  segments <- create_segments(df_sim, n_workers)

  # Log initial segment counts
  total_patients_before <- sum(sapply(segments, nrow))
  print(paste("Initial segments:", length(segments), "Total patients:", total_patients_before))

  # Filter to LOSET segments (working with minimal data)
  segments <- segments[sapply(segments, function(segment) any(segment$loset))]
  total_loset_cases <- sum(sapply(segments, function(segment) sum(segment$loset)))
  print(paste("Filtered to", length(segments), "segments with LOSET cases"))
  print(paste("Total LOSET cases:", total_loset_cases))

  # Process LOSET segments
  print(paste("Creating batches", Sys.time()))
  n_batches <- n_workers * 3
  batched_segments <- create_batches(segments, n_batches)

  print(paste("Starting multisession with", n_workers, "cores.", Sys.time()))
  setup_parallel(n_workers = n_workers)
  on.exit(cleanup_parallel())

  # Initialize progress reporting
  init_progressr()

  results <- future_map(batched_segments, function(batch) {
    map_dfr(batch, ~run_simulation(.x, tte, btte))
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

