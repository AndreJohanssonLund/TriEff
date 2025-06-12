#' Extract Triage Priority Distribution
#'
#' @description
#' Extracts the distribution of triage priorities for LOSET-positive and LOSET-negative
#' patients from emergency department data. This function calculates the proportion of
#' each priority level within each patient group (time-critical vs non-time-critical).
#'
#' @param df Data frame containing patient data with 'priority' and 'loset' columns
#'
#' @return A tibble with columns:
#'   \itemize{
#'     \item priority: Priority level
#'     \item loset_count: Count of LOSET-positive patients with this priority
#'     \item non_loset_count: Count of LOSET-negative patients with this priority
#'     \item loset_distribution: Proportion of all LOSET-positive patients with this priority
#'     \item non_loset_distribution: Proportion of all LOSET-negative patients with this priority
#'   }
#'
#' @keywords internal
extract_triage_distribution <- function(df) {
  total_loset <- sum(df$loset, na.rm = TRUE)
  total_non_loset <- sum(!df$loset, na.rm = TRUE)

  distribution <- df %>%
    group_by(priority) %>%
    summarise(
      loset_count = sum(loset, na.rm = TRUE),
      non_loset_count = sum(!loset, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(
      loset_distribution = loset_count / total_loset,
      non_loset_distribution = non_loset_count / total_non_loset
    )

  return(distribution)
}


#' Apply Triage Priority Distribution
#'
#' @description
#' Applies a triage priority distribution from one ED to another ED's patient data.
#' This function reassigns priorities to patients while preserving all other patient
#' characteristics, allowing assessment of how different triage systems would perform
#' in different settings.
#'
#' @param df Data frame containing patient data
#' @param distribution Priority distribution tibble from \code{extract_triage_distribution()}
#' @param max_distribution_diff Maximum allowed difference between target and achieved
#'   distribution. Default is 0.05 (5%)
#' @param max_tries Maximum number of attempts to achieve target distribution. Default is 10
#'
#' @return Data frame with updated priority assignments
#'
#' @details
#' The function uses an iterative approach to assign priorities:
#' 1. Calculates target counts for each priority based on the input distribution
#' 2. Randomly assigns priorities to patients according to these targets
#' 3. Adjusts assignments to ensure exact totals match
#' 4. Validates that the achieved distribution is within tolerance
#'
#' @keywords internal
apply_triage_distribution <- function(df, distribution, max_distribution_diff = 0.05, max_tries = 10) {
  # Check if we have a valid distribution
  if (nrow(distribution) == 0) {
    warning("Empty distribution provided")
    return(df)
  }

  # If priority column is missing, add it with NA values
  if (!"priority" %in% names(df)) {
    df$priority <- NA_integer_
  }

  # Get all unique priorities from the distribution (the target priorities)
  all_priorities <- distribution$priority
  all_priorities <- all_priorities[!is.na(all_priorities)]
  all_priorities <- sort(all_priorities)

  # First, assign priorities randomly according to distribution proportions
  total_cases <- nrow(df)
  total_loset <- sum(df$loset, na.rm = TRUE)
  total_non_loset <- total_cases - total_loset

  # Calculate target counts for each priority
  priority_counts <- map_dfr(all_priorities, function(p) {
    target_dist <- distribution %>% filter(priority == p)
    if (nrow(target_dist) == 0) {
      return(data.frame(priority = p, loset_count = 0, non_loset_count = 0))
    }

    loset_count <- round(target_dist$loset_distribution * total_loset)
    non_loset_count <- round(target_dist$non_loset_distribution * total_non_loset)

    return(data.frame(
      priority = p,
      loset_count = loset_count,
      non_loset_count = non_loset_count
    ))
  })

  # Adjust counts to ensure exact totals
  while (sum(priority_counts$loset_count) != total_loset) {
    diff <- total_loset - sum(priority_counts$loset_count)
    if (diff > 0) {
      # Add to a random priority
      idx <- sample(1:nrow(priority_counts), 1)
      priority_counts$loset_count[idx] <- priority_counts$loset_count[idx] + 1
    } else {
      # Remove from a random priority with at least 1
      candidates <- which(priority_counts$loset_count > 0)
      if (length(candidates) > 0) {
        idx <- sample(candidates, 1)
        priority_counts$loset_count[idx] <- priority_counts$loset_count[idx] - 1
      }
    }
  }

  while (sum(priority_counts$non_loset_count) != total_non_loset) {
    diff <- total_non_loset - sum(priority_counts$non_loset_count)
    if (diff > 0) {
      # Add to a random priority
      idx <- sample(1:nrow(priority_counts), 1)
      priority_counts$non_loset_count[idx] <- priority_counts$non_loset_count[idx] + 1
    } else {
      # Remove from a random priority with at least 1
      candidates <- which(priority_counts$non_loset_count > 0)
      if (length(candidates) > 0) {
        idx <- sample(candidates, 1)
        priority_counts$non_loset_count[idx] <- priority_counts$non_loset_count[idx] - 1
      }
    }
  }

  # Now assign priorities to the dataframe
  result_df <- df

  # Create assignment vectors
  loset_indices <- which(df$loset)
  non_loset_indices <- which(!df$loset)

  # Shuffle indices
  loset_indices <- sample(loset_indices)
  non_loset_indices <- sample(non_loset_indices)

  # Assign priorities
  loset_start <- 1
  non_loset_start <- 1

  for (i in 1:nrow(priority_counts)) {
    p <- priority_counts$priority[i]
    l_count <- priority_counts$loset_count[i]
    nl_count <- priority_counts$non_loset_count[i]

    # Assign to LOSET cases
    if (l_count > 0) {
      loset_end <- loset_start + l_count - 1
      result_df$priority[loset_indices[loset_start:loset_end]] <- p
      loset_start <- loset_end + 1
    }

    # Assign to non-LOSET cases
    if (nl_count > 0) {
      non_loset_end <- non_loset_start + nl_count - 1
      result_df$priority[non_loset_indices[non_loset_start:non_loset_end]] <- p
      non_loset_start <- non_loset_end + 1
    }
  }

  # Verify the result
  final_dist <- extract_triage_distribution(result_df)

  # Check if we're close enough
  max_diff <- 0
  for (p in all_priorities) {
    target <- distribution %>% filter(priority == p)
    actual <- final_dist %>% filter(priority == p)

    if (nrow(target) > 0 && nrow(actual) > 0) {
      loset_diff <- abs(target$loset_distribution - actual$loset_distribution)
      non_loset_diff <- abs(target$non_loset_distribution - actual$non_loset_distribution)
      max_diff <- max(max_diff, loset_diff, non_loset_diff)
    }
  }

  if (max_diff > max_distribution_diff) {
    warning(paste("Could not achieve distribution within", max_distribution_diff,
                  ". Maximum difference:", round(max_diff, 4)))
  }

  return(result_df)
}


#' Assess Cross-Site Transferability of Triage Effectiveness
#'
#' @description
#' Performs cross-site validation analysis to examine how Triage Effectiveness (TE)
#' measurements transfer between different hospital settings. This analysis extracts
#' triage priority distribution patterns from each site and simulates how each site's
#' triage pattern would perform if applied to other sites' patient populations.
#'
#' @param df Data frame containing patient data that has been initialized with \code{init()}
#' @param calc_method Character. Method to use for calculating TE: "wte" for waiting
#'   time-based TE or "rte" for rank-based TE (default)
#' @param n_iterations Integer. Number of Monte Carlo iterations per site combination
#'   to account for stochastic variation. Default is 100
#' @param n_workers Integer. Number of cores to use for parallel processing. Default
#'   is (detectCores() - 1)
#' @param max_distribution_diff Numeric. Maximum allowed difference between target and
#'   achieved priority distribution. Default is 0.005 (0.5%)
#' @param max_tries Integer. Maximum attempts to achieve target distribution in each
#'   iteration. Default is 10
#' @param seed Integer. Seed for reproducible results. Default is NULL (no seed)
#'
#' @return A list containing:
#'   \itemize{
#'     \item results: Tibble with transferability metrics for each site combination
#'     \item metadata: List with analysis parameters and settings
#'     \item distributions: Tibble with simulation data for each iteration
#'     \item priority_distributions: List of original priority distributions per site
#'   }
#'
#' @details
#' The analysis process:
#' 1. Extracts triage priority distributions for LOSET-positive and LOSET-negative
#'    patients from each site
#' 2. Calculates baseline TE values for each site using their original priority patterns
#' 3. For each combination of original patient data and applied priority distribution:
#'    - Applies the priority distribution while preserving all other patient characteristics
#'    - Simulates theoretical wait times using \code{sim_te()}
#'    - Calculates TE using either WTE or RTE method
#'    - Repeats for n_iterations to capture stochastic variation
#' 4. Summarizes results with means, confidence intervals, and variability measures
#'
#' The difference between original and redistributed TE values quantifies how much
#' triage effectiveness measurements are affected by:
#' - Accuracy of priority assignments across multiple triage levels
#' - Differences in ED throughput and patient flow characteristics
#' - Site-specific case mix and patient populations
#'
#' @examples
#' \dontrun{
#' # Initialize data
#' data <- init(sem_malmo_synth)
#'
#' # Run transferability analysis
#' results <- assess_transferability(
#'   data,
#'   calc_method = "rte",
#'   n_iterations = 100,
#'   n_workers = 4
#' )
#'
#' # View summary results
#' print(results$results)
#'
#' # Create heatmap visualization
#' plot_transferability(results)
#' }
#'
#' @importFrom dplyr filter select group_by summarise mutate bind_rows arrange
#' @importFrom furrr future_map furrr_options
#' @importFrom progressr progressor with_progress
#' @importFrom parallel detectCores
#' @importFrom purrr map map_dfr
#' @importFrom tibble tibble
#'
#' @export
assess_transferability <- function(df,
                                   calc_method = "rte",
                                   n_iterations = 100,
                                   n_workers = parallel::detectCores() - 1,
                                   max_distribution_diff = 0.005,
                                   max_tries = 10,
                                   seed = NULL) {

  # Validate input parameters
  if (!calc_method %in% c("rte", "wte")) {
    stop("calc_method must be either 'rte' or 'wte'")
  }

  # Check that data has been initialized
  if (!all(c("arrival_minute", "resolve_minute") %in% names(df))) {
    stop("Data frame appears uninitialized. Please run init() first.")
  }

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Get units
  units <- unique(df$unit)
  n_units <- length(units)
  total_sims <- n_units * n_units * n_iterations

  # Step 1: Extract distributions and calculate original TEs for each unit
  print(paste("Extracting distributions and calculating original TEs for", n_units, "units..."))
  unit_info <- map(setNames(units, units), ~ {
    unit_data <- df %>% filter(unit == .x)

    # Extract priority distribution
    distribution <- extract_triage_distribution(unit_data)

    # Calculate original TE
    sim_data <- unit_data %>% sim_te()

    if (calc_method == "rte") {
      original_te <- sim_data %>%
        calculate_queue_metrics() %>%
        filter(loset == TRUE, valid == TRUE) %>%
        summarise(mean_rte = mean(theoretical_RTE, na.rm = TRUE)) %>%
        pull(mean_rte)
    } else if (calc_method == "wte") {
      mean_all <- mean(sim_data$theoretical_wait_time, na.rm = TRUE)
      mean_loset <- mean(sim_data$theoretical_wait_time[sim_data$loset], na.rm = TRUE)
      original_te <- 1 - (mean_loset / mean_all)
    }

    list(
      distribution = distribution,
      original_te = original_te,
      n_patients = nrow(unit_data)  # Add patient count
    )
  })

  # Extract patient counts for metadata
  unit_n_patients <- df %>%
    group_by(unit) %>%
    summarise(
      n_patients = n(),
      .groups = 'drop'
    )


  # Step 2: Create all combinations
  combinations <- expand_grid(
    original_unit = units,
    applied_unit = units
  )

  # Step 3: Setup parallel processing
  setup_parallel(n_workers = n_workers)
  on.exit(cleanup_parallel())

  # Step 4: Process combinations in parallel
  print(paste("Running", nrow(combinations), "combinations with", n_iterations, "iterations each..."))

  # Initialize progress reporting
  init_progressr()

  all_results <- with_progress({
    p <- progressor(steps = total_sims)

    # Function to process a single combination
    process_combination <- function(combo_idx) {
      combo <- combinations[combo_idx, ]

      # Get pre-calculated data
      original_info <- unit_info[[combo$original_unit]]
      distribution <- original_info$distribution
      original_te <- original_info$original_te

      # Get applied ED data (remove existing priority column)
      applied_df <- df %>%
        filter(unit == combo$applied_unit) %>%
        select(-priority)

      # Run iterations
      simulation_results <- tibble(
        iteration = 1:n_iterations,
        recreated_te = numeric(n_iterations),
        diff_from_original = numeric(n_iterations)
      )

      for (j in 1:n_iterations) {
        # Apply distribution
        recreated_df <- apply_triage_distribution(
          applied_df,
          distribution,
          max_distribution_diff,
          max_tries
        )

        # Calculate TE
        sim_data <- recreated_df %>% sim_te()

        if (calc_method == "rte") {
          recreated_te <- sim_data %>%
            calculate_queue_metrics() %>%
            filter(loset == TRUE, valid == TRUE) %>%
            summarise(mean_rte = mean(theoretical_RTE, na.rm = TRUE)) %>%
            pull(mean_rte)
        } else if (calc_method == "wte") {
          mean_all <- mean(sim_data$theoretical_wait_time, na.rm = TRUE)
          mean_loset <- mean(sim_data$theoretical_wait_time[sim_data$loset], na.rm = TRUE)
          recreated_te <- 1 - (mean_loset / mean_all)
        }

        # Store results
        simulation_results$recreated_te[j] <- recreated_te
        simulation_results$diff_from_original[j] <- original_te - recreated_te

        # Update progress
        p()
      }

      # Calculate summary statistics
      summary_stats <- tibble(
        original_unit = combo$original_unit,
        applied_unit = combo$applied_unit,
        original_te = original_te,
        mean_recreated_te = mean(simulation_results$recreated_te),
        mean_diff = mean(simulation_results$diff_from_original),
        sd_diff = sd(simulation_results$diff_from_original),
        lower_ci = quantile(simulation_results$diff_from_original, 0.025),
        upper_ci = quantile(simulation_results$diff_from_original, 0.975),
        min_diff = min(simulation_results$diff_from_original),
        max_diff = max(simulation_results$diff_from_original),
        n_iterations = n_iterations
      )

      # Return both summary and simulation data
      return(list(
        summary = summary_stats,
        simulation_data = simulation_results %>%
          mutate(
            original_unit = combo$original_unit,
            applied_unit = combo$applied_unit,
            original_te = original_te
          )
      ))
    }

    # Run all combinations in parallel
    future_map(1:nrow(combinations), process_combination,
               .options = furrr_options(seed = TRUE))
  })

  # Step 5: Extract and organize results
  final_results <- bind_rows(map(all_results, ~ .x$summary))
  all_simulation_data <- bind_rows(map(all_results, ~ .x$simulation_data))

  # Create metadata
  metadata <- list(
    calculation_time = Sys.time(),
    calc_method = calc_method,
    n_iterations = n_iterations,
    n_workers = n_workers,
    max_distribution_diff = max_distribution_diff,
    max_tries = max_tries,
    seed = seed,
    n_units = n_units,
    total_simulations = total_sims,
    unit_n_patients = unit_n_patients
  )

  # Extract priority distributions for reference
  priority_distributions <- map(unit_info, ~ .x$distribution)

  # Return structured results
  return(list(
    results = final_results,
    metadata = metadata,
    distributions = all_simulation_data,
    priority_distributions = priority_distributions
  ))
}
