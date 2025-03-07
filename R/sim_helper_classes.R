#' Create Segment Identifiers in Patient Data
#'
#' This function adds a segment identifier column to the input patient data based on continuous queue periods.
#' Each segment represents a period where the queue starts at 0 and ends at 0, allowing for
#' independent processing of these time periods.
#'
#' @param df A data frame containing patient data, including ID, unit, arrival times, and resolve times.
#' @param n_workers Integer. Number of cores to use for parallel processing. Default is (detectCores() - 1).
#'
#' @return A data frame with an additional 'segment' column that identifies continuous queue segments.
#'
#' @importFrom data.table as.data.table setkey setorder
#' @importFrom dplyr select mutate filter arrange bind_rows group_by ungroup summarize inner_join
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_longer
#' @importFrom purrr map map_dfr map_int
#' @importFrom furrr future_map furrr_options
#'
#' @details
#' The function identifies periods where the queue size goes from 0 to 1, marking the start of a new segment.
#' It then assigns segment IDs to patients based on these boundaries, which can be used for grouped processing
#' in parallel simulations. This segmentation is crucial for maintaining chronological integrity while enabling
#' parallelization.
#'
#' @keywords internal
create_segments <- function(df, n_workers = detectCores() - 1) {
  # Split data by unit for parallel processing
  unit_data <- split(df, df$unit)

  # Set up parallel processing
  setup_parallel(n_workers = n_workers)
  on.exit(cleanup_parallel())

  # Initialize progress reporting
  init_progressr()

  # Process each unit in parallel to identify segments
  segment_results <- future_map(names(unit_data), function(unit_name) {
    unit_df <- unit_data[[unit_name]]

    # Create events table for this unit
    events <- unit_df %>%
      select(id, arrival_minute, resolve_minute) %>%
      pivot_longer(cols = c(arrival_minute, resolve_minute),
                   names_to = "event_type",
                   values_to = "event_minute") %>%
      mutate(event_type = ifelse(event_type == "arrival_minute", "arrival", "resolve")) %>%
      arrange(event_minute, event_type)

    # Calculate queue sizes and segment boundaries
    queue_df <- events %>%
      mutate(
        cum_arrivals = cumsum(event_type == "arrival"),
        cum_resolves = cumsum(event_type == "resolve"),
        queue_size = cum_arrivals - cum_resolves,
        segment_start = queue_size == 1 & lag(queue_size, default = 0) == 0,
        segment_id = cumsum(segment_start)
      )

    # Map segment_id back to original patient data
    id_to_segment <- queue_df %>%
      filter(event_type == "arrival") %>%
      select(id, segment_id) %>%
      distinct()

    # Join segment IDs back to unit data
    unit_df %>%
      left_join(id_to_segment, by = "id") %>%
      # Ensure unique segment IDs across units by prefixing with unit name
      mutate(segment = paste(unit_name, segment_id, sep = "_"))
  }, .options = furrr::furrr_options(seed = TRUE))

  # Combine results from all units
  result_df <- bind_rows(segment_results)

  # Clean up temporary column if it exists
  if ("segment_id" %in% names(result_df)) {
    result_df <- result_df %>% select(-segment_id)
  }

  return(result_df)
}




#' Create Batches from Segments
#'
#' This function organizes segments of patient data into batches of approximately equal size.
#' Batching is used to distribute the workload evenly across available processing cores for
#' efficient parallel computation.
#'
#' @param segments A list of data frame segments created by the create_segments function.
#' @param n_batches Integer. The number of batches to create, typically equal to the number of available cores.
#'
#' @return A list of batches, where each batch contains one or more segments.
#'
#' @details
#' The function aims to create batches of roughly equal total size (in terms of number of patients)
#' to ensure balanced workload distribution in parallel processing. This approach optimizes
#' the use of computational resources during simulation.
#'
#' @keywords internal
create_batches <- function(segments, n_batches) {
  total_size <- sum(map_int(segments, nrow))
  target_size <- ceiling(total_size / n_batches)

  batches <- vector("list", n_batches)
  current_batch <- 1
  current_size <- 0

  for (segment in segments) {
    segment_size <- nrow(segment)
    if (current_size + segment_size > target_size && current_batch < n_batches) {
      current_batch <- current_batch + 1
      current_size <- 0
    }
    batches[[current_batch]] <- c(batches[[current_batch]], list(segment))
    current_size <- current_size + segment_size
  }

  batches
}



#' Run Simulation on a Single Segment
#'
#' This function performs the core simulation logic for calculating Theoretical Triage Efficacy (TTE)
#' and Binary Theoretical Triage Efficacy (BTTE) for a single segment of patient data.
#'
#' @param segment A data.table representing a segment of patient data.
#' @param tte Logical. If TRUE, calculates Theoretical Triage Efficacy (TTE).
#' @param btte Logical. If TRUE, calculates Binary Theoretical Triage Efficacy (BTTE).
#'
#' @return A data.table with calculated wait times and resolve times.
#'
#' @details
#' This function simulates patient flow through a healthcare system in a deterministic manner,
#' strictly adhering to priority order and wait time within priority. It's optimized for performance
#' using data.table operations. It is used both by sim_te and sim_heat.
#'
#' @importFrom data.table := .I
#'
#' @keywords internal
run_simulation <- function(segment, tte = TRUE, btte = FALSE) {

  dt <- data.table::as.data.table(segment)
  data.table::setkey(dt, id)

  resolve_events <- dt[, .(id = NA_integer_, resolve_minute = resolve_minute, priority = NA_integer_)]
  data.table::setorder(resolve_events, resolve_minute)
  resolve_events[, resolve_event_order := .I]

  if(tte) {
    dt[, resolve_order := NA_integer_]
    for (i in 1:nrow(resolve_events)) {
      resolve_time <- resolve_events$resolve_minute[i]
      eligible_index <- which(dt$arrival_minute <= resolve_time & is.na(dt$resolve_order))
      if (length(eligible_index) > 0) {
        order_index <- order(dt$priority[eligible_index], dt$arrival_minute[eligible_index])[1]
        dt$resolve_order[eligible_index[order_index]] <- i
      }
    }
    dt[, theoretical_resolve_time := resolve_events$resolve_minute[resolve_order]]
    dt[, theoretical_wait_time := theoretical_resolve_time - arrival_minute]
  }

  if (btte) {
    dt[, binary_resolve_order := NA_integer_]
    for (i in 1:nrow(resolve_events)) {
      resolve_time <- resolve_events$resolve_minute[i]
      eligible_index <- which(dt$arrival_minute <= resolve_time & is.na(dt$binary_resolve_order))
      if (length(eligible_index) > 0) {
        order_index <- order(dt$priority_binary[eligible_index], dt$arrival_minute[eligible_index])[1]
        dt$binary_resolve_order[eligible_index[order_index]] <- i
      }
    }
    dt[, binary_theoretical_resolve_time := resolve_events$resolve_minute[binary_resolve_order]]
    dt[, binary_theoretical_wait_time := binary_theoretical_resolve_time -arrival_minute]
  }

  # Return as data.table to avoid conversion
  return(dt)
}




#' Validate Data Structure for sim_heat Functions
#'
#' @param df Data frame to validate containing patient data
#'
#' @return TRUE if validation passes, otherwise throws an error with descriptive message
#'
#' @details
#' Performs validation checks for sim_heat and sim_heat_alt functions including:
#' - Required column presence
#' - Column data types
#' - Initialization state
#'
#' Required columns:
#' - arrival_minute (numeric, added by init)
#' - resolve_minute (numeric, added by init)
#' - priority_binary (factor, added by init)
#' - loset (logical)
#' - observed_wait_time (numeric, added by init)
#' - unit (factor)
#'
#' @keywords internal
validate_sim_heat_data <- function(df) {
  # Required columns for sim_heat functions
  required_cols <- c(
    "id",            # Added by init()
    "arrival_minute",     # Added by init()
    "resolve_minute",     # Added by init()
    "loset",             # Required original column
    "observed_wait_time", # Added by init()
    "unit"               # Required original column
  )

  # Check if all required columns exist
  missing_cols <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0) {
    # Check specifically for init-added columns
    init_cols <- c("arrival_minute", "resolve_minute", "observed_wait_time")
    missing_init_cols <- intersect(missing_cols, init_cols)

    if (length(missing_init_cols) > 0) {
      stop("Data frame appears uninitialized. Please run init() first. ",
           "Missing columns: ", paste(missing_cols, collapse = ", "))
    } else {
      stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }
  }

  # Additional validation for column types
  if (!is.numeric(df$arrival_minute)) stop("arrival_minute must be numeric. Please run init() first.")
  if (!is.numeric(df$resolve_minute)) stop("resolve_minute must be numeric. Please run init() first.")
  if (!is.logical(df$loset)) stop("loset must be logical")
  if (!is.numeric(df$observed_wait_time)) stop("observed_wait_time must be numeric. Please run init() first.")

  # Validation passed
  return(TRUE)
}
