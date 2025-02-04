#' Initialize the Data and Settings
#'
#' This function initializes the data frame and settings for the triage efficacy calculations.
#' It ensures that the required columns are present, validates the data, casts variables to their appropriate types,
#' and adds calculated fields. If 'start' and 'stop' are not provided, they default to the earliest arrival time
#' and the latest resolve time in the data frame, respectively. It also logs success and error messages based on the
#' specified logging level.
#'
#' Of note is that we "fix" the datetimes within the function due to an error in R/lubridate that treats hms 00:00:00
#' as non existing, which means that they later are cast to just a date wich returns NA on ymd_hms casting.
#'
#' @param Df Data frame containing the data.
#' @param start Optional. The start datetime for calculating arrival and resolve minutes. If not provided, defaults to the earliest arrival time in the data frame.
#' @param stop Optional. The stop datetime for calculating arrival and resolve minutes. If not provided, defaults to the latest resolve time in the data frame.
#' @param Logg_level An integer indicating the logging level. 0 logs only errors, 1 logs successes as well.
#' @param time_critical_prio 1 or more priorities, specified within c() that will be counted as time critical when calculated sensitivity/specificity.
#' @return A data frame with validated and modified fields.
#' @examples
#' \dontrun{
#' init(Df = your_data_frame, start = start_date, stop = stop_date)
#' init(Df = your_data_frame)  # Automatically uses the earliest and latest times in the data frame
#'}
#' @importFrom dplyr select mutate filter n row_number if_else
#' @importFrom lubridate force_tz ymd_hms
#' @importFrom stats setNames
#' @importFrom utils capture.output
#' @importFrom magrittr %>%
#'
#' @export
#'

init <- function(Df, start = NULL, stop = NULL, Logg_level = 0, time_critical_prio = c(1,2)) {
  # Required columns
  required_columns <- c("arrival", "resolve", "loset", "priority", "unit")

  # Check if all required columns exist
  if (!all(required_columns %in% names(Df))) {
    if (Logg_level >= 0) stop("Error: Missing required columns.")
  }


  # This function is a fix on R/lubridates inclination to treat datetimes
  # that end in 00:00:00 as just dates, interpreting them as Na when
  # then casting them to ymd_hms. It does this by adding a second, that
  # when later looking at the minutes does not matter since we only look
  # at seconds. This is also the reason we recomend initiating the data
  # at the time of running a simulation.
  fix_times <- function(data) {
    # Convert the datetime columns to character first
    data$arrival <- as.character(data$arrival)
    data$resolve <- as.character(data$resolve)

    # Function to process each date-time string
    process_datetime <- function(dt) {
      if (grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", dt)) {
        # Date only, add time
        return(paste(dt, "00:00:01"))
      } else if (grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}T00:00:00Z$", dt)) {
        # ISO 8601 format with 00:00:00, change to 00:00:01
        return(sub("T00:00:00Z$", "T00:00:01Z", dt))
      } else {
        # Already has time or in another format, return as is
        return(dt)
      }
    }

    # Apply the processing function to arrival and resolve columns
    data$arrival <- sapply(data$arrival, process_datetime)
    data$resolve <- sapply(data$resolve, process_datetime)

    # Ensure date-time columns are in UTC using lubridate
    data$arrival <- ymd_hms(data$arrival, tz = "UTC")
    data$resolve <- ymd_hms(data$resolve, tz = "UTC")

    return(data)
  }



  #Fix times before looking at them like times.
  Df <- fix_times(Df)

  # Parse date-time columns using lubridate to ensure they are in UTC
  Df <- Df %>%
    mutate(
      arrival = force_tz(ymd_hms(arrival), "UTC"),
      resolve = force_tz(ymd_hms(resolve), "UTC")
    )

  # Check for NA values in required columns and log if found
  for (col in required_columns) {
    if (any(is.na(Df[[col]]))) {
      if (Logg_level >= 0) stop("Error: NA values found in column", col)
      stop(paste("NA values found in column:", col))
    }
  }

  # Check if any comparison results in NA
  if (any(is.na(Df$resolve < Df$arrival))) {
    if (Logg_level >= 0) {
      na_rows <- which(is.na(Df$resolve < Df$arrival))
      message("Error: NA comparison detected in 'resolve' and 'arrival' times.")
      message("Rows with NA comparison:\n", paste(capture.output(print(Df[na_rows, ])), collapse="\n"))
      message("Variable types:\n", paste(sapply(Df[na_rows, c("resolve", "arrival")], class), collapse=", "))
    }
    stop("NA comparison detected in 'resolve' and 'arrival' times.")
  }

  # Ensure 'resolve' is not before 'arrival'
  if (any(Df$resolve < Df$arrival, na.rm = TRUE)) {
    if (Logg_level >= 0) {
      invalid_rows <- which(Df$resolve < Df$arrival | is.na(Df$resolve) | is.na(Df$arrival))
      message("Error: 'resolve' time is before 'arrival' time for some records.")
      message("Invalid resolve/arrival times:\n", paste(capture.output(print(Df[invalid_rows, c("resolve", "arrival")])), collapse="\n"))
      message("Variable types:\n", paste(sapply(Df[invalid_rows, c("resolve", "arrival")], class), collapse=", "))
    }
    stop("'resolve' time is before 'arrival' time for some records.")
  }

  # Set start and stop if not provided and ensure start and stop times are in UTC
  if (is.null(start)) {
    start <- min(Df$arrival)
  } else {
    start <- force_tz(start, "UTC")
  }

  if (is.null(stop)) {
    stop <- max(Df$resolve)
  } else {
    stop <- force_tz(stop, "UTC")
  }

  # Ensure all 'arrival' and 'resolve' times are within 'start' and 'stop'
  if (any(Df$arrival < start | Df$arrival > stop)) {
    min_arrival <- min(Df$arrival, na.rm = TRUE)
    max_arrival <- max(Df$arrival, na.rm = TRUE)
    if (Logg_level >= 0) {
      message(sprintf("Error: 'arrival' times are out of the specified range. Start: %s, Min Arrival: %s, Max Arrival: %s, Stop: %s",
                      start, min_arrival, max_arrival, stop))
    }
    stop(sprintf("'arrival' times are out of the specified range. Start: %s, Min Arrival: %s, Max Arrival: %s, Stop: %s",
                 start, min_arrival, max_arrival, stop))
  }

  if (any(Df$resolve < start | Df$resolve > stop)) {
    min_resolve <- min(Df$resolve, na.rm = TRUE)
    max_resolve <- max(Df$resolve, na.rm = TRUE)
    if (Logg_level >= 0) {
      message(sprintf("Error: 'resolve' times are out of the specified range. Start: %s, Min Resolve: %s, Max Resolve: %s, Stop: %s",
                      start, min_resolve, max_resolve, stop))
    }
    stop(sprintf("'resolve' times are out of the specified range. Start: %s, Min Resolve: %s, Max Resolve: %s, Stop: %s",
                 start, min_resolve, max_resolve, stop))
  }


  # Validate so no unit totally lacks LOSET
  n_loset <- Df %>%
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


  # Cast variables to their appropriate types
  Df <- Df %>%
    mutate(
      loset = as.logical(loset),
      priority = factor(priority, ordered = TRUE),
      unit = if ("unit" %in% names(Df)) factor(unit) else NULL
    )

  # Add arrival_minute and resolve_minute fields
  Df <- Df %>%
    mutate(
      arrival_minute = round(as.numeric(difftime(arrival, start, units = "mins"))),
      resolve_minute = round(as.numeric(difftime(resolve, start, units = "mins"))),
    ) %>%
    mutate(
      priority_binary = if_else(priority %in% time_critical_prio, 1, 2),
      priority_binary = factor(priority_binary, ordered = TRUE),
      observed_wait_time = resolve_minute - arrival_minute,
      id = row_number()
    ) %>%
    mutate(
      tp = ifelse((priority %in% time_critical_prio & loset == TRUE), TRUE, FALSE),
      fp = ifelse((priority %in% time_critical_prio & loset == FALSE), TRUE, FALSE),
      fn = ifelse((!(priority %in% time_critical_prio) & loset == TRUE), TRUE, FALSE),
      tn = ifelse((!(priority %in% time_critical_prio) & loset == FALSE), TRUE, FALSE)
    )

  # Log successes if Logg_level is 1
  if (Logg_level == 1) {
    message("Data frame columns validated and cast successfully.")
    message("arrival_minute and resolve_minute fields added successfully.")
  }

  return(Df)
}


