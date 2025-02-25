# Helper functions for tests
get_test_data <- function() {
  # Start with creating a smaller dataset for testing
  data <- trieff::sem_malmo_synth[
    trieff::sem_malmo_synth$arrival >= min(trieff::sem_malmo_synth$arrival) &
      trieff::sem_malmo_synth$arrival <= min(trieff::sem_malmo_synth$arrival) + lubridate::days(7),
  ]
  return(data)
}


# Test helpers specific to init testing
setup_test_data <- function() {
  data <- get_test_data()
  # Ensure we have valid test data
  stopifnot(all(c("arrival", "resolve", "loset", "priority", "unit") %in% names(data)))
  data
}

test_that("init accepts valid data and returns expected structure", {
  data <- setup_test_data()
  result <- init(data)

  # Check required output columns exist
  expected_cols <- c("arrival_minute", "resolve_minute", "priority_binary",
                     "observed_wait_time", "tp", "fp", "fn", "tn")
  for(col in expected_cols) {
    expect_true(col %in% names(result),
                info = sprintf("Column '%s' should exist in output", col))
  }

  # Check data types
  expect_type(result$arrival_minute, "double")
  expect_type(result$resolve_minute, "double")
  expect_s3_class(result$priority, "factor")
  expect_type(result$loset, "logical")
  expect_s3_class(result$priority_binary, "factor")
})

test_that("init handles edge case data volumes", {
  data <- setup_test_data()

  # Test single row
  single_row <- data[1,]
  # Make sure this single row has LOSET
  single_row$loset <- TRUE
  expect_error(init(single_row), NA)

  # Test two rows with same arrival time
  same_arrival <- data[1:2,]
  same_arrival$arrival[2] <- same_arrival$arrival[1]
  same_arrival$resolve[2] <- same_arrival$resolve[1] + lubridate::minutes(30)
  # Make sure we have LOSET for each unit
  units <- unique(same_arrival$unit)
  for(i in seq_along(units)) {
    idx <- which(same_arrival$unit == units[i])
    if(length(idx) > 0) {
      # Set LOSET for the first row of each unit
      same_arrival$loset[idx[1]] <- TRUE
    }
  }
  expect_error(init(same_arrival), NA)

  # Test two rows with same resolve time
  same_resolve <- data[1:2,]
  same_resolve$resolve[2] <- same_resolve$resolve[1]
  # Make sure we have LOSET for each unit
  units <- unique(same_resolve$unit)
  for(i in seq_along(units)) {
    idx <- which(same_resolve$unit == units[i])
    if(length(idx) > 0) {
      # Set LOSET for the first row of each unit
      same_resolve$loset[idx[1]] <- TRUE
    }
  }
  expect_error(init(same_resolve), NA)
})

test_that("init enforces correct temporal order", {
  data <- setup_test_data()

  # Test resolve before arrival
  invalid_data <- data[1:10,]
  invalid_data$resolve <- invalid_data$arrival - lubridate::minutes(30)
  # Add LOSET cases for each unit
  invalid_data$loset <- FALSE
  for(unit_val in unique(invalid_data$unit)) {
    idx <- which(invalid_data$unit == unit_val)[1]
    invalid_data$loset[idx] <- TRUE
  }
  expect_error(init(invalid_data),
               pattern = "'resolve' time is before 'arrival' time")

  # Test exact same times
  same_time_data <- data[1:10,]
  same_time_data$resolve <- same_time_data$arrival
  # Add LOSET cases for each unit
  same_time_data$loset <- FALSE
  for(unit_val in unique(same_time_data$unit)) {
    idx <- which(same_time_data$unit == unit_val)[1]
    same_time_data$loset[idx] <- TRUE
  }
  result <- init(same_time_data)
  expect_true(all(result$observed_wait_time == 0))
})

test_that("init handles missing data appropriately", {
  data <- setup_test_data()

  # Test missing arrival
  na_arrival <- data
  na_arrival$arrival[1] <- NA
  expect_error(init(na_arrival))

  # Test missing resolve
  na_resolve <- data
  na_resolve$resolve[1] <- NA
  expect_error(init(na_resolve))

  # Test missing priority
  na_priority <- data
  na_priority$priority[1] <- NA
  expect_error(init(na_priority))
})

test_that("init calculates classification metrics correctly", {
  data <- setup_test_data()
  result <- init(data)

  # Verify tp, fp, fn, tn are mutually exclusive
  expect_true(all(!(result$tp & result$fp)))
  expect_true(all(!(result$tp & result$fn)))
  expect_true(all(!(result$tp & result$tn)))
  expect_true(all(!(result$fp & result$fn)))
  expect_true(all(!(result$fp & result$tn)))
  expect_true(all(!(result$fn & result$tn)))

  # Every row should have exactly one TRUE
  row_sums <- result$tp + result$fp + result$fn + result$tn
  expect_true(all(row_sums == 1))

  # Verify classification logic
  high_priority_cases <- result$priority %in% c(1,2)
  expect_true(all(result$tp == (high_priority_cases & result$loset)))
  expect_true(all(result$fp == (high_priority_cases & !result$loset)))
  expect_true(all(result$fn == (!high_priority_cases & result$loset)))
  expect_true(all(result$tn == (!high_priority_cases & !result$loset)))
})

test_that("init respects custom time_critical_prio parameter", {
  data <- setup_test_data()

  # Test with only priority 1 as time-critical
  result1 <- init(data, time_critical_prio = c(1))
  expect_true(all(result1$tp[result1$priority == 2] == FALSE))

  # Test with priorities 1,2,3 as time-critical
  result2 <- init(data, time_critical_prio = c(1,2,3))
  expect_true(any(result2$tp[result2$priority == 3]))
})

test_that("init handles custom start/stop times correctly", {
  data <- setup_test_data()
  min_time <- min(data$arrival)
  max_time <- max(data$resolve)
  mid_time <- min_time + (max_time - min_time)/2

  # Test custom start time
  filtered_data <- data[data$arrival >= mid_time,]
  result1 <- init(filtered_data, start = mid_time)
  expect_true(all(result1$arrival_minute >= 0))

  # Test with invalid start time (should error)
  expect_error(init(data, start = max_time))

  # Test custom stop time
  filtered_data2 <- data[data$resolve <= mid_time,]
  result2 <- init(filtered_data2, stop = mid_time)
  expect_false(all(result2$resolve_minute <= max(result2$arrival_minute)))

  # Test custom start and stop time together
  middle_data <- data[data$arrival >= mid_time & data$resolve <= max_time,]
  result3 <- init(middle_data, start = mid_time, stop = max_time)
  expect_true(min(result3$arrival_minute) >= 0)

})

test_that("init maintains data consistency", {
  data <- setup_test_data()
  result <- init(data)

  # Check row count preservation
  expect_equal(nrow(data), nrow(result))

  # Check waiting time calculation consistency
  expect_equal(
    result$observed_wait_time,
    result$resolve_minute - result$arrival_minute
  )

  # Verify priority_binary is properly derived from priority
  expect_true(all(result$priority_binary[result$priority %in% c(1,2)] == 1))
  expect_true(all(result$priority_binary[!result$priority %in% c(1,2)] == 2))
})


test_that("init handles datetime edge cases correctly", {
  data <- setup_test_data()

  # Test midnight times (00:00:00)
  midnight_data <- data[1:10,]
  midnight_data$arrival <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  midnight_data$resolve <- as.POSIXct("2024-01-01 00:00:01", tz = "UTC")
  # Add at least one LOSET case per unit to avoid warnings
  midnight_data$loset <- FALSE
  units <- unique(midnight_data$unit)
  for(unit_val in units) {
    idx <- which(midnight_data$unit == unit_val)[1]
    midnight_data$loset[idx] <- TRUE
  }
  expect_no_error(init(midnight_data))

  # Test timezone handling - ensure LOSET cases
  tz_data <- data[1:10,]
  tz_data$arrival <- as.POSIXct("2024-01-01 12:00:00", tz = "America/New_York")
  tz_data$resolve <- as.POSIXct("2024-01-01 13:00:00", tz = "America/New_York")
  tz_data$loset <- FALSE
  for(unit_val in unique(tz_data$unit)) {
    idx <- which(tz_data$unit == unit_val)[1]
    tz_data$loset[idx] <- TRUE
  }
  result <- init(tz_data)
  expect_true(all(!is.na(result$arrival_minute)))
})

test_that("init handles data with no LOSET cases", {
  data <- setup_test_data()
  # Set all LOSET to FALSE to deliberately trigger the warning
  data$loset <- FALSE

  # Expect the specific warning
  expect_warning(
    result <- init(data),
    regexp = "The following units have no LOSET cases"
  )

  # Continue with the checks on the result from the initialization
  # Check that classification metrics are correct
  expect_true(all(!result$tp))  # No true positives
  expect_true(all(!result$fn))  # No false negatives
  expect_true(all(xor(result$fp, result$tn)))  # Every case should be either fp or tn

  # Verify fp/tn based on priority
  high_priority <- result$priority %in% c(1,2)
  expect_true(all(result$fp == high_priority))  # All high priority cases should be false positives
  expect_true(all(result$tn == !high_priority))  # All low priority cases should be true negatives

  # Verify other functionality still works
  expect_true(all(!is.na(result$arrival_minute)))
  expect_true(all(!is.na(result$resolve_minute)))
  expect_true(all(result$resolve_minute >= result$arrival_minute))
})
