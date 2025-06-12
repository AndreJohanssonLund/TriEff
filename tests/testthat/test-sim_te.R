# Helper functions for tests
get_test_data <- function() {
  # Start with creating a smaller dataset for testing
  data <- load_sem_synth() %>%
    dplyr::filter(arrival <= min(arrival) + lubridate::days(7))
  return(data)
}

setup_test_data <- function() {
  data <- get_test_data()
  # Initialize data before simulation
  data <- init(data)
  return(data)
}

test_that("sim_te_alt produces expected output structure", {
  options(trieff.progress = FALSE)
  data <- setup_test_data()
  result <- sim_te_alt(data, n_workers = 1)

  # Check required output columns exist
  expected_cols <- c("theoretical_wait_time", "id", "loset", "priority",
                     "priority_binary", "observed_wait_time")
  for(col in expected_cols) {
    expect_true(col %in% names(result),
                info = sprintf("Column '%s' should exist in output", col))
  }

  # Check data types
  expect_type(result$theoretical_wait_time, "double")
  expect_true(all(result$theoretical_wait_time >= 0))

  # Check row count preservation
  expect_equal(nrow(data), nrow(result))
})

test_that("sim_te and sim_te_alt produce equivalent results", {
  data <- setup_test_data()

  result_norm <- sim_te(data, n_workers = 1)
  result_alt <- sim_te_alt(data, n_workers = 1)

  # Compare theoretical wait times for LOSET cases
  loset_cases <- data$loset
  expect_equal(
    result_norm$theoretical_wait_time[loset_cases],
    result_alt$theoretical_wait_time[loset_cases]
  )
})

test_that("sim_te_alt correctly prioritizes by wait time within priority levels", {
  # Create minimal test case with known arrival/resolve times
  test_data <- data.frame(
    id = 1:4,
    arrival_minute = c(0, 1, 2, 3),
    resolve_minute = c(10, 12, 13, 14),
    priority = factor(c(1, 1, 2, 2)),
    priority_binary = factor(c(1, 1, 2, 2)),
    loset = c(TRUE, TRUE, FALSE, FALSE),
    observed_wait_time = c(10, 10, 10, 10),
    unit = factor(rep("medical", 4)),
    segment = factor(rep("med_1", 4))
  )

  result <- sim_te_alt(test_data, n_workers = 1)

  # First priority 1 case should be resolved first
  expect_true(result$theoretical_wait_time[1] < result$theoretical_wait_time[2])
})

test_that("sim_te_at handles parallel processing correctly", {
  skip_on_cran()  # Skip on CRAN as it might have limited resources
  data <- setup_test_data()

  # Compare results with different numbers of workers
  result_single <- sim_te_alt(data, n_workers = 1)
  result_multi <- sim_te_alt(data, n_workers = 2)

  expect_equal(result_single$theoretical_wait_time,
               result_multi$theoretical_wait_time)
})

test_that("sim_te_alt respects priority order within time windows", {
  # Create test data where we know the expected order
  test_data <- data.frame(
    id = 1:3,
    arrival_minute = c(0, 1, 2),
    resolve_minute = c(10, 11, 12),
    priority = factor(c(1, 2, 1)),  # First and last are high priority
    priority_binary = factor(c(1, 2, 1)),
    loset = c(TRUE, FALSE, TRUE),
    observed_wait_time = c(10, 10, 10),
    unit = factor(rep("medical", 3)),
    segment = factor(rep("med_1", 3))
  )

  result <- sim_te_alt(test_data, n_workers = 1)

  # The last patient (priority 1) should get seen before the middle patient (priority 2)
  # despite arriving later
  expect_true(result$theoretical_wait_time[3] < result$theoretical_wait_time[2])

  # Check that theoretical wait times are consistent with arrival/resolve times
  expect_true(all(result$theoretical_wait_time >= 0))
})

test_that("segmentation does not affect simulation results", {
  # Note: This test deliberately uses the internal run_simulation function to verify
  # that the segmentation logic in sim_te doesn't affect results. While generally
  # we avoid testing internal functions, this is necessary to ensure segmentation
  # maintains simulation accuracy.

  # Create a small dataset where we can process it both ways
  test_data <- data.frame(
    id = 1:6,
    arrival_minute = c(0, 1, 10, 11, 20, 21),  # Three potential segments
    resolve_minute = c(5, 6, 15, 16, 25, 26),
    priority = factor(c(2, 1, 2, 1, 2, 1)),
    priority_binary = factor(c(2, 1, 2, 1, 2, 1)),
    loset = c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE),
    observed_wait_time = rep(5, 6),
    unit = factor(rep("medical", 6)),
    segment = factor(rep("med_1", 6))
  )

  # Run with normal segmentation
  result_segmented <- sim_te_alt(test_data, n_workers = 1)

  # Process all data at once using internal function
  result_single <- run_simulation(test_data)

  # Results should be identical
  expect_equal(result_segmented$theoretical_wait_time,
               result_single$theoretical_wait_time)
})

test_that("sim_te_alt handles data with no LOSET cases", {
  data <- setup_test_data()
  data$loset <- FALSE
  result <- sim_te_alt(data, n_workers = 1)

  # Check that theoretical wait times are still calculated
  expect_false(any(is.na(result$theoretical_wait_time)))

})

test_that("sim_te_alt correctly handles btte parameter", {
  data <- setup_test_data()

  # Test with btte = TRUE
  result <- sim_te_alt(data, tte = FALSE, btte = TRUE, n_workers = 1)
  expect_true("binary_theoretical_wait_time" %in% names(result))
  expect_false("theoretical_wait_time" %in% names(result))

  # Verify binary wait times respect binary priorities
  binary_high_priority <- result$priority_binary == 1
  avg_wait_high <- mean(result$binary_theoretical_wait_time[binary_high_priority])
  avg_wait_low <- mean(result$binary_theoretical_wait_time[!binary_high_priority])
  expect_true(avg_wait_high < avg_wait_low)
})

test_that("segmentation produces identical results to unsegmented processing when done by unit", {
  data <- setup_test_data()

  # Initialize the data
  data <- init(data)

  # Process with normal segmentation using sim_te
  result_segmented <- sim_te(data, n_workers = 1)

  # Process by unit using run_simulation
  df_sim <- data %>%
    select(id, unit, arrival_minute, resolve_minute,
           priority, priority_binary, loset)

  # Split by unit and process each separately
  unit_results <- list()
  for(unit_name in unique(df_sim$unit)) {
    unit_data <- df_sim %>% filter(unit == unit_name)
    unit_results[[unit_name]] <- run_simulation(unit_data)
  }

  # Combine unit results
  result_unsegmented <- bind_rows(unit_results)

  # Sort both results by id
  result_segmented <- result_segmented %>% arrange(id)
  result_unsegmented <- result_unsegmented %>% arrange(id)

  # Compare theoretical wait times
  expect_equal(
    result_segmented$theoretical_wait_time,
    result_unsegmented$theoretical_wait_time,
    info = "Segmented and unsegmented processing should produce identical wait times"
  )

  # Also compare binary theoretical wait times if present
  if ("binary_theoretical_wait_time" %in% names(result_segmented)) {
    expect_equal(
      result_segmented$binary_theoretical_wait_time,
      result_unsegmented$binary_theoretical_wait_time,
      info = "Segmented and unsegmented processing should produce identical binary wait times"
    )
  }

  # Verify order preservation within priority levels for each unit
  for(unit_name in unique(df_sim$unit)) {
    for (prio in unique(data$priority)) {
      prio_cases_seg <- result_segmented %>%
        filter(unit == unit_name, priority == prio) %>%
        arrange(arrival_minute) %>%
        pull(theoretical_wait_time)

      prio_cases_unseg <- result_unsegmented %>%
        filter(unit == unit_name, priority == prio) %>%
        arrange(arrival_minute) %>%
        pull(theoretical_wait_time)

      expect_equal(
        prio_cases_seg,
        prio_cases_unseg,
        info = sprintf("Wait time order within priority level %s for unit %s should be identical",
                       prio, unit_name)
      )
    }
  }
})
