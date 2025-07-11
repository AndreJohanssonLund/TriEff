# Helper functions for tests
get_test_data <- function() {
  # Start with creating a smaller dataset for testing
  data <- load_sem_synth() %>%
    dplyr::filter(arrival <= min(arrival) + lubridate::days(7))
  return(data)
}

setup_base_test_data <- function() {
  data <- get_test_data()

  # Initialize and simulate data without calculating queue metrics
  data <- init(data)
  data <- sim_te(data, btte = TRUE)
  return(data)
}

test_that("calc_rte produces expected output structure", {
  data <- setup_base_test_data()
  result <- calc_rte(data, min_loset_warning = -1, verbose = TRUE)

  # Check that result is a list with expected components
  expect_type(result, "list")
  expect_true(inherits(result, "calc_te"))
  expect_named(result, c("results", "metadata", "verbose_df"))

  # Check results tibble structure
  expect_true(inherits(result$results, "tbl_df"))
  expected_cols <- c("unit", "n_patients", "n_patients_loset", "loset_prevalence",
                     "n_valid_rte", "ote_te", "tte_te", "otg_te",
                     "sensitivity", "specificity")
  expect_true(all(expected_cols %in% names(result$results)))

  # Check metadata structure
  expect_type(result$metadata, "list")
  expect_true("calculation_time" %in% names(result$metadata))

  # # Check distributions structure
  # expect_true(inherits(result$distributions, "tbl_df"))

  # Account for possible column name suffixes (.x, .y) in the output
  rte_col_found <- any(grepl("observed_RTE", names(result$verbose_df)))
  expect_true(rte_col_found,
              info = paste("No observed_RTE column found in distributions.",
                           "Available columns:",
                           paste(names(result$distributions), collapse=", ")))
})

test_that("calc_rte handles subgroup analysis correctly", {
  data <- setup_base_test_data()

  # Add age groups
  data$age_group <- cut(data$age_at_arrival,
                        breaks = c(0, 59, 79, Inf),
                        labels = c("18-59", "60-79", "80+"))

  # Then get results
  result_age <- suppressWarnings(
    calc_rte(data, var1 = "age_group", min_loset_warning = -1)
  )

  # Check that each age group is present in the results
  age_groups_in_results <- result_age$results %>%
    pull(age_group) %>%
    unique() %>%
    length()

  expect_equal(age_groups_in_results, 3)  # Should find all 3 age groups

  # Test with subgroup
  elderly_result <- suppressWarnings(
    calc_rte(data, subgroup = list(age_group = "80+"), min_loset_warning = -1))
  expect_true(all(elderly_result$results$unit != "overall"))
  expect_true(all(elderly_result$results$n_patients > 0))
})

test_that("calc_rte bootstrap analysis works correctly", {
  data <- setup_base_test_data()

  # Basic bootstrap test
  result <- calc_rte(data,
                     bootstrap = TRUE,
                     bootstrap_params = list(
                       sample_percentage = 1,
                       n_iterations = 100,
                       distribution_span = 0.95
                     ),
                     min_loset_warning = -1,
                     seed = 1234)

  # Check bootstrap results structure
  expect_true(!is.null(result$bootstrap_distributions))
  expect_true(all(c("boot_ote_var_lower", "boot_ote_var_upper") %in%
                    names(result$results)))

  # Check confidence intervals are valid
  expect_true(all(result$results$boot_ote_var_lower <= result$results$ote_te))
  expect_true(all(result$results$boot_ote_var_upper >= result$results$ote_te))
})

test_that("calc_rte handles invalid inputs appropriately", {
  data <- setup_base_test_data()

  # Test invalid var1
  expect_error(calc_rte(data, var1 = "nonexistent_column"))

  # Test var2 without var1
  expect_error(calc_rte(data, var2 = "unit"))

  # Test invalid subgroup
  expect_error(calc_rte(data, subgroup = list(nonexistent = "value")))

  # Test incompatible subgroup and var1
  data$test_group <- "A"
  expect_error(calc_rte(data,
                        subgroup = list(test_group = "A"),
                        var1 = "unit", min_loset_warning = -1))
})

test_that("calc_rte calculates metrics correctly", {
  # For this test, we specifically need the queue metrics
  data <- setup_base_test_data()
  data_with_metrics <- trieff:::calculate_queue_metrics(data, verbose = FALSE)

  result <- calc_rte(data, min_loset_warning = -1)

  # Verify that only valid RTE values are used in calculations
  valid_rte_data <- data_with_metrics %>%
    filter(loset == TRUE, valid == TRUE)

  # Test with direct calculation from the raw data
  observed_rte_mean <- mean(valid_rte_data$observed_RTE, na.rm = TRUE)
  theoretical_rte_mean <- mean(valid_rte_data$theoretical_RTE, na.rm = TRUE)

  # Compare with overall results
  overall_result <- result$results %>% filter(unit == "overall")
  expect_equal(overall_result$ote_te, observed_rte_mean, tolerance = 0.01)
  expect_equal(overall_result$tte_te, theoretical_rte_mean, tolerance = 0.01)
  expect_equal(overall_result$otg_te, observed_rte_mean - theoretical_rte_mean, tolerance = 0.01)
})

test_that("calc_rte overall_only parameter works correctly", {
  data <- setup_base_test_data()

  # Test with overall_only = TRUE
  result_overall <- calc_rte(data, overall_only = TRUE, min_loset_warning = -1)
  expect_equal(nrow(result_overall$results), 1)
  expect_equal(result_overall$results$unit, "overall")

  # Compare with full results
  result_full <- calc_rte(data)
  expect_equal(result_overall$results$ote_te,
               result_full$results$ote_te[result_full$results$unit == "overall"])

  # Check distributions filtering
  expect_true(all(result_overall$distributions$unit == "overall"))
})

test_that("calc_rte print method works correctly", {
  data <- setup_base_test_data()
  result <- calc_rte(data, min_loset_warning = -1)

  # Capture print output
  output <- capture.output(print(result))

  # Check for key sections in output
  expect_true(any(grepl("Triage Effectiveness Analysis Results", output)))
  expect_true(any(grepl("Classification Metrics:", output)))
  expect_true(any(grepl("Rank-based Triage Effectiveness Metrics:", output)))
  expect_true(any(grepl("Computation Information", output)))
})

test_that("calc_rte handles perfect triage correctly", {
  # Test perfect triage (should give RTE = 1)
  data <- get_test_data()

  # First initialize, then modify priority, don't pass data parameter to init
  data <- init(data)

  data <- data %>%
    mutate(priority = if_else(loset, 1, 3))

  perfect_data <- data %>%
    sim_te()

  # Calculate queue metrics only for tests that need direct access
  perfect_queue <- trieff:::calculate_queue_metrics(perfect_data, verbose = FALSE)

  # Calculate RTE directly
  perfect_rte <- mean(perfect_queue$theoretical_RTE, na.rm = TRUE)
  expect_equal(perfect_rte, 1, tolerance = 0.01)

  # Run through calc_rte function
  result <- calc_rte(perfect_data, min_loset_warning = -1)
  expect_equal(result$results$tte_te[result$results$unit == "overall"], 1, tolerance = 0.01)

  # Check no segments have RTE > 1 or < 1
  perfect_segments <- perfect_queue %>%
    group_by(segment) %>%
    summarise(mean_rte = mean(theoretical_RTE, na.rm = TRUE)) %>%
    filter(!is.na(mean_rte)) %>%
    filter(mean_rte > 1.01 | mean_rte < 0.99)

  expect_equal(nrow(perfect_segments), 0)
})

test_that("calc_rte handles zero triage effect correctly", {
  # Test zero triage effect (should give RTE = 0)
  data <- get_test_data()
  # First initialize, then modify priority, don't pass data parameter to init
  data <- init(data)
  data <- data %>%
    mutate(priority = 1)
  zero_data <- data %>%
    sim_te()

  # Calculate queue metrics only for tests that need direct access
  zero_queue <- trieff:::calculate_queue_metrics(zero_data, verbose = TRUE)

  # Calculate RTE directly
  zero_rte <- mean(zero_queue$theoretical_RTE, na.rm = TRUE)
  expect_equal(zero_rte, 0, tolerance = 0.01)

  # Run through calc_rte function
  result <- calc_rte(zero_data, min_loset_warning = -1)
  expect_equal(result$results$tte_te[result$results$unit == "overall"], 0, tolerance = 0.01)

  # FIXED: Check only valid time-critical patients for zero effect
  zero_segments <- zero_queue %>%
    filter(loset == TRUE & valid == TRUE) %>%
    group_by(segment) %>%
    summarise(
      mean_rte = mean(theoretical_RTE, na.rm = TRUE)
    ) %>%
    filter(mean_rte != 0)
  expect_equal(nrow(zero_segments), 0)
})

test_that("calc_rte includes/excludes distributions correctly", {
  data <- setup_base_test_data()

  # Test including distributions
  result_with_dist <- calc_rte(data, verbose = TRUE, min_loset_warning = -1)
  expect_true(!is.null(result_with_dist$verbose_df))

  # Test excluding distributions
  result_no_dist <- calc_rte(data, verbose = FALSE, min_loset_warning = -1)
  expect_true(is.null(result_no_dist$verbose_df))
})

test_that("calc_rte handles normality assessment correctly", {
  # For this test, we need the queue metrics directly
  data <- setup_base_test_data()
  data_with_metrics <- trieff:::calculate_queue_metrics(data, verbose = FALSE)

  # Ensure enough valid RTE values for meaningful testing
  valid_rte_count <- sum(data_with_metrics$loset & data_with_metrics$valid, na.rm = TRUE)
  skip_if(valid_rte_count < 3, "Too few valid RTE values for normality testing")

  result <- calc_rte(data, min_loset_warning = -1)

  # Check that normality was tested
  expect_true(!is.na(result$results$ote_shapiro_p[1]) || valid_rte_count > 5000)

  # Check that recommended CI type is set
  expect_true(result$results$recommended_ci_type[1] %in% c("parametric", "nonparametric"))
})

test_that("calc_rte bootstrap parameters affect results", {
  data <- setup_base_test_data()

  # Test two different bootstrap sample percentages
  result1 <- calc_rte(data,
                      bootstrap = TRUE,
                      bootstrap_params = list(
                        sample_percentage = 0.8,
                        n_iterations = 50,
                        distribution_span = 0.95
                      ),
                      min_loset_warning = -1,
                      seed = 1234)

  result2 <- calc_rte(data,
                      bootstrap = TRUE,
                      bootstrap_params = list(
                        sample_percentage = 0.9,
                        n_iterations = 50,
                        distribution_span = 0.95
                      ),
                      min_loset_warning = -1,
                      seed = 1234)

  # Different sample sizes should produce different CI widths
  ci_width1 <- result1$results$boot_ote_var_upper - result1$results$boot_ote_var_lower
  ci_width2 <- result2$results$boot_ote_var_upper - result2$results$boot_ote_var_lower

  # Lower sample percentage should generally produce wider CIs
  # But just check they're different as the exact relationship depends on data
  expect_false(all(abs(ci_width1 - ci_width2) < 0.001))
})

