# Helper functions for tests
get_test_data <- function() {
  # Start with creating a smaller dataset for testing
  data <- trieff::sem_malmo_synth[
    trieff::sem_malmo_synth$arrival >= min(trieff::sem_malmo_synth$arrival) &
      trieff::sem_malmo_synth$arrival <= min(trieff::sem_malmo_synth$arrival) + lubridate::days(7),
  ]
  return(data)
}

setup_test_data <- function() {
  data <- get_test_data()
  # Initialize and simulate data before testing calc_te
  data <- init(data)
  data <- sim_te(data, btte = TRUE)
  return(data)
}

test_that("calc_te produces expected output structure", {
  data <- setup_test_data()
  result <- calc_te(data, min_loset_warning = -1)

  # Check that result is a list with expected components
  expect_type(result, "list")
  expect_true(inherits(result, "calc_te"))
  expect_named(result, c("results", "metadata"))

  # Check results tibble structure
  expect_true(inherits(result$results, "tbl_df"))
  expected_cols <- c("unit", "n_patients", "n_patients_loset", "loset_prevalence",
                     "mean_all", "ote_mean_loset", "ote_te", "tte_mean_loset",
                     "tte_te", "otg_te", "sensitivity", "specificity")
  expect_true(all(expected_cols %in% names(result$results)))

  # Check metadata structure
  expect_type(result$metadata, "list")
  expect_true("calculation_time" %in% names(result$metadata))
})

test_that("calc_te handles subgroup analysis correctly", {
  data <- setup_test_data()

  # Add age groups
  data$age_group <- cut(data$age_at_arrival,
                        breaks = c(0, 59, 79, Inf),
                        labels = c("18-59", "60-79", "80+"))

  # Test with var1
  result_age <- calc_te(data, var1 = "age_group", min_loset_warning = -1)
  expect_true("age_group" %in% names(result_age$results))
  expect_equal(length(unique(result_age$results$age_group)), 3)

  # Test with subgroup
  elderly_result <- calc_te(data,
                            subgroup = list(age_group = "80+"), min_loset_warning = -1)
  expect_true(all(elderly_result$results$unit != "overall"))
  expect_true(all(elderly_result$results$n_patients > 0))
})

test_that("calc_te bootstrap analysis works correctly", {
  data <- setup_test_data()

  # Basic bootstrap test
  result <- calc_te(data,
                    bootstrap = TRUE,
                    bootstrap_params = list(
                      sample_percentage = 0.5,
                      n_iterations = 100,
                      distribution_span = 0.95
                    ),
                    min_loset_warning = -1)

  # Check bootstrap results structure
  expect_true(!is.null(result$bootstrap_distributions))
  expect_true(all(c("boot_ote_var_lower", "boot_ote_var_upper") %in%
                    names(result$results)))

  # Check confidence intervals are valid
  expect_true(all(result$results$boot_ote_var_lower <= result$results$ote_te))
  expect_true(all(result$results$boot_ote_var_upper >= result$results$ote_te))
})

test_that("calc_te handles convergence analysis correctly", {
  data <- setup_test_data()

  result <- calc_te(data,
                    bootstrap = TRUE,
                    bootstrap_params = list(
                      sample_percentage = 0.5,
                      n_iterations = 100,
                      distribution_span = 0.95
                    ),
                    check_convergence = TRUE,
                    min_loset_warning = -1)

  # Check convergence results structure
  expect_true(!is.null(result$convergence))
  expect_true(inherits(result$convergence, "te_convergence"))
  expect_named(result$convergence,
               c("groups", "summary", "metrics_analyzed", "calculation_time"))
})

test_that("calc_te handles invalid inputs appropriately", {
  data <- setup_test_data()

  # Test invalid var1
  expect_error(calc_te(data, var1 = "nonexistent_column"))

  # Test var2 without var1
  expect_error(calc_te(data, var2 = "unit"))

  # Test invalid subgroup
  expect_error(calc_te(data, subgroup = list(nonexistent = "value")))

  # Test incompatible subgroup and var1
  data$test_group <- "A"
  expect_error(calc_te(data,
                       subgroup = list(test_group = "A"),
                       var1 = "unit", min_loset_warning = -1))
})

test_that("calc_te calculates metrics correctly", {
  data <- setup_test_data()
  result <- calc_te(data, min_loset_warning = -1)

  # Test sensitivity calculation
  test_sensitivity <- function(unit_data) {
    tp <- sum(unit_data$tp)
    fn <- sum(unit_data$fn)
    return(tp / (tp + fn))
  }

  # Test specificity calculation
  test_specificity <- function(unit_data) {
    tn <- sum(unit_data$tn)
    fp <- sum(unit_data$fp)
    return(tn / (tn + fp))
  }

  # Verify TE calculation
  verify_te <- function(mean_all, mean_loset) {
    1 - (mean_loset / mean_all)
  }

  # Test for each unit (excluding overall)
  unit_results <- result$results %>%
    filter(unit != "overall")

  expect_equal(unit_results$ote_te,
               verify_te(unit_results$mean_all,
                         unit_results$ote_mean_loset))

  expect_equal(unit_results$tte_te,
               verify_te(unit_results$mean_all,
                         unit_results$tte_mean_loset))
})

test_that("calc_te overall_only parameter works correctly", {
  data <- setup_test_data()

  # Test with overall_only = TRUE
  result_overall <- calc_te(data, overall_only = TRUE, min_loset_warning = -1)
  expect_equal(nrow(result_overall$results), 1)
  expect_equal(result_overall$results$unit, "overall")

  # Compare with full results
  result_full <- calc_te(data)
  expect_equal(result_overall$results$ote_te,
               result_full$results$ote_te[result_full$results$unit == "overall"])
})

test_that("calc_te print method works correctly", {
  data <- setup_test_data()
  result <- calc_te(data, min_loset_warning = -1)

  # Capture print output
  output <- capture.output(print(result))

  # Check for key sections in output
  expect_true(any(grepl("Triage Effectiveness Analysis Results", output)))
  expect_true(any(grepl("Classification Metrics:", output)))
  expect_true(any(grepl("Triage Effectiveness Metrics:", output)))
  expect_true(any(grepl("Computation Information", output)))
})
