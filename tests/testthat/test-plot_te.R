# Helper functions for tests
get_test_data <- function() {
  data <- load_sem_synth() %>%
    dplyr::filter(arrival <= min(arrival) + lubridate::days(7))
  return(data)
}

setup_test_data <- function() {
  data <- get_test_data()
  data <- init(data) %>%
    sim_te()
  return(data)
}

test_that("plot_te creates valid ggplot object", {
  data <- setup_test_data()
  result <- calc_wte(data, min_loset_warning = -1)

  p <- plot_te(result)

  # Check basic plot structure
  expect_s3_class(p, "ggplot")

  # Get layers
  layers <- p$layers

  # Check for essential geoms
  geom_classes <- sapply(layers, function(l) class(l$geom)[1])

  # Should have points for TTE and OTE
  expect_true("GeomPoint" %in% geom_classes)

  # Should have text for labels
  expect_true("GeomText" %in% geom_classes)
})

test_that("plot_te handles bootstrap results correctly", {
  data <- setup_test_data()
  result <- calc_wte(data, bootstrap = TRUE,
                    bootstrap_params = list(
                      sample_percentage = 0.5,
                      n_iterations = 100,
                      distribution_span = 0.95
                    ),
                    min_loset_warning = -1)

  p <- plot_te(result)

  # Get layer classes
  layer_classes <- sapply(p$layers, function(l) class(l$geom)[1])

  # Should have segments for confidence intervals
  expect_true("GeomSegment" %in% layer_classes)
})

test_that("plot_te respects show parameters", {
  data <- setup_test_data()
  result <- calc_wte(data, min_loset_warning = -1)

  # Test hiding TTE
  p1 <- plot_te(result, show_tte = FALSE)
  tte_layers <- sum(sapply(p1$layers, function(l) {
    inherits(l$geom, "GeomPoint") &&
      !is.null(l$aes_params$colour) &&
      l$aes_params$colour == "#00A3A3"
  }))
  expect_equal(tte_layers, 0)

  # Test hiding OTE
  p2 <- plot_te(result, show_ote = FALSE)
  ote_layers <- sum(sapply(p2$layers, function(l) {
    inherits(l$geom, "GeomPoint") &&
      !is.null(l$aes_params$colour) &&
      l$aes_params$colour == "#7D3C98"
  }))
  expect_equal(ote_layers, 0)
})

test_that("plot_te handles label styles correctly", {
  data <- setup_test_data()
  result <- calc_wte(data, min_loset_warning = -1)

  # Test different label styles
  styles <- c("none", "small", "full", "half")

  for(style in styles) {
    p <- plot_te(result, label_style = style)
    expect_s3_class(p, "ggplot")

    if(style == "none") {
      # Should not have additional text layers for legend
      text_layers <- sum(sapply(p$layers, function(l) {
        inherits(l$geom, "GeomText") &&
          !is.null(l$data)
      }))
      expect_true(text_layers <= nrow(result$results) * 2)  # Only value labels
    }
  }

  # Test invalid label style
  expect_error(plot_te(result, label_style = "invalid"))
})

test_that("plot_te handles negative values correctly", {
  data <- setup_test_data()
  result <- calc_wte(data, min_loset_warning = -1)

  # Create test results with negative TE
  result$results$ote_te <- -0.2
  result$results$tte_te <- -0.3

  # Plot with default settings
  p1 <- plot_te(result)

  # Plot with custom min_x
  p2 <- plot_te(result, min_x = -0.5)

  # Check both plots are created successfully
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("plot_te handles missing data appropriately", {
  data <- setup_test_data()
  result <- calc_wte(data,min_loset_warning = -1)

  # Test with missing TTE
  result_no_tte <- result
  result_no_tte$results$tte_te <- NA
  p1 <- plot_te(result_no_tte)
  expect_s3_class(p1, "ggplot")

  # Test with missing bootstrap intervals
  result_no_boot <- result
  result_no_boot$results <- result_no_boot$results %>%
    dplyr::select(-matches("boot_"))
  p2 <- plot_te(result_no_boot)
  expect_s3_class(p2, "ggplot")
})

test_that("plot_te handles aesthetic parameters correctly", {
  data <- setup_test_data()
  result <- calc_wte(data, min_loset_warning = -1)

  # Test aesthetic parameters work
  p1 <- plot_te(result, var_alpha = 0.8)
  p2 <- plot_te(result, dumbell_width = 2)
  p3 <- plot_te(result, var_width = 2)

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
})

test_that("plot_te handles subgroup analysis correctly", {
  data <- setup_test_data()

  # Add age groups
  data$age_group <- cut(data$age_at_arrival,
                        breaks = c(0, 59, 79, Inf),
                        labels = c("18-59", "60-79", "80+"))

  result <- calc_wte(data, var1 = "age_group", min_loset_warning = -1)

  p <- plot_te(result)

  # Check plot is created successfully
  expect_s3_class(p, "ggplot")

  # Should have necessary layers for points and text
  expect_true(any(sapply(p$layers, function(l) inherits(l$geom, "GeomPoint"))))
  expect_true(any(sapply(p$layers, function(l) inherits(l$geom, "GeomText"))))
})
