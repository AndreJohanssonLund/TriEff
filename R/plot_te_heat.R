#' Plot Triage Effectiveness Heatmap
#'
#' @description
#' Creates heatmap visualizations of Triage Effectiveness (TE) data. Can display either
#' a single heatmap or multiple heatmaps showing alternative calculations.
#'
#' @param data Tibble containing TE data with columns:
#'   \itemize{
#'     \item sensitivity: Sensitivity values
#'     \item specificity: Specificity values
#'     \item te: TE values
#'     \item Additional columns for alternative calculations if available
#'   }
#' @param show_labels Logical. If TRUE, displays numeric values on heatmap. Default TRUE.
#' @param show_alt_calc Logical. If TRUE, displays multiple heatmaps with alternative
#'   calculations. Default FALSE. Note that sim_heat_alt needs to be run to have the
#'   relevant variables for this.
#' @param normalize Logical. If TRUE, normalizes values to 0-100 scale. Ignored if
#'   show_alt_calc=TRUE. Default FALSE.
#' @param color_scheme Character. Either "triage", "oceansky", "blueorange" or "thermalsafe". The triage alternative is not color blind safe but is included since it mimics classic triage colors
#' @param title Character. Title for the plot. Default "Triage Effectiveness Heatmap".
#'
#' @return A ggplot object (single heatmap) or a list of patchwork objects (multiple heatmaps)
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradientn scale_fill_manual scale_fill_continuous scale_x_continuous scale_y_continuous labs theme_minimal theme element_blank element_text margin
#' @importFrom scales percent_format rescale
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom dplyr filter pull mutate case_when across where sym
#'
#' @export
plot_te_heatmap <- function(data,
                            show_labels = TRUE,
                            show_alt_calc = FALSE,
                            normalize = FALSE,
                            color_scheme = "oceansky",
                            title = "Triage Effectiveness Heatmap") {

  # Helper function to normalize data
  normalize_data <- function(data, metric = "te", target_100_100 = 1) {
    # Clean names from numeric vectors
    data <- dplyr::mutate(data, dplyr::across(dplyr::where(is.numeric), unname))

    # Get reference value for normalization
    max_value <- data %>%
      dplyr::filter(sensitivity == 1, specificity == 1) %>%
      dplyr::pull(!!rlang::sym(metric))

    # Normalize the specified metric
    data %>%
      dplyr::mutate(
        !!rlang::sym(metric) := dplyr::case_when(
          .data[[metric]] <= 0 ~ .data[[metric]] / abs(max_value) * target_100_100,
          TRUE ~ (.data[[metric]] / max_value) * target_100_100
        )
      )
  }

  # Helper function to create a single heatmap
  create_single_heatmap <- function(data, metric = "te",
                                    show_labels = TRUE,
                                    title = "Triage Efficacy Heatmap") {
    # Format values for display
    if(show_labels) {
      data$label <- dplyr::case_when(
        grepl("^te", metric) ~ sprintf("%.1f", data[[metric]] * 100),
        metric %in% c("cube_all", "log_all") ~ sprintf("%.2f", data[[metric]]),
        grepl("median_all", metric) ~ sprintf("%.0f", data[[metric]]),
        metric == "mean_all" ~ sprintf("%.0f", data[[metric]]),
        TRUE ~ sprintf("%.1f", data[[metric]])
      )
    }

    min_value <- min(data[[metric]], na.rm = TRUE)
    max_value <- max(data[[metric]], na.rm = TRUE)

    # Determine color scale based on data range
    # Validate color_scheme parameter
    valid_schemes <- c("triage", "oceansky", "blueorange", "thermalsafe")

    if (!color_scheme %in% valid_schemes) {
      warning(paste("Invalid color scheme '", color_scheme, "'. Using 'oceansky' instead. Valid options are: ",
                    paste(valid_schemes, collapse = ", ")))
      color_scheme <- "oceansky"
    }

    # Determine color scale based on data range and selected scheme
    if (color_scheme == "triage") {
      # Original triage scheme (not colorblind safe)
      if (min_value < -1) {
        colors <- c("red", "orange", "yellow", "green")
        values <- scales::rescale(c(min_value, -1, 0, max_value))
      } else if (min_value < 0) {
        colors <- c("orange", "yellow", "green")
        values <- scales::rescale(c(min_value, 0, max_value))
      } else {
        colors <- c("yellow", "green")
        values <- scales::rescale(c(0, max_value))
      }
    } else if (color_scheme == "oceansky") {
      # OceanSky scheme - completely colorblind safe
      if (min_value < -1) {
        colors <- c("navy", "lightblue", "lightyellow", "gold")
        values <- scales::rescale(c(min_value, -1, 0, max_value))
      } else if (min_value < 0) {
        colors <- c("lightblue", "lightyellow", "gold")
        values <- scales::rescale(c(min_value, 0, max_value))
      } else {
        colors <- c("lightyellow", "gold")
        values <- scales::rescale(c(0, max_value))
      }
    } else if (color_scheme == "blueorange") {
      # BlueOrange scheme - colorblind safe
      if (min_value < -1) {
        colors <- c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c")
        values <- scales::rescale(c(min_value, -1, 0, max_value))
      } else if (min_value < 0) {
        colors <- c("#abd9e9", "#fdae61", "#d7191c")
        values <- scales::rescale(c(min_value, 0, max_value))
      } else {
        colors <- c("#fdae61", "#d7191c")
        values <- scales::rescale(c(0, max_value))
      }
    } else if (color_scheme == "thermalsafe") {
      # ThermalSafe scheme - colorblind safe
      if (min_value < -1) {
        colors <- c("purple", "lightblue", "yellow", "orange")
        values <- scales::rescale(c(min_value, -1, 0, max_value))
      } else if (min_value < 0) {
        colors <- c("lightblue", "yellow", "orange")
        values <- scales::rescale(c(min_value, 0, max_value))
      } else {
        colors <- c("yellow", "orange")
        values <- scales::rescale(c(0, max_value))
      }
    }


    # Create the base plot
    p <- ggplot2::ggplot(data, ggplot2::aes(x = specificity, y = sensitivity,
                                            fill = .data[[metric]])) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradientn(  # Add this line
        colors = colors,
        values = values,
        limits = c(min_value, max_value)
      )

    # Add labels if requested
    if(show_labels) {
      p <- p + ggplot2::geom_text(ggplot2::aes(label = label), size = 2.5)
    }


    # Add scales and theme
    p <- p +
      ggplot2::scale_x_continuous(
        expand = c(0, 0),
        breaks = seq(0, 1, 0.1),
        labels = scales::percent_format(accuracy = 1)
      ) +
      ggplot2::scale_y_continuous(
        expand = c(0, 0),
        breaks = seq(0, 1, 0.1),
        labels = scales::percent_format(accuracy = 1)
      ) +
      ggplot2::labs(
        title = title,
        x = "Specificity",
        y = "Sensitivity",
        fill = "TE"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, margin = margin(r = 10)),
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
      )

    return(p)
  }

  # If showing alternative calculations
  if(show_alt_calc) {
    # Check if alternative calculations exist
    alt_metrics <- c("te_cube", "te_log", "te_median_local", "te_median_global")
    if(!any(alt_metrics %in% names(data))) {
      stop("Alternative calculations not found in data. Run sim_heat_alt with alt_calc=TRUE.")
    }

    # Create TE metrics batch
    te_plots <- list()
    te_plots[["te"]] <- create_single_heatmap(data, "te", show_labels, "Original TE")

    if("te_cube" %in% names(data))
      te_plots[["te_cube"]] <- create_single_heatmap(data, "te_cube", show_labels,
                                                     "Cube Root Normalized TE")
    if("te_log" %in% names(data))
      te_plots[["te_log"]] <- create_single_heatmap(data, "te_log", show_labels,
                                                    "Log Normalized TE")
    if("te_median_local" %in% names(data))
      te_plots[["te_median_local"]] <- create_single_heatmap(data, "te_median_local",
                                                             show_labels,
                                                             "Median-based TE (Local)")
    if("te_median_global" %in% names(data))
      te_plots[["te_median_global"]] <- create_single_heatmap(data, "te_median_global",
                                                              show_labels,
                                                              "Median-based TE (Global)")

    te_batch <- patchwork::wrap_plots(te_plots, ncol = 3) +
      patchwork::plot_annotation(title = "TE Metrics",
                                 theme = ggplot2::theme(
                                   plot.title = element_text(size = 14, hjust = 0.5)
                                 ))

    # Create normalized TE batch
    norm_data <- normalize_data(data)
    norm_plots <- list()

    norm_plots[["te"]] <- create_single_heatmap(norm_data, "te", show_labels,
                                                "Normalized Original TE")

    if("te_cube" %in% names(data)) {
      norm_data_cube <- normalize_data(data, "te_cube")
      norm_plots[["te_cube"]] <- create_single_heatmap(norm_data_cube,
                                                       "te_cube", show_labels,
                                                       "Normalized Cube Root TE")
    }

    if("te_log" %in% names(data)) {
      norm_data_log <- normalize_data(data, "te_log")
      norm_plots[["te_log"]] <- create_single_heatmap(norm_data_log,
                                                      "te_log", show_labels,
                                                      "Normalized Log TE")
    }

    if("te_median_local" %in% names(data)) {
      norm_data_local <- normalize_data(data, "te_median_local")
      norm_plots[["te_median_local"]] <- create_single_heatmap(
        norm_data_local, "te_median_local", show_labels,
        "Normalized Median-based TE (Local)")
    }

    if("te_median_global" %in% names(data)) {
      norm_data_global <- normalize_data(data, "te_median_global")
      norm_plots[["te_median_global"]] <- create_single_heatmap(
        norm_data_global, "te_median_global", show_labels,
        "Normalized Median-based TE (Global)")
    }

    norm_batch <- patchwork::wrap_plots(norm_plots, ncol = 3) +
      patchwork::plot_annotation(
        title = "Normalized TE Metrics",
        theme = ggplot2::theme(plot.title = element_text(size = 14, hjust = 0.5))
      )

    # Create all metrics batch
    all_plots <- list()

    if("mean_all" %in% names(data))
      all_plots[["mean_all"]] <- create_single_heatmap(data, "mean_all", show_labels,
                                                       "Mean Wait Time (Original)")

    if("cube_all" %in% names(data))
      all_plots[["cube_all"]] <- create_single_heatmap(data, "cube_all", show_labels,
                                                       "Mean Wait Time (Cube Root)")

    if("log_all" %in% names(data))
      all_plots[["log_all"]] <- create_single_heatmap(data, "log_all", show_labels,
                                                      "Mean Wait Time (Log)")

    if("median_all_local" %in% names(data))
      all_plots[["median_all_local"]] <- create_single_heatmap(data, "median_all_local",
                                                               show_labels,
                                                               "Median Wait Time (Local)")

    if("median_all_global" %in% names(data))
      all_plots[["median_all_global"]] <- create_single_heatmap(data, "median_all_global",
                                                                show_labels,
                                                                "Median Wait Time (Global)")

    all_batch <- patchwork::wrap_plots(all_plots, ncol = 3) +
      patchwork::plot_annotation(
        title = "All Patients Wait Time Metrics",
        theme = ggplot2::theme(plot.title = element_text(size = 14, hjust = 0.5))
      )

    return(list(
      te_batch = te_batch,
      all_batch = all_batch,
      normalized_te_batch = norm_batch
    ))
  } else {
    # Single heatmap case
    if(normalize) {
      data <- normalize_data(data)
    }

    return(create_single_heatmap(data, "te", show_labels, title))
  }
}
