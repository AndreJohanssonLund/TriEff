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
#'     \item calc_method: (Optional) Indicates whether data uses "wte" or "rte" calculations
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
#' @param axis_step_size Numeric or NULL. Override for the step size of axis labels (e.g., 0.1 for 10% steps).
#'   If NULL (default), automatically determines optimal step size based on data resolution,
#'   with a minimum of 0.05 (5%). Automatic detection selects 25% steps for data with 25%
#'   resolution, 10% steps for data with 10% resolution, and 5% steps for finer resolutions.
#'
#' @return A ggplot object (single heatmap) or a list of patchwork objects (multiple heatmaps)
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradientn scale_fill_manual scale_fill_continuous scale_x_continuous scale_y_continuous labs theme_minimal theme element_blank element_text element_line margin
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
                            title = "Triage Effectiveness Heatmap",
                            axis_step_size = NULL) {

  # Helper function to determine optimal axis step size
  determine_optimal_step_size <- function(data) {
    # Get unique values of sensitivity
    sens_values <- sort(unique(data$sensitivity))

    # Determine step size based on number of unique values
    n_unique <- length(sens_values)

    if (n_unique <= 5) {
      return(0.25)  # 25% steps (0, 0.25, 0.5, 0.75, 1)
    } else if (n_unique <= 11) {
      return(0.1)   # 10% steps (0, 0.1, 0.2, ..., 1)
    } else if (n_unique <= 21) {
      return(0.05)  # 5% steps (0, 0.05, 0.1, ..., 1)
    } else {
      return(0.05)  # Default to 5% steps for very detailed data
    }
  }

  # Then extract the results data
  results_data <- if (is.data.frame(data)) {
    data  # If data is directly a data frame
  } else if (!is.null(data$results) && is.data.frame(data$results)) {
    data$results  # If data is a list with a results element
  } else {
    stop("Data format not recognized. Expected a data frame or a list with a 'results' element.")
  }

  # Determine calculation method (wte or rte)
  calc_method <- if ("calc_method" %in% names(results_data)) {
    # Get the first value (assuming consistent method across all rows)
    results_data$calc_method[1]
  } else {
    # Default to "wte" for backward compatibility if not specified
    "wte"
  }

  # Determine axis step size if not provided
  if (is.null(axis_step_size)) {
    axis_step_size <- determine_optimal_step_size(results_data)
  }

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
                                    title = "Triage Efficacy Heatmap",
                                    axis_step_size = 0.1) {
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

    # Validate color_scheme parameter
    valid_schemes <- c("triage", "oceansky", "blueorange", "thermalsafe")

    if (!color_scheme %in% valid_schemes) {
      warning(paste("Invalid color scheme '", color_scheme, "'. Using 'oceansky' instead. Valid options are: ",
                    paste(valid_schemes, collapse = ", ")))
      color_scheme <- "oceansky"
    }

    # Define color palettes
    color_palettes <- list(
      oceansky = list(
        extreme_neg = "#2c7bb6",
        neg = "lightblue",
        neutral = "white",
        pos = "gold"
      ),
      thermalsafe = list(
        extreme_neg = "purple",
        neg = "lightblue",
        neutral = "white",
        pos = "orange"
      ),
      blueorange = list(
        extreme_neg = "#2c7bb6",
        neg = "#abd9e9",
        neutral = "white",
        pos = "#d7191c"
      ),
      triage = list(
        extreme_neg = "red",
        neg = "orange",
        neutral = "yellow",
        pos = "green"
      )
    )

    palette <- color_palettes[[color_scheme]]


    # Determine color configuration based on data range
    if (min_value < -1 && max_value > 0) {
      colors <- c(palette$extreme_neg, palette$neg, palette$neutral, palette$pos)
      values <- scales::rescale(c(min_value, -1, 0, max_value))
    } else if (min_value < -1 && max_value > 0) {
      colors <- c(palette$extreme_neg, palette$neg, palette$neutral)
      values <- scales::rescale(c(min_value, -1, max_value))
    } else if (min_value < 0 && max_value > 0) {
      colors <- c(palette$neg, palette$neutral, palette$pos)
      values <- scales::rescale(c(min_value, 0, max_value))
    } else if (min_value < -1) {
      colors <- c(palette$extreme_neg, palette$neg)
      values <- scales::rescale(c(min_value, -1))
    } else if (min_value < 0) {
      colors <- c(palette$neg, palette$neutral)
      values <- scales::rescale(c(min_value, 0))
    } else if (max_value > 0) {
      colors <- c(palette$neutral, palette$pos)
      values <- scales::rescale(c(0, max_value))
    } else {
      colors <- c(palette$neutral, palette$pos)
      values <- scales::rescale(c(0, max_value))
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


    # Add scales and theme with dynamic axis step size
    p <- p +
      ggplot2::scale_x_continuous(
        expand = c(0, 0),
        breaks = seq(0, 1, axis_step_size),
        labels = scales::percent_format(accuracy = 1)
      ) +
      ggplot2::scale_y_continuous(
        expand = c(0, 0),
        breaks = seq(0, 1, axis_step_size),
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
    # Check if any alternative calculations exist based on calculation method
    if (calc_method == "wte") {
      alt_metrics <- c("te_cube", "te_log", "te_median_local", "te_median_global")
      raw_metrics <- c("mean_all", "cube_all", "log_all", "median_all_local", "median_all_global")
      method_title <- "Waiting Time-based TE"
    } else if (calc_method == "rte") {
      alt_metrics <- c("te_cube", "te_log", "te_median", "te_median_cube", "te_median_log")
      raw_metrics <- NULL
      method_title <- "Rank-based TE"
    } else {
      stop("Unknown calculation method: ", calc_method)
    }

    # Check for existence of alternative metrics in the data
    available_alt_metrics <- intersect(alt_metrics, names(results_data))

    if(length(available_alt_metrics) == 0) {
      stop("Alternative calculations not found in data. Run sim_heat_alt with alt_calc=TRUE.")
    }

    # Create TE metrics batch
    te_plots <- list()

    # Base TE plot always included
    method_label <- ifelse(calc_method == "wte", "WTE", "RTE")
    te_plots[["te"]] <- create_single_heatmap(results_data, "te", show_labels,
                                              paste("Original", method_label), axis_step_size)

    # Add alternative calculations based on calculation method and what's available
    if("te_cube" %in% names(results_data))
      te_plots[["te_cube"]] <- create_single_heatmap(results_data, "te_cube", show_labels,
                                                     paste("Cube Root", method_label), axis_step_size)
    # For WTE log transformation (standard log)
    if("te_log" %in% names(results_data))
      te_plots[["te_log"]] <- create_single_heatmap(results_data, "te_log", show_labels,
                                                    paste("Log", method_label), axis_step_size)

    # For RTE signed log transformation
    if("te_signed_log" %in% names(results_data))
      te_plots[["te_signed_log"]] <- create_single_heatmap(results_data, "te_signed_log", show_labels,
                                                           paste("Signed Log", method_label), axis_step_size)
    # Add WTE-specific metrics
    if (calc_method == "wte") {
      if("te_median_local" %in% names(results_data))
        te_plots[["te_median_local"]] <- create_single_heatmap(results_data, "te_median_local",
                                                               show_labels,
                                                               "Median-based WTE (Local)", axis_step_size)
      if("te_median_global" %in% names(results_data))
        te_plots[["te_median_global"]] <- create_single_heatmap(results_data, "te_median_global",
                                                                show_labels,
                                                                "Median-based WTE (Global)", axis_step_size)
    }
    # Add RTE-specific median metrics
    else if (calc_method == "rte") {
      if("te_median" %in% names(results_data))
        te_plots[["te_median"]] <- create_single_heatmap(results_data, "te_median",
                                                         show_labels,
                                                         "Median-based RTE", axis_step_size)
      if("te_median_cube" %in% names(results_data))
        te_plots[["te_median_cube"]] <- create_single_heatmap(results_data, "te_median_cube",
                                                              show_labels,
                                                              "Median Cube Root RTE", axis_step_size)
      if("te_median_log" %in% names(results_data))
        te_plots[["te_median_log"]] <- create_single_heatmap(results_data, "te_median_log",
                                                             show_labels,
                                                             "Median Signed Log RTE", axis_step_size)
    }

    # Organize plots in a 2x2 or 2x3 grid as needed
    num_plots <- length(te_plots)
    if (num_plots <= 4) {
      # Use 2x2 grid for 4 or fewer plots
      ncols <- 2
    } else {
      # Use 2x3 grid for more than 4 plots
      ncols <- 3
    }

    te_batch <- patchwork::wrap_plots(te_plots, ncol = ncols) +
      patchwork::plot_annotation(title = paste(method_title, "Metrics"),
                                 theme = ggplot2::theme(
                                   plot.title = element_text(size = 14, hjust = 0.5)
                                 ))

    # Create normalized TE batch
    norm_data <- normalize_data(results_data)
    norm_plots <- list()

    norm_plots[["te"]] <- create_single_heatmap(norm_data, "te", show_labels,
                                                paste("Normalized Original", method_label), axis_step_size)

    if("te_cube" %in% names(results_data)) {
      norm_data_cube <- normalize_data(results_data, "te_cube")
      norm_plots[["te_cube"]] <- create_single_heatmap(norm_data_cube,
                                                       "te_cube", show_labels,
                                                       paste("Normalized Cube Root", method_label), axis_step_size)
    }

    if("te_log" %in% names(results_data)) {
      norm_data_log <- normalize_data(results_data, "te_log")
      norm_plots[["te_log"]] <- create_single_heatmap(norm_data_log,
                                                      "te_log", show_labels,
                                                      paste("Normalized Signed Log", method_label), axis_step_size)
    }

    # Add WTE-specific normalized metrics
    if (calc_method == "wte") {
      if("te_median_local" %in% names(results_data)) {
        norm_data_local <- normalize_data(results_data, "te_median_local")
        norm_plots[["te_median_local"]] <- create_single_heatmap(
          norm_data_local, "te_median_local", show_labels,
          "Normalized Median-based WTE (Local)", axis_step_size)
      }

      if("te_median_global" %in% names(results_data)) {
        norm_data_global <- normalize_data(results_data, "te_median_global")
        norm_plots[["te_median_global"]] <- create_single_heatmap(
          norm_data_global, "te_median_global", show_labels,
          "Normalized Median-based WTE (Global)", axis_step_size)
      }
    }
    # Add RTE-specific normalized median metrics
    else if (calc_method == "rte") {
      if("te_median" %in% names(results_data)) {
        norm_data_median <- normalize_data(results_data, "te_median")
        norm_plots[["te_median"]] <- create_single_heatmap(
          norm_data_median, "te_median", show_labels,
          "Normalized Median-based RTE", axis_step_size)
      }

      if("te_median_cube" %in% names(results_data)) {
        norm_data_median_cube <- normalize_data(results_data, "te_median_cube")
        norm_plots[["te_median_cube"]] <- create_single_heatmap(
          norm_data_median_cube, "te_median_cube", show_labels,
          "Normalized Median Cube Root RTE", axis_step_size)
      }

      if("te_median_log" %in% names(results_data)) {
        norm_data_median_log <- normalize_data(results_data, "te_median_log")
        norm_plots[["te_median_log"]] <- create_single_heatmap(
          norm_data_median_log, "te_median_log", show_labels,
          "Normalized Median Signed Log RTE", axis_step_size)
      }
    }

    # Organize normalized plots in a 2x2 or 2x3 grid as needed
    num_norm_plots <- length(norm_plots)
    if (num_norm_plots <= 4) {
      # Use 2x2 grid for 4 or fewer plots
      norm_ncols <- 2
    } else {
      # Use 2x3 grid for more than 4 plots
      norm_ncols <- 3
    }

    norm_batch <- patchwork::wrap_plots(norm_plots, ncol = norm_ncols) +
      patchwork::plot_annotation(
        title = paste("Normalized", method_title, "Metrics"),
        theme = ggplot2::theme(plot.title = element_text(size = 14, hjust = 0.5))
      )

    # Create result list with appropriate components
    result_list <- list(
      te_batch = te_batch,
      normalized_te_batch = norm_batch
    )

    # Only include all_batch for WTE calculations
    if (calc_method == "wte" && any(raw_metrics %in% names(results_data))) {
      all_plots <- list()

      if("mean_all" %in% names(results_data))
        all_plots[["mean_all"]] <- create_single_heatmap(results_data, "mean_all", show_labels,
                                                         "Mean Wait Time (Original)", axis_step_size)

      if("cube_all" %in% names(results_data))
        all_plots[["cube_all"]] <- create_single_heatmap(results_data, "cube_all", show_labels,
                                                         "Mean Wait Time (Cube Root)", axis_step_size)

      if("log_all" %in% names(results_data))
        all_plots[["log_all"]] <- create_single_heatmap(results_data, "log_all", show_labels,
                                                        "Mean Wait Time (Log)", axis_step_size)

      if("median_all_local" %in% names(results_data))
        all_plots[["median_all_local"]] <- create_single_heatmap(results_data, "median_all_local",
                                                                 show_labels,
                                                                 "Median Wait Time (Local)", axis_step_size)

      if("median_all_global" %in% names(results_data))
        all_plots[["median_all_global"]] <- create_single_heatmap(results_data, "median_all_global",
                                                                  show_labels,
                                                                  "Median Wait Time (Global)", axis_step_size)

      all_batch <- patchwork::wrap_plots(all_plots, ncol = 3) +
        patchwork::plot_annotation(
          title = "All Patients Wait Time Metrics",
          theme = ggplot2::theme(plot.title = element_text(size = 14, hjust = 0.5))
        )

      result_list$all_batch <- all_batch
    }

    return(result_list)
  } else {
    # Single heatmap case
    if(normalize) {
      results_data <- normalize_data(results_data)
    }

    # Modify title to reflect calculation method if not custom
    if (title == "Triage Effectiveness Heatmap") {
      if (calc_method == "wte") {
        title <- "Waiting Time-based Triage Effectiveness"
      } else if (calc_method == "rte") {
        title <- "Rank-based Triage Effectiveness"
      }
    }

    return(create_single_heatmap(results_data, "te", show_labels, title, axis_step_size))
  }
}


