#' Create a Triage Effectiveness Plot
#'
#' This function generates a plot to visualize Triage Effectiveness (TE) metrics.
#' It can display Theoretical Triage Effectiveness (TTE) and Observed Triage Effectiveness (OTE)
#' for different units, with optional grouping variables and variability intervals.
#'
#' @param data A data frame containing the following columns:
#'   \itemize{
#'     \item unit: Factor or character vector of unit names
#'     \item tte_te: Numeric vector of Theoretical Triage Effectiveness values
#'     \item ote_te: Numeric vector of Observed Triage Effectiveness values
#'     \item boot_tte_var_lower: Optional. Lower variability value for TTE
#'     \item boot_tte_var_upper: Optional. Upper variability value for TTE
#'     \item boot_ote_var_lower: Optional. Lower variability value for OTE
#'     \item boot_ote_var_upper: Optional. Upper variability value for OTE
#'   }
#' @param title Character string for the plot title. Default is "Triage Effectiveness Plot".
#' @param show_tte Logical, whether to display TTE values. Default is NULL. Will show if there is TTE data.
#' @param show_ote Logical, whether to display OTE values. Default is NULL. Will show if there is OTE data.
#' @param show_var Logical, whether to display variability intervals. Default is null. Will show if there is variability data supplied.
#' @param show_percent_tte Logical, whether to display the percentage of TE above the TTE dot or not.
#' @param show_percent_ote Logical, whether to display the percentage of TE above the OTE dor or not.
#' @param var_alpha Numeric between 0 and 1, specifying the transparency of variability interval lines. Default is 0.5.
#' @param dumbell_width Numeric, specifying the width of the line connecting TTE and OTE points. Default is 1. Note that since the line is horizontal, the "width" represents its thickness.
#' @param var_width Numeric, specifying the width of the main variability interval lines - that is the vertical line at the end of the variability span. Default is 1.5.
#' @param var_end_width Numeric, specifying the width (or thickness) of the variability interval end lines. Default is 1.
#' @param var_end_height Numeric, specifying the height of variability interval end lines as offset from center. Default is 0.05.
#' @param label_style categorical, valid inputs - "none", "small", "full" & "half". Decides the style of the label for TTE/OTE on the rights side of the graph. Default is half witch nets the label style as: "Theoretical TE". Small: "TTE", Full: "Theoretical Triage Effectiveness".
#' @param min_x Percentual value that is the minimum value that the plot will stretch to. Default is null which uses the min value in the data rounded to the next 10%. Useful to set if you have several plots with one that has negative values to keep the same aspect ratio.
#' @param ci_vars_override Character string. Explicitly specify which CI variables to use:
#' "parametric", "nonparametric", or "boot". If NULL (default), uses recommended variables
#' from metadata based on calculation method and normality tests.
#' @return A ggplot object representing the Triage Effectiveness plot.
#'
#' @importFrom ggplot2 ggplot aes geom_vline geom_segment geom_point geom_text scale_x_continuous scale_y_continuous scale_color_manual coord_cartesian labs theme_minimal theme element_blank element_text element_line margin
#' @importFrom scales percent_format
#' @importFrom dplyr arrange desc mutate row_number
#' @importFrom grDevices adjustcolor
#' @importFrom rlang sym !!
#'
#' @examples
#' \dontrun{
#' data <- data.frame(
#'   unit = c("Unit A", "Unit B", "Unit C"),
#'   tte_te = c(0.6, 0.7, 0.5),
#'   ote_te = c(0.65, 0.75, 0.55),
#'   boot_tte_var_lower = c(0.55, 0.65, 0.45),
#'   boot_tte_var_upper = c(0.65, 0.75, 0.55),
#'   boot_ote_var_lower = c(0.60, 0.70, 0.50),
#'   boot_ote_var_upper = c(0.70, 0.80, 0.60)
#' )
#' plot_te(data, show_var = TRUE)
#' }
#' @export
plot_te <- function(data, title = "Triage Effectiveness Plot",
                    show_tte = NULL,
                    show_ote = NULL,
                    show_var = NULL,
                    show_percent_tte = TRUE,
                    show_percent_ote = TRUE,
                    var_alpha = 0.5,
                    dumbell_width = 1,
                    var_width = 1.5,
                    var_end_width = 1,
                    var_end_height = 0.1,
                    label_style = "half",
                    min_x = NULL,
                    ci_vars_override = NULL) {

  # First extract the metadata
  metadata <- if (is.list(data) && !is.null(data$metadata)) {
    data$metadata
  } else {
    list()  # Empty list if no metadata available
  }

  # Then extract the results data
  results_data <- if (is.data.frame(data)) {
    data  # If data is directly a data frame
  } else if (!is.null(data$results) && is.data.frame(data$results)) {
    data$results  # If data is a list with a results element
  } else {
    stop("Data format not recognized. Expected a data frame or a list with a 'results' element.")
  }

  # Validate label_style
  if (!label_style %in% c("none", "small", "full", "half")) {
    stop("label_style must be one of: 'none', 'small', 'full', 'half'")
  }

  # Check for metrics availability
  has_tte <- "tte_te" %in% names(results_data) && !all(is.na(results_data$tte_te))
  has_ote <- "ote_te" %in% names(results_data) && !all(is.na(results_data$ote_te))

  # Determine which metrics to show
  show_tte <- if (is.null(show_tte)) has_tte else show_tte
  show_ote <- if (is.null(show_ote)) has_ote else show_ote

  # Stop if no metrics will be shown
  if (!show_tte && !show_ote) {
    stop("No metrics available to display. Function requires either TTE or OTE data.")
  }

  # Determine which CI variables to use
  ci_var_prefix <- if (!is.null(ci_vars_override)) {
    ci_vars_override
  } else if (!is.null(metadata$recommended_ci_vars) && metadata$recommended_ci_vars != "none") {
    metadata$recommended_ci_vars
  } else if (!is.null(metadata$calculation_method)) {
    if (metadata$calculation_method == "rank-based") {
      "parametric" # Default for RTE
    } else {
      "boot" # Default for WTE
    }
  } else {
    "boot" # Fallback default
  }

  # Determine if we should show variability at all
  effective_show_var <- if (is.null(show_var)) {
    if (is.null(metadata$recommended_ci_vars) || metadata$recommended_ci_vars == "none") {
      FALSE # Don't show variability if no CI recommended
    } else {
      TRUE  # Show variability by default if CIs are available
    }
  } else {
    show_var # Use user's explicit preference
  }

  # Generate actual variable names for CI bounds
  tte_lower_var <- paste0(ci_var_prefix, "_tte_var_lower")
  tte_upper_var <- paste0(ci_var_prefix, "_tte_var_upper")
  ote_lower_var <- paste0(ci_var_prefix, "_ote_var_lower")
  ote_upper_var <- paste0(ci_var_prefix, "_ote_var_upper")

  # Check for variability data availability
  can_show_ote_var <- all(c(ote_lower_var, ote_upper_var) %in% names(results_data))
  can_show_tte_var <- all(c(tte_lower_var, tte_upper_var) %in% names(results_data))

  group_var1 = metadata$group_var1
  group_var2 = metadata$group_var2

  # Colors for consistent use throughout
  tte_color <- "#00A3A3"
  ote_color <- "#7D3C98"

  # Set effective show_var based on availability and user preference
  show_ote_var <- can_show_ote_var && effective_show_var
  show_tte_var <- can_show_tte_var && effective_show_var

  # Collect all relevant values for x-axis limits
  x_values <- c(0, 1)
  if (show_ote) x_values <- c(x_values, results_data$ote_te)
  if (show_tte) x_values <- c(x_values, results_data$tte_te)

  if (show_ote_var) {
    x_values <- c(x_values, results_data[[ote_lower_var]], results_data[[ote_upper_var]])
  }

  if (show_tte_var) {
    x_values <- c(x_values, results_data[[tte_lower_var]], results_data[[tte_upper_var]])
  }

  # Calculate actual minimum x value
  actual_min_x <- min(0, x_values, na.rm = TRUE)
  if (actual_min_x < 0) {
    actual_min_x <- floor(actual_min_x * 10) / 10  # Round down to nearest 0.1
  }

  # Handle provided min_x
  if (!is.null(min_x)) {
    if (actual_min_x < min_x) {
      warning(sprintf("Plot contains values (%.2f) below specified min_x (%.2f)",
                      actual_min_x, min_x))
    }
    x_min <- min_x
  } else {
    x_min <- actual_min_x
  }

  # Add y position for each unit
  results_data <- results_data %>%
    arrange(desc(unit)) %>%
    mutate(y_pos = row_number())

  # Handle group labels outside of mutate
  if (!is.null(group_var1) && !is.null(group_var2)) {
    results_data$group_label <- paste(results_data[[group_var1]], results_data[[group_var2]], sep = " / ")
  } else if (!is.null(group_var1)) {
    results_data$group_label <- results_data[[group_var1]]
  } else {
    results_data$group_label <- ""
  }

  # Create the plot
  p <- ggplot(results_data, aes(y = y_pos)) +
    geom_vline(xintercept = seq(x_min, 1, 0.1), color = "lightgray", linewidth = 0.5) +
    geom_vline(xintercept = 1, color = "darkgray", linewidth = 1) +
    geom_vline(xintercept = 0, color = "darkgray", linewidth = 1)

  # Add dumbell plot if both TTE and OTE are shown
  if(show_tte && show_ote) {
    p <- p + geom_segment(aes(x = tte_te, xend = ote_te, yend = y_pos),
                          linewidth = dumbell_width)
  }

  # Add TTE points and variability if requested
  if (show_tte) {
    p <- p + geom_point(aes(x = tte_te), size = 3, color = tte_color)

    if(show_percent_tte) {
      p <- p + geom_text(aes(x = tte_te, label = scales::percent(tte_te, accuracy = 0.1)),
                 hjust = 0.5, vjust = -0.7, color = tte_color, size = 3)
    }
  }

  if (show_tte_var & show_tte) {
    p <- p +
      # Main CI line
      geom_segment(aes(x = .data[[tte_lower_var]], xend = .data[[tte_upper_var]],
                       y = y_pos, yend = y_pos),
                   color = adjustcolor(tte_color, alpha.f = var_alpha),
                   linewidth = var_width) +
      # Left endline
      geom_segment(aes(x = .data[[tte_lower_var]], xend = .data[[tte_lower_var]],
                       y = y_pos - var_end_height, yend = y_pos + var_end_height),
                   color = adjustcolor(tte_color, alpha.f = var_alpha),
                   linewidth = var_end_width) +
      # Right endline
      geom_segment(aes(x = .data[[tte_upper_var]], xend = .data[[tte_upper_var]],
                       y = y_pos - var_end_height, yend = y_pos + var_end_height),
                   color = adjustcolor(tte_color, alpha.f = var_alpha),
                   linewidth = var_end_width)
  }

  # Add OTE points and variability if requested
  if (show_ote) {
    p <- p + geom_point(aes(x = ote_te), size = 3, color = ote_color)

    if(show_percent_ote) {
      p <- p + geom_text(aes(x = ote_te, label = scales::percent(ote_te, accuracy = 0.1)),
                hjust = 0.5, vjust = -0.7, color = ote_color, size = 3)
    }
  }

  if (show_ote_var & show_ote) {
    p <- p +
      # Main CI line
      geom_segment(aes(x = .data[[ote_lower_var]], xend = .data[[ote_upper_var]],
                       y = y_pos, yend = y_pos),
                   color = adjustcolor(ote_color, alpha.f = var_alpha),
                   linewidth = var_width) +
      # Left endline
      geom_segment(aes(x = .data[[ote_lower_var]], xend = .data[[ote_lower_var]],
                       y = y_pos - var_end_height, yend = y_pos + var_end_height),
                   color = adjustcolor(ote_color, alpha.f = var_alpha),
                   linewidth = var_end_width) +
      # Right endline
      geom_segment(aes(x = .data[[ote_upper_var]], xend = .data[[ote_upper_var]],
                       y = y_pos - var_end_height, yend = y_pos + var_end_height),
                   color = adjustcolor(ote_color, alpha.f = var_alpha),
                   linewidth = var_end_width)
  }

  # Finalize the plot
  p <- p +
    geom_text(aes(x = x_min - 0.01, label = unit), hjust = 1, size = 4) +
    geom_text(aes(x = x_min - 0.01, label = group_label), hjust = 1, vjust = 2, size = 3, color = "darkgray") +
    scale_x_continuous(
      limits = c(x_min - 0.3, 1.2),
      labels = scales::percent,
      breaks = seq(x_min, 1, 0.1),
      expand = c(0, 0)
    ) +
    scale_y_continuous(limits = c(0.5, nrow(results_data) + 0.5), expand = c(0, 0)) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_line(linewidth = 0.5),
      plot.title = element_text(hjust = 0.5, size = 16),
      plot.margin = margin(5.5, 40, 5.5, 60, "pt")
    ) +
    labs(title = title)

  # Add legend based on label style
  if (label_style != "none") {
    if (label_style == "small") {
      legend_data <- data.frame(
        x = c(1.05, 1.05),
        y = c(nrow(results_data), nrow(results_data) - 0.15),
        label = c("TTE", "OTE"),
        color = c(tte_color, ote_color)
      )
      p <- p +
        geom_point(data = legend_data, aes(x = x, y = y, color = label), size = 3) +
        geom_text(data = legend_data, aes(x = x + 0.02, y = y, label = label, color = label),
                  hjust = 0, size = 3)
    } else if (label_style == "full") {
      legend_data <- data.frame(
        x = c(1.05, 1.10),
        y = c(0.55, 0.55),
        label = c("Theoretical Triage Effectiveness", "Observed Triage Effectiveness"),
        color = c(tte_color, ote_color)
      )
      p <- p +
        geom_point(data = legend_data, aes(x = x-0.02, y = y, color = label), size = 3) +
        geom_text(data = legend_data, aes(x = x, y = y + 0.15, label = label, color = label),
                  angle = 90, hjust = 0, vjust = 0, size = 3)
    } else if (label_style == "half") {
      legend_data <- data.frame(
        x = c(1.05, 1.10),
        y = c(0.55, 0.55),
        label = c("Theoretical TE", "Observed TE"),
        color = c(tte_color, ote_color)
      )
      p <- p +
        geom_point(data = legend_data, aes(x = x-0.02, y = y, color = label), size = 3) +
        geom_text(data = legend_data, aes(x = x, y = y + 0.15, label = label, color = label),
                  angle = 90, hjust = 0, vjust = 0, size = 3)
    }

    # Add the coord_cartesian after adding the legend elements
    p <- p +
      scale_color_manual(values = c(
        "TTE" = tte_color,
        "OTE" = ote_color,
        "Theoretical Triage Effectiveness" = tte_color,
        "Observed Triage Effectiveness" = ote_color,
        "Theoretical TE" = tte_color,
        "Observed TE" = ote_color
      )) +
      theme(legend.position = "none") +
      coord_cartesian(xlim = c(x_min - 0.3, 1.2),
                      ylim = c(0.5, max(nrow(results_data) + 0.5)),
                      clip = "off")
  }

  return(p)
}

