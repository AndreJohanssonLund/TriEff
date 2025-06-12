#' Plot Transferability Convergence Analysis
#'
#' @description
#' Creates convergence plots for transferability analysis results, showing running means
#' and confidence interval widths for each unit combination. This helps assess whether
#' enough Monte Carlo iterations were used and identify combinations where estimates
#' haven't stabilized.
#'
#' @param results List object returned by \code{assess_transferability()}
#' @param cluster_plot Integer. Number of combinations to include per plot. If 0
#'   (default), all combinations are shown in a single plot. Any positive integer is accepted.
#' @param ncol Integer or NULL. Number of columns in facet grid within each plot.
#'   If NULL (default), automatically calculated as ceiling(sqrt(cluster_plot)) for optimal layout.
#'   If cluster_plot is 0, automatically calculated as ceiling(sqrt(n_combinations)).
#'
#' @return A ggplot object (if cluster_plot is 0) or a list of ggplot objects
#'   (if cluster_plot > 0)
#'
#' @details
#' For each unit combination, the function plots:
#' - Running mean of difference from original TE
#' - Running confidence interval width
#'
#' The plots help identify:
#' - Whether estimates have converged
#' - Which combinations need more iterations
#' - Relative stability across different combinations
#'
#' When \code{cluster_plot} is used, combinations are grouped into separate plots
#' to avoid overcrowding. The layout is automatically optimized to be as square as possible
#' by calculating ncol as ceiling(sqrt(cluster_plot)).
#'
#' @examples
#' \dontrun{
#' # Run transferability analysis
#' results <- assess_transferability(data, n_iterations = 200)
#'
#' # Plot all combinations in one view
#' plot_transferability_convergence(results)
#'
#' # Create clustered plots with 9 combinations each (3x3 layout)
#' convergence_plots <- plot_transferability_convergence(results, cluster_plot = 9)
#'
#' # Create clustered plots with custom column count
#' convergence_plots <- plot_transferability_convergence(results, cluster_plot = 8, ncol = 4)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line labs theme_minimal theme_bw element_text
#' @importFrom ggplot2 facet_wrap scale_color_manual ggtitle
#' @importFrom dplyr arrange group_by mutate ungroup n distinct
#' @importFrom purrr map
#'
#' @export
plot_transferability_convergence <- function(results, cluster_plot = 0, ncol = NULL) {

  # Validate inputs
  if (!is.list(results) || !"distributions" %in% names(results)) {
    stop("Input must be a results list from assess_transferability() with a 'distributions' component")
  }

  if (!is.numeric(cluster_plot) || cluster_plot < 0 || cluster_plot != round(cluster_plot)) {
    stop("cluster_plot must be a non-negative integer (0 for all in one plot)")
  }

  # Extract distributions data
  dist_data <- results$distributions

  # Create combination identifier for easier handling (using ASCII)
  dist_data <- dist_data %>%
    mutate(combination = paste0(original_unit, " -> ", applied_unit)) %>%
    arrange(combination, iteration)

  # Calculate running statistics for each combination
  convergence_data <- dist_data %>%
    group_by(combination) %>%
    arrange(iteration) %>%
    mutate(
      # Running mean
      running_mean = cumsum(diff_from_original) / seq_len(n()),
      # Running variance for CI calculation
      running_var = map_dbl(seq_len(n()), function(i) {
        if (i == 1) 0 else var(diff_from_original[1:i])
      }),
      # Standard error
      running_se = sqrt(running_var / seq_len(n())),
      # 95% CI width (approximation)
      ci_width = 2 * 1.96 * running_se
    ) %>%
    ungroup()

  # Get total number of combinations
  combinations <- unique(convergence_data$combination)
  n_combinations <- length(combinations)

  # Create single plot function
  create_convergence_plot <- function(data, title_suffix = "", auto_ncol = NULL) {
    # Determine ncol for facet_wrap
    if (!is.null(auto_ncol)) {
      facet_ncol <- auto_ncol
    } else if (!is.null(ncol)) {
      facet_ncol <- ncol
    } else {
      # Default to 2 if not specified
      facet_ncol <- 2
    }

    # Create the convergence plot
    p <- ggplot(data, aes(x = iteration)) +
      geom_line(aes(y = running_mean, color = "Running Mean"), linewidth = 0.8) +
      geom_line(aes(y = ci_width, color = "CI Width"), linewidth = 0.8) +
      facet_wrap(~ combination, scales = "free_y", ncol = facet_ncol) +
      scale_color_manual(
        name = "Metric",
        values = c("Running Mean" = "#2E8B57", "CI Width" = "#CD853F"),
        labels = c("Running Mean", "95% CI Width")
      ) +
      labs(
        title = paste("Transferability Convergence Analysis", title_suffix),
        x = "Iteration",
        y = "Value",
        color = "Metric"
      ) +
      theme_bw() +
      theme(
        strip.text = element_text(size = 8),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10)
      )

    return(p)
  }

  # If cluster_plot is 0, show all combinations in one plot
  if (cluster_plot == 0) {
    # Calculate automatic ncol based on total combinations if not specified
    if (is.null(ncol)) {
      auto_ncol <- ceiling(sqrt(n_combinations))
    } else {
      auto_ncol <- ncol
    }

    return(create_convergence_plot(convergence_data, auto_ncol = auto_ncol))
  }

  # If clustering, split data and create multiple plots
  # Calculate automatic ncol for clusters if not specified
  if (is.null(ncol)) {
    cluster_ncol <- ceiling(sqrt(cluster_plot))
  } else {
    cluster_ncol <- ncol
  }

  # Calculate number of plots needed
  n_plots <- ceiling(n_combinations / cluster_plot)

  # Create clusters
  plots <- vector("list", n_plots)

  for (i in 1:n_plots) {
    # Calculate which combinations go in this plot
    start_idx <- (i - 1) * cluster_plot + 1
    end_idx <- min(i * cluster_plot, n_combinations)

    # Select combinations for this cluster
    current_combinations <- combinations[start_idx:end_idx]

    # Filter data for this cluster
    cluster_data <- convergence_data %>%
      filter(combination %in% current_combinations)

    # Create plot for this cluster
    title_suffix <- paste("(Part", i, "of", n_plots, ")")
    plots[[i]] <- create_convergence_plot(cluster_data, title_suffix, cluster_ncol)
  }

  return(plots)
}






#' Plot Transferability Analysis Results
#'
#' @description
#' Creates a heatmap visualization of cross-site transferability analysis results,
#' showing how triage effectiveness measurements transfer between different hospital
#' settings. Each cell displays the mean difference in TE when one site's priority
#' distribution is applied to another site's patient data.
#'
#' @param results List object returned by \code{assess_transferability()}
#' @param show_labels Logical. If TRUE, displays numeric values on heatmap cells. Default TRUE
#' @param show_n_patients Logical. If TRUE, includes patient count in unit labels.
#'   Default TRUE
#' @param color_scheme Character. Color scheme for the heatmap: "oceansky" (default),
#'   "thermalsafe", "blueorange", or "triage". The "oceansky" scheme is colorblind-safe
#' @param scale_limits Numeric vector of length 2. Custom scale limits for color mapping
#'   as c(min, max). If NULL (default), uses data range
#' @param title Character. Custom title for the plot. If NULL (default), automatically
#'   generates title based on calculation method
#' @param subtitle Character. Custom subtitle for the plot. If NULL (default), uses
#'   standard subtitle
#' @param sort_units Character. How to sort units: "char" for alphabetical (A->Z),
#'   "n" for numerical descending (higher numbers first). Default is "char"
#'
#' @return A ggplot object representing the transferability heatmap
#'
#' @details
#' The heatmap displays:
#' - Rows: Original units (source of patient data)
#' - Columns: Applied units (source of priority distribution)
#' - Cell values: Mean difference in TE (original TE - recreated TE)
#' - Cell text: Mean difference and 95% confidence interval
#'
#' Positive values (warmer colors) indicate that applying the column unit's priorities
#' to the row unit's patients results in lower TE than the original. Negative values
#' (cooler colors) indicate higher TE than the original.
#'
#' Diagonal cells represent the baseline (same unit, same distribution) and should
#' ideally be close to zero, representing the stochastic variation around the original TE.
#'
#' @examples
#' \dontrun{
#' # Run transferability analysis
#' results <- assess_transferability(data, calc_method = "rte")
#'
#' # Create basic heatmap
#' plot_transferability(results)
#'
#' # Customize appearance
#' plot_transferability(results,
#'   color_scheme = "thermalsafe",
#'   title = "Custom Title",
#'   scale_limits = c(-0.1, 0.1),
#'   sort_units = "n"
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_text geom_rect scale_fill_gradientn
#' @importFrom ggplot2 scale_fill_manual element_blank element_text labs theme_minimal
#' @importFrom ggplot2 theme guide_colorbar waiver scale_x_discrete scale_y_discrete
#' @importFrom scales rescale
#' @importFrom dplyr mutate select arrange desc
#' @importFrom forcats fct_rev fct_reorder
#'
#' @export
plot_transferability <- function(results,
                                 show_labels = TRUE,
                                 show_n_patients = TRUE,
                                 color_scheme = "oceansky",
                                 scale_limits = NULL,
                                 title = NULL,
                                 subtitle = NULL,
                                 sort_units = "char") {

  # Validate inputs
  if (!is.list(results) || !"results" %in% names(results)) {
    stop("Input must be a results list from assess_transferability() with a 'results' component")
  }

  # Validate sort_units parameter
  if (!sort_units %in% c("char", "n")) {
    stop("sort_units must be either 'char' or 'n'")
  }

  # Extract results data
  heatmap_data <- results$results

  # Extract metadata for default titles and patient counts
  metadata <- results$metadata
  calc_method <- if (!is.null(metadata$calc_method)) metadata$calc_method else "TE"

  # Get patient counts from metadata if available
  unit_n_patients <- if (!is.null(metadata$unit_n_patients)) {
    # Convert tibble to named vector for existing code to work
    counts <- metadata$unit_n_patients$n_patients
    names(counts) <- metadata$unit_n_patients$unit
    counts
  } else {
    NULL
  }

  # Set default titles based on metadata
  if (is.null(title)) {
    method_name <- switch(calc_method,
                          "rte" = "Rank-based TE (RTE)",
                          "wte" = "Waiting Time-based TE (WTE)",
                          "Triage Effectiveness")
    title <- paste("Transferability Analysis:", method_name)
  }

  if (is.null(subtitle)) {
    subtitle <- "Mean difference when applying priority distributions across units"
  }

  # Validate color scheme
  valid_schemes <- c("oceansky", "thermalsafe", "blueorange", "triage")
  if (!color_scheme %in% valid_schemes) {
    warning(paste("Invalid color scheme '", color_scheme, "'. Using 'oceansky' instead. Valid options are: ",
                  paste(valid_schemes, collapse = ", ")))
    color_scheme <- "oceansky"
  }

  # Sort units based on sort_units parameter
  if (sort_units == "char") {
    # Alphabetical sorting (A->Z)
    unit_levels <- sort(unique(c(heatmap_data$original_unit, heatmap_data$applied_unit)))
  } else if (sort_units == "n") {
    # Sort by patient count (higher numbers first)
    all_units <- unique(c(heatmap_data$original_unit, heatmap_data$applied_unit))
    if (!is.null(unit_n_patients)) {
      # Sort units by their patient counts in descending order
      unit_counts_df <- data.frame(
        unit = names(unit_n_patients),
        n_patients = unname(unit_n_patients),
        stringsAsFactors = FALSE
      )
      # Filter to only units present in the data and sort by patient count descending
      unit_levels <- unit_counts_df %>%
        filter(unit %in% all_units) %>%
        arrange(desc(n_patients)) %>%
        pull(unit)
    } else {
      # Fallback to alphabetical if no patient count data
      warning("No patient count data available, falling back to alphabetical sorting")
      unit_levels <- sort(all_units)
    }
  }

  # Convert to factors with proper ordering
  heatmap_data <- heatmap_data %>%
    mutate(
      original_unit = factor(original_unit, levels = unit_levels),
      applied_unit = factor(applied_unit, levels = unit_levels)
    )

  # Get scale limits
  if (is.null(scale_limits)) {
    min_val <- min(heatmap_data$mean_diff, na.rm = TRUE)
    max_val <- max(heatmap_data$mean_diff, na.rm = TRUE)
  } else {
    if (length(scale_limits) != 2) {
      stop("scale_limits must be a numeric vector of length 2 (min, max)")
    }
    min_val <- scale_limits[1]
    max_val <- scale_limits[2]
  }

  # Define color schemes (updated to use extremes + white center for all schemes)
  get_colors <- function(min_value, max_value, scheme) {
    # Define color palettes for each scheme
    color_palettes <- list(
      oceansky = list(
        neg = "#2c7bb6",
        neutral = "white",
        pos = "gold"
      ),
      thermalsafe = list(
        neg = "purple",
        neutral = "white",
        pos = "orange"
      ),
      blueorange = list(
        neg = "#2c7bb6",
        neutral = "white",
        pos = "#d7191c"
      ),
      triage = list(
        neg = "red",
        neutral = "white",
        pos = "green"
      )
    )

    palette <- color_palettes[[scheme]]

    # All schemes now use the same logic: extremes + white center
    # Determine which case we're in
    if (min_value < 0 && max_value > 0) {
      # Values span across zero - use 3 colors: extreme_neg, white, extreme_pos
      colors <- c(palette$neg, palette$neutral, palette$pos)
      values <- scales::rescale(c(min_value, 0, max_value))
    } else if (min_value < 0) {
      # Only negative values - use extreme_neg and white
      colors <- c(palette$neg, palette$neutral)
      values <- scales::rescale(c(min_value, 0))
    } else {
      # Only positive values - use white and extreme_pos
      colors <- c(palette$neutral, palette$pos)
      values <- scales::rescale(c(0, max_value))
    }

    list(colors = colors, values = values)
  }

  # Create unit labels with patient counts if requested and available
  if (show_n_patients && !is.null(unit_n_patients)) {
    # Create a lookup table for unit labels with counts
    unit_labels_x <- sapply(levels(heatmap_data$applied_unit), function(u) {
      count <- unit_n_patients[[u]]
      if (!is.null(count)) {
        paste0(u, "\n(n=", count, ")")
      } else {
        u
      }
    })

    unit_labels_y <- sapply(levels(heatmap_data$original_unit), function(u) {
      count <- unit_n_patients[[u]]
      if (!is.null(count)) {
        paste0(u, "\n(n=", count, ")")
      } else {
        u
      }
    })

    # Apply custom labels
    x_lab_name <- "Applied Unit (Priority Distribution Source)"
    y_lab_name <- "Original Unit (Patient Data Source)"
  } else {
    # Use default labels
    unit_labels_x <- waiver()
    unit_labels_y <- waiver()
    x_lab_name <- "Applied Unit"
    y_lab_name <- "Original Unit"
  }

  color_scheme_data <- get_colors(min_val, max_val, color_scheme)

  # Create the heatmap
  p <- ggplot(heatmap_data, aes(x = applied_unit, y = fct_rev(original_unit), fill = mean_diff)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradientn(
      colors = color_scheme_data$colors,
      values = color_scheme_data$values,
      limits = c(min_val, max_val),
      name = "Mean Difference",
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = 12,
        barheight = 0.8
      )
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_lab_name,
      y = y_lab_name
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
      axis.text.y = element_text(hjust = 1, size = 9),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "bottom",
      panel.grid = element_blank(),
      axis.ticks = element_blank()
    )

  # Apply custom axis labels if needed
  if (show_n_patients && !is.null(unit_n_patients)) {
    p <- p +
      scale_x_discrete(labels = unit_labels_x) +
      scale_y_discrete(labels = unit_labels_y)
  }

  # Add text labels if requested
  if (show_labels) {
    p <- p + geom_text(
      aes(label = sprintf("%.3f\n(%.3f, %.3f)", mean_diff, lower_ci, upper_ci)),
      color = "black",
      size = 2.2,
      fontface = "bold",
      lineheight = 0.8
    )
  }

  # Highlight diagonal cells
  p <- p + geom_rect(
    data = heatmap_data %>% filter(as.character(original_unit) == as.character(applied_unit)),
    aes(xmin = as.numeric(applied_unit) - 0.45,
        xmax = as.numeric(applied_unit) + 0.45,
        ymin = as.numeric(fct_rev(original_unit)) - 0.45,
        ymax = as.numeric(fct_rev(original_unit)) + 0.45),
    fill = NA,
    color = "gray60",
    linewidth = 1,
    inherit.aes = FALSE
  )

  return(p)
}
