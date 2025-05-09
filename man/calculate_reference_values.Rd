% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_wte.R
\name{calculate_reference_values}
\alias{calculate_reference_values}
\title{Calculate Reference Values for TE Calculations}
\usage{
calculate_reference_values(df, subgroup = NULL, var1 = NULL, var2 = NULL)
}
\arguments{
\item{df}{Data frame containing patient data}

\item{subgroup}{Optional list for subgroup analysis}

\item{var1}{Optional string for first comparison variable}

\item{var2}{Optional string for second comparison variable}
}
\value{
A tibble containing:
\itemize{
\item Grouping variables (unit and any comparison variables)
\item mean_all: Mean waiting time for the group
\item n_loset: Number of loset cases in the group
}
}
\description{
Calculate Reference Values for TE Calculations
}
\details{
Calculates essential reference values needed for TE calculations using the
same grouping structure as the main calculations.
}
\keyword{internal}
