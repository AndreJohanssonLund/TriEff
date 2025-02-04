#' Synthetic Emergency Department Data from Malmö Hospital
#'
#' A synthetic dataset generated based on statistical properties of ED visits at Malmö Hospital
#' from the Skåne Emergency Medicine (SEM) database. Generated using synthpop methodology,
#' this dataset mimics real patterns while containing no actual patient data.
#'
#' @format A data frame with 124,311 observations and 16 variables:
#' \describe{
#'   \item{priority}{Factor with 4 levels (1-4) representing triage priority, where 1 is highest priority}
#'   \item{ambulance}{Logical indicating if patient arrived by ambulance}
#'   \item{unit}{Factor with 3 levels: "medical", "orthopedics", "surgical" indicating treating specialty}
#'   \item{chief_complaint}{Factor with 68 levels representing the main reason for ED visit}
#'   \item{age_at_arrival}{Numeric age of patient at time of arrival}
#'   \item{male}{Logical indicating patient gender (TRUE for male)}
#'   \item{planed_return_visit}{Logical indicating if visit was a planned return}
#'   \item{referral}{Logical indicating if patient was referred to the ED}
#'   \item{exit_info}{Factor with 4 levels: "Admitted", "Died", "Home", "Unknown" indicating ED visit outcome}
#'   \item{death_days}{Numeric, days until death if applicable, NA otherwise}
#'   \item{ward}{Factor with 12 levels indicating receiving type of ward for admitted patients}
#'   \item{loset}{Logical indicating if case was time-critical according to LOSET criteria}
#'   \item{arrival_at_hospital}{POSIXct timestamp of arrival at hospital}
#'   \item{exit_ed}{POSIXct timestamp of exit from ED}
#'   \item{resolve}{POSIXct timestamp of first physician contact}
#'   \item{arrival}{POSIXct timestamp of triage assessment}
#' }
#'
#' @details
#' This synthetic dataset is based on ED visits from 2017-2018 at Malmö Hospital. It preserves
#' the statistical relationships and patterns of the original data while containing no real
#' patient information. The dataset includes adult patients (age >= 18) with medical, surgical,
#' or orthopedic complaints. However while the data is similar to SEM it is not identical.
#' Differences exists in waiting times but priorities are similar enough that experimentaiton
#' with the trieff statistical package is feasible.
#'
#' The dataset is used for running tests on the trieff statistical package on build.
#'
#' Time variables follow this sequence:
#' arrival_at_hospital <= arrival (triage) <= resolve (physician) <= exit_ed
#'
#' @source Based on statistical properties of the Skåne Emergency Medicine (SEM) database.
#' Reference: Ekelund, U., et al. (2024). The Skåne Emergency Medicine (SEM) cohort.
#' Scandinavian Journal of Trauma Resuscitation and Emergency Medicine, 32(1).
#' doi: https://doi.org/10.1186/s13049-024-01206-0
#'
"sem_malmo_synth"
"sem_malmo_synth_mini"
