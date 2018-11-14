#' Simulated subtype data
#'
#' A dataset containing 2000 patients: 1200 cases and 800 controls.
#' There are four subtypes, and both numeric and character subtype labels.
#' The subtypes are formed by cross-classification of two binary disease markers,
#' disease marker 1 and disease marker 2.
#' There are three risk factors, two continuous and one binary. One of the
#' continuous risk factors and the binary risk factor are related to the
#' disease subtypes.
#' There are also 30 continuous tumor markers, 20 of which are related to
#' the subtypes and 10 of which represent noise, which could be used in a
#' clustering analysis.
#'
#' @format A data frame with 2000 rows--one row per patient
#' \describe{
#'     \item{case}{Indicator of case control status, 1 for cases and 0 for controls}
#'     \item{subtype}{Numeric subtype label, 0 for control subjects}
#'     \item{subtype_name}{Character subtype label}
#'     \item{marker1}{Disease marker 1}
#'     \item{marker2}{Disease marker 2}
#'     \item{x1}{Continuous risk factor 1}
#'     \item{x2}{Continuous risk factor 2}
#'     \item{x3}{Binary risk factor}
#'     \item{y1}{Continuous tumor marker 1}
#'     \item{y2}{Continuous tumor marker 2}
#'     \item{y3}{Continuous tumor marker 3}
#'     \item{y4}{Continuous tumor marker 4}
#'     \item{y5}{Continuous tumor marker 5}
#'     \item{y6}{Continuous tumor marker 6}
#'     \item{y7}{Continuous tumor marker 7}
#'     \item{y8}{Continuous tumor marker 8}
#'     \item{y9}{Continuous tumor marker 9}
#'     \item{y10}{Continuous tumor marker 10}
#'     \item{y11}{Continuous tumor marker 11}
#'     \item{y12}{Continuous tumor marker 12}
#'     \item{y13}{Continuous tumor marker 13}
#'     \item{y14}{Continuous tumor marker 14}
#'     \item{y15}{Continuous tumor marker 15}
#'     \item{y16}{Continuous tumor marker 16}
#'     \item{y17}{Continuous tumor marker 17}
#'     \item{y18}{Continuous tumor marker 18}
#'     \item{y19}{Continuous tumor marker 19}
#'     \item{y20}{Continuous tumor marker 20}
#'     \item{y21}{Continuous tumor marker 21}
#'     \item{y22}{Continuous tumor marker 22}
#'     \item{y23}{Continuous tumor marker 23}
#'     \item{y24}{Continuous tumor marker 24}
#'     \item{y25}{Continuous tumor marker 25}
#'     \item{y26}{Continuous tumor marker 26}
#'     \item{y27}{Continuous tumor marker 27}
#'     \item{y28}{Continuous tumor marker 28}
#'     \item{y29}{Continuous tumor marker 29}
#'     \item{y30}{Continuous tumor marker 30}
#' }
"subtype_data"
