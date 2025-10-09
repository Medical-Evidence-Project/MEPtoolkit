#' QC Analyte Table (CVA/CVI in proportions)
#'
#' A reference table of clinical analytes with analytical imprecision (CVA) and
#' intraindividual biological variation (CVI), used to compute reliability
#' ceilings and change thresholds.
#'
#' @format A data frame with 38 rows and 8 variables:
#' \describe{
#'   \item{analyte}{character. Analyte name (e.g., "Sodium").}
#'   \item{num_obs}{integer. Number of QC observations.}
#'   \item{mean}{numeric. Mean of QC results (in \code{units}).}
#'   \item{units}{character. Measurement units (e.g., "mmol/L").}
#'   \item{sd}{numeric. SD of QC results (in \code{units}).}
#'   \item{cva}{numeric. Analytical coefficient of variation as a proportion (0–1).}
#'   \item{cvi}{numeric. Intraindividual biological CV as a proportion (0–1).}
#'   \item{cvi_source}{character. Source/citation for CVI (e.g., "EFLM", "Westgard").}
#' }
#'
#' @details
#' Columns \code{cva} and \code{cvi} are stored as proportions (e.g., 0.021, not 2.1%).
#' Use with functions such as \code{analyte_max_cor()} that assume proportional CVs.
#'
#' @source
#' Row-level sources include EFLM, Westgard, and cited journal articles; see \code{cvi_source}.
#'
#' @references
#' McCormack, J. P., & Holmes, D. T. (2020). Your results may vary: the imprecision of medical
#' measurements. \emph{BMJ}, 368:m149.
#'
#' @examples
#' data(analytes)
#' head(analytes)
#'
#' # Select two analytes and show their CVs
#' analytes[analytes$analyte %in% c("Sodium", "Total Cholesterol"),
#'          c("analyte", "cva", "cvi")]
"analytes"
