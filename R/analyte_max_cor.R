# ============================================================
# analyte_max_cor: cohesive data-frame output + one printer
# ============================================================

# -----------------------------
# Internal helper (not exported)
# -----------------------------

#' @keywords internal
#' @noRd
.strip_key <- function(x) {
  tolower(gsub("\\s+", " ", trimws(gsub("[^[:alnum:] ]+"," ", x))))
}

#' R6 Class representing the results from an analyte_max_cor_res call
#'
#' @description
#' Used for pretty printing
Analyte_max_cor_res <- R6::R6Class("Analyte_max_cor_res",
                            public = list(
                              #' @field disp Display table.
                              disp = NULL,
                              #' @field k Number of replicates averaged.
                              k = NULL,
                              #' @field use Which cross-sectional CV was used.
                              use = NULL,
                              #' @description
                              #' Create a new Analyte_max_cor_res object
                              #' @param disp Display table.
                              #' @param k Number of replicates averaged.
                              #' @param use Which cross-sectional CV was used.
                              #' @return A new `Analyte_max_cor_res` object
                              initialize = function(disp, k, use) {
                                self$disp <- disp
                                self$k <- k
                                self$use <- use
                              },
                              #' @description
                              #' Pretty printing of analyte_max_cor results
                              print = function() {
                                cat("\n# Maximum pre-post correlation ceilings\n", sep = "")
                                cat("use = ", self$use, ", k = ", self$k, "\n\n", sep = "")
                                print(self$disp, row.names = FALSE, right = TRUE)
                              }
                            )
)

#' Maximum plausible pre/post correlation among analytes (rounding-aware)
#'
#' Uses CVA/CVI for the given analyte (from \code{MEPtoolkit::analytes}) and the
#' study's pre/post mean & SD to compute the pre–post correlation you would expect
#' if there were no person-to-person heterogeneity in true change; i.e., a
#' maximum expected correlation (a ceiling).
#'
#' You can pass a single case (scalars) or a data frame with columns for multiple
#' cases. For rounding-aware bounds, \strong{m1, s1, m2, s2 must be character strings}
#' such as \code{"64.3"} and \code{"25.4"} so the number of decimals encodes rounding.
#'
#' @details
#' \strong{Goal.} Return the maximum plausible pre–post correlation \eqn{r_{\max}}
#' for an analyte if everyone changes by the same amount (no between-subject
#' heterogeneity). This is the single-measurement ICC ceiling implied by assay plus
#' within-person biology under your study variability.
#'
#' \strong{Inputs used.}
#' \itemize{
#'   \item Loads \code{MEPtoolkit::analytes} and matches \code{analyte}
#'         (case/space/punctuation-insensitive).
#'   \item Gets \eqn{CV_A} (analytic) and \eqn{CV_I} (intraindividual) as
#'         proportions. If \eqn{k} replicates were averaged, uses \eqn{CV_A/\sqrt{k}}.
#'   \item Forms cross-sectional CVs \eqn{CV_{T,\mathrm{pre}} = s_1/m_1} and
#'         \eqn{CV_{T,\mathrm{post}} = s_2/m_2}; uses \code{"pre"}, \code{"post"},
#'         or \code{"average"} per \code{use}.
#' }
#'
#' \strong{Computation (additive \%CV algebra).}
#' \deqn{CV_W^2 = CV_A^2 + CV_I^2,\qquad
#'       CV_G^2 = \mathrm{max}\bigl(CV_T^2 - CV_W^2,\, 0\bigr).}
#' \deqn{r_{\max} = \frac{CV_G^2}{CV_G^2 + CV_W^2}, \qquad
#'       r_{\max}^{(\mathrm{analytic})} = \frac{CV_G^2}{CV_G^2 + CV_A^2}.}
#'
#' \strong{Rounding model.}
#' \itemize{
#'   \item Study inputs: \code{m1}, \code{s1}, \code{m2}, \code{s2} are treated as rounded
#'         to the provided number of decimals; intervals are constructed accordingly via
#'         \code{validate_descriptive()}.
#'   \item Reference table: \code{CVI} is rounded to 3 dp (modeled as \code{cvi - 0.0005} minimum);
#'         \code{CVA} is derived from QC mean/sd, each rounded to 2 dp, then scaled for
#'         replicates by \eqn{1/\sqrt{k}}.
#'   \item The ceiling is maximized by using the largest plausible \eqn{CV_T}
#'         and the smallest plausible within-person variance.
#' }
#'
#' \strong{Assumptions.}
#' \itemize{
#'   \item Same measurand/method as the quality control data.
#'   \item Additive variance algebra on the CV scale:
#'         \eqn{CV_T^2 \approx CV_G^2 + CV_W^2}.
#'   \item Prefer \code{use = "pre"} unless variability post-intervention is a pure
#'         location shift (heterogeneity inflates \eqn{CV_T} at post).
#' }
#'
#' \strong{Units.} The calculation uses SD divided by the mean, so any
#' multiplicative unit change (e.g., mg/dL \eqn{\leftrightarrow} mmol/L) cancels out; you do not
#' need identical units across sources as long as they are linked by a pure rescale.
#' Do not mix scales (e.g., mean on raw scale, SD on log scale). For analytes with
#' affine systems (non-zero intercept) such as HbA1c IFCC vs NGSP, use the matching
#' row or convert both mean and SD to that system before calling the function; affine
#' transforms change CVs.
#'
#' \strong{Edge cases.}
#' \itemize{
#'   \item If \eqn{CV_T^2 \le CV_W^2}, set \eqn{CV_G^2 = 0} and return \eqn{r_{\max} = 0}.
#'   \item If a denominator is numerically zero, return \code{1} for that ratio.
#' }
#'
#' \strong{Available analytes.}
#' The \code{analyte} argument (matching is case/space/punctuation-insensitive)
#' can be any of the following:
#' \itemize{
#'   \item \code{25-hydroxy-Vitamin D}
#'   \item \code{Alanine Aminotransferase}
#'   \item \code{Albumin}
#'   \item \code{Alkaline Phosphatase}
#'   \item \code{Aspartate Aminotransferase}
#'   \item \code{Calcium}
#'   \item \code{Chloride}
#'   \item \code{Creatinine}
#'   \item \code{Gamma Glutamyltransferase}
#'   \item \code{Glucose}
#'   \item \code{HbA1c Diabetes IFCC}
#'   \item \code{HbA1c Diabetes NGSP}
#'   \item \code{HbA1c Healthy IFCC}
#'   \item \code{HbA1c Healthy NGSP}
#'   \item \code{HDL Cholesterol}
#'   \item \code{Hemoglobin}
#'   \item \code{Holotranscobalamin}
#'   \item \code{Iron}
#'   \item \code{Lactate}
#'   \item \code{Lactate Dehydrogenase}
#'   \item \code{LDL Cholesterol}
#'   \item \code{Magnesium}
#'   \item \code{Osmolality}
#'   \item \code{PCO2}
#'   \item \code{Phosphate}
#'   \item \code{Potassium}
#'   \item \code{Rheumatoid Factor}
#'   \item \code{Sodium}
#'   \item \code{Thyroid Stimulating Hormone}
#'   \item \code{Total Bilirubin}
#'   \item \code{Total Cholesterol}
#'   \item \code{Total Protein}
#'   \item \code{Total Testosterone}
#'   \item \code{Transferrin}
#'   \item \code{Triglycerides}
#'   \item \code{Urea}
#'   \item \code{Uric Acid}
#'   \item \code{Vitamin B12}
#' }
#'
#' @param analyte Character scalar or column in \code{data}.
#' @param m1,s1,m2,s2 Character scalars or columns in \code{data} (use strings like \code{"64.3"}).
#' @param k Integer (\eqn{\ge} 1). Number of replicates averaged.
#' @param use Which cross-sectional CV to use: \code{"pre"}, \code{"post"}, or \code{"average"}.
#' @param data Optional data.frame providing columns for analyte, m1, s1, m2, s2.
#' @param append Logical. If \code{TRUE} and \code{data} provided, append ceilings as columns.
#' @param out_cols Length-7 character vector: names for analytic-only and total ceilings.
#'   Defaults to \code{c("r_max_analytic","r_max_total")}.
#' @param output Logical. If \code{TRUE}, also outputs a summary of results in the console.
#'
#' @return
#' \itemize{
#'   \item Single case (no \code{data}, \code{append = FALSE}): an object of class
#'         \code{"analyte_max_cor"} with detailed components and a pretty printer.
#'   \item Multiple cases, \code{append = FALSE}: a data.frame with columns
#'         \code{case}, \code{analyte}, \code{cvi}, \code{cva}, \code{cvt}, \code{cvw}, \code{cvg},
#'         \code{r_max_analytic}, \code{r_max_total}. It is tagged so that
#'         \code{print.analyte_max_cor()} prints a compact table.
#'   \item Multiple cases, \code{append = TRUE}: the input \code{data} with two new columns
#'         named by \code{out_cols}.
#' }
#'
#' @examples
#' analyte_max_cor("25-hydroxy-Vitamin D",
#'                 m1 = "64.3", s1 = "25.4",
#'                 m2 = "88.5", s2 = "23.2",
#'                 k = 1, use = "pre")
#'
#' df <- data.frame(
#'   analyte = c("Sodium","Total Cholesterol"),
#'   m1 = c("140.2","6.60"),
#'   s1 = c("2.10","0.90"),
#'   m2 = c("139.8","6.40"),
#'   s2 = c("2.05","0.90"),
#'   stringsAsFactors = FALSE
#' )
#' analyte_max_cor(analyte = analyte, m1 = m1, s1 = s1, m2 = m2, s2 = s2,
#'                 k = 1, use = "pre", data = df, append = TRUE)
#'
#' @export
analyte_max_cor <- function(analyte,
                            m1, s1, m2, s2,
                            k = 1,
                            use = c("pre", "post", "average"),
                            data = NULL,
                            append = FALSE,
                            out_cols = c("cvi", "cva", "cvt", "cvw", "cvg", "r_max_analytic","r_max_total"),
                            output = FALSE) {
  use <- match.arg(use)

  if(!is.null(data)) {
    if(!is.data.frame(data)) stop("`data` must be a data.frame.")
    s1 <- resolve(enquo(s1), data); s2 <- resolve(enquo(s2), data);
    m1 <- resolve(enquo(m1), data); m2 <- resolve(enquo(m2), data);
    analyte <- resolve(enquo(analyte), data)
  }

  # -- check lengths
  args <- list(s1=s1, s2=s2, m1=m1, m2=m2, analyte=analyte)
  lens <- vapply(args, length, integer(1))
  L <- max(lens)
  ok_len <- lens == 0 | lens == L
  if (!all(ok_len))
    stop("analyte, m1, s1, m2 & s2 must have the same length.")

  a_disp <- character(L)
  cvi <- cva <- cvt <- cvw <- cvg <- rA <- rT <- rep(NA_real_, L)

  r_m1 <- validate_descriptive(m1, allow_na = TRUE)
  r_s1 <- validate_descriptive(s1, allow_na = TRUE)
  r_m2 <- validate_descriptive(m2, allow_na = TRUE)
  r_s2 <- validate_descriptive(s2, allow_na = TRUE)

  for(i in seq_len(L)) {
    tbl <- MEPtoolkit::analytes %>%
                    filter(grepl(.strip_key(.env$analyte[i]), .strip_key(.data$analyte),
                                 fixed = TRUE))
    if(nrow(tbl) > 1) {
      stop("Analyte name ambiguous (several analytes returned), please provide
         an unambiguous name.")
    }
    else if (nrow(tbl)==0){
      stop("Analyte not found.")
    }
    row <- unlist(tbl)

    s1_max <- r_s1$maximum[i]; s2_max <- r_s2$maximum[i]
    m1_min <- r_m1$minimum[i]; m2_min <- r_m2$minimum[i]
    if (m1_min <= 0 || m2_min <= 0) stop("Means must be positive; check inputs/rounding.")

    CVT_pre_max  <- s1_max / m1_min
    CVT_post_max <- s2_max / m2_min
    CVT_use_max  <- switch(use,
                           pre     = CVT_pre_max,
                           post    = CVT_post_max,
                           average = (CVT_pre_max + CVT_post_max)/2)

    CVI_min  <- validate_descriptive(row["cvi"])$minimum
    mean_max <- validate_descriptive(row["mean"])$maximum
    sd_min   <- validate_descriptive(row["sd"])$minimum
    if (mean_max <= 0) stop("QC mean <= 0; cannot form CVA bounds.")
    CVA_min  <- (sd_min / mean_max) / sqrt(k)
    CVW2_min <- CVA_min^2 + CVI_min^2
    CVG2_max <- max(CVT_use_max^2 - CVW2_min, 0)
    denom_total    <- CVG2_max + CVW2_min
    denom_analytic <- CVG2_max + CVA_min^2

    r_total_max    <- if (denom_total    == 0) 1 else CVG2_max / denom_total
    r_analytic_max <- if (denom_analytic == 0) 1 else CVG2_max / denom_analytic

    # display CVs at “reported” values (not edge-case maxima)
    m1o <- r_m1$original[i]; s1o <- r_s1$original[i]
    m2o <- r_m2$original[i]; s2o <- r_s2$original[i]
    CVT_used <- switch(use, pre = s1o/m1o,
                       post = s2o/m2o,
                       average = (s1o/m1o + s2o/m2o)/2)
    CVA_disp <- validate_descriptive(row["cva"])$original / sqrt(k)
    CVI_disp <- validate_descriptive(row["cvi"])$original
    CVW_disp <- sqrt(CVA_disp^2 + CVI_disp^2)
    CVG_disp <- sqrt(pmax(CVT_used^2 - CVW_disp^2, 0))

    ### Prepare outputs
    a_disp[i] <- row["analyte"];
    cvi[i] <- CVI_disp; cva[i] <- CVA_disp; cvt[i] <- CVT_used;
    cvw[i] <- CVW_disp; cvg[i] <- CVG_disp;
    rA[i] <- r_analytic_max; rT[i] <- r_total_max

  }

  # display table that the printer will use (same for single/multi)
  disp <- data.frame(
    case = seq_len(L),
    analyte = a_disp,
    cvi = cvi,
    cva = cva,
    cvt = cvt,
    cvw = cvw,
    cvg = cvg,
    r_max_analytic = rA,
    r_max_total    = rT,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  res <- Analyte_max_cor_res$new(disp, k, use)

  if (append && !is.null(data)) {
    if (!(length(out_cols) == 7 && all(nzchar(out_cols))))
      stop("`out_cols` must be length 7 with non-empty names.")
    out <- data
    out[[out_cols[1]]] <- cvi; out[[out_cols[2]]] <- cva;
    out[[out_cols[3]]] <- cvt; out[[out_cols[4]]] <- cvw;
    out[[out_cols[5]]] <- cvg; out[[out_cols[6]]] <- rA;
    out[[out_cols[7]]] <- rT
  } else {
    # return the res object itself when not appending
    out <- res
  }
  if(output == TRUE)
  {
    print(res)
  }
  return(out)
}
