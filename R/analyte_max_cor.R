# ============================================================
# analyte_max_cor: cohesive data-frame output + one printer
# ============================================================

# -----------------------------
# Internal helpers (not exported)
# -----------------------------

#' @keywords internal
#' @noRd
.strip_key <- function(x) {
  tolower(gsub("\\s+", " ", trimws(gsub("[^[:alnum:] ]+"," ", x))))
}

#' @keywords internal
#' @noRd
.amc_row <- function(analyte, m1, s1, m2, s2, k, use) {
  # Require strings so decimals encode rounding
  need_char <- list(m1 = m1, s1 = s1, m2 = m2, s2 = s2)
  not_char <- names(need_char)[!vapply(need_char, is.character, logical(1))]
  if (length(not_char)) {
    stop("For rounding-aware ceilings, pass m1, s1, m2, s2 as character strings. ",
         "Non-character inputs: ", paste(not_char, collapse = ", "))
  }

  # Load reference table & locate columns
  dat_raw <- MEPtoolkit::analytes
  nm <- tolower(gsub("\\s+", " ", names(dat_raw)))
  i_name <- (which(nm %in% c("test","analyte","name")))[1]
  i_cva  <- (grep("^cva$", nm))[1]; if (is.na(i_cva))  i_cva  <- (grep("^cva",  nm))[1]
  i_cvi  <- (grep("^cvi$", nm))[1]; if (is.na(i_cvi))  i_cvi  <- (grep("^cvi",  nm))[1]
  i_mean <- (grep("^mean$", nm))[1]
  i_sd   <- (grep("^sd$",   nm))[1]
  if (any(is.na(c(i_name, i_cva, i_cvi, i_mean, i_sd)))) {
    stop("MEPtoolkit::analytes must contain analyte name, cva, cvi, mean, and sd columns.")
  }

  tbl <- data.frame(
    analyte_raw = dat_raw[[i_name]],
    analyte_key = .strip_key(dat_raw[[i_name]]),
    cva         = suppressWarnings(as.numeric(dat_raw[[i_cva]])),
    cvi         = suppressWarnings(as.numeric(dat_raw[[i_cvi]])),
    qc_mean_r   = suppressWarnings(as.numeric(dat_raw[[i_mean]])),  # rounded to 2 dp
    qc_sd_r     = suppressWarnings(as.numeric(dat_raw[[i_sd]])),    # rounded to 2 dp
    stringsAsFactors = FALSE
  )
  if (any(!is.finite(tbl$cva)) || any(!is.finite(tbl$cvi)) ||
      any(!is.finite(tbl$qc_mean_r)) || any(!is.finite(tbl$qc_sd_r))) {
    stop("Non-numeric or missing values in analyte CVs or QC mean/sd.")
  }

  # Match analyte
  key <- .strip_key(analyte)
  hits <- which(tbl$analyte_key == key)
  if (!length(hits)) hits <- grep(key, tbl$analyte_key, fixed = TRUE)
  if (!length(hits)) stop("Analyte not found in MEPtoolkit::analytes: '", analyte, "'.")
  if (length(hits) > 1) {
    message("Multiple matches: ", paste(tbl$analyte_raw[hits], collapse = "; "), ". Using first.")
    hits <- hits[1]
  }
  row <- tbl[hits, , drop = FALSE]

  # Rounding ranges for study inputs
  r_m1 <- validate_descriptive(m1, "m1")
  r_s1 <- validate_descriptive(s1, "s1")
  r_m2 <- validate_descriptive(m2, "m2")
  r_s2 <- validate_descriptive(s2, "s2")

  s1_min <- pmax(r_s1$minimum, 0); s1_max <- pmax(r_s1$maximum, 0)
  s2_min <- pmax(r_s2$minimum, 0); s2_max <- pmax(r_s2$maximum, 0)
  m1_min <- r_m1$minimum; m2_min <- r_m2$minimum
  if (m1_min <= 0 || m2_min <= 0) stop("Means must be positive; check inputs/rounding.")

  # Maximize CVT by using max SD and min mean
  CVT_pre_max  <- s1_max / m1_min
  CVT_post_max <- s2_max / m2_min
  CVT_use_max  <- switch(use,
                         pre     = CVT_pre_max,
                         post    = CVT_post_max,
                         average = (CVT_pre_max + CVT_post_max)/2)

  # Reference rounding: CVI at (cvi - 0.0005) min; CVA from QC mean/sd at 2 dp
  CVI_min  <- max(row$cvi - 0.0005, 0)
  mean_min <- row$qc_mean_r - 0.005
  mean_max <- row$qc_mean_r + 0.005
  sd_min   <- max(row$qc_sd_r - 0.005, 0)
  if (mean_min <= 0) stop("QC mean_min <= 0 for analyte; cannot form CVA bounds.")
  CVA_min <- (sd_min / mean_max) / sqrt(k)

  # Algebra on %CV scale
  CVW2_min <- CVA_min^2 + CVI_min^2
  CVG2_max <- max(CVT_use_max^2 - CVW2_min, 0)

  denom_total    <- CVG2_max + CVW2_min
  denom_analytic <- CVG2_max + CVA_min^2
  r_total_max    <- if (denom_total    == 0) 1 else as.numeric(CVG2_max / denom_total)
  r_analytic_max <- if (denom_analytic == 0) 1 else as.numeric(CVG2_max / denom_analytic)

  # Display CVs at original (not max) inputs
  m1_o <- r_m1$original; s1_o <- r_s1$original
  m2_o <- r_m2$original; s2_o <- r_s2$original
  CVT_pre_o  <- s1_o / m1_o
  CVT_post_o <- s2_o / m2_o
  CVT_used_o <- switch(use,
                       pre     = CVT_pre_o,
                       post    = CVT_post_o,
                       average = (CVT_pre_o + CVT_post_o)/2)

  CVA_disp <- row$cva / sqrt(k)
  CVI_disp <- row$cvi
  CVW_disp <- sqrt(CVA_disp^2 + CVI_disp^2)
  CVG_disp <- sqrt(pmax(CVT_used_o^2 - CVW_disp^2, 0))

  # Single-row data.frame (cohesive schema)
  data.frame(
    analyte        = as.character(row$analyte_raw),
    cvi            = as.numeric(CVI_disp),
    cva            = as.numeric(CVA_disp),
    cvt            = as.numeric(CVT_used_o),
    cvw            = as.numeric(CVW_disp),
    cvg            = as.numeric(CVG_disp),
    r_max_analytic = as.numeric(r_analytic_max),
    r_max_total    = as.numeric(r_total_max),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

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
#' @param out_cols Length-2 character vector: names for analytic-only and total ceilings.
#'   Defaults to \code{c("r_max_analytic","r_max_total")}.
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
                            use = c("pre","post","average"),
                            data = NULL,
                            append = FALSE,
                            out_cols = c("r_max_analytic","r_max_total")) {
  use <- match.arg(use)

  # helper to resolve from data or use scalar
  res_arg <- function(expr, nm) {
    if (missing(expr)) return(NULL)
    if (is.null(data)) return(expr)
    e <- substitute(expr)
    if (is.symbol(e)) {
      col <- deparse(e)
      if (!col %in% names(data)) stop("Column '", col, "' not found in `data` for ", nm, ".")
      return(data[[col]])
    }
    if (is.character(expr) && length(expr) == 1 && expr %in% names(data)) return(data[[expr]])
    expr
  }

  .strip_key <- function(x) tolower(gsub("\\s+", " ", trimws(gsub("[^[:alnum:] ]+"," ", x)) ))

  # --- per-row computation (uses your validate_descriptive) -------------------
  .amc_one <- function(analyte, m1, s1, m2, s2, k, use) {
    need_char <- list(m1 = m1, s1 = s1, m2 = m2, s2 = s2)
    not_char <- names(need_char)[!vapply(need_char, is.character, logical(1))]
    if (length(not_char)) {
      stop("For rounding-aware ceilings, pass m1, s1, m2, s2 as character strings: ",
           paste(not_char, collapse = ", "))
    }

    dat_raw <- MEPtoolkit::analytes
    nm <- tolower(gsub("\\s+", " ", names(dat_raw)))
    i_name <- (which(nm %in% c("test","analyte","name")))[1]
    i_cva  <- (grep("^cva$", nm))[1]; if (is.na(i_cva))  i_cva  <- (grep("^cva",  nm))[1]
    i_cvi  <- (grep("^cvi$", nm))[1]; if (is.na(i_cvi))  i_cvi  <- (grep("^cvi",  nm))[1]
    i_mean <- (grep("^mean$", nm))[1]
    i_sd   <- (grep("^sd$",   nm))[1]
    if (any(is.na(c(i_name, i_cva, i_cvi, i_mean, i_sd)))) {
      stop("MEPtoolkit::analytes must contain analyte name, cva, cvi, mean, sd.")
    }

    tbl <- data.frame(
      analyte_raw = dat_raw[[i_name]],
      analyte_key = .strip_key(dat_raw[[i_name]]),
      cva         = suppressWarnings(as.numeric(dat_raw[[i_cva]])),
      cvi         = suppressWarnings(as.numeric(dat_raw[[i_cvi]])),
      qc_mean_r   = suppressWarnings(as.numeric(dat_raw[[i_mean]])),
      qc_sd_r     = suppressWarnings(as.numeric(dat_raw[[i_sd]])),
      stringsAsFactors = FALSE
    )
    key <- .strip_key(analyte)
    hits <- which(tbl$analyte_key == key)
    if (!length(hits)) hits <- grep(key, tbl$analyte_key, fixed = TRUE)
    if (!length(hits)) stop("Analyte not found in MEPtoolkit::analytes: '", analyte, "'.")
    if (length(hits) > 1) hits <- hits[1]
    row <- tbl[hits, , drop = FALSE]

    r_m1 <- validate_descriptive(m1, "m1")
    r_s1 <- validate_descriptive(s1, "s1")
    r_m2 <- validate_descriptive(m2, "m2")
    r_s2 <- validate_descriptive(s2, "s2")

    s1_max <- pmax(r_s1$maximum, 0); s2_max <- pmax(r_s2$maximum, 0)
    m1_min <- r_m1$minimum;         m2_min <- r_m2$minimum
    if (m1_min <= 0 || m2_min <= 0) stop("Means must be positive; check inputs/rounding.")

    CVT_pre_max  <- s1_max / m1_min
    CVT_post_max <- s2_max / m2_min
    CVT_use_max  <- switch(use,
                           pre     = CVT_pre_max,
                           post    = CVT_post_max,
                           average = (CVT_pre_max + CVT_post_max)/2)

    CVI_min  <- max(row$cvi - 0.0005, 0)
    mean_max <- row$qc_mean_r + 0.005
    sd_min   <- max(row$qc_sd_r - 0.005, 0)
    if (mean_max <= 0) stop("QC mean <= 0; cannot form CVA bounds.")
    CVA_min  <- (sd_min / mean_max) / sqrt(k)

    CVW2_min <- CVA_min^2 + CVI_min^2
    CVG2_max <- max(CVT_use_max^2 - CVW2_min, 0)

    denom_total    <- CVG2_max + CVW2_min
    denom_analytic <- CVG2_max + CVA_min^2

    r_total_max    <- if (denom_total    == 0) 1 else as.numeric(CVG2_max / denom_total)
    r_analytic_max <- if (denom_analytic == 0) 1 else as.numeric(CVG2_max / denom_analytic)

    # display CVs at “reported” values (not edge-case maxima)
    m1o <- r_m1$original; s1o <- r_s1$original
    m2o <- r_m2$original; s2o <- r_s2$original
    CVT_used <- switch(use, pre = s1o/m1o, post = s2o/m2o, average = (s1o/m1o + s2o/m2o)/2)
    CVA_disp <- row$cva / sqrt(k); CVI_disp <- row$cvi
    CVW_disp <- sqrt(CVA_disp^2 + CVI_disp^2)
    CVG_disp <- sqrt(pmax(CVT_used^2 - CVW_disp^2, 0))

    list(
      analyte = row$analyte_raw,
      cvi = CVI_disp, cva = CVA_disp, cvt = CVT_used, cvw = CVW_disp, cvg = CVG_disp,
      rA  = r_analytic_max,
      rT  = r_total_max
    )
  }

  # ---- resolve inputs / vectorize ------------------------------------------------
  An <- res_arg(analyte, "analyte")
  M1 <- res_arg(m1, "m1"); S1 <- res_arg(s1, "s1")
  M2 <- res_arg(m2, "m2"); S2 <- res_arg(s2, "s2")

  lens <- vapply(list(An=An,M1=M1,S1=S1,M2=M2,S2=S2),
                 function(x) if (is.null(x)) 0L else length(x), 0L)
  L <- max(lens)
  if (L == 0L) stop("Provide either scalars or columns in `data` for analyte, m1, s1, m2, s2.")
  recyc <- function(x) if (is.null(x) || length(x) == L) x else if (length(x) == 1L) rep(x, L) else x
  An <- recyc(An); M1 <- recyc(M1); S1 <- recyc(S1); M2 <- recyc(M2); S2 <- recyc(S2)

  # compute rows
  a_disp <- character(L)
  cvi <- cva <- cvt <- cvw <- cvg <- rA <- rT <- rep(NA_real_, L)
  for (i in seq_len(L)) {
    res <- .amc_one(An[i], M1[i], S1[i], M2[i], S2[i], k = k, use = use)
    a_disp[i] <- res$analyte
    cvi[i] <- res$cvi; cva[i] <- res$cva; cvt[i] <- res$cvt
    cvw[i] <- res$cvw; cvg[i] <- res$cvg
    rA[i]  <- res$rA;  rT[i]  <- res$rT
  }

  # display table that the printer will use (same for single/multi)
  disp <- data.frame(
    case = seq_len(L),
    analyte = a_disp,
    cvi = as.numeric(cvi),
    cva = as.numeric(cva),
    cvt = as.numeric(cvt),
    cvw = as.numeric(cvw),
    cvg = as.numeric(cvg),
    r_max_analytic = as.numeric(rA),
    r_max_total    = as.numeric(rT),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (append && !is.null(data)) {
    if (!(length(out_cols) == 2 && all(nzchar(out_cols))))
      stop("`out_cols` must be length 2 with non-empty names.")
    out <- data
    out[[out_cols[1]]] <- rA
    out[[out_cols[2]]] <- rT
  } else {
    # return the display table itself when not appending
    out <- disp
  }

  # Attach metadata + the display table so the printer is uniform in all cases
  attr(out, "amc_meta") <- list(use = use, k = k)
  attr(out, "amc_disp") <- disp
  class(out) <- unique(c("analyte_max_cor", class(out)))
  out
}

#' Pretty print for analyte_max_cor (single or multi)
#'
#' Prints a detailed breakdown for a single-case object, or a compact table for
#' multi-case results showing analyte, CV columns (cvi, cva, cvt, cvw, cvg) and
#' the two ceilings.
#'
#' @param x An \code{analyte_max_cor} object (single) or a tagged multi-case
#'          data.frame returned by \code{analyte_max_cor(..., append = FALSE)}.
#' @param digits Digits for numeric columns in the output table (default 3).
#' @param max_rows Maximum number of rows to display for multi-case printing (default 20).
#' @param ... Unused.
#' @export
print.analyte_max_cor <- function(x, digits = 3, max_rows = 20, ...) {
  as_num <- function(v) suppressWarnings(as.numeric(v))
  fmt_num <- function(v, d = digits) {
    v <- as_num(v); ifelse(is.na(v), "NA", formatC(v, format = "f", digits = d))
  }
  fmt_pct <- function(v, d = 2) {
    v <- as_num(v); ifelse(is.na(v), "NA", paste0(formatC(100 * v, format = "f", digits = d), "%"))
  }

  meta <- attr(x, "amc_meta"); if (is.null(meta)) meta <- list(use = NA, k = NA)
  use <- meta$use; k <- meta$k

  disp <- attr(x, "amc_disp")
  # Fallback: try to reconstruct from x if attribute missing
  if (is.null(disp) && inherits(x, "data.frame")) {
    need <- c("case","analyte","cvi","cva","cvt","cvw","cvg","r_max_analytic","r_max_total")
    if (all(need %in% names(x))) disp <- x[, need, drop = FALSE]
  }
  if (is.null(disp)) {
    # Nothing to pretty-print; defer to default
    NextMethod()
    return(invisible(x))
  }

  n_all <- nrow(disp)
  show_idx <- seq_len(min(n_all, max_rows))
  tbl <- disp[show_idx, , drop = FALSE]

  # Coerce numerics and format
  pct_cols <- intersect(c("cvi","cva","cvt","cvw","cvg"), names(tbl))
  num_cols <- intersect(c("r_max_analytic","r_max_total"), names(tbl))
  for (nm in pct_cols) tbl[[nm]] <- fmt_pct(tbl[[nm]])
  for (nm in num_cols) tbl[[nm]] <- fmt_num(tbl[[nm]])

  # IMPORTANT: prevent recursive S3 dispatch
  class(tbl) <- "data.frame"
  attr(tbl, "amc_meta") <- NULL
  attr(tbl, "amc_disp") <- NULL

  cat("\n# Maximum pre-post correlation ceilings\n", sep = "")
  cat("use = ", use, ", k = ", k, "\n\n", sep = "")
  print(tbl, row.names = FALSE, right = TRUE)
  if (n_all > length(show_idx)) {
    cat("\n... ", n_all - length(show_idx), " more rows not shown (increase max_rows to see more).\n", sep = "")
  }
  invisible(x)
}


