#' Compute a paired (pre–post) correlation from summary statistics
#'
#' @description
#' `prepost_r()` computes \emph{bounds} on the paired correlation (pre vs. post)
#' attributable to rounding in the reported summary statistics, and also returns
#' the RIVET point estimate at the reported numbers. It accepts SD of change
#' directly or derives it from a paired \emph{t} or two-sided \emph{p} plus mean change and
#' sample size. Character inputs (e.g., `"37.74"`) allow the function to infer
#' rounding half-ULPs from the number of decimals.
#'
#' `prepost_r_RIVETS()` provides the legacy single-point estimate (no rounding
#' bounds). Both functions share this help page.
#'
#' @details
#' Correlation formula:
#' \deqn{r = \frac{s_1^2 + s_2^2 - s_c^2}{2 s_1 s_2}.}
#' If \eqn{s_c} is missing, it is derived per row (priority \strong{sc > t > p}):
#' \itemize{
#'   \item From paired \eqn{t}: \eqn{s_c = m_c \sqrt{n}/t}.
#'   \item From two-sided paired \eqn{p}: \eqn{t = |qt(p/2, \nu=n-1, lower.tail=FALSE)|}, then \eqn{s_c = m_c \sqrt{n}/t}.
#' }
#'
#' \strong{Rounding model (for `prepost_r()`):}
#' A character value with \eqn{d} decimals is taken to represent an interval
#' \eqn{\pm 0.5 \times 10^{-d}} around the reported number. We propagate these
#' intervals to obtain \code{rho_prepost_min} (worst case: large \code{sc}, corner
#' \code{s1}/\code{s2}) and \code{rho_prepost_max} (best case: small \code{sc}, corner
#' \code{s1}/\code{s2}). We also report the point estimate \code{rho_prepost} at
#' the reported numbers (RIVET).
#'
#' \emph{Vectorization & recycling.} Inputs are evaluated per row. Length-1 inputs
#' recycle; other mismatches error. With `data`, arguments can be bare or quoted
#' column names.
#'
#' \emph{Handling incomplete rows.} With `strict = FALSE` (default), rows lacking
#' sufficient information yield `NA` (and optionally warn). With `strict = TRUE`,
#' the function errors and lists offending rows.
#'
#' \emph{Append behavior.} If `append = TRUE` and `data` is supplied, the function
#' returns `data` with three new columns (defaults: `"rho_prepost_min"`,
#' `"rho_prepost"`, `"rho_prepost_max"`); otherwise it returns a 3-column data frame.
#'
#' @param s1,s2 Character vectors or column names. SD at pre and post. \strong{Must be character} in `prepost_r()` to infer rounding. Numeric allowed in `_RIVET`.
#' @param sc Character or column name. SD of change (if supplied). Character in `prepost_r()`; numeric allowed in `_RIVET`.
#' @param mc Character or column name. Mean change \eqn{(m_2 - m_1)}. If `NULL` and both `m1` and `m2` are supplied, it’s computed row-wise.
#' @param m1,m2 Character or column names. Means at pre and post (used for `mc` and/or deriving `sc`). Character in `prepost_r()`; numeric allowed in `_RIVET`.
#' @param n  Numeric/integer or column name. Paired sample size (treated as exact).
#' @param t  Character or column name. Paired \emph{t} statistic (used when `sc` missing). Character in `prepost_r()`; numeric allowed in `_RIVET`.
#' @param p  Character or column name. Two-sided \emph{p} (used when `sc` and `t` missing). Character in `prepost_r()`; numeric allowed in `_RIVET`.
#' @param data Optional `data.frame`. If provided, arguments are evaluated within it.
#' @param append Logical. If `TRUE` with `data`, append output columns; else return a data frame of results.
#' @param out_cols Character vector. Output column names when `append = TRUE`
#'   (defaults: `"rho_prepost_min"`, `"rho_prepost"`, `"rho_prepost_max"`).
#' @param out_col (legacy, `_RIVET` only) Character. Single output column name (default `"rho_prepost"`).
#' @param strict Logical. If `TRUE`, error out on underspecified rows; else return `NA` (default `FALSE`).
#' @param warn Logical. Warn about underspecified rows when `strict = FALSE` (default `TRUE`).
#'
#' @return
#' - `prepost_r()`: If `append = FALSE`, a data frame with three columns:
#'   `rho_prepost_min`, `rho_prepost`, `rho_prepost_max`. If `append = TRUE` and `data` is given,
#'   returns `data` with those columns appended.
#' - `prepost_r_RIVETS()`: A numeric vector (or appended single column) of point estimates.
#'
#'
#' @examples
#' df <- data.frame(
#'   m1 = c("140.2","6.60"),
#'   s1 = c("2.10","0.90"),
#'   m2 = c("139.8","6.40"),
#'   s2 = c("2.05","0.90"),
#'   t  = c("3.00","2.40"),
#'   n  = c(60, 45),
#'   stringsAsFactors = FALSE
#' )
#'
#' prepost_r(
#'   data = df,
#'   m1 = m1, s1 = s1,
#'   m2 = m2, s2 = s2,
#'   n = n, t = t
#' )
#'
#'
#'
#' @seealso Paired \emph{t} relationships; summary-data meta-analysis for
#' pre–post designs.
#' @importFrom stats qt
#' @import dplyr rlang
#' @rdname prepost_r
#' @export
prepost_r <- function(s1, s2, sc = NULL, mc = NULL,
                      m1 = NULL, m2 = NULL, n = NULL,
                      t = NULL, p = NULL,
                      data = NULL, append = FALSE,
                      out_cols = c("rho_prepost_min",
                                   "rho_prepost",
                                   "rho_prepost_max"),
                      strict = FALSE, warn = TRUE) {

  if(!is.null(data)) {
    if(!is.data.frame(data)) stop("`data` must be a data.frame.")
    sc <- resolve(enquo(sc), data); s1 <- resolve(enquo(s1), data);
    s2 <- resolve(enquo(s2), data); mc <- resolve(enquo(mc), data);
    m1 <- resolve(enquo(m1), data); m2 <- resolve(enquo(m2), data);
    n  <- resolve(enquo(n), data);  t  <- resolve(enquo(t), data);
    p  <- resolve(enquo(p), data)
  }

  if (is.null(n) || any(is.na(n)))
    stop("`n` (paired sample size) is required for every row.")
  n <- as.numeric(n)
  if (any(n <= 1))
    stop("`n` must be >1.")

  if(is.null(s1) || is.null(s2))
    stop("`s1` and `s2` must be non-null.")

  # -- check lengths
  args <- list(s1=s1, s2=s2, sc=sc, mc=mc, m1=m1, m2=m2, n=n, t=t, p=p)
  lens <- vapply(args, length, integer(1))
  L <- max(lens)
  ok_len <- lens == 0 | lens == L
  if (!all(ok_len))
    stop("Non-NULL arguments must have the same length.")

  s1 <- validate_descriptive(s1, allow_na = TRUE)
  s2 <- validate_descriptive(s2, allow_na = TRUE)
  sc <- if(!is.null(sc)) validate_descriptive(sc, allow_na = TRUE) else NULL
  m1 <- if(!is.null(m1)) validate_descriptive(m1, allow_na = TRUE) else NULL
  m2 <- if(!is.null(m2)) validate_descriptive(m2, allow_na = TRUE) else NULL
  mc <- if(!is.null(mc)) validate_descriptive(mc, allow_na = TRUE) else NULL
  t <- if(!is.null(t)) validate_descriptive(t, allow_na = TRUE) else NULL
  p <- if(!is.null(p)) validate_p(p, allow_na = TRUE) else NULL

  rho_min <- rho_mid <- rho_max <- rep(NA_real_, L)
  missing_rows <- integer(0)

  for (i in seq_len(L)) {
    # intervals for s1, s2
    if (is.na(s1[["original"]][i]) || is.na(s2[["original"]][i])) {
      missing_rows <- c(missing_rows, i); next
    }
    else if(s1[["original"]][i]<=0 || s2[["original"]][i]<=0)
      stop(paste0("SDs must be >0 (line ", i, ")."))
    s1min <- s1[["minimum"]][i]
    s1max <- s1[["maximum"]][i]
    s1orig <- s1[["original"]][i]
    s2min <- s2[["minimum"]][i]
    s2max <- s2[["maximum"]][i]
    s2orig <- s2[["original"]][i]

    # derive sc interval & orig by priority sc > t > p
    sc_has <- (!is.null(sc) && !is.na(sc[["original"]][i]))
    t_has  <- (!is.null(t)  && !is.na(t[["original"]][i]))
    p_has  <- (!is.null(p)  && !is.na(p[["original"]][i]))

    sc_min <- sc_max <- sc_orig <- NA_real_

    if(sc_has)
    {
      sc_min <- sc[["minimum"]][i]
      sc_max <- sc[["maximum"]][i]
      sc_orig <- sc[["original"]][i]
    }

    ### Think about the possibility that using m1 and m2 might be more precise
    ### Than using mc?
    if(!sc_has) {
      # mc bounds & orig
      if(is.null(mc) || is.na(mc[["original"]][i])) {
        if (!is.null(m1) && !is.na(m1[["original"]][i]) &&
            !is.null(m2) && !is.na(m2[["original"]][i])) {
          mc_min <- m2[["minimum"]][i] - m1[["maximum"]][i]
          mc_max <- m2[["maximum"]][i] - m1[["minimum"]][i]
          mc_orig <- m2[["original"]][i] - m1[["original"]][i]
        } else {
          missing_rows <- c(missing_rows, i); next
        }
      } else {
        mc_min <- mc[["minimum"]][i]
        mc_max <- mc[["maximum"]][i]
        mc_orig <- mc[["original"]][i]
      }
      # absolute change bounds (use 0 if interval crosses zero)
      mc_abs_min  <- if (mc_min <= 0 && 0 <= mc_max) 0
                     else min(abs(mc_min), abs(mc_max))
      mc_abs_max  <- max(abs(mc_min), abs(mc_max))
      mc_abs_orig <- abs(mc_orig)

      if (t_has) {
        t_abs_min  <- min(abs(t[["minimum"]][i]), abs(t[["maximum"]][i]))
        t_abs_max  <- max(abs(t[["minimum"]][i]), abs(t[["maximum"]][i]))
        t_abs_orig <- abs(t[["original"]][i])
        if(t_abs_orig==0)
          stop(paste0("t must be different from 0 (line ", i, ")."))
        sroot  <- sqrt(n[i])
        # smallest |mc| with largest |t| (tightest)
        sc_min <- mc_abs_min  * sroot / t_abs_max
        # largest  |mc| with smallest |t| (loosest)
        sc_max <- mc_abs_max  * sroot / t_abs_min
        sc_orig <- mc_abs_orig * sroot / t_abs_orig

      } else if (p_has) {
        if(p[["original"]][i]==0 || p[["original"]][i]==1)
          stop(paste0("p must be >0 and <1 (line ", i, ")."))
        t_from_p <- function(pp) abs(stats::qt(pp/2, df = n[i] - 1,
                                               lower.tail = FALSE))
        t_min  <- t_from_p(p[["maximum"]][i])  # largest p -> smallest |t|
        t_max  <- t_from_p(p[["minimum"]][i])  # smallest p -> largest |t|
        t_orig <- t_from_p(p[["original"]][i])
        sroot  <- sqrt(n[i])
        sc_min <- mc_abs_min  * sroot / t_max
        sc_max <- mc_abs_max  * sroot / t_min
        sc_orig <- mc_abs_orig * sroot / t_orig
      } else {
        missing_rows <- c(missing_rows, i); next
      }
    }

    # Calculate rhos
    rho_calculation <-
      function(S1, S2, SC) {(S1^2 + S2^2 - SC^2) / (2 * S1 * S2)}
    rho_cands <- c()
    cases_grid <- expand.grid(s1_range = c(s1min, s1max),
                              s2_range = c(s2min, s2max),
                              sc_range = c(sc_min, sc_max))

    for(j in seq_len(nrow(cases_grid))) {
      rho_cands <- c(rho_cands, rho_calculation(cases_grid[j, "s1_range"],
                                                cases_grid[j, "s2_range"],
                                                cases_grid[j, "sc_range"]))
    }

    rho_max[i] <- max(rho_cands)
    rho_min[i] <- min(rho_cands)
    rho_mid[i] <- rho_calculation(s1orig, s2orig, sc_orig)
  }

  if (length(missing_rows)) {
    if (strict) {
      stop(sprintf("Insufficient information for rows: %s.",
                   paste(missing_rows, collapse = ", ")))
    } else if (warn) {
      warning(sprintf("Insufficient information for rows set to NA: %s.",
                      paste(missing_rows, collapse = ", ")))
    }
  }

  out_df <- data.frame(
    rho_prepost_min = rho_min,
    rho_prepost     = rho_mid,
    rho_prepost_max = rho_max
  )

  if (append && !is.null(data)) {
    data[[out_cols[1]]] <- out_df$rho_prepost_min
    data[[out_cols[2]]] <- out_df$rho_prepost
    data[[out_cols[3]]] <- out_df$rho_prepost_max
    return(data)
  }
  out_df
}

#' @rdname prepost_r
#' @export
prepost_r_RIVETS <- function(s1, s2, sc = NULL, mc = NULL,
                            m1 = NULL, m2 = NULL, n = NULL,
                            t = NULL, p = NULL,
                            data = NULL, append = FALSE, out_col = "rho_prepost",
                            strict = FALSE, warn = TRUE) {

  # capture & resolve
  s1x <- substitute(s1); s2x <- substitute(s2); scx <- substitute(sc)
  mcx <- substitute(mc); m1x <- substitute(m1); m2x <- substitute(m2)
  nx  <- substitute(n);   tx <- substitute(t);   px <- substitute(p)
  resolve <- function(expr) {
    if (is.null(expr)) return(NULL)
    if (!is.null(data)) {
      if (!is.data.frame(data)) stop("`data` must be a data.frame.", call. = FALSE)
      return(eval(expr, envir = data, enclos = parent.frame()))
    }
    eval(expr, envir = parent.frame())
  }
  s1 <- resolve(s1x); s2 <- resolve(s2x); sc <- resolve(scx); mc <- resolve(mcx)
  m1 <- resolve(m1x); m2 <- resolve(m2x); n  <- resolve(nx);  t  <- resolve(tx);  p  <- resolve(px)

  # vectorization & recycling
  args <- list(s1=s1, s2=s2, sc=sc, mc=mc, m1=m1, m2=m2, n=n, t=t, p=p)
  lens <- vapply(Filter(Negate(is.null), args), length, integer(1))
  L <- if (length(lens)) max(lens) else 1L
  ok_len <- function(v) (is.null(v) || length(v) == 1L || length(v) == L)
  bad <- names(args)[!vapply(args, ok_len, logical(1))]
  if (length(bad)) stop("Non-NULL args must have length 1 or ", L, ". Offenders: ",
                        paste(bad, collapse = ", "), call. = FALSE)
  recyc <- function(x) if (is.null(x) || length(x) == L) x else rep(x, L)
  s1 <- recyc(s1); s2 <- recyc(s2); sc <- recyc(sc); mc <- recyc(mc)
  m1 <- recyc(m1); m2 <- recyc(m2); n  <- recyc(n);  t  <- recyc(t);  p  <- recyc(p)

  # coerce to numeric (legacy behavior)
  as_num <- function(x) if (is.null(x)) x else as.numeric(x)
  s1 <- as_num(s1); s2 <- as_num(s2); sc <- as_num(sc); mc <- as_num(mc)
  m1 <- as_num(m1); m2 <- as_num(m2); n  <- as_num(n);  t  <- as_num(t);  p  <- as_num(p)

  if (is.null(n)) stop("`n` (paired sample size) is required.")
  if (any(!is.finite(n)) || any(n <= 1)) stop("`n` must be finite and > 1.", call. = FALSE)

  if (is.null(mc)) mc <- rep(NA_real_, L)
  need_mc <- is.na(mc) & !is.null(m1) & !is.null(m2)
  if (any(need_mc)) mc[need_mc] <- (m2 - m1)[need_mc]

  sc_final <- sc; if (is.null(sc_final)) sc_final <- rep(NA_real_, L)
  use_t <- is.na(sc_final) & !is.null(t)
  if (any(use_t, na.rm = TRUE)) {
    ix <- which(use_t & !is.na(t) & !is.na(n) & !is.na(mc))
    if (length(ix)) sc_final[ix] <- (mc[ix] * sqrt(n[ix])) / t[ix]
  }
  use_p <- is.na(sc_final) & !is.null(p)
  if (any(use_p, na.rm = TRUE)) {
    ix <- which(use_p & !is.na(p) & !is.na(n) & !is.na(mc))
    if (length(ix)) {
      t_from_p <- abs(stats::qt(p[ix]/2, df = n[ix] - 1, lower.tail = FALSE))
      sc_final[ix] <- (mc[ix] * sqrt(n[ix])) / t_from_p
    }
  }

  missing_rows <- which(is.na(sc_final) | is.na(s1) | is.na(s2))
  if (length(missing_rows)) {
    if (strict) {
      stop(sprintf("Insufficient information for rows: %s.", paste(missing_rows, collapse = ", ")), call. = FALSE)
    } else if (warn) {
      warning(sprintf("Insufficient information for rows set to NA: %s.", paste(missing_rows, collapse = ", ")), call. = FALSE)
    }
  }

  rho <- (s1^2 + s2^2 - sc_final^2) / (2 * s1 * s2)

  if (append && !is.null(data)) {
    if (length(rho) == 1L && nrow(data) > 1L) rho <- rep(rho, nrow(data))
    if (length(rho) != nrow(data)) {
      stop(sprintf("Length of '%s' (%d) does not match nrow(data) (%d).", out_col, length(rho), nrow(data)), call. = FALSE)
    }
    out <- data
    out[[out_col]] <- rho
    return(out)
  }

  rho
}
