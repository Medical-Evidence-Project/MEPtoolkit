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
#' @param out_col_min,out_col_mid,out_col_max Character. Output column names when `append = TRUE`
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
#' @seealso Paired \emph{t} relationships; summary-data meta-analysis for pre–post designs.
#' @importFrom stats qt
#' @rdname prepost_r
#' @export
prepost_r <- function(s1, s2, sc = NULL, mc = NULL,
                      m1 = NULL, m2 = NULL, n = NULL,
                      t = NULL, p = NULL,
                      data = NULL, append = FALSE,
                      out_col_min = "rho_prepost_min",
                      out_col_mid = "rho_prepost",
                      out_col_max = "rho_prepost_max",
                      strict = FALSE, warn = TRUE) {

  # -- capture expressions so bare names work with `data`
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

  # -- enforce character inputs (for rounding inference), except n
  need_char <- list(s1=s1, s2=s2, sc=sc, mc=mc, m1=m1, m2=m2, t=t, p=p)
  offenders <- names(need_char)[vapply(need_char, function(z) !is.null(z) && !is.character(z), logical(1))]
  if (length(offenders)) {
    stop("For rounding-aware bounds, these inputs must be character strings: ",
         paste(offenders, collapse = ", "), call. = FALSE)
  }
  if (is.null(n)) stop("`n` (paired sample size) is required.", call. = FALSE)
  n <- as.numeric(n)
  if (any(!is.finite(n)) || any(n <= 1)) stop("`n` must be finite and > 1.", call. = FALSE)

  # -- helpers using your validators (half-ULP bounds from decimals)
  rng_desc <- function(x, name) {
    r <- validate_descriptive(x, name)  # list(minimum, original, maximum) using 0.5 * 10^(-d)
    c(min = r$minimum, max = r$maximum, orig = r$original)
  }
  rng_sd <- function(x, name) { z <- rng_desc(x, name); z[1] <- max(0, z[1]); z }  # SDs ≥ 0
  rng_p  <- function(x) { r <- validate_p(x, "p"); c(min = r$minimum, max = r$maximum, orig = r$original) }

  # -- vectorization & recycling
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

  rho_min <- rho_mid <- rho_max <- rep(NA_real_, L)
  missing_rows <- integer(0)

  # safe corner correlation; returns NA if S1<=0 or S2<=0 or any not finite
  r_corner <- function(S1, S2, SC) {
    if (!is.finite(S1) || !is.finite(S2) || !is.finite(SC) || S1 <= 0 || S2 <= 0) return(NA_real_)
    (S1^2 + S2^2 - SC^2) / (2 * S1 * S2)
  }

  for (i in seq_len(L)) {
    # intervals for s1, s2
    if (is.null(s1) || is.null(s2) || is.na(s1[i]) || is.na(s2[i])) {
      missing_rows <- c(missing_rows, i); next
    }
    s1rng <- rng_sd(s1[i], "s1"); s2rng <- rng_sd(s2[i], "s2")
    s1min <- s1rng["min"]; s1max <- s1rng["max"]; s1orig <- s1rng["orig"]
    s2min <- s2rng["min"]; s2max <- s2rng["max"]; s2orig <- s2rng["orig"]

    # derive sc interval & orig by priority sc > t > p
    sc_has <- (!is.null(sc) && !is.na(sc[i]))
    t_has  <- (!is.null(t)  && !is.na(t[i]))
    p_has  <- (!is.null(p)  && !is.na(p[i]))

    sc_min <- sc_max <- sc_orig <- NA_real_

    if (sc_has) {
      scrng <- rng_sd(sc[i], "sc")
      sc_min <- scrng["min"]; sc_max <- scrng["max"]; sc_orig <- scrng["orig"]

    } else {
      # mc bounds & orig
      if (!is.null(mc) && !is.na(mc[i])) {
        mcrng <- rng_desc(mc[i], "mc")
        mc_min <- mcrng["min"]; mc_max <- mcrng["max"]; mc_orig <- mcrng["orig"]
      } else if (!is.null(m1) && !is.na(m1[i]) && !is.null(m2) && !is.na(m2[i])) {
        m1rng <- rng_desc(m1[i], "m1"); m2rng <- rng_desc(m2[i], "m2")
        mc_min <- m2rng["min"] - m1rng["max"]
        mc_max <- m2rng["max"] - m1rng["min"]
        mc_orig <- m2rng["orig"] - m1rng["orig"]
      } else {
        missing_rows <- c(missing_rows, i); next
      }

      # absolute change bounds (use 0 if interval crosses zero)
      mc_abs_min  <- if (mc_min <= 0 && 0 <= mc_max) 0 else min(abs(mc_min), abs(mc_max))
      mc_abs_max  <- max(abs(mc_min), abs(mc_max))
      mc_abs_orig <- abs(mc_orig)

      if (t_has) {
        trng <- rng_desc(t[i], "t")
        t_abs_min  <- min(abs(trng["min"]), abs(trng["max"]))
        t_abs_max  <- max(abs(trng["min"]), abs(trng["max"]))
        t_abs_orig <- abs(trng["orig"])
        if (!is.finite(t_abs_min) || t_abs_min <= 0) {
          sc_min <- 0; sc_max <- Inf; sc_orig <- NA_real_
        } else {
          sroot  <- sqrt(n[i])
          sc_min <- mc_abs_min  * sroot / t_abs_max  # smallest |mc| with largest |t| (tightest)
          sc_max <- mc_abs_max  * sroot / t_abs_min  # largest  |mc| with smallest |t| (loosest)
          sc_orig <- if (is.finite(t_abs_orig) && t_abs_orig > 0) mc_abs_orig * sroot / t_abs_orig else NA_real_
        }
      } else if (p_has) {
        prng <- rng_p(p[i])
        t_from_p <- function(pp) abs(stats::qt(pp/2, df = n[i] - 1, lower.tail = FALSE))
        t_min  <- t_from_p(prng["max"])  # largest p -> smallest |t|
        t_max  <- t_from_p(prng["min"])  # smallest p -> largest |t|
        t_orig <- t_from_p(prng["orig"])
        if (!is.finite(t_min) || t_min <= 0) {
          sc_min <- 0; sc_max <- Inf; sc_orig <- NA_real_
        } else {
          sroot  <- sqrt(n[i])
          sc_min <- mc_abs_min  * sroot / t_max
          sc_max <- mc_abs_max  * sroot / t_min
          sc_orig <- if (is.finite(t_orig) && t_orig > 0) mc_abs_orig * sroot / t_orig else NA_real_
        }
      } else {
        missing_rows <- c(missing_rows, i); next
      }
    }

    if (!is.finite(sc_min) || sc_min < 0) sc_min <- 0
    if (!is.finite(sc_max) || sc_max < 0) sc_max <- Inf

    # Evaluate corners for bounds; ignore non-finite candidates
    r_max_cands <- c(
      r_corner(s1min, s2min, sc_min),
      r_corner(s1min, s2max, sc_min),
      r_corner(s1max, s2min, sc_min),
      r_corner(s1max, s2max, sc_min)
    )
    r_min_cands <- c(
      r_corner(s1min, s2min, sc_max),
      r_corner(s1min, s2max, sc_max),
      r_corner(s1max, s2min, sc_max),
      r_corner(s1max, s2max, sc_max)
    )

    # clamp candidates then select
    r_max_cands <- pmin(pmax(r_max_cands, -1), 1)
    r_min_cands <- pmin(pmax(r_min_cands, -1), 1)

    rho_max[i] <- if (all(!is.finite(r_max_cands))) NA_real_ else max(r_max_cands, na.rm = TRUE)
    rho_min[i] <- if (all(!is.finite(r_min_cands))) NA_real_ else min(r_min_cands, na.rm = TRUE)

    # Point estimate at reported numbers (RIVET)
    rho_mid[i] <- r_corner(s1orig, s2orig, sc_orig)
    if (is.finite(rho_mid[i])) rho_mid[i] <- pmin(pmax(rho_mid[i], -1), 1) else rho_mid[i] <- NA_real_
  }

  if (length(missing_rows)) {
    if (strict) {
      stop(sprintf("Insufficient information for rows: %s.", paste(missing_rows, collapse = ", ")), call. = FALSE)
    } else if (warn) {
      warning(sprintf("Insufficient information for rows set to NA: %s.", paste(missing_rows, collapse = ", ")), call. = FALSE)
    }
  }

  out_df <- data.frame(
    rho_prepost_min = rho_min,
    rho_prepost     = rho_mid,
    rho_prepost_max = rho_max
  )

  if (append && !is.null(data)) {
    if (nrow(out_df) == 1L && nrow(data) > 1L) out_df <- out_df[rep(1, nrow(data)), , drop = FALSE]
    if (nrow(out_df) != nrow(data)) stop("Output length does not match nrow(data).", call. = FALSE)
    data[[out_col_min]] <- out_df$rho_prepost_min
    data[[out_col_mid]] <- out_df$rho_prepost
    data[[out_col_max]] <- out_df$rho_prepost_max
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




