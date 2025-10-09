#' Flag and quantify pre–post correlation vs analyte ceilings
#'
#' For each row/case this function:
#' (1) computes rounding-aware ceilings with `analyte_max_cor()` (analytic-only and
#'     analytic+biological);
#' (2) obtains the paired correlation via `prepost_r()` and uses its LOWER rounding
#'     bound as the observed correlation `r_obs` used for flagging and p-values
#'     (most generous to the study);
#'     if `prepost_r()` cannot be computed and a reported correlation `r_obs`
#'     (character) is supplied, its LOWER rounding bound is used;
#' (3) computes one-tailed probabilities that a correlation >= `r_obs` would be seen
#'     if the true correlation equaled each ceiling (Fisher z, n > 3).
#'
#' Rounding-aware inputs: `m1, s1, m2, s2` must be character strings (e.g., "64.3")
#' so decimals encode rounding, consistent with `analyte_max_cor()`.
#'
#' @param data Optional data.frame with columns; you may pass bare/quoted column names.
#' @param analyte,m1,s1,m2,s2 Bare/quoted names in `data` or vectors. For rounding-aware
#'        ceilings, these must be character strings (e.g., `m1 = "142.23"`).
#' @param n Column name or numeric vector. Paired sample size (exact).
#' @param r_obs Optional reported correlation as character (e.g., `".980"`). Only used
#'        as a fallback when `prepost_r()` cannot be computed; its lower rounding bound
#'        becomes the numeric `r_obs` used here.
#' @param sc,mc,t,p Optional inputs forwarded to `prepost_r()` (priority `sc > t > p`);
#'        may be columns in `data`. If numeric, coerced to character.
#' @param k Integer (>=1). Replicates averaged (passed to `analyte_max_cor()`).
#' @param use Which study CVT to use in the ceiling calculation: `"pre"`, `"post"`,
#'        or `"average"`.
#' @param append Logical. If TRUE (default) and `data` is supplied, return `data`
#'        with new columns appended; else return only the new columns (and include
#'        `analyte` and `n_used` so printing works).
#'
#' @return A data.frame containing (or appending):
#' \itemize{
#'   \item `n_used` (the n actually used)
#'   \item `r_max_analytic_expected`, `r_max_total_expected`
#'   \item `r_obs` (numeric; lower rounding-bound of the observed correlation)
#'   \item `exceeds_max_analytic`, `exceeds_max_total`
#'   \item `prob_analytic`, `prob_total` (one-tailed Fisher-z at `r_obs`)
#' }
#'
#'
#' @examples
#' # --- Example 1: data.frame input; r derived via paired t (rounding-aware) ---
#' df <- data.frame(
#'   analyte = c("Sodium","Total Cholesterol"),
#'   m1 = c("140.2","6.60"),
#'   s1 = c("2.10","0.90"),
#'   m2 = c("139.8","6.40"),
#'   s2 = c("2.05","0.90"),
#'   t  = c("3.00","2.40"),
#'   n  = c(60, 45),
#'   stringsAsFactors = FALSE
#' )
#' out1 <- analyte_cor_check(
#'   data = df,
#'   analyte = analyte, m1 = m1, s1 = s1, m2 = m2, s2 = s2,
#'   n = n, t = t,
#'   k = 5, use = "pre", append = TRUE
#' )
#' print(out1, max_rows = 10)
#'
#' # --- Example 2: scalar inputs (no data.frame), r derived via t ---
#' out2 <- analyte_cor_check(
#'   analyte = "Sodium",
#'   m1 = "140.2", s1 = "2.10",
#'   m2 = "139.8", s2 = "2.05",
#'   n = 60, t = "2.50",
#'   k = 1, use = "pre", append = FALSE
#' )
#' print(out2)
#'
#' # --- Example 3: reported r only (fallback).
#' # If prepost_r() cannot be computed from the provided stats, the LOWER
#' # rounding bound of r_obs (character) is used for flagging and p-values.
#' out3 <- analyte_cor_check(
#'   analyte = "Total Cholesterol",
#'   m1 = "6.60", s1 = "0.90",
#'   m2 = "6.40", s2 = "0.90",
#'   n = 45,
#'   r_obs = ".980",              # character; uses its lower rounding bound
#'   k = 1, use = "pre", append = FALSE
#' )
#' print(out3)
#'
#' # --- Example 4: r derived via two-sided p (when t/sc not available) ---
#' out4 <- analyte_cor_check(
#'   analyte = "Glucose",
#'   m1 = "92.1",  s1 = "9.8",
#'   m2 = "88.4",  s2 = "9.1",
#'   n = 52,
#'   p = "0.012",                    # two-sided p as character
#'   k = 1, use = "pre", append = FALSE
#' )
#' print(out4)
#'
#'
#'
#' @importFrom stats pnorm qt
#' @export
# analyte_cor_check
analyte_cor_check <- function(
    data = NULL,
    analyte, m1, s1, m2, s2, n,
    r_obs = NULL,            # optional reported correlation (character)
    sc = NULL, mc = NULL,    # optional for prepost_r (only used when r_obs is missing)
    t  = NULL, p  = NULL,    # optional for prepost_r (only used when r_obs is missing)
    k = 1,
    use = c("pre","post","average"),
    append = TRUE
) {
  use <- match.arg(use)

  # --- capture call & resolver so bare names work with `data`
  mcall <- match.call(expand.dots = FALSE)
  arg_present <- function(name) name %in% names(mcall)
  resolve <- function(expr, nm, present) {
    if (!present) return(NULL)
    if (!is.null(data)) {
      e <- substitute(expr)
      if (is.symbol(e)) {
        col <- deparse(e)
        if (!col %in% names(data)) stop("Column '", col, "' not found in `data` for ", nm, ".")
        return(data[[col]])
      }
      if (is.character(expr) && length(expr) == 1L && expr %in% names(data)) return(data[[expr]])
    }
    expr
  }

  An <- resolve(analyte, "analyte", arg_present("analyte"))
  M1 <- resolve(m1,      "m1",      arg_present("m1"))
  S1 <- resolve(s1,      "s1",      arg_present("s1"))
  M2 <- resolve(m2,      "m2",      arg_present("m2"))
  S2 <- resolve(s2,      "s2",      arg_present("s2"))
  N  <- resolve(n,       "n",       arg_present("n"))
  R  <- resolve(r_obs,   "r_obs",   arg_present("r_obs"))
  SC <- resolve(sc,      "sc",      arg_present("sc"))
  MC <- resolve(mc,      "mc",      arg_present("mc"))
  TT <- resolve(t,       "t",       arg_present("t"))
  PP <- resolve(p,       "p",       arg_present("p"))

  # --- lengths + recycling
  lens <- vapply(
    list(An=An,M1=M1,S1=S1,M2=M2,S2=S2,N=N,R=R,SC=SC,MC=MC,TT=TT,PP=PP),
    function(x) if (is.null(x)) 0L else length(x), 0L
  )
  L <- max(lens)
  if (L == 0L) stop("Provide analyte, m1, s1, m2, s2, n (scalars or columns in `data`).")

  recyc <- function(x) if (is.null(x) || length(x) == L) x else if (length(x) == 1L) rep(x, L) else x
  An <- recyc(An); M1 <- recyc(M1); S1 <- recyc(S1); M2 <- recyc(M2); S2 <- recyc(S2)
  N  <- recyc(N);  R  <- recyc(R);  SC <- recyc(SC);  MC <- recyc(MC);  TT <- recyc(TT);  PP <- recyc(PP)

  # --- must-haves / types
  if (is.null(N) || any(is.na(N))) stop("Sample size `n` must be provided.")
  must_char <- c(
    if (!is.character(M1)) "m1",
    if (!is.character(S1)) "s1",
    if (!is.character(M2)) "m2",
    if (!is.character(S2)) "s2"
  )
  if (length(must_char)) {
    stop("m1/s1/m2/s2 must be character per row (e.g., \"64.3\"): ",
         paste(must_char, collapse = ", "))
  }

  # --- generous lower bound for r_obs using validate_descriptive()
  r_obs_lower_from_validate <- function(x) {
    if (is.null(x) || is.na(x)) return(NA_real_)
    s <- trimws(as.character(x))
    s <- sub("^([-+]?)\\.", "\\10.", s)  # ".450" -> "0.450"
    if (!grepl("^[+-]?[0-9]*\\.?[0-9]+$", s))
      stop("r_obs must be a simple number string like \".980\" or \"0.61\".")
    rng <- validate_descriptive(s, "r_obs")
    val <- max(-1, min(1, rng$minimum))
    as.numeric(val)
  }

  # --- compact df for prepost_r() if we need it
  to_chr <- function(x) if (!is.null(x) && !is.character(x)) as.character(x) else x
  SC <- to_chr(SC); MC <- to_chr(MC); TT <- to_chr(TT); PP <- to_chr(PP)
  work_list <- list(s1 = S1, s2 = S2, m1 = M1, m2 = M2, n = N)
  if (!is.null(SC)) work_list$sc <- SC
  if (!is.null(MC)) work_list$mc <- MC
  if (!is.null(TT)) work_list$t  <- TT
  if (!is.null(PP)) work_list$p  <- PP
  work <- do.call(data.frame, c(work_list, stringsAsFactors = FALSE))

  # --- observed r: prefer supplied r_obs; else prepost_r() min bound
  have_R <- !is.null(R)
  if (have_R) {
    r_obs_num <- vapply(R, r_obs_lower_from_validate, numeric(1))
  } else {
    pr <- try(
      prepost_r(
        s1 = s1, s2 = s2, sc = sc, mc = mc,
        m1 = m1, m2 = m2, n = n, t = t, p = p,
        data = work, append = TRUE
      ),
      silent = TRUE
    )
    have_pr <- !inherits(pr, "try-error") && ("rho_prepost_min" %in% names(pr))
    if (!have_pr) {
      stop("Could not compute pre-post correlation via prepost_r(); supply r_obs, ",
           "or provide s1/s2 and either sc or (t or p) with n and mc (or m1 & m2).")
    }
    r_obs_num <- pr$rho_prepost_min
  }

  # --- ceilings from analyte_max_cor() (handle BOTH old-list and new-data.frame shapes)
  rA <- rT <- rep(NA_real_, L)
  for (i in seq_len(L)) {
    res <- try(
      analyte_max_cor(
        analyte = An[i], m1 = M1[i], s1 = S1[i], m2 = M2[i], s2 = S2[i],
        k = k, use = use, append = FALSE
      ),
      silent = TRUE
    )
    if (!inherits(res, "try-error")) {
      # New cohesive version: data.frame with columns
      if (is.data.frame(res) &&
          all(c("r_max_analytic","r_max_total") %in% names(res))) {
        rA[i] <- suppressWarnings(as.numeric(res$r_max_analytic[1]))
        rT[i] <- suppressWarnings(as.numeric(res$r_max_total[1]))
        # Legacy list version with $r_max list
      } else if (is.list(res) && !is.null(res$r_max)) {
        rA[i] <- suppressWarnings(as.numeric(res$r_max$analytic_only))
        rT[i] <- suppressWarnings(as.numeric(res$r_max$total))
      }
    }
  }

  # --- flags + probabilities (use r_obs_num consistently)
  exceeds_analytic <- as.numeric(r_obs_num) > rA
  exceeds_total    <- as.numeric(r_obs_num) > rT

  clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)
  p_ge_r_given_rho <- function(r_obs_vec, rho_true, nvec) {
    ok <- is.finite(r_obs_vec) & is.finite(rho_true) & is.finite(nvec) & (nvec > 3)
    out <- rep(NA_real_, length(r_obs_vec))
    if (!any(ok)) return(out)
    z_obs  <- atanh(clamp(r_obs_vec[ok], -0.999999, 0.999999))
    z_true <- atanh(clamp(rho_true[ok],   -0.999999, 0.999999))
    se <- 1 / sqrt(nvec[ok] - 3)
    out[ok] <- 1 - stats::pnorm((z_obs - z_true) / se)
    pmin(pmax(out, 0), 1)
  }
  probA <- p_ge_r_given_rho(as.numeric(r_obs_num), rA, as.numeric(N))
  probT <- p_ge_r_given_rho(as.numeric(r_obs_num), rT, as.numeric(N))

  # --- build output
  out_core <- data.frame(
    n_used                  = as.numeric(N),
    r_max_analytic_expected = rA,
    r_max_total_expected    = rT,
    r_obs                   = as.numeric(r_obs_num),
    r_obs_min               = as.numeric(r_obs_num),   # keep for plot compatibility
    exceeds_max_analytic    = as.logical(exceeds_analytic),
    exceeds_max_total       = as.logical(exceeds_total),
    prob_analytic           = probA,
    prob_total              = probT,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  res <- if (append && !is.null(data)) {
    cbind(data, out_core)
  } else {
    data.frame(
      analyte = An,
      out_core,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  # --- printer metadata + class (ensure *first* class is ours)
  attr(res, "pcc_meta") <- list(
    analyte_col = if (!is.null(data)) {
      nm <- deparse(substitute(analyte)); if (nm %in% names(res)) nm else "analyte"
    } else "analyte",
    k = k, use = use
  )
  class(res) <- c("analyte_cor_check", setdiff(class(res), "analyte_cor_check"))
  res
}




#' Print method for \code{analyte_cor_check}
#'
#' Produces a compact, human-readable table for one or more cases from
#' \code{analyte_cor_check()}, showing the observed lower-bound correlation used
#' for flagging, each analyte ceiling, flags, and one-tailed probabilities.
#'
#' @param x An object of class \code{analyte_cor_check}, as returned by
#'   \code{analyte_cor_check()}.
#' @param digits Integer. Decimal places for correlations/ceilings. Default 3.
#' @param p_digits Integer. Decimal places for p-values. Default 3.
#' @param max_rows Integer. Max rows to display. Default 20.
#' @param ... Unused; present for S3 method compatibility.
#' @rdname analyte_cor_check
#' @method print analyte_cor_check
#' @exportS3Method print analyte_cor_check
print.analyte_cor_check <- function(x,
                                    digits = 3,
                                    p_digits = 3,
                                    max_rows = 20,
                                    ...) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  meta <- attr(x, "pcc_meta")
  if (is.null(meta)) { NextMethod(); return(invisible(x)) }

  ac <- meta$analyte_col %||% "analyte"

  # Pick an n vector to display (prefer n_used; else a column named 'n'; else NA)
  n_vec <- if ("n_used" %in% names(x)) x[["n_used"]] else if ("n" %in% names(x)) x[["n"]] else rep(NA_real_, nrow(x))

  req <- c("r_obs", "r_max_analytic_expected", "r_max_total_expected",
           "exceeds_max_analytic", "exceeds_max_total",
           "prob_analytic", "prob_total")
  miss <- setdiff(req, names(x))
  if (length(miss)) stop("Missing required columns for printing: ", paste(miss, collapse = ", "))

  fmt_num <- function(v, d = digits) ifelse(is.na(v), "NA", formatC(v, format = "f", digits = d))
  fmt_p   <- function(v) {
    ifelse(is.na(v), "NA",
           ifelse(v < 10^(-p_digits),
                  paste0("<", format(10^(-p_digits), scientific = TRUE)),
                  formatC(v, format = "f", digits = p_digits)))
  }
  flag_str <- function(b) ifelse(is.na(b), "NA", ifelse(b, "[X]", "[OK]"))

  n_all <- nrow(x); show <- seq_len(min(n_all, max_rows))
  tbl <- data.frame(
    case               = show,
    Analyte            = if (ac %in% names(x)) x[[ac]][show] else rep(NA_character_, length(show)),
    n                  = fmt_num(n_vec[show], 0),
    r_obs              = fmt_num(x[["r_obs"]][show]),
    r_max_analytic     = fmt_num(x[["r_max_analytic_expected"]][show]),
    r_max_total        = fmt_num(x[["r_max_total_expected"]][show]),
    flag_analytic      = flag_str(x[["exceeds_max_analytic"]][show]),
    flag_total         = flag_str(x[["exceeds_max_total"]][show]),
    prob_analytic      = fmt_p(x[["prob_analytic"]][show]),
    prob_total         = fmt_p(x[["prob_total"]][show]),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  cat("\n# Pre-post correlation ceiling check\n", sep = "")
  cat("use = ", meta$use, ", k = ", meta$k, "\n\n", sep = "")
  print(tbl, row.names = FALSE, right = TRUE)
  if (n_all > length(show)) {
    cat("\n... ", n_all - length(show), " more rows not shown (increase max_rows to see more).\n", sep = "")
  }
  cat("\nLegend: [OK] = does not exceed r_max; [X] = exceeds r_max (flag)\n\n")
  invisible(x)
}










#' Plot one-tailed exceedance curves for an analyte case
#'
#' @description
#' Given an `analyte_cor_check` result, draws the one-tailed probability
#' curves `Pr(R >= r | rho = r_max)` for both analytic-only and analytic+biological
#' ceilings, shading under each curve across the plotting range. A dotted line marks
#' the observed value (either `r_obs_min` or `r_obs`), and points show the exact
#' one-tailed probabilities at that value.
#'
#' @param x An object returned by `analyte_cor_check()`.
#' @param case Integer row index to plot.
#' @param n_points Number of x-grid points.
#' @param xlim Numeric length-2 bounds for r (default `c(0, 1)`).
#' @param aspect.ratio Plot aspect ratio passed to `theme()`.
#' @param ... Ignored.
#'
#' @return A `ggplot` object (invisibly), and also printed.
#' @method plot analyte_cor_check
#' @export
plot.analyte_cor_check <- function(x,
                                   case = 1,
                                   n_points = 1000,
                                   xlim = c(0, 1),
                                   aspect.ratio = 0.5,
                                   ...) {
  if (!inherits(x, "analyte_cor_check")) stop("x must be an 'analyte_cor_check' object.")
  if (!is.numeric(case) || length(case) != 1 || case < 1 || case > nrow(x))
    stop("`case` must be an integer between 1 and nrow(x).")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("This plot method requires ggplot2.")

  # choose observed column (prefer r_obs_min; else r_obs)
  rcol <- if ("r_obs_min" %in% names(x)) "r_obs_min" else if ("r_obs" %in% names(x)) "r_obs" else NULL
  if (is.null(rcol)) stop("Object lacks 'r_obs_min' and 'r_obs' columns; cannot plot.")
  obs_lab <- if (identical(rcol, "r_obs_min")) "r_obs_min" else "r_obs"

  # pull row values
  n_i  <- if ("n_used" %in% names(x)) as.numeric(x$n_used[case]) else
    if ("n" %in% names(x)) as.numeric(x$n[case]) else NA_real_
  rhoA <- as.numeric(x$r_max_analytic_expected[case])
  rhoT <- as.numeric(x$r_max_total_expected[case])
  rcut <- as.numeric(x[[rcol]][case])

  if (!is.finite(n_i) || n_i <= 3) stop("Selected case lacks a valid n (must be > 3).")
  if (!is.finite(rhoA) || !is.finite(rhoT)) stop("Selected case lacks finite ceiling(s).")
  if (!is.finite(rcut)) stop("Selected case lacks a finite ", obs_lab, ".")

  # clamp & sort xlim
  xlim <- sort(pmax(pmin(xlim, 1), 0))

  # labels
  meta <- attr(x, "pcc_meta")
  analyte_lab <- {
    col <- if (!is.null(meta) && !is.null(meta$analyte_col)) meta$analyte_col else "analyte"
    if (col %in% names(x)) as.character(x[[col]][case]) else NA_character_
  }
  title_txt <- if (!is.na(analyte_lab)) paste0("Case ", case, " - ", analyte_lab) else paste("Case", case)
  subtitle_txt <- if (!is.null(meta)) paste0("use = ", meta$use, ", k = ", meta$k) else NULL

  # one-tailed exceedance via Fisher z
  clamp  <- function(z, lo, hi) pmin(pmax(z, lo), hi)
  one_tail <- function(r, rho, n) {
    z <- (atanh(clamp(r,   -0.999999, 0.999999)) -
            atanh(clamp(rho, -0.999999, 0.999999))) * sqrt(n - 3)
    1 - stats::pnorm(z)
  }

  # grid & curves
  r_grid <- seq(xlim[1], xlim[2], length.out = n_points)
  yA <- one_tail(r_grid, rhoA, n_i)
  yT <- one_tail(r_grid, rhoT, n_i)

  # exact probs at observed
  pA_obs <- one_tail(rcut, rhoA, n_i)
  pT_obs <- one_tail(rcut, rhoT, n_i)

  # data frames
  lvl_names <- c("Analytic ceiling", "Analytic + biological ceiling")
  df_line <- rbind(
    data.frame(r = r_grid, y = yA, type = lvl_names[1]),
    data.frame(r = r_grid, y = yT, type = lvl_names[2])
  )
  df_line$type <- factor(df_line$type, levels = lvl_names)

  # ---- NSE bindings to silence R CMD check (no rlang needed) ----
  r <- y <- type <- NULL

  ann_text <- paste0(
    obs_lab, " = ", formatC(rcut, format = "f", digits = 3), "\n",
    "r_max_analytic = ", formatC(rhoA, format = "f", digits = 3),
    ",  p = ", formatC(pA_obs, format = "f", digits = 3), "\n",
    "r_max_total    = ", formatC(rhoT, format = "f", digits = 3),
    ",  p = ", formatC(pT_obs, format = "f", digits = 3), "\n",
    "n = ", formatC(n_i, format = "f", digits = 0)
  )
  x_ann <- xlim[1] + 0.02 * diff(xlim)
  y_ann <- 0.06

  p <- ggplot2::ggplot(df_line, ggplot2::aes(x = r, y = y, color = type)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = 0, ymax = y, fill = type, group = type),
      alpha = 0.18, colour = NA
    ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_vline(xintercept = rcut, linetype = "dotted", linewidth = 0.7) +
    ggplot2::geom_point(
      data = data.frame(
        r = rcut,
        y = c(pA_obs, pT_obs),
        type = factor(lvl_names, levels = lvl_names)
      ),
      ggplot2::aes(x = r, y = y, color = type),
      size = 2, inherit.aes = FALSE
    ) +
    ggplot2::scale_x_continuous(limits = xlim, breaks = seq(0, 1, by = 0.1),
                                expand = ggplot2::expansion(mult = 0.00)) +
    ggplot2::scale_y_continuous(limits = c(0, 1.012), breaks = seq(0, 1, by = 0.1),
                                expand = ggplot2::expansion(mult = 0.00)) +
    ggplot2::labs(
      title = title_txt,
      subtitle = subtitle_txt,
      x = "Correlation (r)",
      y = "Pr(r \u2265 r_obs | true r = analyte ceiling)",
      color = NULL, fill = NULL
    ) +
    ggplot2::annotate("text", x = x_ann, y = y_ann, label = ann_text,
                      hjust = 0, vjust = 0, size = 4.5) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      aspect.ratio = aspect.ratio,
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::guides(fill = "none")

  print(p)
  invisible(p)
}
