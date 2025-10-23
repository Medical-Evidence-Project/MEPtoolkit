#' R6 Class representing the results from an analyte_max_cor_res call
#'
#' @description
#' Used for pretty printing
Analyte_cor_check_res <- R6::R6Class("Analyte_cor_check_res",
                                   public = list(
                                     #' @field disp Display table.
                                     disp = NULL,
                                     #' @field k Number of replicates averaged.
                                     k = NULL,
                                     #' @field use Which cross-sectional CV was used.
                                     use = NULL,
                                     #' @field analyte The list of analytes
                                     analyte = NULL,
                                     #' @field n Number of participants
                                     n = NULL,
                                     #' @description
                                     #' Create a new Analyte_cor_check_res object
                                     #' @param disp Display table.
                                     #' @param k Number of replicates averaged.
                                     #' @param use Which cross-sectional CV was used.
                                     #' @param analyte List of analytes
                                     #' @param n Number of participants
                                     #' @return A new `Analyte_cor_check_res` object
                                     initialize = function(disp, k, use, analyte, n) {
                                       self$disp <- disp
                                       self$k <- k
                                       self$use <- use
                                       self$analyte <- analyte
                                       self$n <- n
                                     },
                                     #' @description Printing
                                     #' @param digits The number of digits to be displayed
                                     #' @param p_digits The number of digits to be displayed for p
                                     #' @param max_rows The maximum number of rows to be displayed
                                     #' Pretty printing of Analyte_cor_check results
                                     print = function(
                                                      digits = 3,
                                                      p_digits = 3,
                                                      max_rows = 20) {
                                       fmt_num <- function(v, d = digits) ifelse(is.na(v), "NA", formatC(v, format = "f", digits = d))
                                       fmt_p   <- function(v) {
                                         ifelse(is.na(v), "NA",
                                                ifelse(v < 10^(-p_digits),
                                                       paste0("<", format(10^(-p_digits), scientific = TRUE)),
                                                       formatC(v, format = "f", digits = p_digits)))
                                       }
                                       flag_str <- function(b) ifelse(is.na(b), "NA", ifelse(b, "[X]", "[OK]"))

                                       n_all <- nrow(self$disp); show <- seq_len(min(n_all, max_rows))
                                       tbl <- data.frame(
                                         case               = show,
                                         Analyte            = self$analyte,
                                         n                  = fmt_num(self$n[show], 0),
                                         r_obs              = fmt_num(self$disp[["r_obs"]][show]),
                                         r_obs_min          = fmt_num(self$disp[["r_obs_min"]][show]),
                                         r_max_analytic     = fmt_num(self$disp[["r_max_analytic_expected"]][show]),
                                         r_max_total        = fmt_num(self$disp[["r_max_total_expected"]][show]),
                                         flag_analytic      = flag_str(self$disp[["exceeds_max_analytic"]][show]),
                                         flag_total         = flag_str(self$disp[["exceeds_max_total"]][show]),
                                         prob_analytic      = fmt_p(self$disp[["prob_analytic"]][show]),
                                         prob_total         = fmt_p(self$disp[["prob_total"]][show]),
                                         check.names = FALSE,
                                         stringsAsFactors = FALSE
                                       )

                                       cat("\n# Pre-post correlation ceiling check\n", sep = "")
                                       cat("use = ", self$use, ", k = ", self$k, "\n\n", sep = "")
                                       print(tbl, row.names = FALSE, right = TRUE)
                                       if (n_all > length(show)) {
                                         cat("\n... ", n_all - length(show), " more rows not shown (increase max_rows to see more).\n", sep = "")
                                       }
                                       cat("\nLegend: [OK] = does not exceed r_max; [X] = exceeds r_max (flag)\n\n")
                                     },
                                     #' @description Plotting
                                     #' @param case The number of digits to be displayed
                                     #' @param n_points The number of digits to be displayed for p
                                     #' @param xlim The maximum number of rows to be displayed
                                     #' @param aspect.ratio The maximum number of rows to be displayed
                                     #' @param ... The maximum number of rows to be displayed
                                     #' Pretty printing of Analyte_cor_check results
                                     plot = function(
                                                     case = 1,
                                                     n_points = 1000,
                                                     xlim = c(0, 1),
                                                     aspect.ratio = 0.5,
                                                     ...) {
                                       if (!is.numeric(case) || length(case) != 1 || case < 1 || case > nrow(self$disp))
                                         stop("`case` must be an integer between 1 and nrow(x).")
                                       if (!requireNamespace("ggplot2", quietly = TRUE)) stop("This plot method requires ggplot2.")

                                       rcol <- self$disp$r_obs_min
                                       obs_lab <- "r_obs_min"

                                       # pull row values
                                       n_i  <- self$n[case]
                                       rhoA <- self$disp$r_max_analytic_expected[case]
                                       rhoT <- self$disp$r_max_total_expected[case]
                                       rcut <- rcol[case]

                                       if (!is.finite(n_i) || n_i <= 3) stop("Selected case lacks a valid n (must be > 3).")
                                       if (!is.finite(rhoA) || !is.finite(rhoT)) stop("Selected case lacks finite ceiling(s).")
                                       if (!is.finite(rcut)) stop("Selected case lacks a finite ", obs_lab, ".")

                                       # clamp & sort xlim
                                       xlim <- sort(pmax(pmin(xlim, 1), 0))

                                       # labels
                                       analyte_lab <- self$analyte[case]
                                       title_txt <- if (!is.na(analyte_lab)) paste0("Case ", case, " - ", analyte_lab) else paste("Case", case)
                                       subtitle_txt <- paste0("use = ", self$use, ", k = ", self$k)

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
                                   )
)


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
#' @param append Logical. If TRUE and `data` is supplied, return `data`
#'        with new columns appended; else return only the new columns (and include
#'        `analyte` and `n_used` so printing works).
#' @param output Logical. If \code{TRUE}, also outputs a summary of results in the console.
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
    append = FALSE,
    output = FALSE
) {
  use <- match.arg(use)

  if(!is.null(data)) {
    if(!is.data.frame(data)) stop("`data` must be a data.frame.")
    analyte <- resolve(enquo(analyte), data); m1 <- resolve(enquo(m1), data)
    s1 <- resolve(enquo(s1), data); m2 <- resolve(enquo(m2), data)
    s2 <- resolve(enquo(s2), data); n <- resolve(enquo(n), data)
    r_obs <- resolve(enquo(r_obs), data); sc <- resolve(enquo(sc), data)
    mc <- resolve(enquo(mc), data); t <- resolve(enquo(t), data)
    p <- resolve(enquo(p), data)
  }

  # -- check lengths
  args <- list(analyte = analyte, m1 = m1, s1 = s1, m2 = m2, s2 = s2, n = n,
               r_obs = r_obs, sc = sc, mc = mc, t = t, p= p)
  lens <- vapply(args, length, integer(1))
  L <- max(lens)
  ok_len <- lens == 0 | lens == L
  if (!all(ok_len))
    stop("arguments must have the same length.")

  if (is.null(n) || any(is.na(n)))
    stop("Sample size `n` must be provided.")

  work <- NULL
  if(!is.null(r_obs)) {
    r_obs_min <- validate_descriptive(r_obs, allow_na = TRUE)$minimum
    r_obs_orig <- validate_descriptive(r_obs, allow_na = TRUE)$original
    na_r_obs <- is.na(r_obs_min)
    if(sum(na_r_obs>1)) {
      work <- tibble(
        s1 = s1[na_r_obs], s2 = s2[na_r_obs], m1 = m1[na_r_obs],
        m2 = m2[na_r_obs], sc = sc[na_r_obs], mc = mc[na_r_obs],
        n = n[na_r_obs], t = t[na_r_obs], p = p[na_r_obs])
      prepost_r_res <- prepost_r(s1 = s1, s2 = s2, sc = sc, mc = mc,
                                 m1 = m1, m2 = m2, n = n, t = t, p = p,
                                 data = work, append = TRUE)
      r_obs_min[na_r_obs] <- prepost_r_res$rho_prepost_min
      r_obs_orig[na_r_obs] <- prepost_r_res$rho_prepost
    }
  }
  else {
    work <- tibble(
      s1 = s1, s2 = s2, m1 = m1,
      m2 = m2, sc = sc, mc = mc,
      n = n, t = t, p = p)
    prepost_r_res <- prepost_r(s1 = s1, s2 = s2, sc = sc, mc = mc,
                          m1 = m1, m2 = m2, n = n, t = t, p = p,
                          data = work, append = TRUE)
    r_obs_min <- prepost_r_res$rho_prepost_min
    r_obs_orig <- prepost_r_res$rho_prepost
  }

  # --- ceilings from analyte_max_cor()
  res <- analyte_max_cor(analyte = analyte, m1 = m1, s1 = s1, m2 = m2,
                         s2 = s2, k = k, use = use)
  rA <- res$disp$r_max_analytic
  rT <- res$disp$r_max_total

  # --- flags + probabilities
  exceeds_analytic <- r_obs_min > rA
  exceeds_total    <- r_obs_min > rT

  clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)
  p_ge_r_given_rho <- function(r_obs_vec, rho_true, nvec) {
    ok <- is.finite(r_obs_vec) & is.finite(rho_true) & is.finite(nvec) & (nvec > 3)
    out <- rep(NA_real_, length(r_obs_vec))
    if (!any(ok)) return(out)
    z_obs  <- atanh(clamp(r_obs_vec[ok], -0.999999, 0.999999))
    z_true <- atanh(clamp(rho_true[ok],   -0.999999, 0.999999))
    se <- 1 / sqrt(nvec[ok] - 3)
    out[ok] <- 1 - stats::pnorm((z_obs - z_true) / se)
    # See if we need to keep that or not
    # pmin(pmax(out, 0), 1)
  }
  probA <- p_ge_r_given_rho(as.numeric(r_obs_min), rA, n)
  probT <- p_ge_r_given_rho(as.numeric(r_obs_min), rT, n)

  # --- build output
  out_core <- data.frame(
    r_max_analytic_expected = rA,
    r_max_total_expected    = rT,
    r_obs                   = r_obs_orig,
    r_obs_min               = r_obs_min,
    exceeds_max_analytic    = exceeds_analytic,
    exceeds_max_total       = exceeds_total,
    prob_analytic           = probA,
    prob_total              = probT,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  res <- if (append && !is.null(data)) {
    cbind(data, out_core)
  } else {
    data.frame(
      analyte = analyte,
      n = n,
      out_core,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  output_res <- Analyte_cor_check_res$new(disp = res, k = k, use = use,
                                          analyte = analyte, n = n)
  if(output == TRUE) {
    print(output_res)
    plot(output_res)
  }

  if (append == FALSE || is.null(data))
    res <- output_res

  return(res)
}
