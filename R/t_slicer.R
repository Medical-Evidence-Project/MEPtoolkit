#' R6 Class representing the results from a t_slicer call
#'
#' @description
#' Mainly for pretty printing, but all output is available
T_slicer_res <- R6::R6Class("T_slicer_res",
            public = list(
              #' @field Student_RIVETS p value from Student's t test using
              #' RIVETS stats
              Student_RIVETS = NULL,
              #' @field Welch_RIVETS p value from Welch's t test using
              #' RIVETS stats
              Welch_RIVETS = NULL,
              #' @field Student_min p value from Student's t test using
              #' minimum possible values
              Student_min = NULL,
              #' @field Student_max p value from Student's t test using
              #' maximum possible values
              Student_max = NULL,
              #' @field Welch_min p value from Welch's t test using
              #' minimum possible values
              Welch_min = NULL,
              #' @field Welch_max p value from Welch's t test using
              #' maximum possible values
              Welch_max = NULL,
              #' @field Student are statistics consistent with a Student's t
              #' test?
              Student = NULL,
              #' @field Welch are statistics consistent with a Welch's t
              #' test?
              Welch = NULL,
              #' @field Plots the resulting plot
              Plots = NULL,
              #' @description
              #' Create a new T_slicer_res object
              #' @param Student_RIVETS Student_RIVETS
              #' @param Welch_RIVETS Welch_RIVETS
              #' @param Student_min Student_min
              #' @param Student_max Student_max
              #' @param Welch_min Welch_min
              #' @param Welch_max Welch_max
              #' @param Student Student
              #' @param Welch Welch
              #' @param Plots Plots
              #' @param n n
              #' @return A new `T_slicer_res` object
              initialize = function(Student_RIVETS, Welch_RIVETS, Student_min,
                                    Student_max, Welch_min, Welch_max,
                                    Student, Welch, Plots) {
                self$Student_RIVETS <- Student_RIVETS
                self$Welch_RIVETS <- Welch_RIVETS
                self$Student_min <- Student_min
                self$Student_max <- Student_max
                self$Welch_min <- Welch_min
                self$Welch_max <- Welch_max
                self$Student <- Student
                self$Welch <- Welch
                self$Plots <- Plots
              },
              #' @description
              #' Pretty printing of t_slicer results
              #' @param n if NULL (default), outputs results for all
              #' lines. Otherwise, only the one from the specified line.
              print = function(n=NULL) {
                cat("Results from t_slicer\n")
                if(is.null(n))
                  n <- 1:(length(self$Student))
                for(i in n) {
                  cat("\n\n--- Line", i,"---\n")
                  cat("---------------------------------------\n")
                  cat("RIVETS p values:\n")
                  cat("      Student: p =", self$Student_RIVETS[i], "\n")
                  cat("      Welch:   p =", self$Welch_RIVETS[i], "\n\n")
                  cat("Min/max p values:\n")
                  cat("  - Minimum\n")
                  cat("      Student: p =", self$Student_min[i], "\n")
                  cat("      Welch:   p =", self$Welch_min[i], "\n\n")
                  cat("  - Maximum\n")
                  cat("      Student: p =", self$Student_max[i], "\n")
                  cat("      Welch:   p =", self$Welch_max[i], "\n")
                  cat("---------------------------------------\n")
                  cat("Student's t:  ", if(is.na(self$Student[i])) "N/A\n"
                      else if (self$Student[i] == TRUE) "Consistent\n"
                      else "Inconsistent\n")
                  cat("Welch's t:    ", if(is.na(self$Welch[i])) "N/A\n"
                      else if (self$Welch[i] == TRUE) "Consistent\n"
                      else "Inconsistent\n")
                }
              },
              #' @description
              #' Plots the results.
              #' @param n if NULL (default), plots a grid of all graphs for all
              #' lines. Otherwise, plots the graph for the specified line.
              #' @importFrom cowplot plot_grid
              plot = function(n = NULL) {
                if(!is.null(n))
                  plot(self$Plots[[n]])
                else {
                  plot_grid(plotlist = self$Plots, ncol=2)
                }
              }
              )
            )

#' Slice t tests
#'
#' Checks if a t test could be a Welch or a Student test (useful for cases
#' when it is not explicitly stated)
#'
#' t_slicer takes as input 2 sets of descriptive statistics (means and sds)
#' and (optionally) reported p values. It computes the results of the Student's
#' and Welch's t tests, and compares the resulting with the reported p values
#' (if provided) to see which test(s) - if any - match(es)
#'
#' It takes into account the rounding (im)precision, and therefore
#' works with ranges of possible values. It assumes that reported values could
#' have been rounded up or down
#'
#' Plots a graph on top of the return value
#'
#' @param m1 A character vector. The reported mean of the first group
#' @param s1 A character vector. The reported standard deviation of the first group
#' @param n1 A numeric vector. The sample size of the first group
#' @param m2 A character vector. The reported mean of the second group
#' @param s2 A character vector. The reported standard deviation of the second group
#' @param n2 A numeric vector The sample size of the second group
#' @param p (Optional) A character vector. The reported p value
#' @param data A dataframe. If provided, m1, s1, n1, m2, s2, n2 and p are the
#' column names inside this dataframe.
#' @param output A logical. Optional (default = FALSE). Determines if the
#' function should print + plot its output before returning it
#' @return A T_slicer_res object. Contains logicals to represent the consistency
#' of inputs with both tests (if p values were given), min/max/RIVETS p values,
#' pretty printing and a plot
#' @examples
#' t_slicer("1.2", "1.2", 60, "2.1", "2.5", 30, ".08", TRUE)
#' data <- tibble(m1 = c("1.20", "1.25"), s1 = c("1.2", "1.25"), n1 = c(60, 60),
#'                m2 = c("2.1", "2.15"), s2 = c("2.5", "2.55"), n2 = c(30, 30),
#'                p = c(NA, NA))
#' t_slicer(m1 = m1, s1 = s1, n1 = n1, m2 = m2, s2 = s2, n2 = n2,
#'          data = data, output = TRUE)
#'
#' @import dplyr ggplot2
#' @export
t_slicer <- function (m1, s1, n1, m2, s2, n2, p = NULL, data = NULL, output = FALSE)
{
  # Input validation
  if(!is.null(data)) {
    if(!is.data.frame(data)) stop("`data` must be a data.frame.")
    m1 <- resolve(enquo(m1), data); s1 <- resolve(enquo(s1), data);
    n1 <- resolve(enquo(n1), data); m2 <- resolve(enquo(m2), data);
    s2 <- resolve(enquo(s2), data); n2 <- resolve(enquo(n2), data);
    p  <- resolve(enquo(p), data)
  }
  L <- length_check(list(m1=m1, s1=s1, m2=m2, s2=s2, n1=n1, n2=n2, p=p))
  m1 <- validate_descriptive(m1)
  s1 <- validate_descriptive(s1)
  m2 <- validate_descriptive(m2)
  s2 <- validate_descriptive(s2)
  p <- if(!is.null(p)) validate_p(p, allow_na = TRUE)
        else validate_p(rep(NA_character_, L), allow_na = TRUE)

  plots <- list()
  Student_pass <- c()
  Welch_pass <- c()
  Student_min <- c()
  Student_RIVETS <- c()
  Student_max <- c()
  Welch_min <- c()
  Welch_RIVETS <- c()
  Welch_max <- c()

  for (j in seq_len(L)) {
    m1_range <- c(m1$minimum[j], m1$maximum[j])
    m1_num <- m1$original[j]
    s1_range <- c(s1$minimum[j], s1$maximum[j])
    s1_num <- s1$original[j]
    m2_range <- c(m2$minimum[j], m2$maximum[j])
    m2_num <- m2$original[j]
    s2_range <- c(s2$minimum[j], s2$maximum[j])
    s2_num <- s2$original[j]
    p_range <- c(p$minimum[j], p$maximum[j])
    p_num <- p$original[j]

    cases_grid <- expand.grid(m1 = m1_range, s1 = s1_range, m2 = m2_range, s2 = s2_range)

    # RIVETS values
    twelch <- BSDA::tsum.test(mean.x = m1_num, s.x = s1_num, n.x = n1[j], mean.y = m2_num, s.y = s2_num, n.y = n2[j], var.equal = FALSE)$p.value
    tstudent <- BSDA::tsum.test(mean.x = m1_num, s.x = s1_num, n.x = n1[j], mean.y = m2_num, s.y = s2_num, n.y = n2[j], var.equal = TRUE)$p.value

    # Get all possible min and max combinations for means and sds for g1 and g2
    p_twelch_grid <- c()
    p_tstudent_grid <- c()
    for(i in seq_len(nrow(cases_grid)))
    {
      p_twelch_grid <- c(p_twelch_grid, BSDA::tsum.test(mean.x = cases_grid[i, "m1"], s.x = cases_grid[i, "s1"], n.x = n1, mean.y = cases_grid[i, "m2"], s.y = cases_grid[i, "s2"], n.y = n2, var.equal = FALSE)$p.value)
      p_tstudent_grid <- c(p_tstudent_grid, BSDA::tsum.test(mean.x = cases_grid[i, "m1"], s.x = cases_grid[i, "s1"], n.x = n1, mean.y = cases_grid[i, "m2"], s.y = cases_grid[i, "s2"], n.y = n2, var.equal = TRUE)$p.value)
    }

    Student_min <- c(Student_min, min(p_tstudent_grid))
    Student_RIVETS <- c(Student_RIVETS, tstudent)
    Student_max <- c(Student_max, max(p_tstudent_grid))
    Welch_min <- c(Welch_min, min(p_twelch_grid))
    Welch_RIVETS <- c(Welch_RIVETS, twelch)
    Welch_max <- c(Welch_max, max(p_twelch_grid))

    # Plot
    plot_data <- data.frame(
      Student = c(Student_min[j], Student_RIVETS[j], Student_max[j]),
      Welch = c(Welch_min[j], Welch_RIVETS[j], Welch_max[j]),
      Value = c("Minimum", "Mid", "Maximum")
    )

    if(!is.na(p_num)) {
      rect_data <- data.frame(
        Reported = c("Student", "Welch"),
        # Right now rect X coordinates are hard-coded, but would have to be changed
        # if more tests are added
        Xmin = c(0.65, 1.65),
        Xmax = c(1.35, 2.35),
        Ymin = c(p_range[1], p_range[1]),
        Ymax = c(p_range[2], p_range[2])
      )
      if(p_range[1] <= Student_max[j] && p_range[2] >= Student_min[j]) {
        color_rect_student <- "green"
        Student_pass <- c(Student_pass, TRUE)
      }
      else {
        color_rect_student <- "red"
        Student_pass <- c(Student_pass, FALSE)
      }
      if(p_range[1] <= Welch_max[j] && p_range[2] >= Welch_min[j]) {
        color_rect_welch <- "green"
        Welch_pass <- c(Welch_pass, TRUE)
      }
      else {
        color_rect_welch <- "red"
        Welch_pass <- c(Welch_pass, FALSE)
      }
    }
    else {
      Student_pass <- c(Student_pass, NA)
      Welch_pass <- c(Welch_pass, NA)
    }

    plots[[j]] <- ggplot(plot_data) +
      geom_point(aes(x="Welch", y=.data$Welch, shape=.data$Value, fill=.data$Value), size = 5) +
      geom_point(aes(x="Student", y=.data$Student, shape=.data$Value, fill=.data$Value), size = 5) +
      theme_minimal()  +
      labs(title = "T-Test Results", y = "Potential p-values", x = "Tests") +
      theme(axis.text = element_text(size=12),
            axis.title = element_text(size=13, face="bold")) +
      scale_fill_manual(
        name= "Value",
        values = c('Minimum' = "red", 'Mid' = "blue", 'Maximum' = "red")) +
      scale_shape_manual(
        name="Value",
        values = c('Minimum' = 24, 'Mid' = 21, 'Maximum' = 25))
    if(!is.na(p_num)) {
      plots[[j]] <- plots[[j]] + ggnewscale::new_scale("fill") +
        geom_rect(
          data = rect_data,
          aes(xmin=.data$Xmin, xmax=.data$Xmax, ymin=.data$Ymin,
              ymax=.data$Ymax, fill=.data$Reported),
          alpha=.5, color="black") +
        scale_fill_manual(
          values = c('Student' = color_rect_student, 'Welch' = color_rect_welch))
    }

  }

  # Return value
  result <- T_slicer_res$new(Student_RIVETS = Student_RIVETS,
                             Student = Student_pass, Welch = Welch_pass,
                             Welch_max = Welch_max,
                             Welch_RIVETS = Welch_RIVETS,
                             Welch_min = Welch_min,
                             Student_max = Student_max,
                             Student_min = Student_min,
                             Plots = plots)

  if(output == TRUE)
  {
    plot(result)
    print(result)
  }

  return(result)
}
