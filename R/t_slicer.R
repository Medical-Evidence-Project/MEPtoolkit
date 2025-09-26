#' t_slicer
#'
#' Checks if a t test could be a Welch or a Student test (useful for cases
#' when it is not explicitly stated)
#'
#' t_slicer takes as input 2 sets of descriptive statistics (mean and sd)
#' and a reported p value. It computes the results of the Student's and Welch's
#' t tests, and compares the resulting with the reported p values to see which
#' test(s) - if any - match(es)
#'
#' It takes into account the rounding (im)precision, and therefore
#' works with ranges of possible values. It assumes that reported values could
#' have been rounded up or down
#'
#' Plots a graph on top of the return value
#'
#' @param m1 A character string. The reported mean of the first group
#' @param s1 A character string. The reported standard deviation of the first group
#' @param n1 A number. The sample size of the first group
#' @param m2 A character string. The reported mean of the second group
#' @param s2 character string. The reported standard deviation of the second group
#' @param n2 A number. The sample size of the second group
#' @param p A character string. The reported p value
#' @return A list of Booleans for each test (TRUE for a match, FALSE for a mismatch)
#' @examples
#' t_slicer("1.2", "1.2", 60, "2.1", "2.5", 30, ".08")
#' t_slicer("1.2", "1.2", 60, "2.1", "2.5", 60, ".02")
#'
#' @export
t_slicer <- function (m1,s1,n1,m2,s2,n2,p)
{
  # Input validation
  m1 <- validation_descriptive(m1, "m1")
  s1 <- validation_descriptive(s1, "s1")
  m2 <- validation_descriptive(m2, "m2")
  s2 <- validation_descriptive(s2, "s2")
  p <- validation_p(p)

  m1_range <- c(m1$minimum, m1$maximum)
  m1_num <- m1$original
  s1_range <- c(s1$minimum, s1$maximum)
  s1_num <- s1$original
  m2_range <- c(m2$minimum, m2$maximum)
  m2_num <- m2$original
  s2_range <- c(s2$minimum, s2$maximum)
  s2_num <- s2$original
  p_range <- c(p$minimum, p$maximum)
  p_num <- p$original

  cases_grid <- expand.grid(m1 = m1_range, s1 = s1_range, m2 = m2_range, s2 = s2_range)

  # RIVETS values
  twelch <- BSDA::tsum.test(mean.x = m1_num, s.x = s1_num, n.x = n1, mean.y = m2_num, s.y = s2_num, n.y = n2, var.equal = FALSE)
  tstudent <- BSDA::tsum.test(mean.x = m1_num, s.x = s1_num, n.x = n1, mean.y = m2_num, s.y = s2_num, n.y = n2, var.equal = TRUE)

  # Get all possible min and max combinations for means and sds for g1 and g2
  p_twelch_grid <- c()
  p_tstudent_grid <- c()
  for(i in seq_len(nrow(cases_grid)))
  {
    p_twelch_grid <- c(p_twelch_grid, BSDA::tsum.test(mean.x = cases_grid[i, "m1"], s.x = cases_grid[i, "s1"], n.x = n1, mean.y = cases_grid[i, "m2"], s.y = cases_grid[i, "s2"], n.y = n2, var.equal = FALSE)$p.value)
    p_tstudent_grid <- c(p_tstudent_grid, BSDA::tsum.test(mean.x = cases_grid[i, "m1"], s.x = cases_grid[i, "s1"], n.x = n1, mean.y = cases_grid[i, "m2"], s.y = cases_grid[i, "s2"], n.y = n2, var.equal = TRUE)$p.value)
  }
  print(max(p_twelch_grid))
  print(max(p_tstudent_grid))

  # Plot
  plot_data <- data.frame(
    Student = c(min(p_tstudent_grid), tstudent$p.value, max(p_tstudent_grid)),
    Welch = c(min(p_twelch_grid), twelch$p.value, max(p_twelch_grid)),
    Value = c("Minimum", "Mid", "Maximum")
  )

  rect_data <- data.frame(
    Reported = c("Student", "Welch"),
    # Right now rect X coordinates are hard-coded, but would have to be changed
    # if more tests are added
    Xmin = c(0.65, 1.65),
    Xmax = c(1.35, 2.35)
  )

  # Return value
  result <- list(Student = FALSE, Welch = FALSE)

  if(p_range[1] <= max(p_tstudent_grid) && p_range[2] >= min(p_tstudent_grid)) {
    color_rect_student <- "green"
    result$Student <- TRUE
  }
  else
    color_rect_student <- "red"
  if(p_range[1] <= max(p_twelch_grid) && p_range[2] >= min(p_twelch_grid)) {
    color_rect_welch <- "green"
    result$Welch <- TRUE
  }
  else
    color_rect_welch <- "red"

  p <- ggplot(plot_data) +
    geom_point(aes(x="Welch", y=.data$Welch, shape=.data$Value, fill=.data$Value), size = 5) +
    geom_point(aes(x="Student", y=.data$Student, shape=.data$Value, fill=.data$Value), size = 5) +
    theme_minimal() +
    scale_fill_manual(
      name= "Value",
      values = c('Minimum' = "red", 'Mid' = "blue", 'Maximum' = "red")) +
    scale_shape_manual(
      name="Value",
      values = c('Minimum' = 24, 'Mid' = 21, 'Maximum' = 25)) +
    ggnewscale::new_scale("fill") +
    geom_rect(
      data = rect_data,
      aes(xmin=.data$Xmin, xmax=.data$Xmax, ymin=p_range[1],
          ymax=p_range[2], fill=.data$Reported),
      alpha=.5, color="black") +
    scale_fill_manual(
      values = c('Student' = color_rect_student, 'Welch' = color_rect_welch)) +
    labs(title = "T-Test Results", y = "Potential p-values", x = "Tests") +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=13, face="bold"))

  plot(p)
  return(result)
}
