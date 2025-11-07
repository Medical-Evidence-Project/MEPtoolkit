#' R6 Class representing the results from a carlisle() call
#'
#' @description
#' Mainly for pretty printing, but all output is available
Carlisle_res <- R6::R6Class("Carlisle_res",
                            public = list(
                              #' @field Student Carlisle min/RIVETS/max test results using
                              #' the range of Student p values from a T_slicer_res object
                              Student = NULL,
                              #' @field Welch Carlisle min/RIVETS/max test results using
                              #' the range of Welch p values from a T_slicer_res object
                              Welch = NULL,
                              #' @field Provided_p Carlisle min/RIVETS/max test results using
                              #' the range of p values manually provided
                              Provided_p = NULL,
                              #' @description
                              #' Create a new Carlisle_res object
                              #' @param Student Student res
                              #' @param Welch Welch res
                              #' @param Provided_p Provided_p res
                              #' @return A new `Carlisle_res` object
                              initialize = function(Student = NULL,
                                                    Welch = NULL,
                                                    Provided_p = NULL) {
                                self$Student <- Student
                                self$Welch <- Welch
                                self$Provided_p <- Provided_p
                              },
                              #' @description
                              #' Pretty printing of Carlisle results
                              print = function() {
                                cat("Results from Carlisle\n")
                                cat("---------------------------------------\n")
                                if(is.null(self$Provided_p)) {
                                  cat("RIVETS p values:\n")
                                  cat("      Student: p =", self$Student[2], "\n")
                                  cat("      Welch:   p =", self$Welch[2], "\n\n")
                                  cat("Min/max p value combinations:\n")
                                  cat("  - Minimum\n")
                                  cat("      Student: p =", self$Student[1], "\n")
                                  cat("      Welch:   p =", self$Welch[1], "\n\n")
                                  cat("  - Maximum\n")
                                  cat("      Student: p =", self$Student[3], "\n")
                                  cat("      Welch:   p =", self$Welch[3], "\n")
                                }
                                else {
                                  cat("RIVETS p values            :", self$Provided_p[2], "\n")
                                  cat("minimum p value combination:", self$Provided_p[1], "\n")
                                  cat("maximum p value combination:", self$Provided_p[3], "\n")
                                }
                                cat("---------------------------------------\n")
                              }
                            )
)

#' Perform a Carlisle test
#'
#' A Carlisle test is a method for examining the table of participants'
#' baseline characteristics in reports of randomized trials, and checking if
#' groups are too consistently different (or similar) than what would be expected
#' in case of successful randomization. Takes all p values with the assumption
#' of independence of tests and the NULL is true to calculate an aggregated p
#' value using a Fisher-Stouffer method.
#'
#' @param p A character vector. Either the p values from tests of differences
#' between the groups at baseline, or a T_slicer_res object (the results from a
#' t_slicer() call)
#' @param output A logical. Optional (default = FALSE). Determines if the
#' function should print its output before returning it
#' @return A Carlisle_res object. Either contains a single set of min/RIVETS/max
#' results from the Carlisle test if p values were manually provided, or two sets
#' of min/RIVETS/max results (once for Student t tests, one for Welch t tests)
#' if a T_slicer_res object was provided.
#' @importFrom stats qnorm
#' @examples
#' carlisle(c(".99", ".95", ".98", ".90"), output=TRUE)
#' data <- data.frame(m1 = c("1.20", "1.25"), s1 = c("1.2", "1.25"), n1 = c(60, 60),
#'                    m2 = c("2.1", "2.15"), s2 = c("2.5", "2.55"), n2 = c(30, 30),
#'                    p = c(NA, NA))
#' res <- t_slicer(m1 = m1, s1 = s1, n1 = n1, m2 = m2, s2 = s2, n2 = n2,
#'                 data = data, output = TRUE)
#' carlisle(res, output = TRUE)
#'
#' @export
carlisle <- function(p, output = FALSE) {

  if(class(p)[1]!="T_slicer_res") {
    res_list <- list()
    carlisle_res <- c()

    p <- validate_p(p)

    for(i in seq_len(length(p$original))) {
      res_list[[i]] <- c(p$minimum[i], p$maximum[i])
    }
    cases_grid <- expand.grid(res_list)

    for(i in seq_len(nrow(cases_grid))) {
      carlisle_res <- c(carlisle_res, 1 - pnorm(sum(sapply(unlist(cases_grid[i,]), qnorm))/sqrt(length(unlist(cases_grid[i,])))))
    }
    carlisle_RIVETS <- 1 - pnorm(sum(sapply(p$original, qnorm))/sqrt(length(p$original)))

    result <- Carlisle_res$new(Provided_p = c(min(carlisle_res), carlisle_RIVETS, max(carlisle_res)))
  }
  else {
    res_Student <- list()
    res_Welch <- list()
    carlisle_Student_res <- c()
    carlisle_Welch_res <- c()

    for(i in seq_len(length(p$Student_RIVETS))) {
      res_Student[[i]] <- c(p$Student_min[i], p$Student_max[i])
    }

    cases_grid <- expand.grid(res_Student)


    for(i in seq_len(nrow(cases_grid))) {
      carlisle_Student_res <- c(carlisle_Student_res, 1 - pnorm(sum(sapply(unlist(cases_grid[i,]), qnorm))/sqrt(length(unlist(cases_grid[i,])))))
    }
    carlisle_Student_RIVETS <- 1 - pnorm(sum(sapply(p$Student_RIVETS, qnorm))/sqrt(length(p$Student_RIVETS)))

    for(i in seq_len(length(p$Welch_RIVETS))) {
      res_Welch[[i]] <- c(p$Welch_min[i], p$Welch_max[i])
    }

    cases_grid <- expand.grid(res_Welch)
    carlisle_Welch_res <- c()

    for(i in seq_len(nrow(cases_grid))) {
      carlisle_Welch_res <- c(carlisle_Welch_res, 1 - pnorm(sum(sapply(unlist(cases_grid[i,]), qnorm))/sqrt(length(unlist(cases_grid[i,])))))
    }
    carlisle_Welch_RIVETS <- 1 - pnorm(sum(sapply(p$Welch_RIVETS, qnorm))/sqrt(length(p$Welch_RIVETS)))

    result <- Carlisle_res$new(Student = c(min(carlisle_Student_res), carlisle_Student_RIVETS, max(carlisle_Student_res)),
                               Welch = c(min(carlisle_Welch_res), carlisle_Welch_RIVETS, max(carlisle_Welch_res)))

  }

  if(output==TRUE) {
    print(result)
  }
  invisible(result)
}
