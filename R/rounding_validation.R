#' decimal_places
#'
#' Finds and returns the number of decimal places on a statistic
#'
#' @param x Character string. The statistic to examine
#' @return The number of decimal places
#' @examples
#' decimal_places("3")
#' decimal_places(".984")
#' decimal_places("2.00")
decimal_places <- function(x) {
  if (grepl("\\.", x))
    return(nchar(strsplit(x, "\\.")[[1]][2]))
  else
    return(0)
}

#' validate_descriptive
#'
#' Validates a descriptive statistic: checks it is a string that can be
#' converted to a number, and returns the possible range of values if it was
#' rounded
#'
#' @param statistic Character string. The statistic to validate
#' @param name Character string. The name of the variable (for error output)
#' @return A list of numbers: the possible range of values (minimum, original
#' and maximum)
#' @examples
#' validate_descriptive("2.14", "m1")
validate_descriptive <- function (statistic, name)
{
  if(!is.character(statistic)) {
    stop(paste(name, "must be a character string."))
  }
  else {
    suppressWarnings({
      statistic_num <- as.numeric(statistic)
      if(is.na(statistic_num)) {
        stop(paste(name, "must be convertible to a number."))
      }
    }
    )
  }
  statistic_range <- list(
    minimum = statistic_num - 5*10^(-decimal_places(statistic)-1),
    original = statistic_num,
    maximum = statistic_num + 5*10^(-decimal_places(statistic)-1)
  )
  return(statistic_range)
}

#' validate_p
#'
#' Validates a p statistic: checks it is a string that can be
#' converted to a number, and returns the possible range of values if it was
#' rounded
#'
#' @param p Character string. The p value to validate
#' @param name Character string. Optional. The name of the variable
#' (for error output)
#' @return A list of numbers: the possible range of values (minimum, original
#' and maximum)
#' @examples
#' validate_p(".14")
#' p_value <- ".25"
#' validate_p(p_value, "p_value")
#' validate_p("<.001")
validate_p <- function(p, name = "p")
{
  if(!is.character(p)) {
    stop(paste(name, "must be a character string."))
  }
  else {
    suppressWarnings({
      # TODO: Deal with "<.X"
      if((grepl("<", p))) {
        p <- strsplit(p, "<")[[1]][2]
        p_inferior <- TRUE
      }
      else
        p_inferior <- FALSE
      p_num <- as.numeric(p)
      if(is.na(as.numeric(p_num))) {
        stop(paste(name, "must be convertible to a number."))
      }
    }
    )
  }
  if(p_num > 1 || p_num < 0)
  {
    stop(paste(name, "must be between 0 and 1."))
  }
  if(p_inferior == FALSE) {
    p_range <- list(
      minimum = max(0, p_num - 5*10^(-decimal_places(p)-1)),
      original = p_num,
      maximum = min(1, p_num + 5*10^(-decimal_places(p)-1))
    )
  }
  else {
    p_range <- list(
      minimum = .Machine$double.xmin,
      original = p_num,
      maximum = p_num
    )
  }
  return(p_range)
}
