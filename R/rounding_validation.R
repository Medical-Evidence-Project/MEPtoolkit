#' Find and return the number of decimal places on a statistic
#'
#' @param statistic Number or (ideally) character vectors. The statistic to
#' examine
#' @return The number of decimal places
#' @examples
#' decimal_places("3")
#' decimal_places(".984")
#' decimal_places("2.00")
#' decimal_places(c("3.14", NA, "6.8"))
decimal_places <- function(statistic) {
  vapply(statistic, FUN.VALUE = integer(1), function(elem) {
    if (is.na(elem)) return(NA_integer_)
    if (grepl("\\.", elem)) {
      return(nchar(strsplit(elem, "\\.")[[1]][2]))
    }
    else
      return(0L)
  }, USE.NAMES = FALSE)
}

#' Validate a descriptive statistic
#'
#' This function checks if a statistic is a string or vector of
#' strings that can be converted to numbers, and returns the possible range of
#' values before rounding
#'
#' @param statistic Character string or vector of strings.
#' The statistic(s) to validate
#' @param allow_na Logical. Allows the use of NA values
#' (returns NA for min/original/max)
#' @return A list of numbers: the possible range of values (minimum, original
#' and maximum)
#' @examples
#' validate_descriptive("2.14")
#' validate_descriptive(c("2.14", "3.2"))
#' validate_descriptive(c("2.14", NA, "3.2"), allow_na = TRUE)
validate_descriptive <- function (statistic, allow_na = FALSE)
{
  name <- substitute(statistic)
  if(typeof(name) != "symbol")
    name <- "statistic"
  if (!is.character(statistic)) {
    stop(paste(name, "must be a character string or character vector."))
  }
  suppressWarnings({
    statistic_num <- as.numeric(statistic)
  })

  if (any(is.na(statistic)) && !allow_na) {
    stop(paste(name, "contains NA; set allow_na = TRUE to permit NA values."))
  }

  if (any(is.na(statistic_num) & !is.na(statistic))) {
    stop(paste(name, "must be parseable as numeric (or NA)."))
  }

  dp <- decimal_places(statistic)

  statistic_range <- list(
    minimum = statistic_num - 5*10^(-dp-1),
    original = statistic_num,
    maximum = statistic_num + 5*10^(-dp-1)
  )
  return(statistic_range)
}

#' Validate a p statistic
#'
#' Checks it is a string that can be converted to a number, and returns the
#' possible range of values before rounding
#'
#' @param p Character string or vector of strings. The p value(s) to validate
#' @param allow_na Logical. Allows the use of NA values
#' (returns NA for min/original/max)
#' @return A list of numbers: the possible range of values (minimum, original
#' and maximum)
#' @examples
#' validate_p(".14")
#' validate_p(c(".14", "<.001", NA, ".58"), allow_na = TRUE)
validate_p <- function(p, allow_na = FALSE)
{
  name <- substitute(p)
  if(typeof(name) != "symbol")
    name <- "p"

  if (!is.character(p)) {
    stop(paste(name, "must be a character string or character vector."))
  }

  if (any(is.na(p)) && !allow_na) {
    stop(paste(name, "contains NA; set allow_na = TRUE to permit NA values."))
  }

  # Deal with "<.X"
  numeric_parts <- vapply(p, function(elem) {
    if (is.na(elem)) return(NA_character_)
    sub(".*<", "", elem)
  }, FUN.VALUE = character(1L), USE.NAMES = FALSE)

  suppressWarnings({
    numeric_values <- as.numeric(numeric_parts)
  })

  if (any(is.na(numeric_values) & !is.na(numeric_parts))) {
    stop(paste(name, "must be parseable as numeric (or NA)."))
  }

  if(any(numeric_values > 1 | numeric_values < 0, na.rm = TRUE))
  {
    stop(paste(name, "values must be between 0 and 1."))
  }

  dp <- decimal_places(numeric_parts)

  p_range <- list(
    minimum = ifelse(grepl("<", p, fixed = TRUE),
                     .Machine$double.xmin,
                     pmax(0, numeric_values - 5*10^(-dp-1))),
    original = ifelse(grepl("<", p, fixed = TRUE),
                      numeric_values - .Machine$double.xmin,
                      numeric_values),
    maximum = ifelse(grepl("<", p, fixed = TRUE),
                     numeric_values,
                     pmin(1, numeric_values + 5*10^(-dp-1)))
  )

  return(p_range)
}
