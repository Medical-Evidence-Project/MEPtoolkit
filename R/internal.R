#' @keywords internal
#' @noRd
resolve <- function(expr, data) {
  if(quo_is_null(expr)) {
    return(NULL)
  }
  vars <- names(data)
  if(as_string(ensym(expr)) %in% vars) {
    return(data %>% pull(!!expr))
  }
  else if(is.null(eval_tidy(expr))) {
    return(NULL)
  }
  else {
    stop(sprintf("Column `%s` not found in `data`", as_string(ensym(expr))))
  }
}

#' @keywords internal
#' @noRd
.strip_key <- function(x) {
  tolower(gsub("\\s+", " ", trimws(gsub("[^[:alnum:] ]+"," ", x))))
}

#' @keywords internal
#' @noRd
length_check <- function(args) {
  lens <- vapply(args, length, integer(1))
  L <- max(lens)
  ok_len <- lens == 0 | lens == L
  if (!all(ok_len))
    stop("Non-NULL arguments must have the same length.")
  return(L)
}
