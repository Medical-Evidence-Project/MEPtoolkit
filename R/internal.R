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
