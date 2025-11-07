test_that("resolve function can handle usual scenarios", {
  df <- tibble(m=c(2,1))
  df2 <- tibble(n=c(2,1))
  o <- c(1,2)
  p <- NULL
  test_fun <- function(data = NULL, m = NULL)
  {
    if(!is.null(data)) {
      if(!is.data.frame(data)) stop("`data` must be a data.frame.")
      m <- resolve(enquo(m), data)
    }
    return(m)
  }

  expect_null(test_fun(df))
  expect_null(test_fun(df, NULL))
  expect_null(test_fun(df, p))
  expect_error(test_fun(df, o), "Column `o` not found in `data`")
  expect_equal(test_fun(df, m), c(2,1))
  expect_equal(test_fun(df2, n), c(2,1))

})

test_that(".strip_key works on a few scenarios", {
  expect_equal(.strip_key("   analyte"), "analyte")
  expect_equal(.strip_key("   ANALYTE    "), "analyte")
})

test_that("length_check works correctly", {
  expect_equal(length_check(list(a=c(1,2), b=c(3,4))), 2)
  expect_equal(length_check(list(a=c(1,2), b=NULL)), 2)
  expect_error(length_check(list(a=c(1,2), b=c(3))), "Non-NULL arguments must have the same length.")
})
