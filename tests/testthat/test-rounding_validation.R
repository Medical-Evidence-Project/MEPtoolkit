test_that("demical_places returns the right number", {
  expect_equal(decimal_places("1.234"), 3)
  expect_equal(decimal_places("3"), 0)
  expect_equal(decimal_places(".2"), 1)
})

test_that("validate_descriptive correctly stops when not a number+string", {
  mean <- 1.6
  sd <- "1.6c"
  expect_error(validate_descriptive(mean),
               "mean must be a character string.")
  expect_error(validate_descriptive(sd),
               "sd must be convertible to a number.")
  expect_error(validate_descriptive("1.6c"),
               "statistic must be convertible to a number.")
})

test_that("validate_descriptive returns the right numbers", {
  expect_equal(validate_descriptive("1.6"),
               list(minimum=1.55, original=1.6, maximum=1.65))
  expect_equal(validate_descriptive("2"),
               list(minimum=1.5, original=2, maximum=2.5))
})

test_that("validate_p correctly stops when not a number+string", {
  p_value <- .6
  expect_error(validate_p(.6),
               "p must be a character string.")
  expect_error(validate_p(".6c"),
               "p must be convertible to a number.")
  expect_error(validate_p(p_value),
               "p_value must be a character string.")
})

test_that("validate_p correctly stops when 0<p<1", {
  expect_error(validate_p("1.6"),
               "p must be between 0 and 1.")
  expect_error(validate_p("-0.6"),
               "p must be between 0 and 1.")
})

test_that("validate_p returns the right numbers", {
  expect_equal(validate_p("0.6"),
               list(minimum=0.55, original=0.6, maximum=0.65))
  expect_equal(validate_p(".2"),
               list(minimum=.15, original=.2, maximum=.25))
})

test_that("validate_p works with <.X", {
  expect_equal(validate_p("<.001"),
               list(minimum=.Machine$double.xmin, original=.001, maximum=.001))
  expect_equal(validate_p("<.05"),
               list(minimum=.Machine$double.xmin, original=.05, maximum=.05))
})
