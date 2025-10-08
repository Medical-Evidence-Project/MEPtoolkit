test_that("demical_places returns the right number", {
  expect_equal(decimal_places("1.234"), 3)
  expect_equal(decimal_places("3"), 0)
  expect_equal(decimal_places(".2"), 1)
})

test_that("validate_descriptive correctly stops when not a number+string (single
          values)", {
  mean <- 1.6
  sd <- "1.6c"
  expect_error(validate_descriptive(mean),
               "mean must be a character string or character vector.")
  expect_error(validate_descriptive(sd),
               "sd must be parseable as numeric \\(or NA\\).")
  expect_error(validate_descriptive("1.6c"),
               "statistic must be parseable as numeric \\(or NA\\).")
})

test_that("validate_descriptive correctly stops when not a number+string (vector
          of values)", {
            mean <- c(1.6, 3.5)
            sd <- c("1.5", "1.6c")
            expect_error(validate_descriptive(mean),
                         "mean must be a character string or character vector.")
            expect_error(validate_descriptive(sd),
                         "sd must be parseable as numeric \\(or NA\\).")
            expect_error(validate_descriptive(c("1.5", "1.6c")),
                         "statistic must be parseable as numeric \\(or NA\\).")
          })

test_that("validate_descriptive returns the right numbers (single values)", {
  expect_equal(validate_descriptive("1.6"),
               list(minimum=1.55, original=1.6, maximum=1.65))
  expect_equal(validate_descriptive("2"),
               list(minimum=1.5, original=2, maximum=2.5))
})

test_that("validate_descriptive returns the right numbers (vector of values)", {
  expect_equal(validate_descriptive(c("1.6", "2")),
               list(minimum=c(1.55, 1.5), original=c(1.6, 2),
                    maximum=c(1.65, 2.5)))
  expect_equal(validate_descriptive(c("1.6", NA, "2"), allow_na = TRUE),
               list(minimum=c(1.55, NA, 1.5), original=c(1.6, NA, 2),
                    maximum=c(1.65, NA, 2.5)))
})

test_that("validate_descriptive stops for non-allowed NA values", {
  expect_error(validate_descriptive(c("1.6", NA, "2")),
              "statistic contains NA; set allow_na = TRUE to permit NA values.")
  vec <- c("1.6", NA, "2")
  expect_error(validate_descriptive(vec),
               "vec contains NA; set allow_na = TRUE to permit NA values.")
})

test_that("validate_p correctly stops when not a number+string (single
          values)", {
  p_value <- .6
  expect_error(validate_p(.6),
               "p must be a character string or character vector.")
  expect_error(validate_p(".6c"),
               "p must be parseable as numeric \\(or NA\\).")
  expect_error(validate_p(p_value),
               "p_value must be a character string or character vector.")
})

test_that("validate_p correctly stops when not a number+string (vector of
          values)", {
            p_v <- c(.6, .7, .8)
            expect_error(validate_p(p_v),
                         "p_v must be a character string or character vector.")
            p_v <- c(".5", ".5c", ".8")
            expect_error(validate_p(p_v),
                         "p_v must be parseable as numeric \\(or NA\\).")
            expect_error(validate_p(c(.6, .7)),
                         "p must be a character string or character vector.")
})

test_that("validate_p correctly stops when 0<p<1 (single values)", {
  expect_error(validate_p("1.6"),
               "p values must be between 0 and 1.")
  expect_error(validate_p("-0.6"),
               "p values must be between 0 and 1.")
})

test_that("validate_p correctly stops when 0<p<1 (vector of values)", {
  expect_error(validate_p(c(".6", "-.6")),
               "p values must be between 0 and 1.")
  expect_error(validate_p(c(".6", "1.6")),
               "p values must be between 0 and 1.")
})

test_that("validate_p returns the right numbers (single values)", {
  expect_equal(validate_p("0.6"),
               list(minimum=0.55, original=0.6, maximum=0.65))
  expect_equal(validate_p(".2"),
               list(minimum=.15, original=.2, maximum=.25))
})

test_that("validate_p returns the right numbers (vector of values)", {
  expect_equal(validate_p(c("0.6", ".255")),
               list(minimum=c(.55, .2545), original=c(0.6, .255),
                    maximum=c(0.65, .2555)))
})

test_that("validate_p works with <.X (single values)", {
  expect_equal(validate_p("<.001"),
               list(minimum=.Machine$double.xmin, original=.001, maximum=.001))
  expect_equal(validate_p("<.05"),
               list(minimum=.Machine$double.xmin, original=.05, maximum=.05))
})

test_that("validate_p works with <.X (vector of values)", {
  expect_equal(validate_p(c("<.001", ".6", "<.01")),
               list(minimum=c(.Machine$double.xmin, .55, .Machine$double.xmin),
                    original=c(.001 - .Machine$double.xmin, .6,
                               .01 - .Machine$double.xmin),
                    maximum=c(.001, .65, .01)))
})

test_that("validate_p stops for non-allowed NA values", {
  expect_error(validate_p(c(".6", NA, ".1")),
               "p contains NA; set allow_na = TRUE to permit NA values.")
  vec <- c(".6", NA, ".6")
  expect_error(validate_p(vec),
               "vec contains NA; set allow_na = TRUE to permit NA values.")
})

test_that("validate_p works with a variety of values", {
  expect_equal(validate_p(c(".6", NA, "<.001"), allow_na = TRUE),
               list(minimum=c(.55, NA, .Machine$double.xmin),
                    original=c(.6, NA, .001 - .Machine$double.xmin),
                    maximum=c(.65, NA, .001)))
})
