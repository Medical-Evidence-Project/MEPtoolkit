test_that("The 4 combinations of output are possible", {
  res <- t_slicer("1.20", "1.2", 60, "2.1", "2.5", 30, ".02")
  expect_true(res$Student)
  expect_false(res$Welch)

  res <- t_slicer("1.20", "1.2", 60, "2.1", "2.5", 30, ".07")
  expect_false(res$Student)
  expect_true(res$Welch)

  res <- t_slicer("1.20", "1.2", 60, "2.1", "2.5", 60, ".015")
  expect_true(res$Student)
  expect_true(res$Welch)

  res <- t_slicer("1.20", "1.2", 60, "2.1", "2.5", 60, ".03")
  expect_false(res$Student)
  expect_false(res$Welch)
})

test_that("Works with both vectors and dataframes", {
  res <- t_slicer(c("1.20", "1.25"), c("1.2", "1.25"), c(60, 60),
                  c("2.1", "2.15"), c("2.5", "2.55"), c(30, 30), c(".02", ".07"))
  expect_equal(res$Student, c(TRUE, FALSE))
  expect_equal(res$Welch, c(FALSE, TRUE))

  data <- tibble(m1 = c("1.20", "1.25"), s1 = c("1.2", "1.25"), n1 = c(60, 60),
                 m2 = c("2.1", "2.15"), s2 = c("2.5", "2.55"), n2 = c(30, 30),
                 p = c(".02", ".07"))

  res <- t_slicer(m1 = m1, s1 = s1, n1 = n1, m2 = m2, s2 = s2, n2 = n2, p = p,
                  data = data)
  expect_equal(res$Student, c(TRUE, FALSE))
  expect_equal(res$Welch, c(FALSE, TRUE))
})

test_that("Function can deal with missing p values", {
  data <- tibble(m1 = c("1.20", "1.25"), s1 = c("1.2", "1.25"), n1 = c(60, 60),
                 m2 = c("2.1", "2.15"), s2 = c("2.5", "2.55"), n2 = c(30, 30),
                 p = c(".02", NA))

  res <- t_slicer(m1 = m1, s1 = s1, n1 = n1, m2 = m2, s2 = s2, n2 = n2, p = p,
                  data = data)
  expect_equal(res$Student, c(TRUE, NA))
  expect_equal(res$Welch, c(FALSE, NA))

  data <- tibble(m1 = c("1.20", "1.25"), s1 = c("1.2", "1.25"), n1 = c(60, 60),
                 m2 = c("2.1", "2.15"), s2 = c("2.5", "2.55"), n2 = c(30, 30), p = c(NA, NA))

  res <- t_slicer(m1 = m1, s1 = s1, n1 = n1, m2 = m2, s2 = s2, n2 = n2,
                  data = data)
  expect_equal(res$Student, c(NA, NA))
  expect_equal(res$Welch, c(NA, NA))
})

test_that("discrete values work properly", {
  res <- t_slicer("4", "1", 60, "2", "4", 60, ".05")
  expect_true(res$Student)
  expect_true(res$Welch)
})

test_that("A set of results is as expected", {
  res <- t_slicer("1.20", "1.2", 60, "2.1", "2.5", 30, ".02")
  expect_true(res$Student)
  expect_false(res$Welch)
  expect_equal(res$Student_RIVETS,
               BSDA::tsum.test(mean.x = 1.2, s.x = 1.2, n.x = 60, mean.y = 2.1,
                               s.y = 2.5, n.y = 30, var.equal = TRUE)$p.value)
  expect_equal(res$Welch_RIVETS,
               BSDA::tsum.test(mean.x = 1.2, s.x = 1.2, n.x = 60, mean.y = 2.1,
                               s.y = 2.5, n.y = 30, var.equal = FALSE)$p.value)
  expect_equal(res$Student_min,
               BSDA::tsum.test(mean.x = 1.195, s.x = 1.15, n.x = 60, mean.y = 2.15,
                               s.y = 2.45, n.y = 30, var.equal = TRUE)$p.value)
  expect_equal(res$Student_max,
               BSDA::tsum.test(mean.x = 1.205, s.x = 1.25, n.x = 60, mean.y = 2.05,
                               s.y = 2.55, n.y = 30, var.equal = TRUE)$p.value)
  expect_equal(res$Welch_min,
               BSDA::tsum.test(mean.x = 1.195, s.x = 1.15, n.x = 60, mean.y = 2.15,
                               s.y = 2.45, n.y = 30, var.equal = FALSE)$p.value)
  expect_equal(res$Welch_max,
               BSDA::tsum.test(mean.x = 1.205, s.x = 1.25, n.x = 60, mean.y = 2.05,
                               s.y = 2.55, n.y = 30, var.equal = FALSE)$p.value)
})

test_that("Works with edge cases (max/min possible p values within the range of
          the reported p value)", {
  res <- t_slicer("1.20", "1.2", 60, "2.2", "2.5", 30, ".06")
  expect_false(res$Student)
  expect_true(res$Welch)

  res <- t_slicer("1.20", "1.2", 60, "2.1", "2.5", 30, ".05")
  expect_false(res$Student)
  expect_true(res$Welch)

  res <- t_slicer("1.20", "1.2", 60, "2.2", "2.5", 30, ".02")
  expect_true(res$Student)
  expect_false(res$Welch)

  res <- t_slicer("1.20", "1.2", 60, "2.1", "2.5", 30, ".01")
  expect_true(res$Student)
  expect_false(res$Welch)

})

test_that("Outputs behave normally", {
  expect_silent(t_slicer("1.20", "1.2", 60, "2.2", "2.5", 30, ".06"))
  expect_silent(t_slicer("1.20", "1.2", 60, "2.2", "2.5", 30, ".06",
                         output = FALSE))
  expect_output(t_slicer("1.20", "1.2", 60, "2.2", "2.5", 30, ".06",
                         output = TRUE))
})
