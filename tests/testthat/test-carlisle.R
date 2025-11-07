test_that("Carlisle returns some values correctly", {
  res <- carlisle(c(".80", ".80", ".80", ".80"))
  expect_equal(res$Provided_p[2], 0.046164081)

  data <- data.frame(m1 = c("1.20", "1.25"), s1 = c("1.2", "1.25"), n1 = c(60, 60),
                     m2 = c("2.1", "2.15"), s2 = c("2.5", "2.55"), n2 = c(30, 30),
                     p = c(NA, NA))
  res <- t_slicer(m1 = m1, s1 = s1, n1 = n1, m2 = m2, s2 = s2, n2 = n2,
                  data = data) %>%
         carlisle()
  expect_equal(res$Welch[2], 0.9801026)
}
)

test_that("Carlisle output behave normally", {
  expect_silent(carlisle(c(".80", ".80", ".80", ".80")))
  expect_output(carlisle(c(".80", ".80", ".80", ".80"), output = TRUE))

  data <- data.frame(m1 = c("1.20", "1.25"), s1 = c("1.2", "1.25"), n1 = c(60, 60),
                     m2 = c("2.1", "2.15"), s2 = c("2.5", "2.55"), n2 = c(30, 30),
                     p = c(NA, NA))
  expect_output(t_slicer(m1 = m1, s1 = s1, n1 = n1, m2 = m2, s2 = s2, n2 = n2,
                  data = data) %>%
    carlisle(output = TRUE))

}
)
