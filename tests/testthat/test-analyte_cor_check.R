test_that("A result in different scenarios", {
  df <- tibble(analyte = "Sodium",
               m1 = "140.2",
               s1 = "2.10",
               m2 = "139.8",
               s2 = "2.05",
               n = 60,
               t = "2.50")
  out <- analyte_cor_check(
    analyte = "Sodium",
    m1 = "140.2", s1 = "2.10",
    m2 = "139.8", s2 = "2.05",
    n = 60, t = "2.50",
    k = 1, use = "pre"
  )
  expect_equal(out$disp$prob_total, 0.843225136)
  out <- analyte_cor_check(
    data = df,
    analyte = analyte,
    m1 = m1, s1 = s1,
    m2 = m2, s2 = s2,
    n = n, t = t,
    k = 1, use = "pre"
  )
  expect_equal(out$disp$prob_total, 0.843225136)
  out <- analyte_cor_check(
    data = df,
    analyte = analyte,
    m1 = m1, s1 = s1,
    m2 = m2, s2 = s2,
    n = n, t = t,
    k = 1, use = "pre", append = TRUE
  )
  expect_equal(out$prob_total, 0.843225136)
})
