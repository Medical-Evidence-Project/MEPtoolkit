test_that("A result in different scenarios", {
  df <- tibble(analyte = "25-hydroxy-Vitamin D",
               m1 = "64.3",
               s1 = "25.4",
               m2 = "88.5",
               s2 = "23.2")
  res <- analyte_max_cor("25-hydroxy-Vitamin D",
                  m1 = "64.3", s1 = "25.4",
                  m2 = "88.5", s2 = "23.2",
                  k = 1, use = "pre")
  expect_equal(res$k, 1)
  expect_equal(res$use, "pre")
  expect_equal(res$disp$r_max_analytic, 0.987247552)
  expect_equal(res$disp$r_max_total, 0.895884447)
  res <- analyte_max_cor(data = df,
                         analyte = analyte,
                         m1 = m1, s1 = s1,
                         m2 = m2, s2 = s2,
                         k = 1, use = "pre")
  expect_equal(res$k, 1)
  expect_equal(res$use, "pre")
  expect_equal(res$disp$r_max_analytic, 0.987247552)
  expect_equal(res$disp$r_max_total, 0.895884447)
  res <- analyte_max_cor(data = df,
                         analyte = analyte,
                         m1 = m1, s1 = s1,
                         m2 = m2, s2 = s2,
                         k = 1, use = "pre", append = TRUE)
  expect_equal(res$r_max_analytic, 0.987247552)
  expect_equal(res$r_max_total, 0.895884447)
})

