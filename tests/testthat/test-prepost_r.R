test_that("correctly stops when data is given but not a data frame", {
  expect_error(prepost_r(data = c(1, 2), s1 = c("2.10", "0.90"),
                         s2 = c("2.05", "0.90"), n = c(10, 20),
                         t  = c("3.00", "2.40")),
               "`data` must be a data.frame.")
  expect_error(prepost_r(data = matrix(1, 2), s1 = c("2.10","0.90"),
                         s2 = c("2.05", "0.90"), n = c(10, 20),
                         t  = c("3.00", "2.40")),
               "`data` must be a data.frame.")
})

test_that("correctly stops when n is missing", {
  expect_error(prepost_r(s1 = c("2.10", "0.90"), s2 = c("2.05", "0.90"),
                         t  = c("3.00", "2.40")),
               "`n` \\(paired sample size\\) is required for every row.")

  df <- data.frame(m1 = c("140.2","6.60"),
                   s1 = c("2.10","0.90"),
                   m2 = c("139.8","6.40"),
                   s2 = c("2.05","0.90"),
                   n = c(10, NA),
                   t  = c("3.00","2.40"),
                   stringsAsFactors = FALSE)

  expect_error(prepost_r(data = df, s1 = s1, s2 = s2, t=t, n=n),
               "`n` \\(paired sample size\\) is required for every row.")
})

test_that("correctly stops when dataframe argument is wrong", {
  df <- data.frame(m1 = c("140.2","6.60"),
                   s1 = c("2.10","0.90"),
                   m2 = c("139.8","6.40"),
                   s2 = c("2.05","0.90"),
                   n = c(10, NA),
                   t  = c("3.00","2.40"),
                   stringsAsFactors = FALSE)
  tst <- 1

  expect_error(prepost_r(data = df, s1 = s1, s2 = s2, t=tst, n=n),
               "Column `tst` not found in `data`")
})

test_that("correctly stops when n<1", {
  expect_error(prepost_r(s1 = c("2.10", "0.90"), s2 = c("2.05", "0.90"),
                         n = c(10, 1), t  = c("3.00", "2.40")),
               "`n` must be >1.")
})

test_that("Stops when s1 or s2 not present", {
  expect_error(prepost_r(s1 = NULL,
                         s2 = c("2.05", "0.90"),
                         n = c(10, 20),
                         t  = c("3.00", "2.40")),
               "`s1` and `s2` must be non-null.")
})

test_that("Stops when arguments not of the same length", {
  expect_error(prepost_r(s1 = c("2.10", "0.90", "5.30"),
                         s2 = c("2.05", "0.90"),
                         n = c(10, 20),
                         t  = c("3.00", "2.40")),
               "Non-NULL arguments must have the same length.")
})

test_that("Stops when t==0", {
  expect_error(prepost_r(m1 = c("140.2","6.60"),
                         s1 = c("2.10", "0.90"),
                         m2 = c("139.8","6.40"),
                         s2 = c("2.05", "0.90"),
                         n = c(10, 20),
                         t  = c("-3.00", "0.00")),
               "t must be different from 0 \\(line 2\\).")
})

test_that("Stops when p==1 or p==0", {
  expect_error(prepost_r(m1 = c("140.2","6.60"),
                         s1 = c("2.10", "0.90"),
                         m2 = c("139.8","6.40"),
                         s2 = c("2.05", "0.90"),
                         n = c(10, 20),
                         p  = c(".83", "1.00")),
               "p must be >0 and <1 \\(line 2\\).")

  expect_error(prepost_r(m1 = c("140.2","6.60"),
                         s1 = c("2.10", "0.90"),
                         m2 = c("139.8","6.40"),
                         s2 = c("2.05", "0.90"),
                         n = c(10, 20),
                         p  = c("0.00", ".50")),
               "p must be >0 and <1 \\(line 1\\).")
})

test_that("Stops when SD<=0", {
  expect_error(prepost_r(m1 = c("140.2","6.60"),
                         s1 = c("-0.5", "0.90"),
                         m2 = c("139.8","6.40"),
                         s2 = c("2.05", "0.90"),
                         n = c(10, 20),
                         p  = c(".83", ".52")),
               "SDs must be >0 \\(line 1\\).")
})

test_that("Stops/warns when not enough info (strict mode or not)", {
  # Rules are:
  # - s1 & s2 always necessary (see tests above), can't be recalculated
  # - if SC present, we can calculate rho directly, otherwise need t or p
  # - if using t or p, MC is required:
  #   - if MC is present, we can calculate rho using t or p
  #   - if MC is absent, we need m1 and m2 to calculate MC
  expect_error(prepost_r(s1 = c(NA,  "0.90", "2.10", "0.90", "2.10", "0.90"),
                         s2 = c("2.05",  "0.90", "2.05", "0.90", "2.05", "0.90"),
                         sc = c("2.20",  NA,     NA,     NA,     NA,     NA),
                         m1 = c("140.2", "6.60", NA,     "6.60", NA,     NA),
                         m2 = c("139.8", "6.40", NA,     "6.40", NA,     NA),
                         mc = c("1.00",  "2.00", "1.00", NA,     NA,     NA),
                         t = c(NA,       NA,     "1.00", "1.00", "1.00", NA),
                         p = c(NA,       NA,     NA,     NA,     NA,     ".025"),
                         n = c(20, 20, 20, 20, 20, 20),
                         strict = TRUE),
               "Insufficient information for rows: 1, 2, 5, 6.")

  expect_warning(prepost_r(s1 = c(NA,  "0.90", "2.10", "0.90", "2.10", "0.90"),
                         s2 = c("2.05",  "0.90", "2.05", "0.90", "2.05", "0.90"),
                         sc = c("2.20",  NA,     NA,     NA,     NA,     NA),
                         m1 = c("140.2", "6.60", NA,     "6.60", NA,     NA),
                         m2 = c("139.8", "6.40", NA,     "6.40", NA,     NA),
                         mc = c("1.00",  "2.00", "1.00", NA,     NA,     NA),
                         t = c(NA,       NA,     "1.00", "1.00", "1.00", NA),
                         p = c(NA,       NA,     NA,     NA,     NA,     ".025"),
                         n = c(20, 20, 20, 20, 20, 20)),
               "Insufficient information for rows set to NA: 1, 2, 5, 6.")
})

test_that("A few calculations are correct", {
  df <- tibble(s1 = c("1.03"),
               s2 = c("0.83"),
               stt = c("0.94"),
               n = c(100))
  df <- prepost_r(data = df, s1 = s1, s2 = s2, sc = stt,
            n = n, append = TRUE)
  expect_equal(round(df$rho_prepost_min, 3), 0.496)
  expect_equal(round(df$rho_prepost, 3), 0.507)
  expect_equal(round(df$rho_prepost_max, 3), 0.517)

  df <- tibble(s1 = c("0.94"),
               s2 = c("1.02"),
               mc = c("0.07"),
               t = c("0.74"),
               n = c(100))
  df <- prepost_r(data = df, s1 = s1, s2 = s2, mc = mc,
                  n = n, t = t, append = TRUE)
  expect_equal(round(df$rho_prepost_min, 3), 0.455)
  expect_equal(round(df$rho_prepost, 3), 0.537)
  expect_equal(round(df$rho_prepost_max, 3), 0.610)

  df <- data.frame(s1 = c("1.05"),
                   s2 = c("1.03"),
                   m1 = c("-0.09"),
                   m2 = c("-0.01"),
                   t = c("0.68"),
                   n = c(100))
  df <- prepost_r(data = df, s1 = s1, s2 = s2, m1 = m1, m2 = m2,
                  t = t, n = n, append = TRUE)
  expect_equal(round(df$rho_prepost_min, 3), 0.170)
  expect_equal(round(df$rho_prepost, 3), 0.360)
  expect_equal(round(df$rho_prepost_max, 3), 0.522)

  df <- data.frame(s1 = c("1.05"),
                   s2 = c("1.03"),
                   mc = c("0.08"),
                   p = c("0.498"),
                   n = c(100))
  df <- prepost_r(data = df, s1 = s1, s2 = s2, mc = mc,
                  n = n, p = p, append = TRUE)
  expect_equal(round(df$rho_prepost_min, 3), 0.269)
  expect_equal(round(df$rho_prepost, 3), 0.361)
  expect_equal(round(df$rho_prepost_max, 3), 0.445)

  df <- data.frame(s1 = c("1.05"),
                   s2 = c("1.03"),
                   m1 = c("-0.09"),
                   m2 = c("-0.01"),
                   p = c("0.498"),
                   n = c(100))
  df <- prepost_r(data = df, s1 = s1, s2 = s2, m1 = m1, m2 = m2,
                  n = n, p = p, append = TRUE)
  expect_equal(round(df$rho_prepost_min, 3), 0.181)
  expect_equal(round(df$rho_prepost, 3), 0.361)
  expect_equal(round(df$rho_prepost_max, 3), 0.516)
})
