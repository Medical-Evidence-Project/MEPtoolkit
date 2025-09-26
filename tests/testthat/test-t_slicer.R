test_that("The 4 combinations of output are possible", {
  expect_identical(t_slicer("1.20", "1.2", 60, "2.1", "2.5", 30, ".02"),
                   list(Student=TRUE, Welch=FALSE))
  expect_identical(t_slicer("1.20", "1.2", 60, "2.1", "2.5", 30, ".07"),
                   list(Student=FALSE, Welch=TRUE))
  expect_identical(t_slicer("1.20", "1.2", 60, "2.1", "2.5", 60, ".015"),
                   list(Student=TRUE, Welch=TRUE))
  expect_identical(t_slicer("1.20", "1.2", 60, "2.1", "2.5", 60, ".03"),
                   list(Student=FALSE, Welch=FALSE))
})

test_that("dicrete values work properly", {
  expect_identical(t_slicer("4", "1", 60, "2", "4", 60, ".05"),
                   list(Student=TRUE, Welch=TRUE))
})

test_that("Works with edge cases (max/min possible p values within the range of
          the reported p value)", {
  expect_identical(t_slicer("1.20", "1.2", 60, "2.2", "2.5", 30, ".06"),
                   list(Student=FALSE, Welch=TRUE))
  expect_identical(t_slicer("1.20", "1.2", 60, "2.1", "2.5", 30, ".05"),
                   list(Student=FALSE, Welch=TRUE))
  expect_identical(t_slicer("1.20", "1.2", 60, "2.2", "2.5", 30, ".02"),
                   list(Student=TRUE, Welch=FALSE))
  expect_identical(t_slicer("1.20", "1.2", 60, "2.1", "2.5", 30, ".01"),
                   list(Student=TRUE, Welch=FALSE))
})
