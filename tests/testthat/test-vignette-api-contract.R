testthat::test_that("Vignette book API remains exported", {
  expected <- c("atlastable", "basisprod", "coat", "multicompplot", "multirec", "ptest", "rbfunc", "rec", "simbrain", "sizechange")
  testthat::expect_true(all(expected %in% getNamespaceExports("mand")))
})
