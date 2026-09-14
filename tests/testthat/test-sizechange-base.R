testthat::test_that("sizechange keeps the published formal arguments", {
  testthat::expect_identical(names(formals(sizechange)),
                             c("img1", "simscale", "refsize", "..."))
})

testthat::test_that("sizechange resizes arrays without imager", {
  x <- array(seq_len(3 * 4 * 5), c(3, 4, 5))
  linear <- sizechange(x, refsize = c(6, 8, 10))
  nearest <- sizechange(x, refsize = c(6, 8, 10),
                        interpolation = "nearest")
  testthat::expect_identical(dim(linear), c(6L, 8L, 10L))
  testthat::expect_identical(dim(nearest), c(6L, 8L, 10L))
  testthat::expect_true(all(is.finite(linear)))
  testthat::expect_true(all(nearest %in% x))
})

testthat::test_that("simscale remains supported", {
  x <- array(0, c(3, 4, 5))
  testthat::expect_identical(dim(sizechange(x, simscale = 2)),
                             c(6L, 8L, 10L))
})
