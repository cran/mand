testthat::test_that("published optional-dependency API signatures remain unchanged", {
  testthat::expect_identical(names(formals(ptest)),
    c("object", "Z", "newdata", "testZ", "regmethod", "methods1",
      "metric", "number1", "repeats1", "params"))
  testthat::expect_identical(names(formals(simbrain)),
    c("baseimg", "diffimg", "sdevimg", "mask", "n0", "c1", "sd1",
      "rho", "zeromask", "reduce", "output", "seed"))
})

testthat::test_that("simbrain works with rho NULL", {
  baseimg <- array(1, c(2, 2, 2))
  diffimg <- array(0, c(2, 2, 2))
  out <- simbrain(baseimg, diffimg, n0 = 2, rho = NULL, seed = 1)
  testthat::expect_identical(dim(out$S), c(4L, 8L))
})
