testthat::test_that("internal heat palette has expected structure", {
  heat_palette <- getFromNamespace(".mand_hotmetal", "mand")
  palette64 <- heat_palette(64L)
  testthat::expect_length(palette64, 64L)
  testthat::expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", palette64)))
  rgb <- grDevices::col2rgb(palette64)
  testthat::expect_lt(rgb["red", 1L], rgb["red", 64L])
  testthat::expect_lt(rgb["green", 1L], rgb["green", 64L])
  testthat::expect_lt(rgb["blue", 1L], rgb["blue", 64L])
})

testthat::test_that("internal heat palette supports arbitrary size", {
  heat_palette <- getFromNamespace(".mand_hotmetal", "mand")
  testthat::expect_length(heat_palette(1L), 1L)
  testthat::expect_length(heat_palette(10L), 10L)
  testthat::expect_length(heat_palette(128L), 128L)
  testthat::expect_error(heat_palette(0L), "positive integer")
  testthat::expect_error(heat_palette(1.5), "positive integer")
})

testthat::test_that("internal heat palette remains internal", {
  testthat::expect_true(exists(
    ".mand_hotmetal", envir = asNamespace("mand"), inherits = FALSE
  ))
  testthat::expect_false(
    ".mand_hotmetal" %in% getNamespaceExports("mand")
  )
})

testthat::test_that("coat retains its published formal arguments", {
  testthat::expect_identical(
    names(formals(mand::coat)),
    c(
      "x", "y", "pseq", "xyz", "col.x", "col.y", "breaks.y",
      "zlim.x", "zlim.y", "rownum", "colnum", "plane", "xlab",
      "ylab", "axes", "oma", "mar", "bg", "paron", "cross.hair",
      "chxy", "color.bar", "regionplot", "atlasdataset",
      "regionname", "regionlegend", "atlasname", "ROIids", "..."
    )
  )
})
