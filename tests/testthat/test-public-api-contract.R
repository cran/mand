testthat::test_that("book-facing API is exported", {
  expected <- c(
    "atlastable", "basisprod", "coat", "imgdatamat", "multicoat",
    "multicompplot", "multirec", "ptest", "rbfunc", "rec", "simbrain",
    "sizechange"
  )
  testthat::expect_true(all(expected %in% getNamespaceExports("mand")))
})

testthat::test_that("book-facing formal arguments are unchanged", {
  contract <- list(
    atlastable = c("x", "y", "atlasdataset", "ROIids", "..."),
    basisprod = c("A", "B"),
    coat = c("x", "y", "pseq", "xyz", "col.x", "col.y", "breaks.y",
      "zlim.x", "zlim.y", "rownum", "colnum", "plane", "xlab", "ylab",
      "axes", "oma", "mar", "bg", "paron", "cross.hair", "chxy",
      "color.bar", "regionplot", "atlasdataset", "regionname",
      "regionlegend", "atlasname", "ROIids", "..."),
    imgdatamat = c("imgfnames", "mask", "ROI", "atlas", "atlasdataset",
      "ROIids", "zeromask", "schange", "..."),
    multicoat = c("imgs", "y", "row4imp", "col4imp", "trm", "..."),
    multicompplot = c("object", "x", "comps", "row4comp", "col4comp",
      "pseq4comp", "..."),
    multirec = c("object", "imagedim", "B", "mask", "midx", "comps",
      "XY", "signflip"),
    ptest = c("object", "Z", "newdata", "testZ", "regmethod", "methods1",
      "metric", "number1", "repeats1", "params"),
    rbfunc = c("imagedim", "seppix", "hispec", "mask", "brainpos"),
    rec = c("Q", "imagedim", "B", "mask"),
    simbrain = c("baseimg", "diffimg", "sdevimg", "mask", "n0", "c1",
      "sd1", "rho", "zeromask", "reduce", "output", "seed"),
    sizechange = c("img1", "simscale", "refsize", "...")
  )
  for (name in names(contract)) {
    testthat::expect_identical(
      names(formals(getExportedValue("mand", name))), contract[[name]],
      info = name
    )
  }
})
