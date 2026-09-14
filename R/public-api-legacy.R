#' Result report with atlas data 
#'
#' This function refers to the results obtained by the analysis in an atlas image, and reports a summary of the results for each anatomical region.
#' 
#' \code{atlastable} requires the atlas image and data frame including the ROI id and the name.
#'
#' @name atlastable
#' @aliases atlastable
#' @rdname atlastable
#' @docType methods
#' @export
#'
#' @param x A three-dimensional atlas image.
#' @param y an array for the result image.
#' @param atlasdataset a matrix or data.frame. The colnames should include "ROIid" and "ROIname".
#' @param ROIids a vector indicating ROI id shown in the result.
#' @param ... further arguments passed to or from other methods.
#' 
#' @examples
#' data(diffimg)
#' data(atlasdatasets)
#' data(atlas)
#' atlasname = "aal3"
#' atlasdataset = atlasdatasets[[atlasname]]
#' tmpatlas = atlas[[atlasname]]
#' atlastable(tmpatlas, diffimg, atlasdataset=atlasdataset, ROIids = c(1:2, 41:44))
atlastable <- function(x, y, atlasdataset = NULL, ROIids = NULL, ...) {
    if (is.null(atlasdataset)) 
        stop("atlasdataset should be required")
    if (is.null(ROIids)) 
        ROIids = atlasdataset$ROIid
    atlasdataset1 = atlasdataset[as.character(atlasdataset$ROIid) %in% as.character(ROIids), ]
    stats = do.call(rbind, lapply(ROIids, function(ROIid1) {
        x1 = (x == ROIid1)
        tmp = y * x1
        sizenum = sum(tmp != 0)
        sumvalue = sum(tmp)
        sizepct = sizenum/sum(x1)
        c(sizenum = sizenum, sizepct = sizepct, sumvalue = sumvalue, summary(tmp))
    }))
    out = list(table = cbind(atlasdataset1, stats))
    class(out) = "atlastable"
    out
}

#' Product Radial Basis Function
#'
#' This is a function to product the output for the rbfunc function with data matrix for a dimension reduction.
#' 
#' \code{basisprod} requires one list and one matrix.
#'
#' @name basisprod
#' @aliases basisprod
#' @rdname basisprod
#' @docType methods
#' @export
#'
#' @param A a list or a matrix correponding to the output for the \code{rbfunc} function with the argument hispec=FALSE or data matrix, respectivey. 
#' @param B a list or a matrix.
#' 
#' @examples
#' 
#' imagedim1=c(10,10,10)
#' 
#' B1 = rbfunc(imagedim=imagedim1, seppix=4, hispec=TRUE)
#' B2 = rbfunc(imagedim=imagedim1, seppix=4, hispec=FALSE)
#' 
#' n = 50
#' S = matrix(rnorm(n*prod(imagedim1)), nrow = n, ncol = prod(imagedim1))
#' 
#' SB1 = S %*% B1
#' SB12 = tcrossprod(S, t(B1))
#' all(SB1-SB12 == 0)
#' 
#' SB2 = basisprod(S, B2)
#' all(SB1-SB2 == 0)
#' 
#' BS1 = t(B1) %*% t(S)
#' BS2 = basisprod(B2, S)
#' all(BS1-t(BS2) == 0)
#' 
basisprod <- function(A, B) {
    if (inherits(B, "list")) {
        do.call(cbind, lapply(B, function(b2) {
            A[, b2[, 1]] %*% b2[, 2]
        }))
    }
    else if (inherits(A, "list")) {
        do.call(cbind, lapply(A, function(b2) {
            B[, b2[, 1]] %*% b2[, 2]
        }))
    }
    else {
        stop("Either A or B should be output of rbfunc")
    }
}

#' Coat Function
#'
#' This is a function for plotting an image. The analysis result can be overcoated on the template.
#' 
#' \code{coat} requires a image array.
#'
#' @name coat
#' @aliases coat
#' @rdname coat
#' @docType methods
#' @export
#'
#' @param x A three-dimensional background image.
#' @param y image2 to be overcoated.
#' @param pseq a vector plot sequence.
#' @param xyz a vector position to be plotted.
#' @param col.x a color vector for image1.
#' @param col.y a color vector for image2.
#' @param breaks.y a vector breaks value for y.
#' @param zlim.x a vector plot limitation values for z of x.
#' @param zlim.y a vector plot limitation values for z of y.
#' @param rownum a numeric, the number of row for the plot.
#' @param colnum a numeric, the number of colnum for the plot.
#' @param plane a vector plot sequence.
#' @param xlab a character for a label in the x axis.
#' @param ylab a character for a label in the y axis.
#' @param axes a logical. TRUE presents the axes.
#' @param oma a vector for outer margin area.
#' @param mar a vector for margin.
#' @param bg a character for color of background.
#' @param paron a logical. TRUE means par is used.
#' @param cross.hair a logical.
#' @param chxy a vector cross hair position to be plotted.
#' @param color.bar a logical.
#' @param regionplot a logical.
#' @param atlasdataset a matrix or data.frame. colnames shold include "ROIid" and "ROIname".
#' @param regionname a character.
#' @param regionlegend a logical.
#' @param atlasname a character.
#' @param ROIids a vector
#' @param ... further arguments passed to or from other methods.
#' 
#' @examples
#' 
#' data(exbrain)
#' coat(exbrain)
#' 
coat <- function(x, y = NULL, pseq = NULL, xyz = NULL, col.x = gray(0:64/64), col.y = NULL, breaks.y = NULL, zlim.x = NULL, 
    zlim.y = NULL, rownum = 5, colnum = NULL, plane = c("axial", "coronal", "sagittal", "all")[1], xlab = "", ylab = "", 
    axes = FALSE, oma = rep(0, 4), mar = rep(0, 4), bg = "black", paron = TRUE, cross.hair = FALSE, chxy = NULL, color.bar = TRUE, 
    regionplot = FALSE, atlasdataset = NULL, regionname = c("atlas", "stat")[1], regionlegend = FALSE, atlasname = "", ROIids = 1:9, 
    ...) {
    if (paron) {
        oldpar <- par(no.readonly = TRUE)
        on.exit(par(oldpar))
    }
    else {
        regionlegend = FALSE
    }
    if (is.null(zlim.x)) 
        zlim.x = range(x, na.rm = TRUE)
    if (!is.null(y)) {
        if (is.null(zlim.y)) 
            zlim.y = range(y, na.rm = TRUE)
    }
    dim.x = dim(x)
    X <- dim.x[1]
    Y <- dim.x[2]
    Z <- dim.x[3]
    if (plane != "all") {
        pseqmax = switch(plane, axial = Z, coronal = Y, sagittal = X)
        if (is.null(pseq)) 
            pseq = rev(seq(1, pseqmax, length = 30))
    }
    rownum = min(c(rownum, length(pseq)))
    if (is.null(rownum) & !is.null(colnum)) 
        rownum = ceiling(length(pseq)/colnum)
    if (is.null(colnum)) 
        colnum = ceiling(length(pseq)/rownum)
    plotnum = rownum * colnum
    breaks.x <- seq(min(x, zlim.x, na.rm = TRUE), max(x, zlim.x, na.rm = TRUE), length = length(col.x) + 1)
    if (!is.null(y)) {
        if (is.null(col.y)) {
            col.y = .mand_hotmetal()
            prescol.y = FALSE
        }
        else {
            prescol.y = TRUE
        }
        if (!is.na(col.y[1])) 
            col.y = c(NA, col.y)
        col.y2 = col.y
        if (is.null(breaks.y)) 
            breaks.y <- seq(zlim.y[1], zlim.y[2], length = length(col.y) + 1)
        if (any(sign(breaks.y) < 0)) {
            negidxmax = max(which(sign(breaks.y) < 0))
            if (negidxmax < length(breaks.y)) {
                if (!prescol.y) {
                  breaks.y = c(breaks.y[1:negidxmax], 0, breaks.y[(negidxmax + 1):length(breaks.y)])
                  col.y = c(rev(tim.colors2(negidxmax - 1)), NA, .mand_hotmetal(length(breaks.y) - (negidxmax + 1)))
                  col.y2 = c(rev(tim.colors2(negidxmax)), .mand_hotmetal(length(breaks.y) - (negidxmax + 1)))
                }
            }
        }
    }
    if (is.null(y)) 
        color.bar = FALSE
    if (regionplot) {
        if (is.null(atlasdataset)) 
            stop("atlasdataset should be required")
        atlasdataset1 = atlasdataset[as.character(atlasdataset$ROIid) %in% as.character(ROIids), ]
        if (regionname == "stat") {
            breaks.y = c(0, sort(ROIids) - 0.5, max(ROIids) + 1)
        }
        else {
            y[!(y %in% ROIids)] = 0
            breaks.y = c(0, sort(ROIids) - 0.5, max(ROIids) + 1)
        }
        colors2 = c("red", "green3", "orange", "magenta", "blue", "deeppink", "darkgreen", "yellow", "brown", "cyan")
        col.y = c(NA, colors2[1:(length(ROIids))])
        color.bar = FALSE
    }
    if (plane == "axial") {
        if (paron) {
            par(oma = oma, mar = mar, bg = bg)
            if (color.bar | regionlegend) {
                layout(rbind(matrix(1:plotnum, ncol = colnum, byrow = TRUE), rep(plotnum + 1, colnum)))
            }
            else {
                layout(matrix(1:plotnum, ncol = colnum, byrow = TRUE))
            }
        }
        for (z2 in pseq) {
            graphics::image(1:X, 1:Y, x[, , z2], col = col.x, breaks = breaks.x, zlim = zlim.x, axes = axes, xlab = xlab, 
                ylab = ylab, ...)
            if (!is.null(y)) {
                graphics::image(1:X, 1:Y, y[, , z2], col = col.y, breaks = breaks.y, zlim = zlim.y, add = TRUE)
            }
        }
        if (color.bar) {
            colbar(zlim = zlim.y, col.y = col.y2, nticks = 4, horizontal = TRUE)
        }
        else if (regionlegend) {
            par(oma = c(1, 0, 4, 0), bg = bg, new = TRUE)
            plot(1, type = "n", axes = FALSE, ann = FALSE)
            mtext(atlasname, col = "white", cex = 2, line = -3)
            graphics::legend("left", legend = abbreviate(atlasdataset1$ROIname, 10), col = col.y[-1], ncol = 5, pch = 15, 
                text.font = 1, text.col = col.y[-1])
        }
    }
    else if (plane == "sagittal") {
        if (paron) 
            par(mfrow = c(rownum, colnum), oma = oma, mar = mar, bg = bg)
        for (x2 in pseq) {
            graphics::image(1:Y, 1:Z, x[x2, , ], col = col.x, breaks = breaks.x, zlim = zlim.x, axes = axes, xlab = xlab, 
                ylab = ylab)
            if (!is.null(y)) {
                graphics::image(1:Y, 1:Z, y[x2, , ], col = col.y, breaks = breaks.y, zlim = zlim.y, add = TRUE)
            }
        }
    }
    else if (plane == "coronal") {
        if (paron) 
            par(mfrow = c(rownum, colnum), oma = oma, mar = mar, bg = bg)
        for (y2 in pseq) {
            graphics::image(1:X, 1:Z, x[, y2, ], col = col.x, breaks = breaks.x, zlim = zlim.x, axes = axes, xlab = xlab, 
                ylab = ylab)
            if (!is.null(y)) {
                graphics::image(1:X, 1:Z, y[, y2, ], col = col.y, breaks = breaks.y, zlim = zlim.y, add = TRUE)
            }
        }
    }
    else if (plane == "all") {
        if (is.null(xyz)) 
            xyz = c(round(X/2), round(Y/2), round(Z/2))
        if (paron) 
            par(mfrow = c(1, 3), oma = oma, mar = c(0, 0, 3, 0), bg = "black")
        graphics::image(1:X, 1:Y, x[, , xyz[3]], col = col.x, breaks = breaks.x, zlim = zlim.x, axes = axes, xlab = xlab, 
            ylab = ylab, main = "axial", col.main = "white")
        if (!is.null(y)) {
            graphics::image(1:X, 1:Y, y[, , xyz[3]], col = col.y, breaks = breaks.y, zlim = zlim.y, add = TRUE)
            if (cross.hair) {
                if (is.null(chxy)) 
                  chxy = which(y[, , xyz[3]] == max(y[, , xyz[3]]), arr.ind = TRUE)[1, ]
                abline(v = chxy[1], col = "green")
                abline(h = chxy[2], col = "green")
            }
        }
        graphics::image(1:Y, 1:Z, x[xyz[1], , ], col = col.x, breaks = breaks.x, zlim = zlim.x, axes = axes, xlab = xlab, 
            ylab = ylab, main = "sagittal", col.main = "white")
        if (!is.null(y)) {
            graphics::image(1:Y, 1:Z, y[xyz[1], , ], col = col.y, breaks = breaks.y, zlim = zlim.y, add = TRUE)
            if (cross.hair) {
                if (is.null(chxy)) 
                  chxy = which(y[xyz[1], , ] == max(y[xyz[1], , ]), arr.ind = TRUE)[1, ]
                abline(v = chxy[1], col = "green")
                abline(h = chxy[2], col = "green")
            }
        }
        graphics::image(1:X, 1:Z, x[, xyz[2], ], col = col.x, breaks = breaks.x, zlim = zlim.x, axes = axes, xlab = xlab, 
            ylab = ylab, main = "coronal", col.main = "white")
        if (!is.null(y)) {
            graphics::image(1:X, 1:Z, y[, xyz[2], ], col = col.y, breaks = breaks.y, zlim = zlim.y, add = TRUE)
            if (cross.hair) {
                if (is.null(chxy)) 
                  chxy = which(y[, xyz[2], ] == max(y[, xyz[2], ]), arr.ind = TRUE)[1, ]
                abline(v = chxy[1], col = "green")
                abline(h = chxy[2], col = "green")
            }
        }
    }
}

#' Creat Data Matrix Function
#'
#' This is a function that creates a data matrix for analysis from a file saved in image format.
#' 
#' \code{imgdatamat} requires image file names.
#'
#' @name imgdatamat
#' @aliases imgdatamat
#' @rdname imgdatamat
#' @docType methods
#' @export
#'
#' @param atlas An optional three-dimensional atlas image.
#' @param imgfnames a vector for (nifti) file names to be used.
#' @param mask a vector for brain mask data.
#' @param ROI a logical for roi data set.
#' @param atlasdataset a matrix or data.frame. colnames shold include "ROIid" and "ROIname".
#' @param ROIids a vector
#' @param zeromask a logical for masking voxel with all zeros.
#' @param schange a logical for change dimension.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return \item{S}{data matrix}
#' @return \item{brainpos}{binary brain position.}
#' @return \item{imagedim}{three dimensional vector for image dimension}
#' 
#' @examples
#' # imgfnames1 = c("img1.nii", "img2.nii")
#' # imgdata = imgdatamat(imgfnames1)
#' 
imgdatamat <- function(imgfnames, mask = NULL, ROI = FALSE, atlas = NULL, atlasdataset = NULL, ROIids = NULL, zeromask = FALSE, 
    schange = FALSE, ...) {
    if (length(imgfnames) > 1) {
        x0 = oro.nifti::readNIfTI(imgfnames[1], reorient = FALSE)
    }
    else {
        x00 = oro.nifti::readNIfTI(imgfnames, reorient = FALSE)
        if (length(dim(x00)) == 4) 
            x0 = x00[, , , 1]
    }
    if (schange) 
        x0 = sizechange(x0, ...)
    imagedim = dim(x0)
    if (is.null(mask)) {
        brainpos = array(TRUE, imagedim)
    }
    else {
        if (schange) 
            mask = sizechange(mask, ...)
        brainpos = (mask != 0)
    }
    n = ifelse(length(imgfnames) > 1, length(imgfnames), dim(x00)[4])
    S = roi = NULL
    for (i in 1:n) {
        if (length(imgfnames) > 1) {
            x = oro.nifti::readNIfTI(imgfnames[i], reorient = FALSE)
        }
        else {
            x = x00[, , , i]
        }
        if (schange) 
            x = sizechange(x, ...)
        if (ROI) {
            roi0 = atlastable(atlas, x, atlasdataset = atlasdataset, ROIids = ROIids)$table
            rownames(roi0) = abbreviate(roi0$ROIname, 10)
            roi = rbind(roi, t(roi0[, "sumvalue", drop = FALSE]))
        }
        x = c(x)
        x = ifelse(is.na(x), 0, x)
        x = x[brainpos]
        S = rbind(S, x)
    }
    if (zeromask) {
        zeropos = apply(S, 2, function(x) all(x == 0))
        S = S[, !zeropos]
        brainpos[brainpos] = !zeropos
    }
    if (ROI) {
        rownames(roi) = NULL
    }
    else {
        roi = NULL
    }
    list(S = S, brainpos = brainpos, imagedim = imagedim, roi = roi)
}

#' Multi Coat Function
#'
#' This is a function for plotting an image. The analysis result can be overcoated on the template.
#' 
#' \code{multicoat} requires a image array.
#'
#' @name multicoat
#' @aliases multicoat
#' @rdname multicoat
#' @docType methods
#' @export
#'
#' @param imgs A list of three-dimensional image arrays to be displayed.
#' @param y image(s) to be overcoated. In the case of a single image, the same y images are overlaid on each x. In the case of a list of images equal to the number of x, each image is overlaid on each x.
#' @param row4imp the number of rows per a image
#' @param col4imp the number of columns per a image
#' @param trm the index to trim the top and bottom of the slice
#' @param ... further arguments passed to or from other methods.
#' 
#' @examples
#' 
#' data(exbrain)
#' 
multicoat <- function(imgs, y = NULL, row4imp = 6, col4imp = 1, trm = NULL, ...) {
    nimg = length(imgs)
    imagnames = names(imgs)
    if (is.null(imagnames)) 
        imagnames = paste("Image", 1:length(imgs))
    x = imgs[[1]]
    pseq4imps = lapply(imgs, function(img) {
        nslice = dim(img)[3]
        if (is.null(trm)) 
            trm = ifelse(nslice > 100, 15, ifelse(nslice > 50, 5, 1))
        seq(trm, nslice - trm + 1, length = row4imp * col4imp)
    })
    laymat = do.call(cbind, lapply(1:nimg, function(com) rbind(rep((com - 1) * (row4imp * col4imp + 1) + 1, col4imp), matrix((com - 
        1) * (row4imp * col4imp + 1) + (2:(row4imp * col4imp + 1)), ncol = col4imp, byrow = FALSE))))
    par(oma = c(0, 0, 2, 0), mar = c(1, 1, 1, 1), bg = "black")
    layout(laymat, heights = c(1.5, rep(3, row4imp)))
    for (c1 in 1:nimg) {
        plot.new()
        text(0.5, 0.5, paste(imagnames[c1]), cex = 1, font = 1, col = "white")
        for (p1 in rev(pseq4imps[[c1]])) {
            coat(imgs[[c1]], pseq = p1, color.bar = FALSE, paron = F, ...)
        }
    }
}

#' Multi components plot 
#'
#' This is a function that plots the vectorized image returned to its original dimensions by the multirec function.
#' 
#' \code{multicompplot} requires the output result of \code{msma} function.
#'
#' @name multicompplot
#' @aliases multicompplot
#' @rdname multicompplot
#' @docType methods
#' @export
#'
#' @param x A three-dimensional template image.
#' @param object an object of class "\code{multirec}." Usually, a result of a call to \code{multirec}
#' @param comps a component sequence to be plotted.
#' @param row4comp the number of rows per a component
#' @param col4comp the number of columns per a component
#' @param pseq4comp the number of images per a component
#' @param ... further arguments passed to or from other methods.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("msma", quietly = TRUE)) {
#' data(baseimg)
#' data(diffimg)
#' data(mask)
#' data(template)
#' img1 = simbrain(baseimg = baseimg, diffimg = diffimg, mask=mask)
#' B1 = rbfunc(imagedim=img1$imagedim, seppix=2, hispec=FALSE, mask=img1$brainpos)
#' SB1 = basisprod(img1$S, B1)
#' fit111 = msma::msma(SB1, comp=2)
#' ws = multirec(fit111, imagedim=img1$imagedim, B=B1, mask=img1$brainpos)
#' multicompplot(ws, template)
#' }
#' }
multicompplot <- function(object, x, comps = NULL, row4comp = 6, col4comp = 1, pseq4comp = NULL, ...) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    if (is.null(pseq4comp)) 
        pseq4comp = seq(5, dim(x)[3] - 5, length = row4comp * col4comp)
    if (is.null(comps)) 
        comps = object$comps
    ncomp = length(comps)
    laymat = do.call(cbind, lapply(1:ncomp, function(com) rbind(rep((com - 1) * (row4comp * col4comp + 1) + 1, col4comp), 
        matrix((com - 1) * (row4comp * col4comp + 1) + (2:(row4comp * col4comp + 1)), ncol = col4comp, byrow = FALSE))))
    par(oma = c(0, 0, 2, 0), mar = c(1, 1, 1, 1), bg = "black")
    layout(laymat, heights = c(1.5, rep(3, row4comp)))
    for (c1 in 1:ncomp) {
        plot.new()
        text(0.5, 0.5, paste("Comp", comps[c1]), cex = 1, font = 1, col = "white")
        for (p1 in rev(pseq4comp)) {
            coat(x, object$outstat[[comps[c1]]], pseq = p1, color.bar = FALSE, paron = F, ...)
        }
    }
}

#' Multi components reconstruction 
#'
#' This is a function that returns the weight vector of multiple components obtained by the \code{msma} function applied after dimension reduction by the radial basis function to the same dimension as the original image.
#' 
#' \code{multirec} requires the output result of \code{msma} function.
#'
#' @name multirec
#' @aliases multirec
#' @rdname multirec
#' @docType methods
#' @export
#'
#' @param object an object of class \code{msma}. Usually, a result of a call to \code{msma}
#' @param imagedim a vector for original dimension.
#' @param B a list or a matrix.
#' @param mask a list or a matrix.
#' @param midx a block number to be plotted.
#' @param comps a component sequence to be plotted.
#' @param XY a character, indicating "X" or "Y". "XY" for the scatter plots using X and Y scores from \code{msma}.
#' @param signflip a logical if the sign in the block is flipped to pose the super as possitive.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("msma", quietly = TRUE)) {
#' data(baseimg)
#' data(diffimg)
#' data(mask)
#' img1 = simbrain(baseimg = baseimg, diffimg = diffimg, mask=mask)
#' B1 = rbfunc(imagedim=img1$imagedim, seppix=2, hispec=FALSE, mask=img1$brainpos)
#' SB1 = basisprod(img1$S, B1)
#' fit111 = msma::msma(SB1, comp=2)
#' ws = multirec(fit111, imagedim=img1$imagedim, B=B1, mask=img1$brainpos)
#' }
#' }
multirec <- function(object, imagedim, B = NULL, mask = NULL, midx = 1, comps = NULL, XY = c("X", "Y", "XY")[1], signflip = FALSE) {
    if (is.null(comps)) 
        comps = 1:(object$comp[1])
    outstat3 = lapply(comps, function(vidx) {
        Q = object[[paste0("wb", XY)]][[midx]][, vidx]
        outstat1 = rec(Q, imagedim, B = B, mask = mask)
        outstat2 = outstat1 * ifelse(signflip, -1, 1)
    })
    list(outstat = outstat3, comps = comps)
}

#' Prediction Model Function
#'
#' This is the function that creates and evaluates the predictive model.
#' 
#' \code{ptest} requires the output result of \code{msma} function.
#'
#' @name ptest
#' @aliases ptest
#' @rdname ptest
#' @docType methods
#' @export
#'
#' @param object a matrix indicating the explanatory variable(s), or an object of class \code{msma}, which is a result of a call to \code{msma} .
#' @param Z a vector, response variable(s) for the construction of the prediction model. The length of Z is the number of subjects for the training. 
#' @param newdata a matrix for the prediction.
#' @param testZ a vector, response variable(s) for the prediction evaluation. The length of testZ is the number of subjects for the validation. 
#' @param regmethod a character for the name of the prediction model. This corresponds to the \code{method} argument of the \code{train} function in the \code{caret} package.
#' @param methods1 a character for the name of the evaluation method.
#' @param metric a character for the name of summary metric to select the optimal model. 
#' @param number1 a number of folds or number of resampling iterations
#' @param repeats1 a number of repeats for the repeated cross-validation
#' @param params a data frame with possible tuning values. 
#' 
#' @return \item{object}{{an object of class "\code{msma}", usually, a result of a call to \code{msma}}}
#' @return \item{trainout}{a predictive model output from the train function in the caret package with scores computed by the msma function as predictors}
#' @return \item{scorecvroc}{the training evaluation measure and values of the tuning parameters}
#' @return \item{evalmeasure}{evaluation measures and information criterion for the msma model}
#' @return \item{traincnfmat}{a confusion matrix in training data}
#' @return \item{predcnfmat}{a confusion matrix in test data}
#' 
#' @examples
#' 
#' \donttest{
#' if (requireNamespace("msma", quietly = TRUE)) {
#' data(baseimg)
#' data(diffimg)
#' data(mask)
#' data(template)
#' img1 = simbrain(baseimg = baseimg, diffimg = diffimg, mask=mask)
#' B1 = rbfunc(imagedim=img1$imagedim, seppix=2, hispec=FALSE, mask=img1$brainpos)
#' SB1 = basisprod(img1$S, B1)
#' fit111 = msma::msma(SB1, comp=2)
#' predmodel = ptest(fit111, Z=img1$Z)
#' }
#' }
#' 
ptest <- function (object, Z = Z, newdata = NULL, testZ = NULL, regmethod = "glm", methods1 = c("boot", "boot632", "cv", 
    "repeatedcv", "LOOCV", "LGOCV")[4], metric = "ROC", number1 = 10, repeats1 = 5, params = NULL) 
{
    if (!requireNamespace("caret", quietly = TRUE)) {
        stop("ptest() requires the optional package 'caret'.")
    }
    glmFuncs <- caret::lrFuncs
    fivestats <- function(...) c(caret::twoClassSummary(...), caret::defaultSummary(...), prSummary2(...))
    glmFuncs$summary <- fivestats
    ctrl <- caret::trainControl(method = methods1, number = number1, repeats = repeats1, classProbs = TRUE, summaryFunction = fivestats)
    if (inherits(object, "msma")) {
        Ss = object$ssX
        if (inherits(Ss, "list")) {
            Ss = do.call(cbind, Ss)
        }
    }
    else {
        Ss = object
    }
    colnames(Ss) = paste("c", c(1:ncol(Ss)), sep = "")
    swdata = data.frame(Z = as.factor(ifelse(Z == 1, "Y", "N")), Ss)
    n = nrow(swdata)
    if (regmethod == "rf") 
        swdata = swdata[, c(1, 1 + which(apply(swdata[, -1], 2, function(x) length(unique(x)) > 2)))]
    set.seed(721)
    if (regmethod %in% c("glmnet", "rpart")) {
        trainout = caret::train(Z ~ ., data = swdata, method = regmethod, metric = metric, trControl = ctrl, tuneGrid = params)
    }
    else if (regmethod == "mxnet") {
        trainout = caret::train(Z ~ ., data = swdata, method = regmethod, metric = metric, trControl = ctrl, tuneGrid = params, 
            verbose = FALSE)
    }
    else if (regmethod == "nnet") {
        trainout = caret::train(Z ~ ., data = swdata, method = regmethod, metric = metric, trControl = ctrl, tuneGrid = params, 
            trace = FALSE, MaxNWts = 5000)
    }
    else {
        trainout = caret::train(Z ~ ., data = swdata, method = regmethod, metric = metric, trControl = ctrl, tuneGrid = params, 
            trace = FALSE)
    }
    cvout = trainout$results[which.max(trainout$results$ROC), ]
    traincnfmat = caret::confusionMatrix(predict(trainout), swdata$Z, positive = "Y", mode = "everything")
    if (inherits(object, "msma")) {
        evalmeasure = c(scorecvauc = cvout$ROC, bic1 = object$bic[[object$comp[1]]][1], bic1 = object$bic[[object$comp[1]]][2])
    }
    else {
        evalmeasure = cvout
    }
    if (!is.null(newdata)) {
        if (inherits(object, "msma")) {
            pred0 = predict(object, newX = newdata)
            pred1 = pred0$sbX[[1]]
        }
        else {
            pred1 = newdata
        }
        colnames(pred1) = paste("c", c(1:ncol(pred1)), sep = "")
        predZ = predict(trainout, newdata = pred1)
        predcnfmat = caret::confusionMatrix(predZ, testZ, positive = "Y", mode = "everything")
    }
    else {
        predcnfmat = NA
    }
    list(object = object, trainout = trainout, scorecvroc = cvout, evalmeasure = evalmeasure, traincnfmat = traincnfmat, 
        predcnfmat = predcnfmat)
}

#' Radial Basis Function
#'
#' This makes a radial basis function. 
#' 
#' \code{rbfunc} requires the dimensions of the original image to be applied and the knot interval. The output is obtained as a matrix, with the number of rows corresponding to the number of voxels in the original image and the number of columns determined by the knot spacing. By setting hispec = TRUE, you can get the output in list format with a smaller memory.
#'
#' @name rbfunc
#' @aliases rbfunc
#' @rdname rbfunc
#' @docType methods
#' @export
#'
#' @param imagedim a vector indicating image three dimension.
#' @param seppix a numeric. distance between knots.
#' @param hispec a logical. TRUE produces a matrix output. FALSE produces a list output to reduce the data memorry. 
#' @param mask a vector. 
#' @param brainpos a logical vector. 
#' 
#' @examples
#' 
#' imagedim1=c(10,10,10)
#' 
#' B1 = rbfunc(imagedim=imagedim1, seppix=4, hispec=TRUE)
#' B2 = rbfunc(imagedim=imagedim1, seppix=4, hispec=FALSE)
#' 
rbfunc <- function(imagedim, seppix, hispec = FALSE, mask = NULL, brainpos = NULL) {
    h = sqrt(sum(rep(seppix, length(imagedim))^2))
    g = function(r) ifelse(r <= h, h^3 + 3 * h^2 * (h - r) + 3 * h * (h - r)^2 - 3 * (h - r)^3, ifelse(r <= 2 * h, (2 * h - 
        r)^3, 0))/(4 * h^2)
    k0 = as.matrix(expand.grid(lapply(imagedim, function(x) seq(1, x, by = seppix))))
    v0 = as.matrix(expand.grid(lapply(imagedim, function(x) 1:x)))
    if (!is.null(mask)) {
        v = v0[mask[v0] != 0, ]
        k = k0[mask[k0] != 0, ]
    }
    else if (!is.null(brainpos)) {
        v = v0[brainpos[v0], ]
        k = k0[brainpos[k0], ]
    }
    else {
        v = v0
        k = k0
    }
    if (hispec) {
        B = apply(k, 1, function(y) apply(v, 1, function(x) g(sqrt(sum((x - y)^2)))))
    }
    else {
        h2 = ceiling(2 * h)
        vc = apply(v, 1, function(x) paste(x, collapse = "_"))
        B = lapply(1:nrow(k), function(ki) {
            k1 = k[ki, ]
            v1 = as.matrix(expand.grid(lapply(1:length(k1), function(x) {
                max(c(1, k1[x] - h2)):min(c(k1[x] + h2, imagedim[x]))
            })))
            if (!is.null(mask)) 
                v1 = v1[mask[v1] != 0, ]
            if (!is.null(brainpos)) 
                v1 = v1[brainpos[v1], ]
            b1 = apply(v1, 1, function(x) g(sqrt(sum((x - k1)^2))))
            nzidx = b1 > 0
            v2 = v1[nzidx, ]
            b2 = b1[nzidx]
            v2c = apply(v2, 1, function(x) paste(x, collapse = "_"))
            nzidx2 = which(vc %in% v2c)
            B0 = cbind(nzidx2, b2)
            B0
        })
    }
    B
}

#' Reconstruction
#'
#' This is a function that restores the vectorized image to its original dimensions, reduced in dimension by the radial basis function.
#' 
#' \code{rec} requires a vector to be converted to a array.
#'
#' @name rec
#' @aliases rec
#' @rdname rec
#' @docType methods
#' @export
#'
#' @param Q a vector for reduced data. 
#' @param imagedim a vector for original dimension.
#' @param B a list or a matrix indicating the basis function used in the dimension reduction.
#' @param mask a list or a matrix indicating the mask image used in the dimension reduction.
#' 
#' @examples
#' 
#' imagedim1=c(10,10,10)
#' recvec = rec(rnorm(prod(imagedim1)), imagedim1)
#' 
rec <- function(Q, imagedim, B = NULL, mask = NULL) {
    Q = c(Q)
    if (!is.null(B)) {
        QB0 = do.call(rbind, lapply(1:length(B), function(b) {
            tmp = B[[b]]
            tmp[, 2] = tmp[, 2] * Q[b]
            cbind(b, tmp)
        }))
        QB = aggregate(QB0[, 3], by = list(QB0[, 2]), FUN = sum)
        outstat0 = QB[, 2]
    }
    else {
        outstat0 = Q
    }
    outstat1 = array(0, imagedim)
    if (!is.null(mask)) {
        outstat1[mask] = outstat0
    }
    else {
        outstat1 = outstat0
    }
    outstat1
}

#' Generate simulation data Function
#'
#' This is a function for simulation data based on the real base brain image data and difference in brain between healty and disease groups.
#' 
#' \code{simbrain} requires a base brain image data and mean difference image data.
#'
#' @name simbrain
#' @aliases simbrain
#' @rdname simbrain
#' @docType methods
#' @export
#'
#' @param baseimg an array for the basis image.
#' @param diffimg an array for the difference image.
#' @param sdevimg an array for the standard deviation image.
#' @param mask an array for the mask image.
#' @param n0 a numeric, which is a sample size per group.
#' @param c1 a numeric, the strength of the difference
#' @param sd1 a numeric, standard deviation for the individual variation.
#' @param rho a numeric, correlation coefficient in the noize
#' @param zeromask a logical, whether mask the position with zero values for all subjects.
#' @param reduce a vector.
#' @param output a vector.
#' @param seed a numeric for seed for random variables.
#' 
#' @return \item{S}{data matrix}
#' @return \item{Z}{binary group variable}
#' @return \item{brainpos}{binary brain position.}
#' @return \item{imagedim}{three dimensional vector for image dimension}
#' 
#' @examples
#' data(baseimg)
#' data(diffimg)
#' sim1 = simbrain(baseimg = baseimg, diffimg = diffimg)
#' 
simbrain <- function (baseimg, diffimg, sdevimg = NULL, mask = NULL, n0 = 10, c1 = 0.5, sd1 = 0.01, rho = NULL, zeromask = FALSE, 
    reduce = c("no", "rd1", "rd2")[1], output = c("rdata", "nifti")[1], seed = 1) 
{
    if (!is.null(rho) && !requireNamespace("mnormt", quietly = TRUE)) {
        stop("simbrain(rho != NULL) requires the optional package 'mnormt'.")
    }
    imagedim = dim(baseimg)
    if (is.null(mask)) {
        brainpos = array(TRUE, imagedim)
    }
    else {
        brainpos = (mask != 0)
    }
    n = (2 * n0)
    Z = rep(c(0, 1), each = n0)
    if (is.null(sdevimg)) 
        sdevimg = array(1, dim = imagedim)
    if (!is.null(rho)) {
        cordmat = expand.grid(lapply(imagedim, seq))
        sigma1 = rho^as.matrix(dist(cordmat))
        sigma1 = diag(sd1 * c(sdevimg)) %*% sigma1 %*% diag(sd1 * c(sdevimg))
    }
    set.seed(seed)
    S = do.call(rbind, lapply(1:n, function(i) {
        if (is.null(rho)) {
            noise = array(rnorm(prod(imagedim), sd = sd1 * c(sdevimg)), dim = imagedim)
        }
        else {
            noise = array(c(mnormt::rmnorm(1, sigma = sigma1)), dim = imagedim)
        }
        data1 = baseimg + c1 * Z[i] * diffimg
        data2 = data1
        data2[data2 != 0] = data2[data2 != 0] + noise[data2 != 0]
        data2[data2 < 0] = 0
        data2 = c(data2)
        data2[brainpos]
    }))
    if (zeromask) {
        zeropos = apply(S, 2, function(x) all(x == 0))
        S = S[, !zeropos]
        brainpos[brainpos] = !zeropos
    }
    list(S = S, Z = Z, brainpos = brainpos, imagedim = imagedim)
}

.mand_resize_axis_linear <- function(x, new_n, axis) {
    old_n <- dim(x)[axis]
    if (old_n == new_n) 
        return(x)
    old_pos <- seq_len(old_n)
    new_pos <- seq(1, old_n, length.out = new_n)
    perm <- c(axis, setdiff(seq_along(dim(x)), axis))
    inv_perm <- order(perm)
    xp <- aperm(x, perm)
    d <- dim(xp)
    mat <- matrix(xp, nrow = d[1L])
    out <- vapply(seq_len(ncol(mat)), function(j) stats::approx(old_pos, mat[, j], xout = new_pos, method = "linear", rule = 2)$y, 
        numeric(new_n))
    out <- array(out, dim = c(new_n, d[-1L]))
    aperm(out, inv_perm)
}

.mand_resize_axis_nearest <- function(x, new_n, axis) {
    old_n <- dim(x)[axis]
    if (old_n == new_n) 
        return(x)
    index <- round(seq(1, old_n, length.out = new_n))
    index <- pmax(1L, pmin(old_n, index))
    subs <- rep(list(quote(expr = )), length(dim(x)))
    subs[[axis]] <- index
    do.call(`[`, c(list(x), subs, list(drop = FALSE)))
}

.mand_resize_array3d <- function(x, refsize, interpolation = c("linear", "nearest")) {
    interpolation <- match.arg(interpolation)
    if (length(dim(x)) != 3L) 
        stop("x must be a 3-D array.")
    refsize <- as.integer(refsize)
    if (length(refsize) != 3L || anyNA(refsize) || any(refsize < 1L)) {
        stop("refsize must contain three positive integers.")
    }
    storage.mode(x) <- "double"
    out <- x
    for (axis in 1:3) {
        out <- if (interpolation == "linear") {
            .mand_resize_axis_linear(out, refsize[axis], axis)
        }
        else {
            .mand_resize_axis_nearest(out, refsize[axis], axis)
        }
    }
    dim(out) <- refsize
    out
}

.mand_resize_axis_linear <- function(x, new_n, axis) {
    old_n <- dim(x)[axis]
    if (old_n == new_n) 
        return(x)
    old_pos <- seq_len(old_n)
    new_pos <- seq(1, old_n, length.out = new_n)
    perm <- c(axis, setdiff(seq_along(dim(x)), axis))
    inv_perm <- order(perm)
    xp <- aperm(x, perm)
    d <- dim(xp)
    mat <- matrix(xp, nrow = d[1L])
    out <- vapply(seq_len(ncol(mat)), function(j) stats::approx(old_pos, mat[, j], xout = new_pos, method = "linear", rule = 2)$y, 
        numeric(new_n))
    out <- array(out, dim = c(new_n, d[-1L]))
    aperm(out, inv_perm)
}

.mand_resize_axis_nearest <- function(x, new_n, axis) {
    old_n <- dim(x)[axis]
    if (old_n == new_n) 
        return(x)
    index <- round(seq(1, old_n, length.out = new_n))
    index <- pmax(1L, pmin(old_n, index))
    subs <- rep(list(quote(expr = )), length(dim(x)))
    subs[[axis]] <- index
    do.call(`[`, c(list(x), subs, list(drop = FALSE)))
}

.mand_resize_array3d <- function(x, refsize, interpolation = c("linear", "nearest")) {
    interpolation <- match.arg(interpolation)
    if (length(dim(x)) != 3L) 
        stop("x must be a 3-D array.")
    refsize <- as.integer(refsize)
    if (length(refsize) != 3L || anyNA(refsize) || any(refsize < 1L)) {
        stop("refsize must contain three positive integers.")
    }
    storage.mode(x) <- "double"
    out <- x
    for (axis in 1:3) {
        out <- if (interpolation == "linear") {
            .mand_resize_axis_linear(out, refsize[axis], axis)
        }
        else {
            .mand_resize_axis_nearest(out, refsize[axis], axis)
        }
    }
    dim(out) <- refsize
    out
}

#' Size change Function
#'
#' TThis is a function that changes the resolution of the image.
#' 
#' \code{sizechange} requires the array data.
#'
#' @name sizechange
#' @aliases sizechange
#' @rdname sizechange
#' @docType methods
#' @export
#'
#' @param img1 a array or nifti class, which is a image data to be changed the size.
#' @param simscale a numeric.
#' @param refsize a vector with length 3, which is a size to be changed.
#' @param ... further arguments passed to or from other methods.
#' 
#' @examples
#' 
#' data(exbrain)
#' exbrain2 = sizechange(exbrain, simscale=1/2)
#' 
sizechange <- function(img1, simscale = NULL, refsize = NULL, ...) {
    dots <- list(...)
    interpolation <- if (!is.null(dots$interpolation)) {
        match.arg(dots$interpolation, c("linear", "nearest"))
    }
    else {
        "linear"
    }
    orgsize <- dim(img1)
    if (length(orgsize) < 3L) 
        stop("img1 must be at least three-dimensional.")
    orgsize <- as.integer(orgsize[1:3])
    if (!is.null(simscale) && is.null(refsize)) {
        simscale <- as.numeric(simscale)
        if (!length(simscale) %in% c(1L, 3L) || any(!is.finite(simscale)) || any(simscale <= 0)) 
            stop("simscale must be positive and length 1 or 3.")
        if (length(simscale) == 1L) 
            simscale <- rep(simscale, 3L)
        refsize <- ceiling(orgsize * simscale)
    }
    if (is.null(refsize)) 
        stop("Specify simscale or refsize.")
    refsize <- as.integer(refsize)
    if (length(refsize) != 3L || anyNA(refsize) || any(refsize < 1L)) {
        stop("refsize must contain three positive integers.")
    }
    is_nifti <- methods::is(img1, "nifti")
    source_array <- if (is_nifti) 
        img1@.Data
    else as.array(img1)
    if (length(dim(source_array)) == 4L && dim(source_array)[4L] == 1L) {
        source_array <- source_array[, , , 1L, drop = TRUE]
    }
    if (length(dim(source_array)) != 3L) 
        stop("img1 must contain one 3-D image.")
    resized <- .mand_resize_array3d(source_array, refsize, interpolation)
    if (!is_nifti) 
        return(resized)
    result <- img1
    result@.Data <- resized
    result@dim_[1L] <- 3L
    result@dim_[2:4] <- refsize
    result@dim_[5:8] <- 1L
    result@pixdim[2:4] <- result@pixdim[2:4] * orgsize/refsize
    result
}

