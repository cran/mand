# Internal compatibility helpers restored from public mand 2.0.
# These functions are intentionally not exported.

tim.colors2 <- function (n = 64) 
{
    orig <- c("#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF", "#0000DF", "#0000EF", "#0000FF", "#0010FF", "#0020FF", "#0030FF", "#0040FF", "#0050FF", "#0060FF", "#0070FF", "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF", "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF", "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F", "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50", "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00", "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", "#FFAF00", 
        "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", "#FF6000", "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000", "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000", "#AF0000", "#9F0000", "#8F0000", "#800000")[1:32]
    orig = rep(orig, each = 2)
    if (n == 64) {
        return(orig)
    }
    rgb.tim <- t(col2rgb(orig))
    temp <- matrix(NA, ncol = 3, nrow = n)
    x <- seq(0, 1, length.out = 64)
    xg <- seq(0, 1, length.out = n)
    for (k in 1:3) {
        hold <- splines::interpSpline(x, rgb.tim[, k])
        hold <- predict(hold, xg)$y
        hold[hold < 0] <- 0
        hold[hold > 255] <- 255
        temp[, k] <- round(hold)
    }
    rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
}

colbar <- function (zlim, col.y = .mand_hotmetal(), nticks = 4, horizontal = FALSE, axiscol = "white", bg = "black", ...) 
{
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    min = zlim[1]
    max = zlim[2]
    ticks = round(seq(min, max, len = nticks), 2)
    scale1 = (length(col.y) - 1)/(max - min)
    rectcords = sapply(1:(length(col.y) - 1), function(i) {
        y = (i - 1)/scale1 + min
        c(y, 0, y + 1/scale1, 10)
    })
    if (horizontal) {
        x1 = c(min, max)
        y1 = c(0, 10)
        axispos = 1
        rdix = 1:4
    }
    else {
        x1 = c(0, 10)
        y1 = c(min, max)
        axispos = 2
        rdix = c(2, 1, 4, 3)
    }
    par(omi = c(dev.size()[2] * 0.1, 0, dev.size()[2] * 0.3, 0), bg = bg, new = TRUE)
    plot(x1, y1, type = "n", bty = "n", xaxt = "n", xlab = "", yaxt = "n", ylab = "", ylim = c(0, 20), ...)
    axis(axispos, ticks, las = 1, col = axiscol, col.axis = axiscol)
    for (i in 1:ncol(rectcords)) rect(rectcords[rdix[1], i], rectcords[rdix[2], i], rectcords[rdix[3], i], rectcords[rdix[4], i], col = col.y[i], border = NA)
}

prSummary2 <- function (data, lev = NULL, model = NULL) 
{
    if (length(levels(data$obs)) > 2) 
        stop(paste("Your outcome has", length(levels(data$obs)), "levels. `prSummary`` function isn't appropriate.", call. = FALSE))
    if (!all(levels(data[, "pred"]) == levels(data[, "obs"]))) 
        stop("Levels of observed and predicted data do not match.", call. = FALSE)
    if (!lev[1] %in% colnames(data)) 
        stop(paste("Class probabilities are needed to score models using the", "area under the PR curve. Set `classProbs = TRUE`", "in the trainControl() function."), call. = FALSE)
    c(Precision = caret::precision(data = data$pred, reference = data$obs, relevant = lev[1]), Recall = caret::recall(data = data$pred, reference = data$obs, relevant = lev[1]), F = caret::F_meas(data = data$pred, reference = data$obs, relevant = lev[1]))
}

