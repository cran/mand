#' Multivariate Analysis for Neuroimaging Data Package
#'
#' A Package for implementation of multivariate data analysis for neuroimaging data.
#'
#' @name mand-package
#' @aliases mand-package
#' @rdname mand-package
#' @docType package
#' @keywords documentation
#'
#' @author Atsushi Kawaguchi. \email{kawa_a24@@yahoo.co.jp}
#' @references 
#' Kawaguchi, A. (2021). Multivariate Analysis for Neuroimaging Data. CRC Press.
#' @import msma oro.nifti oro.dicom imager caret
#' @importFrom graphics abline matplot plot par axis layout rect legend title mtext plot.new text
#' @importFrom grDevices gray rainbow col2rgb rgb dev.size
#' @importFrom stats cor predict quantile rbinom rnorm runif aggregate dist integrate na.omit splinefun
#' @importFrom utils data
NULL

#' Example Brain Data
#'
#' The data are from a MRI gray matter brain data for one subject.
#' 
#' @docType data
#' @keywords datasets
#' @name exbrain
#' @usage data(exbrain)
#' @format A array
NULL

#' Brain Template
#'
#' The data is the brain tempalte. This is an average brain image, and is mainly used for overlaying analysis results.
#' 
#' @docType data
#' @keywords datasets
#' @name template
#' @usage data(template)
#' @format A array
NULL

#' Base Brain Data
#'
#' The data is the base brain data. This is an average image of a healthy person, and is used when generating artificial data.
#' 
#' @docType data
#' @keywords datasets
#' @name baseimg
#' @usage data(baseimg)
#' @format A array
NULL

#' Difference Brain Data
#'
#' The data is the difference brain data. This represents the difference between the average images of healthy subjects and patients with Alzheimer's disease, and is used when generating artificial data.
#' 
#' @docType data
#' @keywords datasets
#' @name diffimg
#' @usage data(diffimg)
#' @format A array
NULL

#' Standard Deviation Brain Data
#'
#' The data is the standard deviation brain data. This represents the common standard deviation between the average images of healthy subjects and patients with Alzheimer's disease, and is used when generating artificial data.
#' 
#' @docType data
#' @keywords datasets
#' @name sdevimg
#' @usage data(sdevimg)
#' @format A array
NULL

#' Brain Mask
#'
#' The data is the brain mask. This is used to exclude extra-brain regions from the analysis.
#' 
#' @docType data
#' @keywords datasets
#' @name mask
#' @usage data(mask)
#' @format A array
NULL

#' Atlas data set
#'
#' The data is the atlas data. Various atlases are stored. Each matrix has "ROIid" and "ROIname" as column names.
#' 
#' @docType data
#' @keywords datasets
#' @name atlasdatasets
#' @usage data(atlasdatasets)
#' @format A list of matrix
NULL

#' Atlas set
#'
#' The data is the atlas image data. An image whose element is "ROIid" is stored for each atlas.
#' 
#' @docType data
#' @keywords datasets
#' @name atlas
#' @usage data(atlas)
#' @format A list of array
NULL


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
rbfunc = function(imagedim, seppix, hispec=FALSE, mask=NULL, brainpos=NULL)
{
h = sqrt(sum(rep(seppix, length(imagedim))^2))
g = function(r) ifelse(r <= h, h^3 + 3 * h^2 * (h - r) + 3 * h * (h - r)^2 - 3 * (h - r)^3, ifelse(r <= 2*h, (2*h - r)^3, 0)) / (4 * h^2)
k0 = as.matrix(expand.grid(lapply(imagedim, function(x) seq(1, x, by = seppix))))
v0 = as.matrix(expand.grid(lapply(imagedim, function(x) 1:x)))
if(!is.null(mask)){
v = v0[mask[v0] != 0, ]
k = k0[mask[k0] != 0, ]
}else if(!is.null(brainpos)){
v = v0[brainpos[v0], ]
k = k0[brainpos[k0], ]
}else{
v = v0
k = k0
}
if(hispec){
B = apply(k, 1, function(y) apply(v, 1, function(x) g(sqrt(sum((x - y)^2)))))
}else{
h2 = ceiling(2*h)
vc = apply(v, 1, function(x) paste(x, collapse="_"))
# neighbourhood search
B = lapply(1:nrow(k), function(ki){
k1 = k[ki, ]
v1 = as.matrix(expand.grid(lapply(1:length(k1), function(x){ max(c(1, k1[x]-h2)) : min(c(k1[x]+h2,imagedim[x]))} )))
if(!is.null(mask)) v1 = v1[mask[v1] != 0, ]
if(!is.null(brainpos)) v1 = v1[brainpos[v1], ]
b1 = apply(v1, 1, function(x) g(sqrt(sum((x - k1)^2))))
nzidx = b1 > 0
v2 = v1[nzidx,]
b2 = b1[nzidx]
v2c = apply(v2, 1, function(x) paste(x, collapse="_"))
nzidx2 = which(vc %in% v2c)
B0 = cbind(nzidx2, b2)
B0
})

}

B
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
basisprod = function(A, B){
if(inherits(B, "list")){
do.call(cbind, lapply(B, function(b2){A[,b2[,1]] %*% b2[,2]}))
}else if(inherits(A, "list")){
do.call(cbind, lapply(A, function(b2){B[,b2[,1]] %*% b2[,2]}))
}else{
stop("Either A or B should be output of rbfunc")
}

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
rec = function(Q, imagedim, B=NULL, mask=NULL){
Q = c(Q)
if(!is.null(B)){
QB0 = do.call(rbind, lapply(1:length(B), function(b){
tmp = B[[b]]
tmp[,2] = tmp[,2]*Q[b]
cbind(b, tmp)
}))
QB = aggregate(QB0[,3], by=list(QB0[,2]), FUN=sum)
outstat0 = QB[,2]
}else{
outstat0 = Q
}
outstat1 = array(0, imagedim)
if(!is.null(mask)){outstat1[mask] = outstat0}else{outstat1 = outstat0}
outstat1
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
#' @param x image1. Base image.
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
coat = function(x, y=NULL, pseq=NULL, xyz=NULL, col.x=gray(0:64/64), col.y=NULL, breaks.y=NULL, zlim.x=NULL, zlim.y=NULL, rownum=5, colnum=NULL, plane=c("axial", "coronal", "sagittal", "all")[1], xlab="", ylab="", axes=FALSE, oma=rep(0,4), mar=rep(0,4), bg="black", paron=TRUE, cross.hair=FALSE, chxy=NULL, color.bar=TRUE, regionplot=FALSE, atlasdataset=NULL, regionname=c("atlas", "stat")[1], regionlegend=FALSE, atlasname = "", ROIids=1:9, ...) {

if(paron){
oldpar <- par(no.readonly = TRUE)
on.exit(par(oldpar))
}else{
regionlegend = FALSE
}

#####
if(is.null(zlim.x)) zlim.x = range(x, na.rm=TRUE)
if(!is.null(y)){ if(is.null(zlim.y)) zlim.y = range(y, na.rm=TRUE)}

## set dimensions
dim.x = dim(x); X <- dim.x[1]; Y <- dim.x[2]; Z <- dim.x[3]

#####
if(plane != "all"){
pseqmax = switch(plane, axial=Z, coronal=Y, sagittal=X)
if(is.null(pseq)) pseq = rev(seq(1, pseqmax, length=30))
}

#
rownum = min(c(rownum, length(pseq)))
if(is.null(rownum) & !is.null(colnum)) rownum = ceiling(length(pseq)/colnum)
if(is.null(colnum)) colnum = ceiling(length(pseq)/rownum)

plotnum = rownum*colnum

breaks.x <- seq(min(x, zlim.x, na.rm=TRUE), max(x, zlim.x, na.rm=TRUE), length=length(col.x)+1)

if(!is.null(y)){
if(is.null(col.y)){ col.y = hotmetal(); prescol.y=FALSE}else{prescol.y=TRUE}
if(!is.na(col.y[1])) col.y = c(NA, col.y)
col.y2 = col.y
if(is.null(breaks.y)) breaks.y <- seq(zlim.y[1], zlim.y[2], length=length(col.y)+1)
if(any(sign( breaks.y )<0)){ 
negidxmax = max(which(sign( breaks.y )<0))
if(negidxmax < length(breaks.y)){
if(!prescol.y) {
breaks.y = c(breaks.y[1:negidxmax], 0, breaks.y[(negidxmax+1):length(breaks.y)])
col.y = c(rev(tim.colors2(negidxmax-1)), NA, hotmetal(length(breaks.y)-(negidxmax+1)))
col.y2 = c(rev(tim.colors2(negidxmax)), hotmetal(length(breaks.y)-(negidxmax+1)))
#col.y[negidxmax] = NA 
}}}}

if(is.null(y)) color.bar = FALSE

#####
if(regionplot){
if(is.null(atlasdataset)) stop("atlasdataset should be required")
atlasdataset1 = atlasdataset[as.character(atlasdataset$ROIid) %in% as.character(ROIids),]
if(regionname=="stat"){
breaks.y = c(0, sort(ROIids)-0.5, max(ROIids)+1)
#breaks.y = c(0, (1:length(ROIids))-0.5, length(ROIids)+1);
}else{
y[!(y %in% ROIids)] = 0
breaks.y = c(0, sort(ROIids)-0.5, max(ROIids)+1)
}
#colors2 = c("red", "green3", "cyan", "magenta", "yellow", "brown", "antiquewhite", "aquamarine", "darkgoldenrod")
colors2 = c("red", "green3", "orange", "magenta", "blue", "deeppink", "darkgreen", "yellow", "brown", "cyan")
#col.y = c(NA, colors2[1:(length(breaks.y)-2)])
col.y = c(NA, colors2[1:(length(ROIids))])
color.bar=FALSE
}

########## Plot ##########

##### axial #####
if(plane=="axial"){
if(paron){
#par(mfrow=c(rownum, colnum), oma=oma, mar=mar, bg=bg)
par(oma=oma, mar=mar, bg=bg)
if(color.bar | regionlegend){
layout(rbind(matrix(1:plotnum, ncol=colnum, byrow=TRUE), rep(plotnum+1, colnum)))
}else{
layout(matrix(1:plotnum, ncol=colnum, byrow=TRUE))
}
}
for (z2 in pseq) {
graphics::image(1:X, 1:Y, x[,,z2], col=col.x, breaks=breaks.x, zlim=zlim.x, axes=axes, xlab=xlab, ylab=ylab, ...)
if(!is.null(y)){
graphics::image(1:X, 1:Y, y[,,z2], col=col.y, breaks=breaks.y, zlim=zlim.y, add=TRUE)
}}
if(color.bar){ colbar(zlim=zlim.y, col.y=col.y2, nticks=4, horizontal=TRUE)
}else if(regionlegend){
#par(oma= c(1, 0, 4, 0), mar= c(4, 2, 4, 2), bg=bg, new=TRUE)
par(oma= c(1, 0, 4, 0), bg=bg, new=TRUE)
plot(1, type = "n", axes = FALSE, ann = FALSE); mtext(atlasname, col="white", cex = 2, line=-3)
graphics::legend("left", legend = abbreviate(atlasdataset1$ROIname, 10), col = col.y[-1], ncol = 5, pch = 15, text.font = 1, text.col = col.y[-1])
}
}else if(plane=="sagittal"){
##### sagittal #####
if(paron) par(mfrow=c(rownum, colnum), oma=oma, mar=mar, bg=bg)
for (x2 in pseq) {
graphics::image(1:Y, 1:Z, x[x2,,], col=col.x, breaks=breaks.x, zlim=zlim.x, axes=axes, xlab=xlab, ylab=ylab)
if(!is.null(y)){ 
graphics::image(1:Y, 1:Z, y[x2,,], col=col.y, breaks=breaks.y, zlim=zlim.y, add=TRUE)
}}
}else if(plane=="coronal"){
##### coronal #####
if(paron) par(mfrow=c(rownum, colnum), oma=oma, mar=mar, bg=bg)
for (y2 in pseq) {
graphics::image(1:X, 1:Z, x[,y2,], col=col.x, breaks=breaks.x, zlim=zlim.x, axes=axes, xlab=xlab, ylab=ylab)
if(!is.null(y)){
graphics::image(1:X, 1:Z, y[,y2,], col=col.y, breaks=breaks.y, zlim=zlim.y, add=TRUE)
}}
}else if(plane=="all"){
##### axial, sagittal, coronal #####
if(is.null(xyz)) xyz = c(round(X/2), round(Y/2), round(Z/2))
if(paron) par(mfrow=c(1, 3), oma=oma, mar=c(0,0,3,0), bg="black")
graphics::image(1:X, 1:Y, x[,,xyz[3]], col=col.x, breaks=breaks.x, zlim=zlim.x, axes=axes, xlab=xlab, ylab=ylab, main="axial", col.main="white")
if(!is.null(y)){
graphics::image(1:X, 1:Y, y[,,xyz[3]], col=col.y, breaks=breaks.y, zlim=zlim.y, add=TRUE)
if(cross.hair){
if(is.null(chxy)) chxy = which(y[,,xyz[3]]==max(y[,,xyz[3]]), arr.ind=TRUE)[1,]
abline(v=chxy[1], col="green")
abline(h=chxy[2], col="green")
}}
graphics::image(1:Y, 1:Z, x[xyz[1],,], col=col.x, breaks=breaks.x, zlim=zlim.x, axes=axes, xlab=xlab, ylab=ylab, main="sagittal", col.main="white")
if(!is.null(y)){
graphics::image(1:Y, 1:Z, y[xyz[1],,], col=col.y, breaks=breaks.y, zlim=zlim.y, add=TRUE)
if(cross.hair){
if(is.null(chxy)) chxy = which(y[xyz[1],,]==max(y[xyz[1],,]), arr.ind=TRUE)[1,]
abline(v=chxy[1], col="green")
abline(h=chxy[2], col="green")
}}
graphics::image(1:X, 1:Z, x[,xyz[2],], col=col.x, breaks=breaks.x, zlim=zlim.x, axes=axes, xlab=xlab, ylab=ylab, main="coronal", col.main="white")
if(!is.null(y)){
graphics::image(1:X, 1:Z, y[,xyz[2],], col=col.y, breaks=breaks.y, zlim=zlim.y, add=TRUE)
if(cross.hair){
if(is.null(chxy)) chxy = which(y[,xyz[2],]==max(y[,xyz[2],]), arr.ind=TRUE)[1,]
abline(v=chxy[1], col="green")
abline(h=chxy[2], col="green")
}}
}


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
#' @param imgs list of images. Base images.
#' @param y image2 to be overcoated.
#' @param row4imp the number of rows per a image
#' @param col4imp the number of columns per a image
#' @param trm the index to trim the top and bottom of the slice
#' @param ... further arguments passed to or from other methods.
#' 
#' @examples
#' 
#' data(exbrain)
#' 
multicoat = function(imgs, y=NULL, row4imp=6, col4imp=1, trm=NULL,...) {

nimg = length(imgs)
imagnames = names(imgs)
if(is.null(imagnames)) imagnames = paste("Image", 1:length(imgs))
x = imgs[[1]]
pseq4imps = lapply(imgs, function(img){
nslice = dim(img)[3]
if(is.null(trm)) trm = ifelse(nslice>100, 15, ifelse(nslice>50, 5, 1))
seq(trm, nslice-trm+1, length=row4imp*col4imp)
})
  
laymat = do.call(cbind, lapply(1:nimg, function(com) rbind(rep((com-1)*(row4imp*col4imp+1)+1, col4imp), 
                                                     matrix((com-1)*(row4imp*col4imp+1) + (2:(row4imp*col4imp+1)), ncol=col4imp, byrow=FALSE)))
)
  
  par(oma=c(0,0,2,0), mar=c(1,1,1,1), bg="black")
  layout(laymat, heights=c(1.5, rep(3, row4imp)))
  for(c1 in 1:nimg){
    plot.new(); text(0.5,0.5,paste(imagnames[c1]), cex=1, font=1, col="white")
    for(p1 in rev(pseq4imps[[c1]])){
      coat(imgs[[c1]], pseq=p1, color.bar=FALSE, paron=F,...)
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
#' @param midx a block number.
#' @param comps a component sequence to be plotted.
#' @param XY a character, indicating "X" or "Y". "XY" for the scatter plots using X and Y scores from \code{msma}.
#' @param signflip a logical if the sign in the block is flipped to pose the super as possitive.
#'
#' @examples
#' \donttest{
#' data(baseimg)
#' data(diffimg)
#' data(mask)
#' img1 = simbrain(baseimg = baseimg, diffimg = diffimg, mask=mask)
#' B1 = rbfunc(imagedim=img1$imagedim, seppix=2, hispec=FALSE, mask=img1$brainpos)
#' SB1 = basisprod(img1$S, B1)
#' fit111 = msma(SB1, comp=2)
#' ws = multirec(fit111, imagedim=img1$imagedim, B=B1, mask=img1$brainpos)
#' }
multirec = function(object, imagedim, B=NULL, mask=NULL, midx=1, comps=NULL, XY=c("X", "Y", "XY")[1], signflip=FALSE){

if(is.null(comps)) comps=1:(object$comp[1])

outstat3 = lapply(comps, function(vidx){
Q = object[[paste0("wb", XY)]][[midx]][,vidx]
outstat1 = rec(Q, imagedim, B=B, mask=mask)
outstat2 = outstat1 * ifelse(signflip,-1,1)
})

list(outstat=outstat3, comps=comps)
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
#' @param object an object of class "\code{multirec}." Usually, a result of a call to \code{multirec}
#' @param x template image
#' @param comps a component sequence to be plotted.
#' @param row4comp the number of rows per a component
#' @param col4comp the number of columns per a component
#' @param pseq4comp the number of images per a component
#' @param ... further arguments passed to or from other methods.
#'
#' @examples
#' \donttest{
#' data(baseimg)
#' data(diffimg)
#' data(mask)
#' data(template)
#' img1 = simbrain(baseimg = baseimg, diffimg = diffimg, mask=mask)
#' B1 = rbfunc(imagedim=img1$imagedim, seppix=2, hispec=FALSE, mask=img1$brainpos)
#' SB1 = basisprod(img1$S, B1)
#' fit111 = msma(SB1, comp=2)
#' ws = multirec(fit111, imagedim=img1$imagedim, B=B1, mask=img1$brainpos)
#' multicompplot(ws, template)
#' }
multicompplot = function(object, x, comps=NULL, row4comp=6, col4comp=1, pseq4comp=NULL,...){

oldpar <- par(no.readonly = TRUE)
on.exit(par(oldpar))

if(is.null(pseq4comp)) pseq4comp=seq(5, dim(x)[3]-5, length=row4comp*col4comp)
if(is.null(comps)) comps=object$comps

ncomp=length(comps)

laymat = do.call(cbind, lapply(1:ncomp, function(com) rbind(rep((com-1)*(row4comp*col4comp+1)+1, col4comp), 
matrix((com-1)*(row4comp*col4comp+1) + (2:(row4comp*col4comp+1)), ncol=col4comp, byrow=FALSE)))
)

par(oma=c(0,0,2,0), mar=c(1,1,1,1), bg="black")
layout(laymat, heights=c(1.5, rep(3, row4comp)))
for(c1 in 1:ncomp){
plot.new(); text(0.5,0.5,paste("Comp",comps[c1]), cex=1, font=1, col="white")
for(p1 in rev(pseq4comp)){
coat(x, object$outstat[[comps[c1]]], pseq=p1, color.bar=FALSE, paron=F,...)
}}

}

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
#' @param x an array for the atlas image.
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
atlastable = function(x, y, atlasdataset=NULL, ROIids=NULL, ...) {

if(is.null(atlasdataset)) stop("atlasdataset should be required")
if(is.null(ROIids)) ROIids = atlasdataset$ROIid
atlasdataset1 = atlasdataset[as.character(atlasdataset$ROIid) %in% as.character(ROIids),]

stats = do.call(rbind, lapply(ROIids, function(ROIid1){
x1 = (x == ROIid1)
tmp = y * x1
sizenum = sum(tmp!=0)
sumvalue = sum(tmp)
sizepct = sizenum / sum(x1)
c(sizenum=sizenum, sizepct=sizepct, sumvalue=sumvalue, summary(tmp))
}))

out=list(table=cbind(atlasdataset1, stats))
class(out) = "atlastable"
out
}


#' @rdname atlastable
#' @method print atlastable
#' @family print
#' @export
print.atlastable = function(x, ...)
{
t1 = x$table
absminmax = apply(abs(t1[,c("Min.", "Max.")]), 1, max)
disnames = c("sizepct", "sumvalue", "Min.", "Mean", "Max.")
idx = which(colnames(t1) %in% disnames)
t1[, idx] = round(t1[, idx], 3)
print(t1[order(-absminmax)[1:min(c(10,nrow(t1)))], c("ROIid", "ROIname", disnames)])
cat("\n")
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
#' @param imgfnames a vector for (nifti) file names to be used.
#' @param mask a vector for brain mask data.
#' @param ROI a logical for roi data set.
#' @param atlas an array for the atlas.
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
imgdatamat = function(imgfnames, mask=NULL, ROI=FALSE, atlas=NULL, atlasdataset=NULL, ROIids=NULL, zeromask=FALSE, schange=FALSE, ...)
{
if(length(imgfnames)>1){ x0 = readNIfTI(imgfnames[1], reorient = FALSE)}else{
x00 = readNIfTI(imgfnames, reorient = FALSE)
if(length(dim(x00)) == 4) x0 = x00[,,,1]
}

if(schange) x0 = sizechange(x0, ...)
imagedim = dim(x0)

if(is.null(mask)){brainpos = array(TRUE, imagedim)}else{
if(schange) mask = sizechange(mask, ...)
 brainpos = (mask != 0)
}

n = ifelse(length(imgfnames)>1, length(imgfnames), dim(x00)[4])
S = roi = NULL
for(i in 1:n){
if(length(imgfnames)>1){ x = readNIfTI(imgfnames[i], reorient = FALSE)}else{
x = x00[,,,i]
}
if(schange) x = sizechange(x, ...)
if(ROI){
roi0 = atlastable(atlas, x, atlasdataset=atlasdataset, ROIids=ROIids)$table
rownames(roi0) = abbreviate(roi0$ROIname, 10)
roi = rbind(roi, t(roi0[,"sumvalue",drop=FALSE]))
}
x = c(x)
x = ifelse(is.na(x), 0, x) 
x = x[brainpos]
S = rbind(S, x)
}

if(zeromask){
zeropos = apply(S, 2, function(x) all(x == 0))
S = S[,!zeropos] 
brainpos[brainpos] = !zeropos 
}

if(ROI){rownames(roi) = NULL}else{roi = NULL}

list(S=S, brainpos=brainpos, imagedim=imagedim, roi=roi)
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
simbrain = function(baseimg, diffimg, sdevimg=NULL, mask=NULL, n0 = 10, c1 = 0.5, sd1 = 0.01, rho=NULL, zeromask=FALSE, reduce = c("no", "rd1", "rd2")[1], output = c("rdata", "nifti")[1], seed=1){

imagedim = dim(baseimg)
if(is.null(mask)){brainpos = array(TRUE, imagedim)}else{ brainpos = (mask != 0)}
n = (2*n0)
Z = rep(c(0,1), each=n0)

if(is.null(sdevimg)) sdevimg = array(1, dim=imagedim)

if(!is.null(rho)) {
cordmat = expand.grid(lapply(imagedim, seq))
sigma1 = rho^as.matrix(dist(cordmat))
sigma1 = diag(sd1*c(sdevimg)) %*% sigma1 %*% diag(sd1*c(sdevimg))
}

set.seed(seed)

S = do.call(rbind, lapply(1:n, function(i){
if(is.null(rho)){
noise = array(rnorm(prod(imagedim), sd=sd1*c(sdevimg)), dim=imagedim)
}else{
noise = array(c(rmnorm(1, sigma=sigma1)), dim=imagedim)
}
data1 = baseimg + c1 * Z[i] * diffimg
data2 = data1
data2[data2 != 0] = data2[data2 != 0] + noise[data2 != 0]
data2[data2 < 0] = 0
data2 = c(data2)
data2[brainpos]
}))

if(zeromask){
zeropos = apply(S, 2, function(x) all(x == 0))
S = S[,!zeropos] 
brainpos[brainpos] = !zeropos 
}

list(S=S, Z=Z, brainpos=brainpos, imagedim=imagedim)
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
sizechange = function(img1, simscale=NULL, refsize=NULL, ...){

orgsize = dim(img1)
if(!is.null(simscale) & is.null(refsize)) refsize = ceiling(dim(img1)*simscale)
res = refsize / orgsize

if(inherits(img1, "nifti")){cimg1 = img1@.Data}else{cimg1 = img1}
if(length(dim(cimg1))==3) dim(cimg1) = c(dim(cimg1), 1)
cimg1 = as.cimg(cimg1)

cimg1r = resize(cimg1, refsize[1], refsize[2], refsize[3],...)

if(inherits(img1, "nifti")){
img1r = img1
img1r@.Data = cimg1r[,,,1]
img1r@dim_[2:4] = refsize[1:3]
img1r@pixdim[2:4] = img1r@pixdim[2:4] / res[1:3]
}else{
img1r = cimg1r[,,,1]
}

img1r
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
#' data(baseimg)
#' data(diffimg)
#' data(mask)
#' data(template)
#' img1 = simbrain(baseimg = baseimg, diffimg = diffimg, mask=mask)
#' B1 = rbfunc(imagedim=img1$imagedim, seppix=2, hispec=FALSE, mask=img1$brainpos)
#' SB1 = basisprod(img1$S, B1)
#' fit111 = msma(SB1, comp=2)
#' predmodel = ptest(fit111, Z=img1$Z)
#' }
#' 
ptest = function(object, Z=Z, newdata=NULL, testZ=NULL, regmethod = "glm", methods1 = c("boot", "boot632", "cv", "repeatedcv", "LOOCV", "LGOCV")[4], metric = "ROC", number1 = 10, repeats1 = 5, params=NULL){

## Compute the area under the ROC curve, sensitivity, specificity, accuracy and Kappa
glmFuncs <- caret::lrFuncs
fivestats <- function(...) c(caret::twoClassSummary(...), caret::defaultSummary(...), prSummary2(...))
glmFuncs$summary <- fivestats
ctrl <- caret::trainControl(method = methods1, number = number1, repeats = repeats1, classProbs = TRUE, summaryFunction = fivestats)

if(inherits(object, "msma")){ Ss = object$ssX; if(inherits(Ss, "list")){Ss=do.call(cbind,Ss)} }else{Ss = object}
colnames(Ss) = paste("c", c(1:ncol(Ss)), sep="")

swdata = data.frame(Z = as.factor(ifelse(Z == 1, "Y", "N")), Ss)
n = nrow(swdata)

if(regmethod == "rf") swdata = swdata[, c(1, 1+which(apply(swdata[,-1], 2, function(x) length(unique(x)) > 2)))]

set.seed(721)
if(regmethod %in% c("glmnet", "rpart")){
trainout = caret::train(Z~., data=swdata, method=regmethod, metric = metric, trControl = ctrl, tuneGrid=params)
}else if(regmethod == "mxnet"){
trainout = caret::train(Z~., data=swdata, method=regmethod, metric = metric, trControl = ctrl, tuneGrid=params, verbose = FALSE)
}else if(regmethod == "nnet"){
trainout = caret::train(Z~., data=swdata, method=regmethod, metric = metric, trControl = ctrl, tuneGrid=params, trace = FALSE, MaxNWts = 5000)
}else{
trainout = caret::train(Z~., data=swdata, method=regmethod, metric = metric, trControl = ctrl, tuneGrid=params, trace = FALSE)
}
cvout = trainout$results[which.max(trainout$results$ROC),]

traincnfmat = caret::confusionMatrix(predict(trainout), swdata$Z, positive="Y", mode="everything")

if(inherits(object, "msma")){ 
evalmeasure = c(scorecvauc=cvout$ROC, bic1=object$bic[[object$comp[1]]][1], bic1=object$bic[[object$comp[1]]][2])
}else{
evalmeasure = cvout
}

if(!is.null(newdata)){
if(inherits(object, "msma")){ 
pred0 = predict(object, newX=newdata)
pred1 = pred0$sbX[[1]]
}else{
pred1 = newdata
}
colnames(pred1) = paste("c", c(1:ncol(pred1)), sep="")
predZ = predict( trainout, newdata=pred1)
predcnfmat = caret::confusionMatrix(predZ, testZ, positive="Y", mode="everything")
}else{
predcnfmat = NA
}

list(object=object, trainout=trainout, scorecvroc=cvout, evalmeasure=evalmeasure, traincnfmat=traincnfmat, predcnfmat=predcnfmat)
}


#' @rdname mand-internal
prSummary2 = function (data, lev = NULL, model = NULL) 
{
    if (length(levels(data$obs)) > 2) 
        stop(paste("Your outcome has", length(levels(data$obs)), 
            "levels. `prSummary`` function isn't appropriate.", 
            call. = FALSE))
    if (!all(levels(data[, "pred"]) == levels(data[, "obs"]))) 
        stop("Levels of observed and predicted data do not match.", 
            call. = FALSE)
    if (!lev[1] %in% colnames(data)) 
        stop(paste("Class probabilities are needed to score models using the", 
            "area under the PR curve. Set `classProbs = TRUE`", 
            "in the trainControl() function."), call. = FALSE)
    c(Precision = caret::precision(data = data$pred, reference = data$obs, relevant = lev[1]), 
Recall = caret::recall(data = data$pred, reference = data$obs, relevant = lev[1]),
F = caret::F_meas(data = data$pred, reference = data$obs, relevant = lev[1]))
}


#' @rdname mand-internal
colbar = function(zlim, col.y=hotmetal(), nticks=4, horizontal=FALSE, axiscol="white", bg="black", ...) {
oldpar <- par(no.readonly = TRUE)
on.exit(par(oldpar))

min=zlim[1]; max=zlim[2]
ticks=round(seq(min, max, len=nticks),2)
scale1 = (length(col.y)-1)/(max-min)
rectcords = sapply(1:(length(col.y)-1), function(i) {
y = (i-1)/scale1 + min
c(y, 0, y+1/scale1, 10)
})

if(horizontal){
x1 = c(min,max); y1 = c(0,10)
axispos = 1
rdix = 1:4
}else{
x1 = c(0,10); y1 = c(min,max)
axispos = 2
rdix = c(2,1,4,3)
}

#par(oma= c(1, 0, 4, 0), mar= c(4, 2, 4, 2), bg=bg, new=TRUE)
#par(oma= c(4, 0, 4, 0), bg=bg, new=TRUE)
par(omi= c(dev.size()[2]*0.1, 0, dev.size()[2]*0.3, 0), bg=bg, new=TRUE)
plot(x1, y1, type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', ylim=c(0,20), ...)
axis(axispos, ticks, las=1, col=axiscol, col.axis=axiscol)
for (i in 1:ncol(rectcords)) rect(rectcords[rdix[1],i], rectcords[rdix[2],i], rectcords[rdix[3],i], rectcords[rdix[4],i], col=col.y[i], border=NA)
}


#' @rdname mand-internal
tim.colors2 = function (n = 64) 
{
    orig <- c("#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF", 
        "#0000DF", "#0000EF", "#0000FF", "#0010FF", "#0020FF", 
        "#0030FF", "#0040FF", "#0050FF", "#0060FF", "#0070FF", 
        "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF", 
        "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF", 
        "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F", 
        "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50", 
        "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00", 
        "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", "#FFAF00", 
        "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", "#FF6000", 
        "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000", 
        "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000", 
        "#AF0000", "#9F0000", "#8F0000", "#800000")[1:32]
    orig = rep(orig, each=2)
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


#' @rdname mand-internal
rmnorm <- function(n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), method=c("eigen", "svd", "chol"), pre0.9_9994 = FALSE, checkSymmetry = TRUE)
{
if (checkSymmetry && !isSymmetric(sigma, tol = sqrt(.Machine$double.eps), check.attributes = FALSE)) {
stop("sigma must be a symmetric matrix")
}
if (length(mean) != nrow(sigma))
stop("mean and sigma have non-conforming size")

method <- match.arg(method)

R <- if(method == "eigen") {
ev <- eigen(sigma, symmetric = TRUE)
if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))){
    warning("sigma is numerically not positive semidefinite")
}
## ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% t(ev$vectors)
## faster for large  nrow(sigma):
t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))
}
else if(method == "svd"){
s. <- svd(sigma)
if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))){
    warning("sigma is numerically not positive semidefinite")
}
t(s.$v %*% (t(s.$u) * sqrt(pmax(s.$d, 0))))
}
else if(method == "chol"){
R <- chol(sigma, pivot = TRUE)
R[, order(attr(R, "pivot"))]
}

retval <- matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = !pre0.9_9994) %*%  R
retval <- sweep(retval, 2, mean, "+")
colnames(retval) <- names(mean)
retval
}

