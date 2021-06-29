## -----------------------------------------------------------------------------
library(mand)

## -----------------------------------------------------------------------------
data(template)

## -----------------------------------------------------------------------------
dim(template)

## -----------------------------------------------------------------------------
coat(template)

## ----fig.width = 5, fig.height = 2--------------------------------------------
coat(x=template, plane="all")

## ----eval=FALSE---------------------------------------------------------------
#  fnames1 = c("data1.nii", "data2.nii")
#  imgmat = imgdatamat(fnames1, simscale=1/4)

## ----eval=FALSE---------------------------------------------------------------
#  imgmat = imgdatamat(fnames1, simscale=1/4, ROI=TRUE)

## ----echo=FALSE, eval=FALSE---------------------------------------------------
#  '

## -----------------------------------------------------------------------------
data(diffimg)
coat(template, diffimg)

## -----------------------------------------------------------------------------
data(atlasdatasets)
atlasname = "aal"
atlasdataset = atlasdatasets[[atlasname]]
head(atlasdataset)

## -----------------------------------------------------------------------------
data(atlas)
tmpatlas = atlas[[atlasname]]
dim(tmpatlas)

## -----------------------------------------------------------------------------
tmpatlas[11:15,11:15,10]

## -----------------------------------------------------------------------------
coat(template, tmpatlas, regionplot=TRUE, 
atlasdataset=atlasdataset, ROIids = c(1:2, 37:40), 
regionlegend=TRUE)

## -----------------------------------------------------------------------------
atlastable(tmpatlas, diffimg, atlasdataset, ROIids = c(1:2, 
37:40))

## -----------------------------------------------------------------------------
hipmask = (tmpatlas == 37) + (tmpatlas == 38)
template2 = template * hipmask

## ----fig.width = 5, fig.height = 2--------------------------------------------
par(mfrow=c(1,3), mar=rep(1,4))
coat(template, pseq=11, paron=FALSE)
coat(hipmask, pseq=11, paron=FALSE)
coat(template2, pseq=11, paron=FALSE)

## -----------------------------------------------------------------------------
sum(template[which(hipmask==1, arr.ind = TRUE)])/1000

## -----------------------------------------------------------------------------
data(baseimg)
data(diffimg)
data(mask)
data(sdevimg)

## -----------------------------------------------------------------------------
dim(baseimg)

## -----------------------------------------------------------------------------
diffimg2 = diffimg * (tmpatlas %in% 37:40)

## -----------------------------------------------------------------------------
img1 = simbrain(baseimg = baseimg, diffimg = diffimg2, 
sdevimg=sdevimg, mask=mask, n0=20, c1=0.01, sd1=0.05)

## -----------------------------------------------------------------------------
dim(img1$S)

## -----------------------------------------------------------------------------
coat(rec(img1$S[1,], img1$imagedim, mask=img1$brainpos))

## -----------------------------------------------------------------------------
sdimg = apply(img1$S, 2, sd)
coat(template, rec(sdimg, img1$imagedim, mask=img1$brainpos))

## -----------------------------------------------------------------------------
(fit111 = msma(img1$S, comp=2))

## -----------------------------------------------------------------------------
plot(fit111, v="score", axes = 1:2, plottype="scatter")

## -----------------------------------------------------------------------------
midx = 1 ## the index for the modality
vidx = 1 ## the index for the component
Q = fit111$wbX[[midx]][,vidx]
outstat1 = rec(Q, img1$imagedim, mask=img1$brainpos)

## -----------------------------------------------------------------------------
coat(template, outstat1)

## -----------------------------------------------------------------------------
B1 = rbfunc(imagedim=img1$imagedim, seppix=2, hispec=FALSE, 
mask=img1$brainpos)

## -----------------------------------------------------------------------------
SB1 = basisprod(img1$S, B1)

## -----------------------------------------------------------------------------
dim(img1$S)

## -----------------------------------------------------------------------------
dim(SB1)

## -----------------------------------------------------------------------------
(fit211 = msma(SB1, comp=2))

## -----------------------------------------------------------------------------
Q = fit211$wbX[[1]][,1]
outstat1 = rec(Q, img1$imagedim, B=B1, mask=img1$brainpos)

## -----------------------------------------------------------------------------
outstat2 = -outstat1
coat(template, outstat2)

## -----------------------------------------------------------------------------
(fit112 = msma(SB1, comp=2, lambdaX=0.075))

## -----------------------------------------------------------------------------
Q = fit112$wbX[[midx]][,vidx]
outstat1 = rec(Q, img1$imagedim, B=B1, mask=img1$brainpos)
outstat2 = outstat1
coat(template, outstat2)

## -----------------------------------------------------------------------------
atlastable(tmpatlas, outstat2, atlasdataset)

## -----------------------------------------------------------------------------
Z = img1$Z

## -----------------------------------------------------------------------------
(fit113 = msma(SB1, Z=Z, comp=2, lambdaX=0.075, muX=0.5))

## -----------------------------------------------------------------------------
Q = fit113$wbX[[1]][,1]
outstat1 = rec(Q, img1$imagedim, B=B1, mask=img1$brainpos)
outstat2 = -outstat1
coat(template, outstat2)

## -----------------------------------------------------------------------------
atlastable(tmpatlas, outstat2, atlasdataset)

## -----------------------------------------------------------------------------
Q = fit113$wbX[[1]][,2]
outstat1 = rec(Q, img1$imagedim, B=B1, mask=img1$brainpos)
outstat2 = -outstat1
coat(template, outstat2)

## -----------------------------------------------------------------------------
atlastable(tmpatlas, outstat2, atlasdataset)

## ----fig.width = 7, fig.height = 6.5------------------------------------------
ws = multirec(fit113, imagedim=img1$imagedim, B=B1, 
mask=img1$brainpos)
multicompplot(ws, template, col4comp=4)

