## ----echo=FALSE, error=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
options(warn=-1)
knitr::opts_chunk$set(eval = FALSE)

## ----eval=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  if(!require("mand")) install.packages("mand")

## ----eval=TRUE, echo=TRUE, message=FALSE, error=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(mand)

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  library(oro.nifti)
#  library(oro.dicom)
#  library("imager")

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  data(template)

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  dim(template)

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  coat(template)

## ----fig.width = 5, fig.height = 2------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  coat(x=template, plane="all")

## ----eval=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  fnames1 = c("data1.nii", "data2.nii")
#  imgmat = imgdatamat(fnames1, simscale=1/4)

## ----eval=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  imgmat = imgdatamat(fnames1, simscale=1/4, ROI=TRUE)

## ----echo=FALSE, eval=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  '

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  data(diffimg)
#  coat(template, diffimg)

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  data(atlasdatasets)
#  atlasname = "aal"
#  atlasdataset = atlasdatasets[[atlasname]]
#  head(atlasdataset)

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  data(atlas)
#  tmpatlas = atlas[[atlasname]]
#  dim(tmpatlas)

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  tmpatlas[11:15,11:15,10]

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  coat(template, tmpatlas, regionplot=TRUE,
#  atlasdataset=atlasdataset, ROIids = c(1:2, 37:40),
#  regionlegend=TRUE)

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  atlastable(tmpatlas, diffimg, atlasdataset, ROIids = c(1:2,
#  37:40))

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  hipmask = (tmpatlas == 37) + (tmpatlas == 38)
#  template2 = template * hipmask

## ----fig.width = 5, fig.height = 2------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  par(mfrow=c(1,3), mar=rep(1,4))
#  coat(template, pseq=11, paron=FALSE)
#  coat(hipmask, pseq=11, paron=FALSE)
#  coat(template2, pseq=11, paron=FALSE)

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  sum(template[which(hipmask==1, arr.ind = TRUE)])/1000

