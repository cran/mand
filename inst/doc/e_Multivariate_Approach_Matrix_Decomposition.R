## ----echo=FALSE, error=FALSE--------------------------------------------------
options(warn=-1)
knitr::opts_chunk$set(eval = FALSE)

## -----------------------------------------------------------------------------
#  if(!require("mand")) install.packages("mand")

## -----------------------------------------------------------------------------
#  library(mand)

## -----------------------------------------------------------------------------
#  data(baseimg)
#  data(diffimg)
#  data(mask)
#  data(sdevimg)

## -----------------------------------------------------------------------------
#  dim(baseimg)

## -----------------------------------------------------------------------------
#  diffimg2 = diffimg * (tmpatlas %in% 37:40)

## -----------------------------------------------------------------------------
#  img1 = simbrain(baseimg = baseimg, diffimg = diffimg2, sdevimg=sdevimg, mask=mask, n0=20, c1=0.01, sd1=0.05)

## -----------------------------------------------------------------------------
#  dim(img1$S)

## -----------------------------------------------------------------------------
#  coat(rec(img1$S[1,], img1$imagedim, mask=img1$brainpos))

## -----------------------------------------------------------------------------
#  sdimg = apply(img1$S, 2, sd)
#  coat(template, rec(sdimg, img1$imagedim, mask=img1$brainpos))

## -----------------------------------------------------------------------------
#  (fit111 = msma(img1$S, comp=2))

## ----fig.width = 4, fig.height = 3--------------------------------------------
#  plot(fit111, v="score", axes = 1:2, plottype="scatter")

## -----------------------------------------------------------------------------
#  midx = 1 ## the index for the modality
#  vidx = 1 ## the index for the component
#  Q = fit111$wbX[[midx]][,vidx]
#  outstat1 = rec(Q, img1$imagedim, mask=img1$brainpos)

## -----------------------------------------------------------------------------
#  coat(template, outstat1)

## -----------------------------------------------------------------------------
#  B1 = rbfunc(imagedim=img1$imagedim, seppix=2, hispec=FALSE,
#  mask=img1$brainpos)

## -----------------------------------------------------------------------------
#  SB1 = basisprod(img1$S, B1)

## -----------------------------------------------------------------------------
#  dim(img1$S)

## -----------------------------------------------------------------------------
#  dim(SB1)

## -----------------------------------------------------------------------------
#  (fit211 = msma(SB1, comp=2))

## -----------------------------------------------------------------------------
#  Q = fit211$wbX[[1]][,1]
#  outstat1 = rec(Q, img1$imagedim, B=B1, mask=img1$brainpos)

## -----------------------------------------------------------------------------
#  outstat2 = -outstat1
#  coat(template, outstat2)

## -----------------------------------------------------------------------------
#  (fit112 = msma(SB1, comp=2, lambdaX=0.075))

## -----------------------------------------------------------------------------
#  Q = fit112$wbX[[midx]][,vidx]
#  outstat1 = rec(Q, img1$imagedim, B=B1, mask=img1$brainpos)
#  outstat2 = outstat1
#  coat(template, outstat2)

## -----------------------------------------------------------------------------
#  atlastable(tmpatlas, outstat2, atlasdataset)

## -----------------------------------------------------------------------------
#  Z = img1$Z

## -----------------------------------------------------------------------------
#  (fit113 = msma(SB1, Z=Z, comp=2, lambdaX=0.075, muX=0.5))

## -----------------------------------------------------------------------------
#  Q = fit113$wbX[[1]][,1]
#  outstat1 = rec(Q, img1$imagedim, B=B1, mask=img1$brainpos)
#  outstat2 = -outstat1
#  coat(template, outstat2)

## -----------------------------------------------------------------------------
#  atlastable(tmpatlas, outstat2, atlasdataset)

## -----------------------------------------------------------------------------
#  Q = fit113$wbX[[1]][,2]
#  outstat1 = rec(Q, img1$imagedim, B=B1, mask=img1$brainpos)
#  outstat2 = -outstat1
#  coat(template, outstat2)

## -----------------------------------------------------------------------------
#  atlastable(tmpatlas, outstat2, atlasdataset)

## ----fig.width = 7, fig.height = 6.5------------------------------------------
#  ws = multirec(fit113, imagedim=img1$imagedim, B=B1,
#  mask=img1$brainpos)
#  multicompplot(ws, template, col4comp=4)

## -----------------------------------------------------------------------------
#  seppixs = 2:7
#  fit115s = lapply(seppixs, function(sp){
#  B1 = rbfunc(imagedim=img1$imagedim, seppix=sp,
#   hispec=FALSE, mask=img1$brainpos)
#  SB1 = basisprod(img1$S, B1)
#  fit=msma(SB1, Z=Z, comp=2, lambdaX=0.075, muX=0.5)
#  list(fit=fit, B1=B1)
#  })

## ----fig.width = 5, fig.height = 3--------------------------------------------
#  par(mfrow=c(2,3), mar=c(1,2,1,2))
#  for(i in 1:length(seppixs)){
#  Q = fit115s[[i]]$fit$wbX[[midx]][,vidx]
#  outstat1 = rec(Q, img1$imagedim, B=fit115s[[i]]$B1,
#  mask=img1$brainpos)
#  coat(template, -outstat1, pseq=10,color.bar=FALSE,
#  paron=FALSE, main=paste("seppix =", seppixs[i]))
#  }

## -----------------------------------------------------------------------------
#  lambdaXs = round(seq(0, 0.2, by=0.005), 3)
#  fit114s = lapply(lambdaXs, function(lam)
#  msma(SB1, Z=Z, comp=2, lambdaX=lam, muX=0.5, type="lasso") )

## ----fig.width = 5, fig.height = 3--------------------------------------------
#  lambdaXs2 = c(0, 0.025, 0.05, 0.075, 0.1, 0.15)
#  par(mfrow=c(2,3), mar=c(1,2,1,2))
#  for(i in which(lambdaXs %in% lambdaXs2)){
#  Q = fit114s[[i]]$wbX[[1]][,1]
#  outstat1 = rec(Q, img1$imagedim, B=B1, mask=img1$brainpos)
#  coat(template, -outstat1, pseq=10,color.bar=FALSE,
#   paron=FALSE, main=paste("lambda =", lambdaXs[i]))
#  }

## -----------------------------------------------------------------------------
#  nzwbXs = unlist(lapply(fit114s, function(x) x$nzwbX[2]))
#  BICs = unlist(lapply(fit114s, function(x) x$bic[2]))
#  (optlam = lambdaXs[which.min(BICs)])
#  (optnzw = nzwbXs[which.min(BICs)])

## ----fig.width = 7, fig.height = 4--------------------------------------------
#  par(mfrow=c(1,2))
#  plot(lambdaXs, BICs, ylab="BIC")
#  abline(v=optlam, col="red", lty=2)
#  plot(nzwbXs, BICs, ylab="", log="x")
#  abline(v=optnzw, col="red", lty=2)

## -----------------------------------------------------------------------------
#  penalties2 = c("lasso", "hard", "scad", "mcp")
#  etas = list(lasso=1, hard=1, scad=c(1, 3.7), mcp=c(2, 3))

## ----fig.width = 6, fig.height = 3.5------------------------------------------
#  xs = seq(-6, 6, by=0.1)
#  par(mfrow=c(2,3), mar=c(2,2,3,2))
#  for(p1 in penalties2){
#  eta1 = etas[[p1]]
#  for(e1 in eta1){
#  sout1 = sparse(xs, 2, type=p1, eta=e1)
#  plot(xs, sout1, xlab="", ylab="",
#  main=paste(p1, "(eta =", e1, ")"), type="b")
#  }}

## ----fig.width = 5, fig.height = 3--------------------------------------------
#  par(mfrow=c(2,3), mar=c(1,2,1,2))
#  for(p1 in penalties2){
#  eta1 = etas[[p1]]
#  for(e1 in eta1){
#  fit = msma(SB1, Z=Z, comp=2, lambdaX=0.025, muX=0.5,
#   type=p1, eta=e1)
#  Q = fit$wbX[[midx]][,vidx]
#  outstat1 = rec(Q, img1$imagedim, B=B1, mask=img1$brainpos)
#  outstat2 = -outstat1
#  coat(template, outstat2, pseq=10, color.bar=FALSE,
#   paron=FALSE, main=paste(p1, "(eta =", e1, ")"))
#  }}

## -----------------------------------------------------------------------------
#  fit114 = msma(SB1, Z=Z, comp=30, muX=0.5)

## ----fig.width = 4, fig.height = 3.5------------------------------------------
#  plot(fit114, v="cpev")
#  abline(h=0.8, lty=2)

## -----------------------------------------------------------------------------
#  (ncomp1 = ncompsearch(SB1, Z=Z, muX=0.5,
#  comps = 50, criterion="BIC"))

## -----------------------------------------------------------------------------
#  (ncomp2 = ncompsearch(SB1, Z=Z, muX=0.5,
#  comps = c(1, seq(5, 30, by=5)), criterion="CV"))
#  

## ----fig.width = 7, fig.height = 4--------------------------------------------
#  par(mfrow=c(1,2))
#  plot(ncomp1)
#  plot(ncomp2)

## -----------------------------------------------------------------------------
#  maxncomp = 5
#  opts = sapply(1:maxncomp, function(c1){
#  opt=regparasearch(SB1, Z=Z, comp=c1, muX=0.5)$optlambdaX
#  fit = msma(SB1, Z=Z, comp=c1, lambdaX=opt, muX=0.5)
#  nz = rep(NA, maxncomp)
#  nz[1:c1] = fit$nzwbX
#  c(c1, round(opt,3), nz)
#  })

## ----results='asis'-----------------------------------------------------------
#  opts1=t(opts)
#  colnames(opts1) = c("#comp", "lambda",
#  paste("comp",1:maxncomp))
#  kable(opts1, "latex", booktabs = T)

## ----fig.width = 4, fig.height = 3.5------------------------------------------
#  (ncomp3 = ncompsearch(SB1, Z=Z, muX=0.5,
#  lambdaX=0.075, comps = 30, criterion="BIC"))
#  plot(ncomp3)

## -----------------------------------------------------------------------------
#  (opt11 = optparasearch(SB1, Z=Z, muX=0.5, comp=5,
#  search.method = "regparaonly", criterion="BIC"))

## -----------------------------------------------------------------------------
#  (fit311 = msma(SB1, Z=Z, muX=0.5,
#  comp=opt11$optncomp, lambdaX=opt11$optlambdaX))

## -----------------------------------------------------------------------------
#  (opt12 = optparasearch(SB1, Z=Z, muX=0.5,
#  search.method = "regpara1st", criterion="BIC"))
#  fit312 = msma(SB1, Z=Z, muX=0.5,
#  comp=opt12$optncomp, lambdaX=opt12$optlambdaX)

## -----------------------------------------------------------------------------
#  (opt13 = optparasearch(SB1, Z=Z, muX=0.5,
#  search.method = "ncomp1st", criterion="BIC"))
#  fit313 = msma(SB1, Z=Z, muX=0.5,
#  comp=opt13$optncomp, lambdaX=opt13$optlambdaX)

## -----------------------------------------------------------------------------
#  (opt14 = optparasearch(SB1, Z=Z, muX=0.5,
#  search.method = "simultaneous", criterion="BIC"))
#  fit314 = msma(SB1, Z=Z, muX=0.5,
#  comp=opt14$optncomp, lambdaX=opt14$optlambdaX)

## ----error=FALSE, message=FALSE-----------------------------------------------
#  if(!require("NMF")) install.packages("NMF")
#  library(NMF)

## -----------------------------------------------------------------------------
#  res = nmf(SB1, 2)

## -----------------------------------------------------------------------------
#  Q = t(coef(res))[,1]
#  outstat1 = rec(Q, img1$imagedim, B=B1, mask=img1$brainpos)
#  coat(template, outstat1)

## ----error=FALSE, message=FALSE-----------------------------------------------
#  if(!require("ica")) install.packages("ica")
#  library(ica)

## -----------------------------------------------------------------------------
#  imod = icaimax(SB1,2)

## -----------------------------------------------------------------------------
#  Q = imod$M[,1]
#  outstat1 = rec(Q, img1$imagedim, B=B1, mask=img1$brainpos)
#  coat(template, outstat1)

## ----echo=FALSE, message=FALSE------------------------------------------------
#  if(!require("dendextend"))install.packages(c("dendextend"))

## -----------------------------------------------------------------------------
#  library(dendextend)

## -----------------------------------------------------------------------------
#  hcmsma111 = hcmsma(fit111)

## ----fig.width = 6, fig.height = 4--------------------------------------------
#  dend = as.dendrogram(hcmsma111$hcout)
#  d1 = color_branches(dend, k=4, groupLabels=TRUE)
#  labels_colors(d1) = Z[as.numeric(labels(d1))]+1
#  plot(d1)

## -----------------------------------------------------------------------------
#  clus=cutree(d1, 4, order_clusters_as_data = FALSE)
#  clus=clus[as.character(1:length(clus))]
#  table(Z, clus)

## ----fig.width = 6, fig.height = 4--------------------------------------------
#  hcmsma112 = hcmsma(fit112)
#  dend = as.dendrogram(hcmsma112$hcout)
#  d1 = color_branches(dend, k=4, groupLabels=TRUE)
#  labels_colors(d1) = Z[as.numeric(labels(d1))]+1
#  plot(d1)
#  clus=cutree(d1, 4, order_clusters_as_data = FALSE)
#  clus=clus[as.character(1:length(clus))]
#  table(Z, clus)

## ----fig.width = 6, fig.height = 4--------------------------------------------
#  hcmsma113 = hcmsma(fit113)
#  dend = as.dendrogram(hcmsma113$hcout)
#  d1 = color_branches(dend, k=4, groupLabels=TRUE)
#  labels_colors(d1) = Z[as.numeric(labels(d1))]+1
#  plot(d1)
#  clus=cutree(d1, 4, order_clusters_as_data = FALSE)
#  clus=clus[as.character(1:length(clus))]
#  table(Z, clus)

