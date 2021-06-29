## ----echo=FALSE, error=FALSE--------------------------------------------------
options(warn=-1)
knitr::opts_chunk$set(eval = FALSE)

## -----------------------------------------------------------------------------
#  n = 100; seed = 2
#  dataset1 = strsimdata(n = n, ncomp=2,
#  Xps=c(4,4), Ztype="binary", seed=seed)

## -----------------------------------------------------------------------------
#  X2 = dataset1$X;
#  Z = dataset1$Z
#  str(dataset1[c("X","Z")])

## -----------------------------------------------------------------------------
#  dataset1$WX

## -----------------------------------------------------------------------------
#  dataset1$nZeroX

## -----------------------------------------------------------------------------
#  dataset1$ZcoefX

## -----------------------------------------------------------------------------
#  (opt212 = optparasearch(X2, Z=Z, muX=0.5,
#  search.method = "ncomp1st", criterion="BIC", whichselect="X"))

## -----------------------------------------------------------------------------
#  (fit212 = msma(X2, Z=Z, muX=0.5, comp=opt212$optncomp,
#  lambdaX=opt212$optlambdaX))

## ----fig.width = 7, fig.height = 6--------------------------------------------
#  par(mfrow=c(2,2), oma = c(0, 0, 2, 0))
#  plot(fit212, axes = 1, plottype="bar",
#  block="super", XY="X")
#  plot(fit212, axes = 2, plottype="bar",
#  block="super", XY="X")
#  plot(fit212, axes = 1, plottype="bar",
#  block="block", XY="X")
#  plot(fit212, axes = 2, plottype="bar",
#  block="block", XY="X")

## ----fig.width = 7, fig.height = 3.5------------------------------------------
#  par(mfrow=c(1,2))
#  for(i in 1:2){
#  t1=t.test(fit212$ssX[,i]~Z)
#  boxplot(fit212$ssX[,i]~Z,
#  main=paste("Comp", i),
#  sub=paste("t-test p =", round(t1$p.value,4)))
#  }

## -----------------------------------------------------------------------------
#  dataset2 = strsimdata(n = n, ncomp=2, Xps=c(4,4),
#  Yps=c(3,5), Ztype="binary", cz=c(10,10), seed=seed)

## -----------------------------------------------------------------------------
#  X2 = dataset2$X; Y2 = dataset2$Y
#  Z = dataset2$Z
#  str(dataset2[c("X","Y","Z")])

## -----------------------------------------------------------------------------
#  dataset2$WX

## -----------------------------------------------------------------------------
#  dataset2$nZeroX

## -----------------------------------------------------------------------------
#  dataset2$WY
#  dataset2$nZeroY

## -----------------------------------------------------------------------------
#  dataset2$ZcoefX
#  dataset2$ZcoefY

## -----------------------------------------------------------------------------
#  (opt222 = optparasearch(X2, Y2, Z, muX=0.3, muY=0.3,
#  search.method = "ncomp1st", criterion="BIC",
#  criterion4ncomp="BIC", whichselect=c("X","Y")))
#  (fit222 = msma(X2, Y2, Z,
#  muX=0.3, muY=0.3, comp=opt222$optncomp,
#  lambdaX=opt222$optlambdaX, lambdaY=opt222$optlambdaY))

## ----fig.width = 7, fig.height = 6--------------------------------------------
#  par(mfrow=c(2,2), oma = c(0, 0, 2, 0))
#  plot(fit222, axes = 1, plottype="bar",
#  block="super", XY="X")
#  plot(fit222, axes = 2, plottype="bar",
#  block="super", XY="X")
#  plot(fit222, axes = 1, plottype="bar",
#  block="block", XY="X")
#  plot(fit222, axes = 2, plottype="bar",
#  block="block", XY="X")

## ----fig.width = 7, fig.height = 6--------------------------------------------
#  par(mfrow=c(2,2), oma = c(0, 0, 2, 0))
#  plot(fit222, axes = 1, plottype="bar",
#  block="super", XY="Y")
#  plot(fit222, axes = 2, plottype="bar",
#  block="super", XY="Y")
#  plot(fit222, axes = 1, plottype="bar",
#  block="block", XY="Y")
#  plot(fit222, axes = 2, plottype="bar",
#  block="block", XY="Y")

## ----fig.width = 7, fig.height = 3.5------------------------------------------
#  par(mfrow=c(1,2))
#  for(i in 1:2) plot(fit222, axes = i, XY="XY")

## ----fig.width = 7, fig.height = 6--------------------------------------------
#  par(mfrow=c(2,2))
#  for(xy in c("X","Y")){for(i in 1:2){
#  ss = fit222[[paste0("ss", xy)]][,i]
#  t1=t.test(ss~Z)
#  boxplot(ss~Z, main=paste(xy, "Comp", i),
#  sub=paste("t-test p =", round(t1$p.value,4)))
#  }}

