## ----echo=FALSE, error=FALSE--------------------------------------------------
options(warn=-1)
knitr::opts_chunk$set(eval = FALSE)

## -----------------------------------------------------------------------------
#  n = 20; x = 1:n
#  set.seed(1); y = abs(rnorm(n)); names(y) = x

## -----------------------------------------------------------------------------
#  cuty = function(y1, cdt){
#  rbind(ifelse(y1-cdt<0, y1, cdt),
#  ifelse(y1-cdt<0, 0, y1-cdt))
#  }

## ----fig.width = 7, fig.height = 4--------------------------------------------
#  par(mfrow=c(1, 2), mar=c(4,3,4,3))
#  for(cdt in c(0.5, 1.5)){
#  y2 = cuty(y, cdt)
#  barplot(y2, col=c(8,2), main=paste("CDT =", cdt))
#  abline(h = cdt, lty=2, col=2)
#  }

## -----------------------------------------------------------------------------
#  cdt = 0.75

## ----fig.width = 7, fig.height = 4--------------------------------------------
#  par(mfrow=c(1,2), mar=c(4,3,4,3))
#  for(f1 in c(2,8)){
#  sy = gsmooth(x, y, f1); y2 = cuty(sy, cdt)
#  barplot(y2, col=c(8,2), main=paste("FWHM =", f1), ylim=c(0,2))
#  abline(h = cdt, lty=2, col=2)
#  }

## -----------------------------------------------------------------------------
#  exp(-0.53*(96/325*3.3*(3.3^2-1)^(2/3)))

## -----------------------------------------------------------------------------
#  s=96; h=3.3; r=325;
#  Es=(4*log(2))^(3/2)*h*(h^2-1)/(2*pi)^(3/2);
#  d=(gamma(5/2)*Es)^(2/3);
#  exp(-d*(s/r)^(2/3))

## -----------------------------------------------------------------------------
#  (325*(-log(0.05)/0.53)^(3/2))/(3.3*(3.3^2-1))

## -----------------------------------------------------------------------------
#  ss = c(1, 10, 50)
#  h = 4
#  b = (4*log(2)/(2*pi))*(gamma(5/2))^(2/3)
#  fs = 2:10

## ----fig.width = 5, fig.height = 3.5------------------------------------------
#  par(mar=c(2,2,2,6))
#  for(sidx in 1:length(fs)){
#  p = exp(-b*(ss[sidx]/fs^3*h*((h)^2-1))^(2/3))
#  if(sidx == 1){
#  plot(fs, p, lty = sidx, type="l", ylim=c(0,1),
#  xlab="FWHM", ylab="p-value", main="")
#  }else{
#  points(fs, p, lty = sidx, type="l")
#  }}
#  abline(h=0.05, lty=5, lwd=2)
#  par(xpd=TRUE)
#  legend(par()$usr[2], par()$usr[4],
#  legend=paste("s =", ss[order(ss)]), lty=order(ss))

## -----------------------------------------------------------------------------
#  TFCE0 = function(y1, cdt, E=0.5, H=2){
#  cy = 1*(y1 >= cdt)
#  
#  clsta0 = c(1, as.numeric(names(which(diff(cy)!=0))))
#  clend0 = c(clsta0[-1]-1, length(cy))
#  clsta = clsta0[cy[clsta0]==1]
#  clend = clend0[cy[clsta0]==1]
#  
#  clustidxs = lapply(1:length(clsta),
#  function(i) clsta[i]:clend[i])
#  clustsize = unlist(lapply(clustidxs, length))
#  
#  x = 1:length(y1)
#  clust = unlist(lapply(x, function(x1){
#  a = which(unlist(lapply(clustidxs,
#  function(cidx) x1 %in% cidx)))
#  ifelse(length(a)>0, a, NA)
#  }))
#  
#  a = clustsize[clust]^E * cdt^H
#  ifelse(is.na(a),0, a)
#  }

## -----------------------------------------------------------------------------
#  TFCE1 = function(f1){
#  for(cdt1 in cdts2){
#  sy = gsmooth(x, y, f1)
#  y2 = cuty(sy, cdt1)
#  barplot(y2, col=c(8,2), main=paste("CDT =", cdt1))
#  abline(h = cdt1, lty=2, col=2)
#  tfce = TFCE0(sy,cdt1)
#  barplot(tfce, main=paste("TFCE CDT =", cdt1), ylim=c(0,10))
#  }
#  }

## ----fig.width = 6, fig.height = 5--------------------------------------------
#  cdts2 = c(0.75, 0.5, 0.25)
#  par(mfrow=c(length(cdts2), 2), mar=c(3,4,2,4))
#  TFCE1(f1=2)

## ----fig.width = 7, fig.height = 6--------------------------------------------
#  f1s = c(0, 2, 4, 8)
#  par(mfrow=c(length(f1s),2), mar=c(3,4,2,4))
#  for(f1 in f1s){
#  if(f1 == 0){barplot(y, main="Original")}else{
#  sy = gsmooth(x, y, f1)
#  barplot(sy, main=paste("FWHM =", f1), ylim=c(0,2))
#  }
#  if(f1 > 0){ sy = gsmooth(x, y, f1) }else{ sy = y }
#  cdts = c(0, sort(unique(sy)))
#  tfce = colSums(do.call(rbind,lapply(cdts, function(cdt1){
#  TFCE0(sy,cdt1)
#  })))
#  barplot(tfce, main=paste("TFCE FWHM =", f1))
#  }

## -----------------------------------------------------------------------------
#  n = 5; x0 = rep(c(0,1), each = n)
#  set.seed(1); y = ifelse(x0==0, rnorm(n, 0), rnorm(n, 1))

## -----------------------------------------------------------------------------
#  fitfunc = function(idx){
#  x = x0[idx]
#  fit = summary(lm(y~x))
#  coef(fit)[2,1:3]
#  }

## -----------------------------------------------------------------------------
#  (fit0 = fitfunc(1:(2*n)))

## -----------------------------------------------------------------------------
#  set.seed(1); idx = sample(2*n); round(fitfunc(idx), 3)
#  
#  set.seed(2); idx = sample(2*n); round(fitfunc(idx), 3)
#  
#  set.seed(1000); idx = sample(2*n); round(fitfunc(idx), 3)

## -----------------------------------------------------------------------------
#  tststs = sapply(1:1000, function(s1){
#  set.seed(s1); idx = sample(2*n); fitfunc(idx)[3]
#  })

## -----------------------------------------------------------------------------
#  (q1 = quantile(tststs, prob=0.975))
#  ifelse(abs(fit0[3])>q1, "reject", "not reject")

## -----------------------------------------------------------------------------
#  (q2 = qt(0.975, 2*n-2) )
#  ifelse(abs(fit0[3])>q2, "reject", "not reject")

## ----fig.width = 4, fig.height = 3.5------------------------------------------
#  hist(tststs, xlab="T stat", main="", freq=FALSE,
#  ylim=c(0, 0.4))
#  abline(v=fit0[3], col=2, lty=2)
#  abline(v=q1, col=4, lty=3); abline(v=q2, col=3, lty=4)
#  f1 = function(x)dt(x, 2*n-2)
#  curve(f1, -5, 5, add=TRUE, col=3)

## -----------------------------------------------------------------------------
#  (pp = 2*(mean(tststs > fit0[3])))

## -----------------------------------------------------------------------------
#  (tp = 2*(1 - pt(fit0[3], 2*n-2)))

## -----------------------------------------------------------------------------
#  nreps = round(c(seq(10, 100, length=10),
#  seq(200, length(tststs)-100, length=10)))
#  q1s = sapply(nreps, function(x) {
#  sapply(1:1000, function(y) {
#  set.seed(y); quantile(tststs[sample(1:length(tststs), x)],
#  0.975)})})
#  colnames(q1s)=nreps

## ----fig.width = 4, fig.height = 3.5------------------------------------------
#  boxplot(q1s, main="", xlab="Number of permutations",
#  ylab="Two-sided 95% point", type="b")
#  abline(h=q2, col=3, lty=2)

## -----------------------------------------------------------------------------
#  clustsize = function(y1, cdt){
#  cy = 1*(y1 >= cdt)
#  
#  if(all(cy==0)){
#  clustsize = 0
#  }else{
#  clsta0 = c(1, as.numeric(names(which(diff(cy)!=0))))
#  clend0 = c(clsta0[-1]-1, length(cy))
#  clsta = clsta0[cy[clsta0]==1]
#  clend = clend0[cy[clsta0]==1]
#  
#  clustidxs = lapply(1:length(clsta),
#  function(i) clsta[i]:clend[i])
#  clustsize = unlist(lapply(clustidxs, length))
#  }
#  clustsize
#  }

## -----------------------------------------------------------------------------
#  gsmooth = function(x, y, FWHM){
#  sigma = FWHM / 2*sqrt(2*log(2))
#  sy = sapply(x, function(x1)
#  weighted.mean(y,
#  dnorm(x, x1, sigma)/sum(dnorm(x, x1,sigma))) )
#  names(sy) = x
#  sy
#  }

## -----------------------------------------------------------------------------
#  cuty = function(y1, cdt){
#  rbind(ifelse(y1-cdt<0, y1, cdt),
#  ifelse(y1-cdt<0, 0, y1-cdt))
#  }

## -----------------------------------------------------------------------------
#  n = 20; x = 1:n; f1 = 2; cdt1 = 0.75

## ----fig.width = 7, fig.height = 6--------------------------------------------
#  sidxs = 1:4
#  par(mfrow=c(length(sidxs),1), mar=c(3,4,2,4))
#  for(sidx in sidxs){
#  sigma1 = ifelse(sidx==1, 1, 0.9)
#  set.seed(sidx); y = abs(rnorm(n,0,sigma1)); names(y) = x
#  sy = gsmooth(x, y, f1); y2 = cuty(sy, cdt1)
#  csize = clustsize(sy, cdt1)
#  barplot(y2, col=c(8,2),
#  main=paste("Max Cluster Size =", max(csize)), ylim=c(0,1.2))
#  abline(h = cdt1, lty=2, col=2)
#  }

## -----------------------------------------------------------------------------
#  sidxs2 = 1:1001
#  maxcsizes = sapply(sidxs2, function(sidx){
#  sigma1 = ifelse(sidx==1, 1, 0.9)
#  set.seed(sidx); y = abs(rnorm(n,0,sigma1)); names(y) = x
#  sy = gsmooth(x, y, f1)
#  csize = clustsize(sy, cdt1)
#  max(csize)
#  })
#  (obs = maxcsizes[1])
#  (q1 = quantile(maxcsizes[-1], 0.95))

## ----fig.width = 4, fig.height = 3.5------------------------------------------
#  hist(maxcsizes, main="", xlab="Max Cluster Size")
#  abline(v = q1, col=2, lty=2); abline(v = obs, col=3, lty=2)

## ----fig.width = 4, fig.height = 3.5------------------------------------------
#  nreps = seq(10, length(maxcsizes), length=100)
#  q1s = sapply(nreps, function(x)
#  quantile(maxcsizes[2:x], 0.95))
#  par(mfrow=c(1,1))
#  plot(nreps, q1s, main="", xlab="Number of permutations",
#  ylab="95% point", type="b")

