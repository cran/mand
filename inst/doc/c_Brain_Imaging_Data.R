## ----echo=FALSE, error=FALSE--------------------------------------------------
options(warn=-1)
knitr::opts_chunk$set(eval = FALSE)

## -----------------------------------------------------------------------------
#  X = matrix(c(
#  -1, -1, 1, 1,
#  -1, 1, -1, 1
#  ),ncol=4, byrow=TRUE)

## -----------------------------------------------------------------------------
#  par(mfrow=c(1,1), cex=0.4)
#  plot(NULL, xlim=c(-8, 8), ylim=c(-8, 8), xlab="x", ylab="y")
#  
#  rect(xleft=min(X[1,]), ybottom=min(X[2,]),
#  xright=max(X[1,]), ytop=max(X[2,]))
#  
#  vs = sapply(1:4, function(x)
#  paste("(", paste(X[,x], collapse=", "), ")"))
#  text(X[1,], X[2,], paste(LETTERS[1:4], vs), pos=c(1,3,1,3))

## -----------------------------------------------------------------------------
#  Y = X + c(5,5)

## ----fig.width = 2.5, fig.height = 2.5----------------------------------------
#  par(mfrow=c(1,1), cex=0.4)
#  plot(NULL, xlim=c(-8, 8), ylim=c(-8, 8), xlab="x", ylab="y")
#  rect(xleft=min(X[1,]), ybottom=min(X[2,]),
#  xright=max(X[1,]), ytop=max(X[2,]))
#  rect(xleft=min(Y[1,]), ybottom=min(Y[2,]),
#  xright=max(Y[1,]), ytop=max(Y[2,]))
#  vs = sapply(1:4, function(x)
#  paste("(", paste(Y[,x], collapse=", "), ")"))
#  text(Y[1,], Y[2,], paste(LETTERS[1:4], vs), pos=c(1,3,1,3))

## -----------------------------------------------------------------------------
#  T = matrix(c(
#  2, 0,
#  0, 1
#  ),ncol=2,byrow=TRUE)
#  s = c(-5,-5)
#  Y = T %*% X + s

## ----fig.width = 2.5, fig.height = 2.5----------------------------------------
#  par(mfrow=c(1,1), cex=0.4)
#  plot(NULL, xlim=c(-8, 8), ylim=c(-8, 8), xlab="x", ylab="y")
#  rect(xleft=min(X[1,]), ybottom=min(X[2,]),
#  xright=max(X[1,]), ytop=max(X[2,]))
#  rect(xleft=min(Y[1,]), ybottom=min(Y[2,]),
#  xright=max(Y[1,]), ytop=max(Y[2,]))
#  vs = sapply(1:4, function(x)
#  paste("(", paste(Y[,x], collapse=", "), ")"))
#  text(Y[1,], Y[2,], paste(LETTERS[1:4], vs), pos=c(1,3,1,3))

## -----------------------------------------------------------------------------
#  theta = pi/6
#  T = matrix(c(
#  cos(theta), -sin(theta),
#  sin(theta), cos(theta)
#  ),ncol=2,byrow=TRUE)
#  s = c(-5, 5)
#  Y = T %*% X + s

## ----fig.width = 2.5, fig.height = 2.5----------------------------------------
#  par(mfrow=c(1,1), cex=0.4)
#  plot(NULL, xlim=c(-8, 8), ylim=c(-8, 8), xlab="x", ylab="y")
#  rect(xleft=min(X[1,]), ybottom=min(X[2,]),
#  xright=max(X[1,]), ytop=max(X[2,]))
#  lines(Y[1, c(1, 2)], Y[2, c(1, 2)])
#  lines(Y[1, c(1, 3)], Y[2, c(1, 3)])
#  lines(Y[1, c(2, 4)], Y[2, c(2, 4)])
#  lines(Y[1, c(3, 4)], Y[2, c(3, 4)])
#  vs = sapply(1:4, function(x)
#  paste("(", paste(round(Y[,x],1), collapse=", "), ")"))
#  text(Y[1,], Y[2,], paste(LETTERS[1:4], vs), pos=c(1,3,1,3))

## -----------------------------------------------------------------------------
#  T = matrix(c(
#  1, 0,
#  tan(pi/6), 1
#  ),ncol=2,byrow=TRUE)
#  s = c(5, -5)
#  Y = T %*% X + s

## ----fig.width = 2.5, fig.height = 2.5----------------------------------------
#  par(mfrow=c(1,1), cex=0.4)
#  plot(NULL, xlim=c(-8, 8), ylim=c(-8, 8), xlab="x", ylab="y")
#  rect(xleft=min(X[1,]), ybottom=min(X[2,]),
#  xright=max(X[1,]), ytop=max(X[2,]))
#  lines(Y[1, c(1, 2)], Y[2, c(1, 2)])
#  lines(Y[1, c(1, 3)], Y[2, c(1, 3)])
#  lines(Y[1, c(2, 4)], Y[2, c(2, 4)])
#  lines(Y[1, c(3, 4)], Y[2, c(3, 4)])
#  vs = sapply(1:4, function(x)
#  paste("(", paste(round(Y[,x],1), collapse=", "), ")"))
#  text(Y[1,], Y[2,], paste(LETTERS[1:4], vs), pos=c(1,3,1,3))

## ----fig.width = 5, fig.height = 2--------------------------------------------
#  simscale1 = 4
#  img1r = sizechange(template, simscale=simscale1)
#  dim(img1r)
#  coat(img1r, plane="all")

## ----fig.width = 5, fig.height = 2--------------------------------------------
#  simscale1 = 1/2
#  img1r = sizechange(template, simscale=simscale1)
#  dim(img1r)
#  coat(img1r, plane="all")

## -----------------------------------------------------------------------------
#  (a = sqrt(2*log(2)))

## ----fig.width = 4, fig.height = 3--------------------------------------------
#  plot(dnorm, -4, 4, xaxt = "n", ylab="f(x)")
#  axis(1, -4:4, -4:4)
#  lines(c(0,0), c(-1,dnorm(0)), lty=2, col=1)
#  lines(c(-a,a), c(dnorm(0)/2,dnorm(0)/2), lty=2, col=1)
#  lines(c(a,a), c(-1,dnorm(a)), lty=2, col=2)
#  lines(-c(a,a), c(-1,dnorm(a)), lty=2, col=2)

## ----fig.width = 7, fig.height = 3.5------------------------------------------
#  set.seed(1)
#  n = 20
#  x = 1:n
#  y = abs(rnorm(n)); names(y) = x
#  barplot(y)

## -----------------------------------------------------------------------------
#  g1 = function(x, x1, FWHM){
#  sigma = FWHM / 2*sqrt(2*log(2))
#  d1 = dnorm(x, x1, sigma)
#  names(d1) = x
#  d1
#  }

## -----------------------------------------------------------------------------
#  gsmooth = function(x, y, FWHM){
#  sigma = FWHM / 2*sqrt(2*log(2))
#  sy = sapply(x, function(x1)
#  weighted.mean(y,
#  dnorm(x, x1, sigma)/sum(dnorm(x, x1, sigma))) )
#  names(sy) = x; sy; }

## ----fig.width = 5, fig.height = 5--------------------------------------------
#  x1 = 5
#  col1 = rep(1, n); col1[x1] = 2
#  par(mfrow=c(4,1), mar=c(3,4,1,4))
#  barplot(y, col=col1)
#  barplot(g1(x,x1,2), col=col1, ylim=c(0, 0.4))
#  barplot(y*g1(x,x1,2), col=col1)
#  sy = rep(0, n); sy[x1] = sum(y*g1(x,x1,2)); names(sy) = x
#  barplot(sy, col=col1, ylim=c(0,2))

## ----fig.width = 6, fig.height = 5--------------------------------------------
#  col2 = rep(1, n)
#  par(mfrow=c(2,1))
#  barplot(y, col=col2, main="Original")
#  f1 = 2
#  sy = gsmooth(x, y, f1)
#  barplot(sy, main=paste("FWHM =", f1), col=col2, ylim=c(0,2))

## ----fig.width = 5, fig.height = 5--------------------------------------------
#  f1s = c(2, 4, 8)
#  layout(cbind(c(8,1,2,3), c(4:7)))
#  par(mar=c(3,4,2,4))
#  for(f1 in f1s){
#  barplot(g1(x,8,f1), main=paste("FWHM =", f1), ylim=c(0,0.4))
#  }
#  barplot(y, main="Original")
#  for(f1 in f1s){
#  sy = gsmooth(x, y, f1)
#  barplot(sy, main=paste("FWHM =", f1), ylim=c(0,2))
#  }

