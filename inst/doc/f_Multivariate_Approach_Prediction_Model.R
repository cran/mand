## ----echo=FALSE, error=FALSE--------------------------------------------------
options(warn=-1)
knitr::opts_chunk$set(eval = FALSE)

## -----------------------------------------------------------------------------
#  if(!require("mand")) install.packages("mand")

## -----------------------------------------------------------------------------
#  library(mand)

## -----------------------------------------------------------------------------
#  fit113 = msma(SB1, Z=Z, comp=2, lambdaX=0.075, muX=0.5)
#  Ss = fit113$ssX
#  colnames(Ss) = paste("c", c(1:ncol(Ss)), sep="")
#  swdata113 = data.frame(
#  Z = as.factor(ifelse(Z == 1, "Y", "N")), Ss)

## -----------------------------------------------------------------------------
#  glmfit = glm(Z~., data=swdata113, family=binomial)
#  summary(glmfit)

## -----------------------------------------------------------------------------
#  test = predict(glmfit, type="response")>=0.5

## -----------------------------------------------------------------------------
#  (err.table = table(swdata113$Z, test))

## -----------------------------------------------------------------------------
#  1 - sum(diag(err.table)) / sum(err.table)

## -----------------------------------------------------------------------------
#  t(apply(err.table, 1, function(x) x / sum(x)))

## -----------------------------------------------------------------------------
#  matrix(
#  c("specificity", "false positive rate",
#  "false negative rate","sensitivity")
#  , ncol=2)

## ----fig.width = 5, fig.height = 3.5------------------------------------------
#  x = seq(min(swdata113$c1), max(swdata113$c1), length = 30)
#  y = seq(min(swdata113$c2), max(swdata113$c2),
#  length = length(x))
#  prob = function(x, y) 1/(1+exp(-predict(glmfit,
#  newdata=data.frame(c1=x, c2=y))))
#  z = outer(x, y, prob)
#  filled.contour(x,y,z, xlab="Component 1", ylab="Component 2")

## ---- message=FALSE-----------------------------------------------------------
#  library(e1071)

## -----------------------------------------------------------------------------
#  set.seed(1)
#  tuneSVM = tune(svm, Z~., data=swdata113,
#  ranges = list(gamma = 2^(0:2), cost = c(4, 6, 8)),
#  tunecontrol = tune.control(cross = nrow(swdata113)))

## -----------------------------------------------------------------------------
#  summary(tuneSVM)

## ----fig.width = 5, fig.height = 3.5------------------------------------------
#  plot(tuneSVM, color.palette = heat.colors)

## -----------------------------------------------------------------------------
#  bestGamma = tuneSVM$best.parameters$gamma
#  bestC = tuneSVM$best.parameters$cost

## -----------------------------------------------------------------------------
#  set.seed(1)
#  svmfit = svm(Z~., data=swdata113,
#  cost = bestC, gamma = bestGamma,
#  probability=TRUE, kernel="radial", cross=nrow(swdata113))
#  summary(svmfit)

## -----------------------------------------------------------------------------
#  pred = predict(svmfit, newdata=swdata113, probability=TRUE,
#  decision.values=TRUE)
#  (err.table = table(swdata113$Z, pred))
#  1 - sum(diag(err.table)) / sum(err.table)

## ----fig.width = 5, fig.height = 3.5------------------------------------------
#  plot(svmfit, swdata113, c2~c1)

## -----------------------------------------------------------------------------
#  opt11 = optparasearch(SB1, Z=Z, comp=20,
#  search.method = "regparaonly", criterion="BIC")
#  (fit311 = msma(SB1, Z=Z,
#  comp=opt11$optncomp, lambdaX=opt11$optlambdaX))
#  Ss = fit311$ssX
#  colnames(Ss) = paste("c", c(1:ncol(Ss)), sep="")
#  swdata311 = data.frame(
#  Z = as.factor(ifelse(Z == 1, "Y", "N")), Ss)

## -----------------------------------------------------------------------------
#  if(!require("rpart.plot")) install.packages("rpart.plot")
#  library(rpart)
#  library(rpart.plot)

## -----------------------------------------------------------------------------
#  set.seed(1)
#  (treefit = rpart(Z~., data=swdata311,
#  control = rpart.control(minsplit = 4)))

## -----------------------------------------------------------------------------
#  prp(treefit, type=4, extra=1, faclen=0, nn=TRUE)

## -----------------------------------------------------------------------------
#  printcp(treefit)

## -----------------------------------------------------------------------------
#  plotcp(treefit)

## -----------------------------------------------------------------------------
#  (treefit1 = prune(treefit, cp=0.05))

## -----------------------------------------------------------------------------
#  (treefit2 = snip.rpart(treefit, 7))

## -----------------------------------------------------------------------------
#  prp(treefit2, type=4, extra=1, faclen=0, nn=TRUE)

## -----------------------------------------------------------------------------
#  pred = predict(treefit2, type="class")

## -----------------------------------------------------------------------------
#  (err.table = table(swdata311$Z, pred))

## -----------------------------------------------------------------------------
#  1 - sum(diag(err.table))/sum(err.table)

## -----------------------------------------------------------------------------
#  swdata3112 = head(swdata311)

## -----------------------------------------------------------------------------
#  set.seed(1)
#  (idrand = sample(1:6, replace=TRUE))

## -----------------------------------------------------------------------------
#  (bsample = swdata3112[idrand, 1:4])

## -----------------------------------------------------------------------------
#  (oobsample = swdata3112[!(1:nrow(swdata3112) %in%
#  unique(idrand)), 1:4])

## -----------------------------------------------------------------------------
#  library(randomForest)
#  library(e1071)

## -----------------------------------------------------------------------------
#  set.seed(1)
#  tuneRF = tune(randomForest, Z~., data=swdata311,
#  ranges = list(mtry = c(4,6,8), ntree = c(300, 500, 1000),
#  nodesize= c(1,2,3)),
#  tunecontrol = tune.control(cross = nrow(swdata311)))

## -----------------------------------------------------------------------------
#  summary(tuneRF)

## -----------------------------------------------------------------------------
#  bestmtry = tuneRF$best.parameters$mtry
#  bestntree = tuneRF$best.parameters$ntree
#  bestnodesize = tuneRF$best.parameters$nodesize

## -----------------------------------------------------------------------------
#  set.seed(1)
#  (rffit = randomForest(Z~., data=swdata311, proximity=TRUE,
#  mtry = bestmtry, ntree = bestntree, nodesize=bestnodesize))

## ----fig.width = 5, fig.height = 5--------------------------------------------
#  varImpPlot(rffit)

## -----------------------------------------------------------------------------
#  Q = fit311$wbX[[1]][,6]
#  outstat1 = rec(Q, img1$imagedim, B=B1, mask=img1$brainpos)
#  outstat2 = -outstat1
#  coat(template, outstat2)

## -----------------------------------------------------------------------------
#  atlastable(tmpatlas, outstat2, atlasdataset)

## -----------------------------------------------------------------------------
#  pred = predict(rffit, type="class")

## -----------------------------------------------------------------------------
#  (err.table = table(swdata311$Z, pred))

## -----------------------------------------------------------------------------
#  1 - sum(diag(err.table))/sum(err.table)

## -----------------------------------------------------------------------------
#  t(apply(err.table, 1, function(x) x / sum(x)))

## ----fig.width = 6, fig.height = 3.5------------------------------------------
#  par(mfrow=c(1,2))
#  partialPlot(rffit, swdata311, c1)
#  partialPlot(rffit, swdata311, c2)

## ----fig.width = 4.5, fig.height = 3------------------------------------------
#  par(mfrow=c(1,1), mar=c(3,3,3,8))
#  z = rffit$proximity
#  n = nrow(swdata311)
#  filled.contour(x=1:n, y=1:n, z=z, color = terrain.colors)

## ----fig.width = 4.5, fig.height = 3------------------------------------------
#  par(mfrow=c(1,1), mar=c(4,3,2,2))
#  MDSplot(rffit, factor(swdata311$Z),
#  pch=as.numeric(swdata311$Z)-1)

## ----fig.width = 3.5, fig.height = 3------------------------------------------
#  par(mfrow=c(1,1), mar=c(4,3,2,2))
#  plot(randomForest::outlier(rffit), type="h",
#  col= as.numeric(swdata311$Z))

## -----------------------------------------------------------------------------
#  img2 = simbrain(baseimg = baseimg, diffimg = diffimg2,
#  sdevimg=sdevimg, mask=mask, n0=500, c1=0.01, sd1=0.1,
#  zeromask=FALSE, seed=2)

## -----------------------------------------------------------------------------
#  testZ = as.factor(ifelse(img2$Z == 1, "Y", "N"))

## -----------------------------------------------------------------------------
#  SB2 = basisprod(img2$S, B1)

## -----------------------------------------------------------------------------
#  ptest113 = ptest(object=fit113, Z=Z, newdata=SB2,
#  testZ=testZ, regmethod = "glm")

## -----------------------------------------------------------------------------
#  summary(ptest113$trainout$finalModel)

## -----------------------------------------------------------------------------
#  ptest113$trainout

## -----------------------------------------------------------------------------
#  ptest113$predcnfmat

## -----------------------------------------------------------------------------
#  ptest112 = ptest(object=fit112, Z=Z, newdata=SB2, testZ=testZ,
#  regmethod = "glm")
#  ptest112$predcnfmat$overall["Accuracy"]

## -----------------------------------------------------------------------------
#  ptest311 = ptest(object=img1$S, Z=Z, newdata=img2$S,
#  testZ=testZ, regmethod = "glm")
#  ptest312 = ptest(object=SB1, Z=Z, newdata=SB2, testZ=testZ,
#  regmethod = "glm")

## -----------------------------------------------------------------------------
#  ptest311$predcnfmat$overall["Accuracy"]
#  ptest312$predcnfmat$overall["Accuracy"]

## -----------------------------------------------------------------------------
#  mus = seq(0, 1, by=0.25)
#  comps = c(1,2,5,10)
#  
#  paramtest1 = lapply(comps, function(c1){lapply(mus,
#  function(mu1){
#  tmpfit = msma(SB1, Z=Z, comp=c1, lambdaX=0.075, muX=mu1)
#  tmpptest = ptest(object=tmpfit, Z=Z, newdata=SB2,
#  testZ=testZ, regmethod = "glm")
#  tmpptest$predcnfmat
#  })})
#  
#  out1 = do.call(rbind, lapply(paramtest1, function(x)
#  do.call(cbind, lapply(x,
#  function(y)y$overall["Accuracy"]))))
#  rownames(out1)=comps; colnames(out1)=mus
#  

## ----results='asis'-----------------------------------------------------------
#  kable(out1, "latex", booktabs = T)

## -----------------------------------------------------------------------------
#  svmgrid = expand.grid(sigma = c(0.001, 0.025, 0.05),
#  C = c(0.5, 0.75, 1))
#  ptest211 = ptest(object=fit311, Z=Z, newdata=SB2, testZ=testZ,
#  regmethod = "svmRadial", metric="ROC",
#  param=svmgrid)
#  ptest211$trainout

## -----------------------------------------------------------------------------
#  ptest211$predcnfmat

## -----------------------------------------------------------------------------
#  treegrid = NULL
#  ptest212 = ptest(object=fit311, Z=Z, newdata=SB2, testZ=testZ,
#  regmethod = "rpart", metric="ROC",
#  param=treegrid)
#  ptest212$trainout

## -----------------------------------------------------------------------------
#  ptest212$predcnfmat

## -----------------------------------------------------------------------------
#  rfgrid  =  NULL
#  ptest213 = ptest(object=fit311, Z=Z, newdata=SB2, testZ=testZ,
#  regmethod = "rf", metric="ROC",
#  param=rfgrid)
#  ptest213$trainout

## -----------------------------------------------------------------------------
#  ptest213$predcnfmat

## -----------------------------------------------------------------------------
#  layers0 = c(1, 5, 10); layers1 = c(0, 1, 5, 10)
#  rate0 = c(0, 0.25, 0.5, 0.75)
#  activation=c("relu", "sigmoid", "tanh", "softrelu")
#  mxnet.params = expand.grid(layer1=layers0, layer2=layers1,
#  layer3=0, learning.rate=0.1, momentum=0.9, dropout=0,
#  activation=activation[3])

## ----eval=FALSE---------------------------------------------------------------
#  ptest215 = ptest(object=fit113, Z=Z, newdata=SB2, testZ=testZ,
#  regmethod = "mxnet", metric="Accuracy",
#  param=mxnet.params)
#  ptest215$trainout

