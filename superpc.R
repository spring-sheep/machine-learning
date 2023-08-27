rm(list = ls())
library(stringr)
library(survival)
library(superpc)
library(dplyr)
library(tibble)
library(BART)

#superpc----

data <- list(x=t(est_dd[,-c(1,2)]),y=est_dd$OS.time,censoring.status=est_dd$OS,featurenames=colnames(est_dd)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5)
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs <- lapply(val_dd_list,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$OS.time,censoring.status=w$OS,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})