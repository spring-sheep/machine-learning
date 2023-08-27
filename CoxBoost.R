rm(list = ls())
library(stringr)
library(survival)
library(CoxBoost)
library(dplyr)
library(tibble)
library(BART)

#CoxBoost----

set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost')
result <- rbind(result,cc)