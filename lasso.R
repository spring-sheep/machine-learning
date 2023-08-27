rm(list = ls())
library(stringr)
library(survival)
library(dplyr)
library(tibble)
library(BART)

#lasso----
library(glmnet)
set.seed(6)
x = as.matrix(est_dd[,-c(1,2)])
y = est_dd$OS 
cvfit = cv.glmnet(x, y)
fit <- glmnet(x=x, y=y,lambda=cvfit$lambda.min)

rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, newx= as.matrix(x[,-c(1,2)]), s=cvfit$lambda.min)))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso')
result <- rbind(result,cc)