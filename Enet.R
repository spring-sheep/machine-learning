rm(list = ls())
library(stringr)
library(survival)
library(glmnet)
library(dplyr)
library(tibble)
library(BART)

#Enet----

x1 <- as.matrix(est_dd[,var])
x2 <- as.matrix(Surv(est_dd$OS.time,est_dd$OS))

for (alpha in seq(0,1,0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
  rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),
                                                                           s=fit$lambda.min)))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('Enet','[α=',alpha,']')
  result <- rbind(result,cc)
}