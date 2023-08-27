rm(list = ls())
library(stringr)
library(survival)
library(randomForestSRC)
library(dplyr)
library(tibble)
library(BART)

#RSF----
set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd,
             ntree = 1000,nodesize = rf_nodesize,
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- 'RSF'
result <- rbind(result,cc)