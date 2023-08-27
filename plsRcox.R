rm(list = ls())
library(stringr)
library(survival)
library(plsRcox)
library(dplyr)
library(tibble)
library(BART)

#plsRcox----

set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=est_dd[,var],time=est_dd$OS.time,status=est_dd$OS),nt=10,verbose = FALSE)
fit <- plsRcox(est_dd[,var],time=est_dd$OS.time,event=est_dd$OS,nt=as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('plsRcox')
result <- rbind(result,cc)