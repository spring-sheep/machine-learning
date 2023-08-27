rm(list = ls())
library(stringr)
library(survival)
library(CoxBoost)
library(dplyr)
library(tibble)
library(BART)

#StepCox----

for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd),direction = direction)
  rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']')
  result <- rbind(result,cc)
}
