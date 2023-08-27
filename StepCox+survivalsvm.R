rm(list = ls())
library(stringr)
library(survival)
library(CoxBoost)
library(dplyr)
library(tibble)
library(BART)

#StepCox+survivalsvm----

for (direction in c("both", "backward","forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd),direction = direction)
  rid <- names(coef(fit))
  rid <- str_replace_all(rid,"`","")
  est_dd2 <- est_data[,c('OS.time','OS',rid)]
  val_dd_list2 <- lapply(val_data_list,function(x){x[,c('OS.time','OS',rid)]})
  
  fit = survivalsvm(Surv(OS.time,OS)~., data= est_dd2, gamma.mu = 1)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + survival-SVM')
  result <- rbind(result,cc)
}