if(!require('glmnet')) {install.packages('glmnet')} #
if(!require('survC1')) {install.packages('survC1')} #
if(!require('survAUC')) {install.packages('survAUC')}
if(!require('tidyverse')) {install.packages('tidyverse')} #
if(!require('survival')) {install.packages('survival')} #

setwd("~/gene_signature/4_gene_selection")

load("../data/1.data.processed.rda")
#data_icgc <- data[data$dataset =="ICGC", ]


cols <- setdiff(colnames(data), c("neoplasm", "metastasis", "ajcc.stage", "dataset"))

#data_icgc <- data_icgc[, cols]

col_sign <-  c("LIMCH1", "MT1G", "FOXJ1", "LINC01732", "SEMA3G", "AR", "HHLA2", "SAA1", "ZIC2", "AL353637.1")
#col_sign <-  c("AL353637.1", "AR", "FOXJ1", "HHLA2", "LINC01732", "MT1G", "SEMA3G", "ZIC2")


library(data.table)
comb_sign <- as.data.frame(do.call(CJ, replicate(10, 0:1, FALSE)))
names(comb_sign) <- col_sign
comb_sign <- comb_sign[-1,]  
comb_sign <- (comb_sign==1)  

df_results <- as.data.frame(matrix(nrow = nrow(comb_sign), ncol = 4))
names(df_results) <- c("genes", "pval", "auc", "ci")

for(i in c(1:nrow(comb_sign))){
  cols <- names(comb_sign[i,])[comb_sign[i,]]
  if(length(cols) > 5){
  #print(col_sign)
  set.seed(12) ## for reproducibility
  
  fold1 <- c(1:nrow(data))[data$dataset == "ICGC"]
  
  x <- (data[, names(data) %in% cols] +1)
  y <- Surv(data$obs.time, data$status)
  y_tr <- y[-fold1]
  x_tr <- x[-fold1,]
  y_te <- y[fold1]
  x_te <- x[fold1,]
  
  set.seed(12) 
  cvfit_tr <- cv.glmnet(data.matrix(x_tr),y_tr, nfolds=5, gamma = 1, relax = T, family="cox", type.measure="C")
  preds <- predict(cvfit_tr,data.matrix(x_te), s="lambda.min")
  levs <- cut_number(preds,3)
  levs <- as.character(as.numeric(levs))
  x_te$risk <- fct_relevel(fct_recode(levs, Low = "1", Moderate = "2", High = "3"), "High", "Moderate", "Low")
  
  fit <- survfit(y_te~levs)
  out <- survdiff(y_te~levs)
  p.val <- 1 - pchisq(out$chisq, length(out$n) - 1)
  p.val
  times <- seq(10, 2300, 10)   # 4000               
  AUC_Uno <- AUC.uno(y_tr, y_te, preds, times)
  auc_uno <- AUC_Uno$iauc
  
  mydata <- data.frame(as.matrix(y_te),preds)
  out <- Est.Cval(mydata, 2300, nofit=TRUE)
  cind <- out$Dhat
  df_results[i, ] <- list(paste0(cols, collapse = ", "), p.val, auc_uno, cind)
  print(paste0(i," - ", "pval:", p.val ," / auc: ",auc_uno," / ci: ",  cind, " / genes: ", paste0(cols, collapse = ", ")) )
  }
}

save(df_results, file = "genes.fs.mrmr.combinations.rda", compress = T)
