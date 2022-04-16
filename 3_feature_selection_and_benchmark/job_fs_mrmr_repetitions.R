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
col_sign <-  c("age", "AL353637.1", "DPP6", "FOXJ1", "GNB3", "HHLA2", "IL4", "LIMCH1", "OTX1", "SAA1", "ZIC2")
#col_sign <- c("FOXJ1", "LINC01732", "SEMA3G", "AR", "HHLA2", "AL353637.1")

library(data.table)

df_results <- as.data.frame(matrix(nrow = 0, ncol = 5))

names(df_results) <- c("iter","genes", "pval", "auc", "ci")

for(i in c(1:100)){
    set.seed(i)  ## for random values
    
    fold1 <- c(1:nrow(data))[data$dataset == "ICGC"]
    
    x <- (data[, names(data) %in% col_sign] +1)
    y <- Surv(data$obs.time, data$status)
    y_tr <- y[-fold1]
    x_tr <- x[-fold1,]
    y_te <- y[fold1]
    x_te <- x[fold1,]
    
    set.seed(i)  ## for random values
    cvfit_tr <- cv.glmnet(data.matrix(x_tr),y_tr, nfolds=5, gamma = 1, relax = T, family="cox", type.measure="C")
    preds <- predict(cvfit_tr,data.matrix(x_te), s="lambda.min")
    levs <- cut_number(preds,3)
    levs <- as.character(as.numeric(levs))
    x_te$risk <- fct_relevel(fct_recode(levs, Low = "1", Moderate = "2", High = "3"), "High", "Moderate", "Low")
    
    fit <- survfit(y_te~levs)
    out <- survdiff(y_te~levs)
    p.val <- 1 - pchisq(out$chisq, length(out$n) - 1)
    times <- seq(10, 2000, 10)         
    AUC_Uno <- AUC.uno(y_tr, y_te, preds, times)
    auc_uno <- AUC_Uno$iauc
    
    mydata <- data.frame(as.matrix(y_te),preds)
    out <- Est.Cval(mydata, 2000, nofit=TRUE)
    cind <- out$Dhat
    df_results <- rbind(df_results, list(i, paste0(col_sign, collapse = ", "), p.val, auc_uno, cind))
    print(paste0(i," - ", "pval:", round(p.val, digits = 5) ," / auc: ",round(auc_uno, digits = 3)," / ci: ",  round(cind, digits = 3), " / genes: ", paste0(col_sign, collapse = ", ")) )
}

names(df_results) <- c("iter","genes", "pval", "auc", "ci")

save(df_results, file = "genes.fs.mrmr.rep.rda", compress = T)
