if(!require('tidyverse')) {install.packages('tidyverse')} #
if(!require('splitstackshape')) {install.packages('splitstackshape')} #
if(!require('glmnet')) {install.packages('glmnet')} #
if(!require('survC1')) {install.packages('survC1')} #
if(!require('survAUC')) {install.packages('survAUC')}
if(!require('survival')) {install.packages('survival')} 
if(!require('survminer')) {install.packages('survminer')} 
if(!require('survivalROC')) {install.packages('survivalROC')} 
if(!require('ggpmisc')) {install.packages('ggpmisc')} 
if(!require('glmpca')) {install.packages('glmpca')} 
if(!require('glmnet')) {install.packages('glmnet')} 
if(!require('factoextra')) {install.packages('factoextra')}
if(!require('FactoMineR')) {install.packages('FactoMineR')} 

setwd(".")
load(file = "../data/2.gsign_kidney.rda")
load("../data/1.data.processed.rda")
#load("../data/2.data.full.count.rda")
col_clini <- c("obs.time", "status", "age", "gender", "neoplasm", "metastasis", "ajcc.stage", "dataset")
col_genes <- setdiff(colnames(data), col_clini)


data_tcga <- data[data$dataset =="TCGA", ]
data_icgc <- data[data$dataset =="ICGC", ]

gsign_kidney <- gsign_kidney %>%
  dplyr::filter(!id %in% c("sign.3", "sign.4",  "sign.11", "sign.13", "sign.14", "fs.genetic", "fs.mRMR.status", "fs.mRMR.time")) %>%
  mutate(id =  fct_recode(id,
                          Ha.2014.CaInfo = "sign.1",
                          Zhan.2015.CMMM = "sign.2",
                          Dai.2016.Oncotarget = "sign.5",
                          Wan.2017.IJCancer = "sign.6",
                          Chang.2018.Medicine = "sign.7",
                          YuanLeiChen.2018.JCP = "sign.8",
                          Hu.2019.IJMS = "sign.9",
                          LiangChen.2018.JCP = "sign.10",
                          #Pan.2019.MSM.5genes = "sign.11", # removed
                          Pan.2019.MSM = "sign.12",
                          Wu.2019.FrontiersOnco = "sign.15",
                          Jiang.2020.ACS = "sign.16",
                          Zou.2020.PeerJ = "sign.17",
                          LingChen.2020.Hereditas = "sign.18",
                          DCosta.2020.SciReports = "sign.19",
                          wrap.boruta = "fs.boruta",
                          wrap.rfe = "fs.mlr.rfe",
                          wrap.mRMR = "fs.mRMR"
  ))



rownames(gsign_kidney) <- gsign_kidney$id

gsign_kidney["wrap.mRMR", "signature"][[1]] <- list(c("age", "AL353637.1", "DPP6", "FOXJ1", "GNB3", "HHLA2", "IL4", "LIMCH1", "OTX1", "SAA1", "ZIC2"))
                                                    
gsign_kidney["mRMR.search", "signature"][[1]] <- list(c("AR", "AL353637.1", "FOXJ1", "LINC01732", "SEMA3G", "HHLA2"))


bmk_icgc <- setNames(data.frame(matrix(nrow=0, ncol=6)), c("sign", "iteration", "p.val", "cind", "auc_uno"))
for(i in c(1:100)){
for(sign in rownames(gsign_kidney)) {
  #print(sign)
  col_sign <- intersect(gsign_kidney[sign, "signature"][[1]], names(data_tcga))
    print(paste0(i, " - ", sign))
    #set.seed(1) ## for reproducibility
    set.seed(i)
    ix <- c(1:nrow(data))[data$dataset == "TCGA"]
    
    x <- (data[, col_sign] +1)
    y <- Surv(data$obs.time, data$status)
    y_tr <- y[ix]
    x_tr <- x[ix,]
    y_te <- y[-ix]
    x_te <- x[-ix,]
    
    # nFolds_inner <- 3
    # foldid_inner <- rep(0, nrow(x_tr))
    # dat <- data[ix, c("obs.time","status")] %>%  dplyr::mutate(id = c(1:nrow(x_tr)))
    # for(k in 1:nFolds_inner){
    #   j=nFolds_inner-(k-1)
    #   if(j>1){
    #     a=stratified(dat, group = "status", size = 1/j)
    #     foldid_inner[a$id] <- k
    #     dat=dat%>%filter(id %in% setdiff(dat$id,a$id))
    #   } else{
    #     foldid_inner[dat$id] <- k
    #   }
    # }
    # 
    #                  
    # cvfit_tr <- cv.glmnet(data.matrix(x_tr),
    #                       y_tr, 
    #                       nfolds=nFolds_inner,
    #                       foldid = foldid_inner,
    #                       standardize = FALSE,
    #                       gamma = 1, relax = T, family="cox", type.measure="C")
    cvfit_tr <- cv.glmnet(data.matrix(x_tr),y_tr, nfolds=5, gamma = 1, relax = T, family="cox", type.measure="C")
    
    preds <- predict(cvfit_tr,data.matrix(x_te), s="lambda.min")
    
    levs <- cut_number(preds,3)
    
    fit <- survfit(y_te~levs)
    out <- survdiff(y_te~levs)
    p.val <- 1 - pchisq(out$chisq, length(out$n) - 1)
    print(paste0("p-val: ", p.val))
    
    
    mydata <- data.frame(as.matrix(y_te),preds)
    out <- Est.Cval(mydata, 2555, nofit=TRUE)
    cind <- out$Dhat
    #print(paste0("C-index: ", cind))
    
    times <- seq(20, 2555, 10)                  
    AUC_Uno <- AUC.uno(y_tr, y_te, preds, times)
    auc_uno <- AUC_Uno$iauc
    print(paste0("AUC_Uno: ", auc_uno)) 
    df <- setNames(data.frame(sign, i, p.val, cind, auc_uno), c("sign", "iteration", "p.val", "cind", "auc_uno"))
    
    bmk_icgc <- rbind(bmk_icgc,  df)
  }
} 
bmk_icgc$padj <- 1

for(sign in unique(bmk_icgc$sign)) {
  bmk_icgc[bmk_icgc$sign == sign, "padj"] <- p.adjust(bmk_icgc[bmk_icgc$sign == sign, "p.val"], 
                                                      n = length(bmk_icgc[bmk_icgc$sign == sign,"p.val"]),
                                                      method = "fdr")
}

save(bmk_icgc, file="validation_tcga_icgc_2555.mrmr0.rda")
