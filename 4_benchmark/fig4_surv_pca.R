if(!require('tidyverse')) {install.packages('tidyverse')} 
if(!require('glmnet')) {install.packages('glmnet')} #
if(!require('survC1')) {install.packages('survC1')} #
if(!require('survAUC')) {install.packages('survAUC')}
if(!require('survival')) {install.packages('survival')} 
if(!require('survminer')) {install.packages('survminer')} 
if(!require('survivalROC')) {install.packages('survivalROC')} 

# Plotting survival curves


load("../data/1.data.processed.rda")
load(file = "../data/2.gsign_kidney.rda")

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
                          Boruta = "fs.boruta",
                          RFE = "fs.mlr.rfe",
                          mRMR = "fs.mRMR"
  ))

rownames(gsign_kidney) <- gsign_kidney$id


data.filt <- data[data$obs.time > 20 & data$obs.time < 2555, ]
data_tcga <- data[data$dataset =="TCGA", ]
data_icgc <- data[data$dataset =="ICGC", ]
data_icgc.filt <- data.filt[data.filt$dataset =="ICGC", ]

col_sign <- c("age", "AR", "AL353637.1", "FOXJ1", "LINC01732", "SEMA3G",  "HHLA2")
col_sign <- c("AR", "FOXJ1", "LINC01732", "SEMA3G",  "MT1G")
set.seed(12) ## for reproducibility

fold1 <- sample(1:nrow(data_tcga),size=floor(nrow(data_tcga)/3))
fold2 <- sample(setdiff(1:nrow(data_tcga), fold1), size=floor(nrow(data_tcga)/3))
fold3 <- setdiff(1:nrow(data_tcga), union(fold1,fold2))

x <- (data_tcga[, col_sign] +1)
y <- Surv(data_tcga$obs.time, data_tcga$status)
y_tr <- y[-fold1]
x_tr <- x[-fold1,]
y_te <- y[fold1]
x_te <- x[fold1,]

set.seed(12) 
cvfit_tr <- cv.glmnet(data.matrix(x_tr),y_tr, nfolds=5, gamma = 1, relax = T, family="cox", type.measure="C")
preds <- predict(cvfit_tr,data.matrix(x_te), s="lambda.min")
levs_tcga <- cut_number(preds,3)
levs_tcga <- as.character(as.numeric(levs_tcga))
x_te$risk <- fct_relevel(fct_recode(levs_tcga, Low = "1", Moderate = "2", High = "3"), "High", "Moderate", "Low")

fit <- survfit(y_te~levs_tcga)
out <- survdiff(y_te~levs_tcga)
p.val <- 1 - pchisq(out$chisq, length(out$n) - 1)

times <- seq(10, 4000, 10)   # 2555               
AUC_Uno <- AUC.uno(y_tr, y_te, preds, times)
auc_uno <- AUC_Uno$iauc

km_data <- surv_fit(y_te ~ risk, data = cbind(x_te,y_te))


#makes the scientific notation using "AeB" explicitly write out Ax10^B
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e(.*)$", "\\1e\\2", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", " \U00D7 10^", l)
  # return this as a string
  l
}
library(ggtext)
library(glue)
psurv_mRMR <-  ggsurvplot(km_data, conf.int = FALSE, data =  cbind(x_te,y_te), fun = NULL,
                          #fontsize = 11,
                          ggtheme = theme_grey(), 
                          palette = c("red3", "orange1", "green4"),
                          title = paste0("Survival prediction with mRMR genes on TCGA-KIRC\nTraining with 66% / Testing with 33% / AUC-Uno: ", signif(auc_uno, 2)),
                          
                          risk.table.y.text = FALSE,
                          risk.table.title = "Nº at risk (nº of deceased)",
                          tables.height = 0.25,
                          risk.table = "nrisk_cumevents",
                          pval = format(2^31-1, scientific =T))
                          #"p-val \U00D7 10^{-05}")
                          #pval = bquote('p-val:' 1.171723 \U00D7 10^-05')
  # bquote('x axis'~(Å^2))
psurv_mRMR
                            #paste0("italic(p-val:)~1.171723 × 10^(-05)", fancy_scientific(p.val)) )
                          #pval = paste0("p-val: ", signif(p.val, 3)) )

# Compute PCA  -----
if(!require('factoextra')) {install.packages('factoextra')} 
if(!require('FactoMineR')) {install.packages('FactoMineR')} 

set.seed(12) 
res.pca <- PCA(data_tcga[fold1, col_sign], graph = FALSE)

status <- as.factor(ifelse(data_tcga$status[fold1]==1, "Deceased", "Alive"))
risk <- fct_relevel(fct_recode(levs_tcga, Low = "1", Moderate = "2", High = "3"), "High", "Moderate", "Low")

tb <- as.data.frame.matrix(table(risk,status))
data.tb <- tibble(x = Inf, y = -Inf, 
                  tb = list(tb))


tt <- ttheme_default(core=list(fg_params=list(hjust=0, x=0.1)),
                     rowhead=list(fg_params=list(col="orange", fontface=3L)),
                     rowhead=list(fg_params=list(hjust=0, x=0)))

p_pca <- fviz_pca_ind(res.pca,
                      label = "none", # hide individual labels
                      habillage = risk, # color by groups
                      #habillage = status, # color by groups
                      show.clust.cent = FALSE,
                      palette = c("red3", "orange1", "green4"),
                      addEllipses = TRUE, 
                      ellipse.level = 0.95,
                      ellipse.type ="confidence",
                      ggtheme = theme_grey(),
                      title = "PCA with mRMR genes on TCGA-KIRC predicted samples\n") + 
  geom_rug(sides="b", aes(color = risk)) +
  geom_table(data = data.tb, aes(x, y, label = tb), 
             table.theme = ttheme_gtlight,
             table.rownames = T, 
             size = 4,
             show.legend = T,
             hjust = 1, vjust = -0.25
  ) + 
  theme(legend.position="top")



p_mRMR <- plot_grid(psurv_mRMR$plot, psurv_mRMR$table,  ncol = 1, rel_heights = c(1, 0.4))

pcow <- cowplot::plot_grid(
  p_mRMR, NULL,  p_pca,
  ncol = 3,
  rel_widths = c(0.8,0.1, 1),
  labels = c("a","", "b")
)

# pcow



splots <- list()
splots_auc <- list()
splots_gpca <- list()

i <- 1
sign <- "mRMR"

col_sign <- c("age", "AL353637.1", "DPP6", "FOXJ1", "GNB3", "HHLA2", "IL4", "LIMCH1", "OTX1", "SAA1", "ZIC2")
col_sign <- c("AR", "FOXJ1", "LINC01732", "SEMA3G",  "MT1G")

set.seed(12) ## for reproducibility

fold1 <- c(1:nrow(data))[data$dataset == "ICGC"]

x <- (data[, col_sign] +1)
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
times <- seq(10, 2000, 10)   # 4000               
AUC_Uno <- AUC.uno(y_tr, y_te, preds, times)
auc_uno <- AUC_Uno$iauc
auc_uno  
km_data <- surv_fit(y_te ~ risk, data = cbind(x_te,y_te))

psurv_mRMR_icgc  <-  ggsurvplot(km_data, conf.int = FALSE, data =  x_te, fun = NULL,
                                # fontsize = 12,
                                ggtheme = theme_grey(), 
                                palette = c("red3", "orange1", "green4"),
                                #title = "Survival prediction with mRMR genes\n(train TCGA-KIRC / test ICGC-RECA / auc_uno: round(auc_uno, digits = 2))",
                                title = paste0("Survival prediction with mRMR genes\nTraining with TCGA-KIRC / Testing with ICGC-RECA / AUC-Uno: ", round(auc_uno, digits = 2)),
                                #title = paste0("Signature ", sign ," / p-val: ", signif(p.val, 3), " / auc_uno: ", round(auc_uno, digits = 2)),
                                risk.table.y.text = FALSE,
                                risk.table.title = "Nº at risk (nº of deceased)",
                                tables.height = 0.25,
                                risk.table = "nrisk_cumevents",
                                pval = paste0("p-val: ", signif(p.val, 3)))


# ##AUC
survivalROC_helper <- function(t) {
  
  survivalROC(Stime    = data$obs.time[fold1],
              status   = data$status[fold1],
              marker   = preds,
              predict.time = t,
              method       = "NNE",
              span = 0.25 * length(fold1)^(-0.20))
}
labels_obs.time <- c("2555" = "7 years", "3650" = "10 years", "4000" = "12 years" )
## Evaluate every 3650 days
survivalROC_data <- tibble(t = c(2555)) %>%
  dplyr::mutate(survivalROC = purrr::map(t, survivalROC_helper),
                ## Extract scalar AUC
                auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
                ## Put cut off dependent values in a data_frame
                df_survivalROC = purrr::map(survivalROC, function(obj) {
                  as_tibble(obj[c("cut.values","TP","FP")])
                })) %>%
  dplyr::select(-survivalROC) %>%
  unnest(cols = c(df_survivalROC)) %>%
  arrange(t, FP, TP)

## Plot AUC
splots_auc[[i]] <- survivalROC_data %>%
  ggplot(mapping = aes(x = FP, y = TP)) +
  ggtitle(paste0("AUC of ",  sign, " (train TCGA-KIRC / test ICGC-RECA)")) +
  geom_point() +
  geom_line() +
  geom_label(data = survivalROC_data %>% dplyr::select(t,auc) %>% unique,
             mapping = aes(label = sprintf("%.2f", auc)), x = 0.5, y = 0.5) +
  facet_wrap( ~ t, labeller = as_labeller(labels_obs.time),  ncol = 4) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.key = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank())

# Compute PCA  -----
library(factoextra)
library(FactoMineR)


res.pca_icgc <- PCA(data[fold1, col_sign], graph = FALSE)

status <- as.factor(ifelse(data$status[fold1]==1, "Deceased", "Alive"))
risk <- fct_relevel(fct_recode(levs, Low = "1", Moderate = "2", High = "3"), "High", "Moderate", "Low")

tb <- as.data.frame.matrix(table(risk,status))
data.tb <- tibble(x = Inf, y = -Inf, 
                  tb = list(tb))


tt <- ttheme_default(core=list(fg_params=list(hjust=0, x=0.1)),
                     rowhead=list(fg_params=list(col="orange", fontface=3L)),
                     rowhead=list(fg_params=list(hjust=0, x=0)))

p_pca_icgc <- fviz_pca_ind(res.pca_icgc,
                           label = "none", # hide individual labels
                           habillage = risk, # color by groups
                           #habillage = status, # color by groups
                           show.clust.cent = FALSE,
                           palette = c("red3", "orange1", "green4"),
                           addEllipses = TRUE, 
                           ellipse.level = 0.95,
                           ellipse.type ="confidence",
                           ggtheme = theme_grey(),
                           title = "PCA with mRMR genes on ICGC-RECA\n") + 
  geom_rug(sides="b", aes(color = risk)) +
  geom_table(data = data.tb, aes(x, y, label = tb), 
             table.theme = ttheme_gtlight,
             table.rownames = T, 
             size = 4,
             show.legend = T,
             hjust = 1, vjust = -0.25
  ) + 
  theme(legend.position="top")



p_mRMR_icgc <- plot_grid(psurv_mRMR_icgc$plot, psurv_mRMR_icgc$table,  ncol = 1, rel_heights = c(1, 0.4))



pcow_icgc <- cowplot::plot_grid(
  p_mRMR_icgc, NULL,  p_pca_icgc,
  ncol = 3,
  rel_widths = c(0.8,0.1, 1),
  labels = c("c","", "d")
)



pdf("../plots/fig4_surv_pca.pdf", width = 16, height = 11)
cowplot::plot_grid(
  pcow,
  pcow_icgc,
  ncol = 1
)
dev.off()

```

