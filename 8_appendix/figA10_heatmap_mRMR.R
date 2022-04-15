if(!require('tidyverse')) {install.packages('tidyverse')} 
if(!require('pheatmap')) {install.packages('pheatmap')} 

setwd("~/bio/gene_signature/8_appendix")

load(file = "../data/2.gsign_kidney.rda")
load("../data/1.data.processed.rda")

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

sign <- c("AR", "AL353637.1", "DPP6", "FOXJ1", "GNB3", "HHLA2", "IL4", "LIMCH1", "LINC01732", "OTX1", "SAA1", "SEMA3G", "ZIC2")

df0 <- data[rownames(data)[order(data$status)] ,]
df0 <-df0[df0$metastasis !="MX" & !is.na(df0$metastasis), ]
df <- df0[, sign]


ann_heatmap <- data.frame(status = ifelse(df0$status == 1, "deceased", "alive"), 
                          metastasis = df0$metastasis,
                          dataset = df0$dataset, 
                          # age = data$age,
                          ajcc.stage = df0$ajcc.stage)

row.names(ann_heatmap) <- rownames(df)

ann_colors <- list(  
  status = c("#FF1F5B", "#009ADE"), 
  metastasis = c("#F2ACCA", "#06592A"),
  dataset = c("#045275", "#FFF3B2"),
  ajcc.stage = c("#FFF3B2", "#FFC61E", "#FF1F5B", "#009392")
)


names(ann_colors$status) <- c("deceased", "alive")
names(ann_colors$metastasis) <- c("M0", "M1")
names(ann_colors$dataset) <- c("TCGA", "ICGC")
names(ann_colors$ajcc.stage) <- c("T1", "T2", "T3", "T4")

# method 	
#the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).

pdf("../figs/figA10_heatmap_mRMR.pdf", width = 12, height = 9)
pheatmap(df, 
         colorRampPalette(c("navy", "white", "firebrick1"))(100),
         annotation_row = ann_heatmap,
         annotation_colors = ann_colors,
         legend = TRUE, 
         #cluster_rows = FALSE,
         #cluster_cols = FALSE,
         method = "average",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         #kmeans_k = 2,
         show_rownames = FALSE,
         angle_col = 45,
         cutree_cols = 2,
         cutree_rows = 2,
         scale = 'column')
dev.off()
