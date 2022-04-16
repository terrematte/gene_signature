# Installing and Loading Libraries            

if(!require('tidyverse')){install.packages('tidyverse')}
if(!require('grid')){install.packages('grid')}
if(!require('gridExtra')) {install.packages('gridExtra')} #  Provides a number of user-level functions to work with "grid" graphics
if(!require('ggrepel')) {install.packages('ggrepel')} 
if(!require('cowplot')) {install.packages('cowplot')} # Provides various features that help with creating publication-quality figures
if(!require('GOplot')) {install.packages('GOplot')}
if(!require('VennDiagram')) {install.packages('VennDiagram')}
if(!require('futile.logger')) {install.packages('futile.logger')} # For VennDiagramLogger

# Load Data

load(file = "../data/2.gsign_kidney.rda")

gsign_kidney["fs.mRMR", "signature"][[1]] <- list(c("AL353637.1", "AR",  "DPP6", "FOXJ1", "FOXJ1", "GNB3", "HHLA2", "IL4", "LIMCH1", "LINC01732", "OTX1", "SAA1", "SEMA3G", "ZIC2"))

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
                          GBM = "filt.gbm",
                          Rpart = "filt.rpart",
                          XGBoost = "filt.xgboost",
                          Boruta = "fs.boruta",
                          RFE = "fs.mlr.rfe",
                          mRMR = "fs.mRMR"
  ))

rownames(gsign_kidney) <- gsign_kidney$id


# Matrix of Intersections

cols <- c("gene", "genes_DEA_M0", "genes_DEA_NT", "genes_eqtls", "genes_papers",  as.character(gsign_kidney$id))
genes_matrix <- as.data.frame(matrix(nrow = length(genes_all), ncol = length(cols)))
colnames(genes_matrix) <- cols
genes_matrix$gene <- genes_all
rownames(genes_matrix) <- genes_all

for(sign_id in gsign_kidney$id) {
  genes_matrix[,sign_id]  <- as.numeric(genes_matrix$gene %in% gsign_kidney[sign_id, "signature"][[1]])
}

gene_list <- list(genes_DEA_NT, genes_DEA_M0,  genes_eqtls, genes_papers) # genes_risk,
names(gene_list) <- c("genes_DEA_NT", "genes_DEA_M0", "genes_eqtls", "genes_papers") # "genes_risk", 

for(l in names(gene_list)) {
  genes_matrix[, l] <- as.numeric(genes_matrix$gene %in% gene_list[[l]]) 
}

rownames(dea.M0.M1) <- dea.M0.M1$symbol
rownames(dea.NT.TP) <- dea.NT.TP$symbol

genes_matrix$logFC.M0.M1 <-  dea.M0.M1[genes_matrix$gene,"logFC"]
genes_matrix$logFC.NT.TP <-  dea.NT.TP[genes_matrix$gene,"logFC"]

genes_matrix$logFC.M0.M1[is.na(genes_matrix$logFC.M0.M1)] <- 0
genes_matrix$logFC.NT.TP[is.na(genes_matrix$logFC.NT.TP)] <- 0

genes_matrix$nlog10FDR.M0.M1 <- -log10(dea.M0.M1[genes_matrix$gene,"FDR"])
genes_matrix$nlog10FDR.NT.TP <- -log10(dea.NT.TP[genes_matrix$gene,"FDR"])

genes_matrix$nlog10FDR.NT.TP[is.na(genes_matrix$nlog10FDR.NT.TP)] <- 0
genes_matrix$nlog10FDR.M0.M1[is.na(genes_matrix$nlog10FDR.M0.M1)] <- 0

# sets_list <- c("Ha.2014.CaInfo", "Zhan.2015.CMMM", "Dai.2016.Oncotarget", "Wan.2017.IJCancer", "Chang.2018.Medicine", 
#           "YuanLeiChen.2018.JCP", "Hu.2019.IJMS", "LiangChen.2018.JCP", "Pan.2019.MSM", "Wu.2019.FrontiersOnco", "Jiang.2020.ACS", 
#           "Zou.2020.PeerJ", "LingChen.2020.Hereditas", "DCosta.2020.SciReports", "filt.gbm", "Rpart", "XGBoost", "Boruta",
#           "RFE", "mRMR", "genes_DEA_M0", "genes_DEA_NT", "genes_eqtls") 

mRMR <- gsign_kidney["mRMR","signature"][[1]]
Boruta <- gsign_kidney["Boruta","signature"][[1]]
RFE <- gsign_kidney["RFE","signature"][[1]]
GBM <- gsign_kidney["GBM","signature"][[1]]
Rpart <- gsign_kidney["Rpart","signature"][[1]]
XGBoost <- gsign_kidney["XGBoost","signature"][[1]]


Pan.2019.MSM <- gsign_kidney["Pan.2019.MSM","signature"][[1]]
Jiang.2020.ACS <- gsign_kidney["Jiang.2020.ACS","signature"][[1]]
Chang.2018.Medicine  <- gsign_kidney["Chang.2018.Medicine","signature"][[1]]
Zou.2020.PeerJ   <- gsign_kidney["Zou.2020.PeerJ","signature"][[1]]
Wan.2017.IJCancer   <- gsign_kidney["Wan.2017.IJCancer","signature"][[1]]
Zhan.2015.CMMM <- gsign_kidney["Zhan.2015.CMMM","signature"][[1]]
Dai.2016.Oncotarget <- gsign_kidney["Dai.2016.Oncotarget","signature"][[1]]
YuanLeiChen.2018.JCP <- gsign_kidney["YuanLeiChen.2018.JCP","signature"][[1]]
Hu.2019.IJMS <- gsign_kidney["Hu.2019.IJMS","signature"][[1]]
LiangChen.2018.JCP <- gsign_kidney["LiangChen.2018.JCP","signature"][[1]]
Wu.2019.FrontiersOnco <- gsign_kidney["Wu.2019.FrontiersOnco","signature"][[1]]
LingChen.2020.Hereditas <- gsign_kidney["LingChen.2020.Hereditas","signature"][[1]]
DCosta.2020.SciReports <- gsign_kidney["DCosta.2020.SciReports","signature"][[1]]

# Chord data

# use this function to make row or column names bold
# parameters:
#   mat: the matrix passed to pheatmap
#   rc_fun: either rownames or colnames
#   rc_names: vector of names that should appear in boldface
make_bold_names <- function(mat, rc_fun, rc_names) {
  bold_names <- rc_fun(mat)
  ids <- rc_names %>% match(rc_fun(mat))
  ids %>%
    walk(
      function(i)
        bold_names[i] <<-
        bquote(bold(.(rc_fun(mat)[i]))) %>%
        as.expression()
    )
  bold_names
}

source("codes/GOChord.R")


# Chord mRMR

sets_list <- c("Rpart", "XGBoost", "Boruta", "RFE", "GBM", "mRMR") 

#Check the set genes intersected with mRMR
d <- genes_matrix
target_cols <- c("mRMR")
other_cols <- colnames(d)[!(colnames(d) %in% c(sets_list, "gene", "genes_papers", "logFC.M0.M1", "logFC.NT.TP", "nlog10FDR.M0.M1", "nlog10FDR.NT.TP"))]
cols_sets <- d %>% filter_at(vars(all_of(target_cols)), ~.==1) %>% select_at(vars(any_of(other_cols))) %>%  colSums(.) > 0
cols_sets <- names(cols_sets)[cols_sets]

circ_dcols <- c("Category","ID","Term","Genes","adj_pval")

circ_data <- as.data.frame(matrix(nrow = length(cols_sets), ncol = length(circ_dcols)))

for(i in c(1:length(cols_sets)) ){
  circ_data[i,] <- list("X", cols_sets[i],cols_sets[i], paste0( get(cols_sets[i]), collapse=", "), 0 )
}

names(circ_data) <- circ_dcols
mRMR <- setdiff(mRMR, c("age"))


circ_genelist <- genes_matrix[mRMR, c("gene", "nlog10FDR.M0.M1")]
rownames(circ_genelist) <- circ_genelist$gene
colnames(circ_genelist) <- c("ID", "logFC")

# Generate the plotting object
circ <- circle_dat(circ_data, circ_genelist)

# Create the plot
chord <- chord_dat(data = circ, genes = circ_genelist, process = circ_data$ID)

chord  <- cbind(chord, genes_matrix[rownames(chord),"logFC.NT.TP"])

p_chord <- GOChordnew(chord,
                      gene.order = 'logFC',
                      nlfc = 2,
                      gene.size = 4.5, 
                      gene.space = 0.22,
                      space = 0.01,  
                      border.size= 0.000001, 
                      process.label = 12,
                      ribbon.col = c('#FF1F5B','#009ADE', '#00CD6C', '#AF58BA' , '#FFC61E', '#F28522', '#A0B1BA'),  #, '#A6761D', DC3977
                      lfc.col=c("#CF597E", "ivory", "#009392")
)

pdf(file="figA07_Chord_mRMR.pdf", width = 11, height = 12)
cowplot::plot_grid(
  p_chord,
  nrow = 1,
  rel_widths = c(0.1, 1, 0.1)
)
dev.off()

