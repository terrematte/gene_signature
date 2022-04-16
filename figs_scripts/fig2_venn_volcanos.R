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

# Venn diagrams of signatures

genes_DEA_M0 <- readLines("../2_gene_expression_differentiation/genes.DEA.M0.vs.M1.lst")
genes_DEA_NT <- readLines("../2_gene_expression_differentiation/genes.DEA.NT.vs.TP.lst")
genes_eqtls <- readLines("../data/genes_kidney_cortex_eqtls.lst")
genes_papers <- unique(unlist(gsign_kidney$signature[c(1:14)])) # readLines("../data/genes_papers.lst")

genes_DEA_M0 <- gsub("-", ".", genes_DEA_M0)
genes_DEA_NT <- gsub("-", ".", genes_DEA_NT)
genes_eqtls <- gsub("-", ".", genes_eqtls)
genes_papers <- gsub("-", ".", genes_papers)
genes_mRMR <- c("AL353637.1", "AR",  "DPP6", "FOXJ1", "GNB3", "HHLA2", "IL4", "LIMCH1", "LINC01732", "OTX1", "SAA1", "SEMA3G", "ZIC2")

genes_selection <-  unique(union(union(union(genes_DEA_M0, genes_DEA_NT), genes_papers), genes_eqtls))
genes_all <- union(union(union(union(genes_DEA_M0, genes_DEA_NT), genes_eqtls), genes_papers), genes_mRMR)

named.list.for.venn <-
   list(`DEA_M0` = genes_DEA_M0,
        `DEA_NT` = genes_DEA_NT,
        genes_eqtls = genes_eqtls,
        genes_papers = genes_papers)

library(VennDiagram)

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
pvenn <- venn.diagram(
    x = named.list.for.venn,
    #filename = "../figs/fig2_venn.png",
    filename = NULL,
    category.names = c("Genes DEA M0 vs M1" , "Genes DEA NT vs TP" , "Genes Kidney eqtls", "Genes of Selected Papers"),
    imagetype = "png",
    #col = "transparent",,
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"), #"#ed85b0"
    #alpha = 0.50,
    # label.col = c("orange", "white", "orange", "white", "darkblue", "white", "darkgreen", "white"),
    cex = .9,
    fontfamily = "sans",
    #fontface = "bold",
    #cat.col = c("darkblue", "darkgreen", "orange"),
    cat.fontfamily = "sans",
    cat.cex = .75,
    #cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.dist = c(0.055, 0.055, 0.1, 0.1)
    )

rm(list = setdiff(ls(), c("pvenn", "gsign_kidney", "genes_DEA_M0", "genes_DEA_NT", "genes_eqtls", "genes_papers", "genes_mRMR", "genes_all"))) 



# Load data DEA


load("../data/1.data.processed.rda")
load("../data/dea.M0.M1.rda")
load("../data/dea.NT.TP.rda")






mydata <- dea.NT.TP %>% 
  dplyr::filter(dea.NT.TP$symbol %in% genes_all) 

mydata <- mydata %>% 
    mutate(gene_label=ifelse(symbol %in% genes_mRMR, symbol, "")) %>%    # For only genes on mRMR
    mutate(sig=ifelse((logFC >= 3 & pvalue < 0.05), "up", ifelse((logFC <= -3 & pvalue < 0.05), "down", "notsig")))

colours <- setNames(c( "#009392", "darkgrey", "#CF597E"), c("down", "notsig", "up"))
# colours <- setNames(c("cornflowerblue", "grey", "firebrick"), c("down", "notsig", "up"))
# c('#CF597E','#E9E29C','#009392')
p1 <- ggplot(mydata, aes(x = logFC, y = -log10(FDR), col=sig)) +
  geom_point(alpha = 0.5)  +
  scale_color_manual(values=colours) +
  geom_vline(xintercept=c(-3,3), linetype="dotted") +
  geom_hline(yintercept=c(-log10(0.01)), linetype="dotted") +
  ggtitle("Normal vs Tumor") +
  xlab("Gene expression change\n log2(FC)") +
  ylab("Significance\n -log10(FDR)") +
  xlim(c(-50,50)) +
  ylim(c(-2,550)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), legend.position="none") +
  geom_label_repel(data = mydata[mydata$logFC > -3 & mydata$logFC < 3 , ],
           aes(label=gene_label),
           seed = 123,
           xlim = c(NA, 25),  # reduce overlaps
           ylim = c(430, NA),
           direction = "both",
           max.time = 5,
           max.iter = 1e5,
           force = 1,
           size = 3,
           max.overlaps = Inf) +
  geom_label_repel(data = mydata[mydata$logFC > 3, ],
            aes(label=gene_label),
            seed = 123,
            xlim = c(30, NA),  # reduce overlaps
            ylim = c(-1,500),
            direction = "y",
            max.time = 5,
            max.iter = 1e5,
            force = 2,
            size = 3,
            max.overlaps = Inf) +
  geom_label_repel(data = mydata[mydata$logFC < -3, ],
            aes(label=gene_label),
            seed = 123,
            xlim = c(NA, -30),  # reduce overlaps
            ylim = c(0, 450),
            #nudge_x = -30  - subset(mydata, logFC < -3)$logFC,
            hjust = 1,
            direction = "y",
            max.time = 5,
            max.iter = 1e5,
            force = 1,
            size = 3,
            max.overlaps = Inf)


mydata <- dea.M0.M1 %>% 
  dplyr::filter(symbol %in% genes_all) 


mydata <- mydata %>% 
    mutate(gene_label=ifelse(symbol %in% genes_mRMR, symbol, "")) %>%    # For only genes on mRMR
    mutate(sig=ifelse((logFC >= 2 & pvalue < 0.05), "up", ifelse((logFC <= -2 & pvalue < 0.05), "down", "notsig")))

p2 <- ggplot(mydata, aes(x = logFC, y = -log10(FDR), col=sig)) +
  geom_point(alpha = 0.5)  +
  scale_color_manual(values=colours) +
        geom_vline(xintercept=c(-2,2), linetype="dotted") +
        geom_hline(yintercept=c(-log10(0.05)), linetype="dotted") +
        ggtitle("Non-metastatic vs Metastatic") +
        xlab("Gene expression change\n log2(FC)") +
        ylab("Significance\n -log10(FDR)") +
        xlim(-15,15) +
        ylim(-1,70) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), legend.position="none") +
        geom_label_repel(data = mydata[mydata$logFC < -2, ],
                         aes(label=gene_label),
                         seed = 123,
                         xlim = c(NA, -7.5),  # reduce overlaps
                         ylim = c(0, 60),
                         nudge_x = -7.5  - subset(mydata, logFC < -2)$logFC,
                         hjust = 1,
                         direction = "y",
                         max.time = 5,
                         max.iter = 1e5,
                         force = 0.2,
                         size = 3,
                         max.overlaps = Inf) +
        geom_label_repel(data = mydata[mydata$logFC > 2, ],
                         aes(label=gene_label),
                         seed = 123,
                         xlim = c(7.5, NA),  # reduce overlaps
                         ylim = c(-1, 70),
                         nudge_x = 7.5  - subset(mydata, logFC > 2)$logFC,
                         direction    = "y",
                         hjust = 0,
                         max.time = 5,
                         max.iter = 1e5,
                         force = 1,
                         size = 3,
                         max.overlaps = Inf) +
          geom_label_repel(data = mydata[mydata$logFC > -2 & mydata$logFC < 2 , ],
                         aes(label=gene_label),
                         seed = 123,
                         xlim = c(NA, NA),  # reduce overlaps
                         ylim = c(40, NA),
                         direction = "both",
                         max.time = 5,
                         max.iter = 1e5,
                         force = 20,
                         size = 3,
                         max.overlaps = Inf)

    
pdf("../figs/fig2_venn_volcanos.pdf", width = 15, height = 4)
  cowplot::plot_grid(NULL, grobTree(pvenn), NULL, p1, NULL,  p2, NULL, 
                     label_size = 20,
                     labels = c("", "a", "", "b", "", "c", ""), 
                     nrow=1, rel_widths = c(0.1, 1.1, 0.2, 1, 0.2, 1, 0.1))
dev.off()


