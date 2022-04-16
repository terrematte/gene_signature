

if(!require("tidyverse")){install.packages("tidyverse")}
if(!require("survival")){install.packages("survival")}
if(!require("glmnet")){install.packages("glmnet")}
if(!require("rms")){install.packages("rms")} # vif - One possible source for a `vif`-function


# Load Data ---- 

setwd("~/gene_signature/plots")
load("1.data.processed.rda")

# Selecting the feature
sign <- c("AL353637.1", "AR",  "DPP6", "FOXJ1", "GNB3", "HHLA2", "IL4", "LIMCH1", "LINC01732", "OTX1", "SAA1", "SEMA3G", "ZIC2")

genes_vif <- setNames(data.frame(matrix(nrow=0, ncol=3)), c("sign", "gene", "vif"))

sign <- intersect(colnames(data), sign)

kirc_sign <- data[, c("obs.time", "status", sign) ]

surv_obj <- Surv(time = kirc_sign$obs.time, event = kirc_sign$status)

mdl1 <- coxph(as.formula(paste0("surv_obj ~ ", paste(sign, collapse= " + "))), data = kirc_sign)

cvif <- vif(mdl1) 
df_genes_vif <- as.data.frame(cvif)
df_genes_vif$gene  <- rownames(df_genes_vif)
rownames(df_genes_vif) <- NULL
genes_vif <- rbind(genes_vif, df_genes_vif)
print(cvif)


pdf("figA04_Variance_Inflation_Factor.pdf", width = 12, height = 7)
ggplot(data=genes_vif, aes(x=gene, y=cvif)) +
  geom_bar(stat="identity", fill="steelblue") +
  geom_text(aes(label=round(cvif, digits = 2)), vjust = -0.6, hjust = - 0.3, angle = 45, size=4) +
  theme_minimal() +
  theme(plot.title=element_text(hjust=0.5, size=15), 
        axis.text.x = element_text(angle = 45, hjust=1, size=15),
        axis.text.y = element_text(angle = 0, hjust=1, size=13)) +
  ylim(0, 2.5) +
  xlab("Genes") +
  ylab("Variance Inflation Factors") +
  theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
  theme(axis.title.x = element_text(size = rel(1.8), angle = 00))
dev.off()
