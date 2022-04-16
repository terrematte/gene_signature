if(!require(tidyverse)){install.packages("tidyverse")}
if(!require(caret)){install.packages("caret")}
if(!require(ggpubr)){install.packages("ggpubr")}
if(!require(BBmisc)){install.packages("BBmisc")}

load("../data/2.data.full.count.rda")
load("../data/1.data.processed.rda")

col_clini <- c("obs.time", "status", "age", "gender", "neoplasm", "metastasis", "ajcc.stage", "dataset")
col_genes <- setdiff(colnames(data), col_clini)

data_long <- full.data %>% 
  dplyr::select(!c("obs.time", "status", "age", "gender", "neoplasm", "metastasis", "ajcc.stage")) %>% 
  pivot_longer(-dataset,names_to = "genes", values_to = "cts")  

df <- data_long %>% 
  # for each gene
  group_by(genes, dataset) %>% 
  # get mean and variance
  summarise(median = median(cts)) %>% 
  pivot_wider(names_from = "dataset", values_from = "median")


sp_raw <- ggscatter(df, x = "TCGA", y = "ICGC",
                add = "reg.line",  # Add regressin line
                xlab = "TCGA-KIRC", ylab = "ICGC-RECA",
                title = "Raw counts",
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
                ) +
                scale_x_continuous(labels = scales::comma_format(big.mark = ',',
                                decimal.mark = '.')) +
                scale_y_continuous(labels = scales::comma_format(big.mark = ',',
                                                   decimal.mark = '.'))

#makes the scientific notation using "AeB" explicitly write out Ax10^B
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e(.*)$", "\\1e\\2", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", '~`\U00D7 10`^', l)
  # return this as a string
  l
}

# Add correlation coefficient
p1 <-sp_raw + stat_cor(method = "spearman", 
                       label.x = 40000, 
                       label.y = 70000,
                       aes(label = paste(..rr.label.., fancy_scientific(..p.label..), sep = "~`,`~") ) # 
                       ) 

data_long <- data %>% 
  dplyr::select(!c("obs.time", "status", "age", "gender", "neoplasm", "metastasis", "ajcc.stage")) %>% 
  pivot_longer(-dataset,names_to = "genes", values_to = "cts") 

df <- data_long %>% 
  # for each gene
  group_by(genes, dataset) %>% 
  # get mean and variance
  summarise(median = median(cts)) %>% 
  pivot_wider(names_from = "dataset", values_from = "median")

sp_bc <- ggscatter(df, x = "TCGA", y = "ICGC",
                   xlab = "TCGA-KIRC", ylab = "ICGC-RECA",
                add = "reg.line",  # Add regressin line
                title = "Box-Cox transformation",
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
p2 <- sp_bc + stat_cor(method = "pearson", label.x = 3, label.y = 12,
                       aes(label = paste(..rr.label.., fancy_scientific(..p.label..), sep = "~`,`~") ) ) 


pdf("../plots/fig_supl_scatter.pdf", width = 12, height = 4.5)
  cowplot::plot_grid(p1, p2, labels = c("a", "b"), nrow=1)
dev.off()



# DESeq2+vst ===========================================================================



# DESeq2+vst ===========================================================================
library(DESeq2)

vst2 <- function(countdata, condition){
  colData <- data.frame(matrix( nrow=nrow(countdata) ))
  colData$condition <- as.factor(condition)
  rownames(colData) <- rownames(countdata)
  countdata <-  DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(t(countdata)),
                                               colData = colData, design = ~ condition)
  countdata <- DESeq2::estimateSizeFactors( countdata )
  cdsBlind <- DESeq2::estimateDispersions( countdata, fitType="mean")
  vstdata <- DESeq2::varianceStabilizingTransformation( cdsBlind , fitType="mean")
  return( t(assay(vstdata)))
}

data.vst2 <- vst2(full.data[, col_genes], full.data$dataset)
data.vst <- cbind(full.data[, col_clini], data.vst2)
data.vst <- as.data.frame(data.vst)

data_long <- data.vst %>% 
  dplyr::select(!c("obs.time", "status", "age", "gender", "neoplasm", "metastasis", "ajcc.stage")) %>% 
  pivot_longer(-dataset,names_to = "genes", values_to = "cts") 

df <- data_long %>% 
  # for each gene
  group_by(genes, dataset) %>% 
  # get mean and variance
  summarise(median = median(cts)) %>% 
  pivot_wider(names_from = "dataset", values_from = "median")  %>% 
  data.frame()

sp_devst <- ggscatter(df, x = "TCGA", y = "ICGC",
                xlab = "TCGA-KIRC", ylab = "ICGC-RECA",
                add = "reg.line",  # Add regressin line
                title = "Variance-stabilizing transformation with DESeq2",
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE) # Add confidence interval

# Add correlation coefficient
p_devst <- sp_devst +
  stat_cor(method = "pearson", label.x = 3, label.y = 16,
                               aes(label = paste(..rr.label.., fancy_scientific(..p.label..), sep = "~`,`~")) ) 

# caret+rang ===========================================================================
data.range <- full.data
data.range[, col_genes] <- full.data[, col_genes]+1
set.seed(1)
# calculate the pre-process parameters from the dataset
#preprocessParams <- preProcess(data[data$dataset == "TCGA",col_genes], method=c("BoxCox"), fudge=1)
#cols <- setdiff(col_genes, c("GDF1", "RPL36P4", "RPL23AP1", "AC005363.1"))
preprocessParams <- caret::preProcess(full.data[full.data$dataset == "TCGA", col_genes], rangeBounds = c(0,1), method="range")

# transform the dataset using the parameters
data.range[data.range$dataset == "TCGA",col_genes] <- stats::predict(preprocessParams, full.data[full.data$dataset == "TCGA", col_genes])
data.range[data.range$dataset == "ICGC",col_genes] <- stats::predict(preprocessParams, full.data[full.data$dataset == "ICGC", col_genes])
#data.range[,col_genes] <- predict(preprocessParams, full.data[, col_genes])

data_long <- data.range %>% 
  dplyr::select(-c("obs.time", "status", "age", "gender", "neoplasm", "metastasis", "ajcc.stage")) %>% 
  pivot_longer(-dataset,names_to = "genes", values_to = "cts") 

df <- data_long %>% 
  # for each gene
  group_by(genes, dataset) %>% 
  # get mean and variance
  summarise(median = median(cts)) %>% 
  pivot_wider(names_from = "dataset", values_from = "median")

#df <- df[df$ICGC < 1,]
sp_range <- ggscatter(df, x = "TCGA", y = "ICGC",
                      xlab = "TCGA-KIRC", ylab = "ICGC-RECA",
                      add = "reg.line",  # Add regressin line
                      title = "Scaling between zero and one \n(with Caret R package and 'range' method)", 
                      add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                      conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
p_range <- sp_range + 
  stat_cor(method = "pearson", label.x = 0.3, label.y = 7,
           aes(label = paste(..rr.label.., fancy_scientific(..p.label..), sep = "~`,`~") ) ) + 
  xlim(0, 1)

p_range

# BBmisc+rang ===========================================================================


data.rscale <- full.data
set.seed(123)
data.rscale[, col_genes] <- normalize(full.data[, col_genes], method="range", range=c(0,1))

data_long <- data.rscale %>% 
  dplyr::select(-c("obs.time", "status", "age", "gender", "neoplasm", "metastasis", "ajcc.stage")) %>% 
  pivot_longer(-dataset,names_to = "genes", values_to = "cts") 

df <- data_long %>% 
  # for each gene
  group_by(genes, dataset) %>% 
  # get mean and variance
  summarise(median = median(cts)) %>% 
  pivot_wider(names_from = "dataset", values_from = "median")

#df <- df[df$ICGC < 1,]
sp_rscale <- ggscatter(df, x = "TCGA", y = "ICGC",
                      xlab = "TCGA-KIRC", ylab = "ICGC-RECA",
                      add = "reg.line",  # Add regressin line
                      title = "Scaling between zero and one \n(with BBmisc R package and 'range' method)", 
                      add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                      conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
p_rscale <- sp_rscale + 
  stat_cor(method = "pearson", label.x = 0.3, label.y = 1,
           aes(label = paste(..rr.label.., fancy_scientific(..p.label..), sep = "~`,`~") ) ) + 
  xlim(0, 1)

p_rscale

# log2(count+1) ===========================================================================

data.log <- full.data
data.log[, col_genes] <- log2(data.log[, col_genes]+1)

data_long <- data.log %>%
  dplyr::select(!c("obs.time", "status", "age", "gender", "neoplasm", "metastasis", "ajcc.stage")) %>%
  pivot_longer(-dataset,names_to = "genes", values_to = "cts")

df <- data_long %>%
  # for each gene
  group_by(genes, dataset) %>%
  # get mean and variance
  summarise(median = median(cts)) %>%
  pivot_wider(names_from = "dataset", values_from = "median")

#df <- df[df$ICGC < 1,]
sp_log <- ggscatter(df, x = "TCGA", y = "ICGC",
                      xlab = "TCGA-KIRC", ylab = "ICGC-RECA",
                      add = "reg.line",  # Add regressin line
                      title = "log2(count+1) normalisation",
                      add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                      conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
p_log <- sp_log + stat_cor(method = "pearson", 
                           label.x = 7, 
                           label.y = 16,
                           aes(label = paste(..rr.label.., fancy_scientific(..p.label..), sep = "~`,`~") ) )
p_log

# # center+scale+nzv ===========================================================================
# 
# data.csn <- full.data
# data.csn[, col_genes] <- data.csn[, col_genes]+1
# 
# # calculate the pre-process parameters from the dataset
# preprocessParams <- preProcess(full.data[full.data$dataset == "TCGA",col_genes], method=c("center","scale","nzv"))
# 
# # transform the dataset using the parameters
# data.csn[data.csn$dataset == "TCGA",col_genes] <- predict(preprocessParams, full.data[full.data$dataset == "TCGA", col_genes])
# data.csn[data.csn$dataset == "ICGC",col_genes] <- predict(preprocessParams, full.data[full.data$dataset == "ICGC", col_genes])
# #data.csn$pp  <- "center+scale+nzv"
# 
# 
# data_long <- data.csn %>% 
#   dplyr::select(!c("obs.time", "status", "age", "gender", "neoplasm", "metastasis", "ajcc.stage")) %>% 
#   pivot_longer(-dataset,names_to = "genes", values_to = "cts") 
# 
# df <- data_long %>% 
#   # for each gene
#   group_by(genes, dataset) %>% 
#   # get mean and variance
#   summarise(median = median(cts)) %>% 
#   pivot_wider(names_from = "dataset", values_from = "median")
# 
# #df <- df[df$ICGC < 1,]
# sp_csn <- ggscatter(df, x = "TCGA", y = "ICGC",
#                       xlab = "TCGA-KIRC", ylab = "ICGC-RECA",
#                       add = "reg.line",  # Add regressin line
#                       title = "Median gene expression by center+scale+nzv \n near-zero-variance (caret+range)", 
#                       add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#                       conf.int = TRUE # Add confidence interval
# )
# # Add correlation coefficient
# p_csn <- sp_csn + stat_cor(method = "pearson", label.x = -0.2, label.y = 40)
# p_csn

pdf("../figs/figA02_normalization.pdf", width = 12, height = 15)
cowplot::plot_grid(p1, p_log,  p_devst, p2, p_range, p_rscale, labels = "auto", nrow=3)
dev.off()
