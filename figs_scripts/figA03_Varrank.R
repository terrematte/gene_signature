if(!require("mRMRe")){BiocManager::install("mRMRe")}
if(!require("survival")){BiocManager::install("survival")}
if(!require("varrank")){BiocManager::install("varrank")}

load("../data/1.data.processed.rda")
load("../data/1.gsign_kidney.processed.rda")

list_mRMR_basic <- c("HHLA2", "FOXJ1", "DPP6", "LIMCH1", "SAA1", "ZIC2", "IL4", "OTX1", "AL353637.1", "GNB3")

varrank_df <- data[, c("ajcc.stage", list_mRMR_basic)]

varrank <- varrank(data.df = varrank_df, 
                   method = "estevez", 
                   #method = "peng",
                   variable.important = "ajcc.stage", 
                   #variable.important = "status", 
                   discretization.method = "sturges", 
                   n.var = 10,
                   algorithm = "forward", 
                   scheme="mid", 
                   verbose = FALSE)

pdf("../figs/figA03_Varrank.pdf", width = 7, height = 6)
  plot(x = varrank, colsep = F, rowsep = F,cellnote = F, labelscex =1)
dev.off()
