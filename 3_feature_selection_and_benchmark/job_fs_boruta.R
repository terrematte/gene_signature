# Installing and Loading Libraries      -----      

setwd("~/gene_signature/4_gene_selection")

# devtools::install_github("collectivemedia/tictoc")
packages_cran = c("tictoc","tidyverse", "caret", "Boruta", "xgboost", "survival")

package.check <- lapply(packages_cran, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

rm(packages_cran, package.check)

load("../data/1.data.processed.rda")

cols <- setdiff(colnames(data), c("neoplasm", "metastasis", "ajcc.stage", "dataset"))

data <- data[, cols]

#  Feature Selection Boruta for Survival ---

tic("Total - Feature Selection Boruta for Survival")

b_vit <- list()
b <-  list()
genes.fs.boruta <- NULL

x <- data[,-c(1,2)]
y <- Surv(time = data$obs.time, event = data$status)

for(i in c(1:10)){
  print(i)
  b[[i]] <-  Boruta(x = x , y = y, maxRuns = 150, getImp=getImpRfZ)
  b_vit[[i]] <- getSelectedAttributes(b[[i]])
  genes.fs.boruta <- union(genes.fs.boruta, getSelectedAttributes(b[[i]])) 
}

toc()

print(length(genes.fs.boruta))

print(genes.fs.boruta)

news <- NULL

for(i in c(1:10)){
  print(paste0(i, ": ", length(b_vit[[i]]), " added ", sum(!b_vit[[i]] %in% news)))
  news <- union(news, b_vit[[i]])
}

save(genes.fs.boruta, b, file = "genes.fs.boruta.rda", compress = T)


# Total - Feature Selection Boruta for Survival: 44704.924 sec elapsed
# [1] 43
# [1] "ZIC2"       "CHAT"       "AMH"        "OTX1"       "BARX1"     
# [6] "TROAP"      "CKAP4"      "ITPKA"      "NUF2"       "KRT75"     
# [11] "KIF18B"     "SLC18A3"    "AL355796.1" "RPL10P19"   "LINC02154" 
# [16] "LINC00973"  "IL4"        "HOTAIRM1"   "Z84485.1"   "LINC02362" 
# [21] "CASP9"      "CCNF"       "RTL1"       "BID"        "CHGA"      
# [26] "RANBP3L"    "ZIC5"       "SLC16A12"   "SPATC1L"    "age"       
# [31] "CD44"       "KRI1"       "RUFY4"      "AC073324.1" "AC091812.1"
# [36] "AC156455.1" "AGAP6"      "AC128685.1" "SEMA3G"     "IGFN1"     
# [41] "KLRC2"      "ANXA8"      "AURKB"     
# [1] "1: 16 added 16"
# [1] "2: 15 added 4"
# [1] "3: 15 added 3"
# [1] "4: 20 added 3"
# [1] "5: 20 added 3"
# [1] "6: 23 added 5"
# [1] "7: 21 added 2"
# [1] "8: 17 added 2"
# [1] "9: 24 added 3"
# [1] "10: 17 added 2"
