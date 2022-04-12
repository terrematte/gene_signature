
setwd("~/gene_signature")

if(!require('tidyverse')){install.packages('tidyverse')} 
if(!require('clusterProfiler')){BiocManager::install("clusterProfiler")} 
if(!require('enrichplot')){BiocManager::install("enrichplot")} 
if(!require('org.Hs.eg.db')){BiocManager::install("org.Hs.eg.db")} 
if(!require('annotables')){BiocManager::install("stephenturner/annotables")} 
if(!require('ReactomePA')){BiocManager::install("ReactomePA")} 
if(!require('DOSE')){BiocManager::install("DOSE")} 
if(!require('enrichplot')){BiocManager::install("enrichplot")} 
if(!require(gprofiler2)){install.packages("gprofiler2")}

geneslist_mRMR <- c("CLEC18C", "AL121987.1", "MT1G", "MYADML2", "OPN4", "NFE4", "CASR", "L138830.2", "AL391095.1", "ANXA8", "TPSD1", "AL357507.1", "ORM2", "LINC00865", "AC106047.1", "IDO1", "NR0B1", "CRH", "GNB3", "HHLA2", "FOXJ1", "GOLGA6L9", "DACH2", "LINC01732", "FLG", "LINC00896", "LINC01956", "EEF1DP5", "IGFN1", "DPP6", "HRH2", "SULT4A1", "AC078950.1", "AP000697.1", "PLG", "KLRC2", "KISS1", "C3orf85", "SAA1", "LINC00973", "MICG", "HOTAIRM1", "SEMA3G", "VSX1", "LINC01436", "ANKRD33", "AC002451.1", "CHAT", "AC093001.1", " ZIC2", "SLC16A12", "AR", "AL353637.1", "AC008663.1", "AC105384.2", "TNNT1", "AL355796.1", "BECN2", "AC091812.1", "IL4", "OTX1", "LINC00524", "WFDC3", "LIMCH1")


sign.entrez1 <- mapIds(org.Hs.eg.db, keys = geneslist_mRMR, keytype = "SYMBOL", column="ENTREZID")
sign.miss <- names(sign.entrez1[is.na(sign.entrez1)])
sign.entrez2 <- grch37 %>%
  dplyr::filter(symbol %in% sign.miss) %>%
  pull(entrez)
sign.miss <- sign.miss[is.na(sign.entrez2)]
sign.entrez3 <- grch38 %>%
  dplyr::filter(symbol %in% sign.miss ) %>%
  pull(entrez)
sign.miss <- sign.miss[is.na(sign.entrez3)]
sign_entrez <- union(union(sign.entrez1, sign.entrez2),sign.entrez3)
sign_entrez <- sign_entrez[!is.na(sign_entrez)]

geneslist <- list()

geneslist$mRMR  <- as.vector(sign_entrez)

# https://rdrr.io/bioc/clusterProfiler/man/compareCluster.html
ck <- clusterProfiler::compareCluster(geneCluster = geneslist,
                                      fun="enrichKEGG",
                                      organism="hsa", pvalueCutoff=0.05)

dotplot(ck)

# enrichNCG

edo_DGN <- DOSE::enrichDGN(geneslist$mRMR)
edox_DGN <- setReadable(edo_DGN, 'org.Hs.eg.db', 'ENTREZID')



edox_DGN@result <- edox_DGN@result %>% 
  #  dplyr::filter(!Description %in% c("Acne", "Dermatitis, Irritant", "Weight decreased", "Hypogonadotropic hypogonadism","Helicobacter pylori (H. pylori) infection in conditions classified elsewhere and of unspecified site","EPIDERMAL DIFFERENTIATION COMPLEX", "Skin Diseases, Infectious", "Chronic idiopathic urticaria", "Bullous pemphigoid", "Dermatitis, Occupational", "Keratoconjunctivitis, Vernal", "Testicular hypogonadism", "Gonadal Dysgenesis, 46,XY", "46, XY female", "Candidiasis, Vulvovaginal", "Depressive Symptoms", "EPIDERMAL DIFFERENTIATION COMPLEX", "Substance Withdrawal Syndrome")) %>%
  dplyr::filter(Description %in% c("Anemia of chronic disease", "Pre−Eclampsia", "Eclampsia", "Chronic kidney disease stage 5", "NEPHROLITHIASIS, CALCIUM OXALATE", "Klinefelter Syndrome" ,"Neoplasm Recurrence, Local", "Hyperplasia", "Carcinoma, Transitional Cell", "Reticulosarcoma")) %>% 
  dplyr::mutate(Description =  fct_recode(Description, 
                                          "Nephrolithiasis, Calcium Oxalate" = "NEPHROLITHIASIS, CALCIUM OXALATE"
  ))

edox2 <- pairwise_termsim(edox_DGN)
p1 <- treeplot(edox2)

ck <- clusterProfiler::compareCluster(geneCluster = geneslist,
                                      fun="enrichNCG")

dotplot(ck)


ck <- clusterProfiler::compareCluster(geneCluster = geneslist,
                                      fun="enrichDO")

dotplot(ck)


ck <- clusterProfiler::compareCluster(geneCluster = geneslist,
                                      fun="enrichDGN")

dotplot(ck)




ck <- clusterProfiler::compareCluster(geneCluster = geneslist,
                                      fun="enrichPathway")

dotplot(ck)



edo <- enrichDGN(geneslist$mRMR)
edo2 <- gseDO(geneslist$mRMR)

dotplot(edo, showCategory=30) 
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
edox2 <- pairwise_termsim(edox)
p1 <- treeplot(edox2)


# https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html#visualization
if(!require(gprofiler2)){install.packages("gprofiler2")}

# enrichment analysis
gp = gost(geneslist_mRMR, 
          ordered_query = FALSE, 
          user_threshold = 0.05,
          exclude_iea = TRUE,
          correction_method = "fdr",
             organism = "hsapiens")

gostplot(gp, capped = F, interactive = TRUE)

p1 = gostplot(gp, interactive = FALSE)

publish_gostplot(p1) #highlight_terms =c("CORUM:1514", "CORUM:2160", "CORUM:6087", "GO:1990405", 
#                                        "REAC:R-HSA-500792", "REAC:R-HSA-4164", "WP:WP176"))


# enrichment analysis
gp_EA = gost(geneslist_EA, 
             user_threshold = 0.05,
             organism = "hsapiens")

multi_gostres2 <- gost(query = list("EA" = geneslist_EA,
                                    "AA" = geneslist_AA), 
                       multi_query = TRUE)

gostplot(multi_gostres2, capped = TRUE, interactive = TRUE)


# enrichment analysis

multi_gost <- gost(query = list("PTC" = genes_mutated_PTC,
                                "notPTC" = genes_mutated_notPTC), 
                   significant = TRUE,
                   ordered_query = TRUE,
                   correction_method = "fdr",
                   #correction_method = "g_SCS",
                   user_threshold = 0.05,
                   multi_query = TRUE,
                   organism = "hsapiens")

p <- gostplot(multi_gost, capped = F, interactive = F)

pp <- publish_gostplot(p, highlight_terms =  c("REAC:R-HSA-975956", "REAC:R-HSA-927802","REAC:R-HSA-975957", "GO:0006955", "GO:0043069"), 
                       width = NA, height = NA, filename = NULL )

```

