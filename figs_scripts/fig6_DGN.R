if(!require(disgenet2r)){devtools::install_bitbucket("ibi_group/disgenet2r")}


disgenet_api_key <- get_disgenet_api_key(
  email = "patrick.terrematte.013@ufrn.edu.br", 
  password = "6BPE9aL444dStAY" )

Sys.setenv(DISGENET_API_KEY = disgenet_api_key)

genes_mRMR <- c("AL353637.1", "AR",  "DPP6", "FOXJ1", "GNB3", "HHLA2","IL4", "LIMCH1", "LINC01732", "OTX1", "SAA1", "SEMA3G", "ZIC2")

data <- gene2disease(
  gene     = genes_mRMR,
  #vocabulary = "ENTREZ",
  score = c(0, 1),
  database = "ALL",
  #score = c(0.2, 1),
  verbose  = TRUE
)

plot( data,
      class  ="Heatmap",
      limit  = 50, nchars = 10 )

plot( data,
      class="DiseaseClass", nchars=60)


# C0007134 - Renal cell carcinoma
# C0740457 - Malignant tumor of kidney
# C1306837 - Papillary renal cell carcinoma
# C0686172 - Kidney carcinoma in situ
# C1378703 - Kidney Carcinoma
# C1868672 - Nephrotic syndrome, type 2
# C1336839 - Renal cell carcinoma, papillary, 1
# CN074294 - Renal cell carcinoma, nonpapillary
# C1857453 - Renal hypoplasia/aplasia
# CN227682 - Non-syndromic renal or urinary tract malformation
# C0431693 - Renal cysts and diabetes syndrome
# C1408258 - Kidney damage
# C0542518 - Enlarged kidney
# CN238732 - Lethal polycystic kidney disease  
# C2609414 - Acute kidney injury

list_of_genes <- disease2gene(disease = "C0004352", database = "GENOMICS_ENGLAND")

data_GDA <- disease2evidence( disease  = c("C0007134", "C0740457", "C1306837", "C0686172", 
                                           "C1868672",  "C1378703",  "C1336839", "CN074294",
                                           "C1857453", "CN227682", "C0431693", "C1408258", 
                                           "C0542518", "CN238732", "C2609414"),
                              gene = genes_mRMR,
                              type = "GDA",
                              database = "ALL",
                              #score    = c( 0.4,1 ),
                              warnings = TRUE)

results <- extract(data_GDA)

# Enrichment ----


genes_mRMR <- c("AL353637.1", "AR",  "DPP6", "FOXJ1", "GNB3", "HHLA2","IL4", "LIMCH1", "LINC01732", "OTX1", "SAA1", "SEMA3G", "ZIC2")
res_enrich <-disease_enrichment( entities =genes_mRMR, vocabulary = "HGNC",
                                 database = "ALL" )

table1 <- res_enrich@qresult[res_enrich@qresult$FDR < 0.05, c("ID", "Description", "disease_semantic_type",  "FDR", "Ratio",  "BgRatio", "shared_symbol")]

# C0022658 - Kidney Diseases
# C1333015 - Childhood Kidney Wilms Tumor
# C2931852 - Clear-cell metastatic renal cell carcinoma
# C1333001 - Childhood Renal Cell Carcinoma
# 
# C0006826 - Malignant Neoplasms
# C0027627 - Neoplasm Metastasis
# C0027651 - Neoplasms
# C1269955 - Tumor Cell Invasion
# 
# C0027726 -  Nephrotic Syndrome
# C0279680 - Transitional cell carcinoma of bladder
# c(C0022658, C1333015, C2931852, C1333001, C0006826, C0027627, C0027651, C1269955, C0027726, C0279680)

res_enrich@qresult <- res_enrich@qresult[res_enrich@qresult$ID %in% c("C0022658", "C1333015", "C2931852", "C1333001", "C0006826", 
                                                                      "C0027627", "C0027651", "C1269955", "C0027726", "C0279680"), ]
p_enrich <- plot(res_enrich, class = "Enrichment", cutoff= 0.05) +
  ggtitle("") +  
  theme(plot.title = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +
  theme_bw() 


# --- 

##' Get genes
##' 
##' @method geneID enrichResult
##' @export
geneID <- function(x) as.character(x@qresult$shared_symbol)


##' Get genes for each Categories
##'
##' @method geneInCategory enrichResult
##' @export
##' @importFrom stats setNames
geneInCategory.enrichResult <- function(x)
  setNames(strsplit(geneID(x), ";", fixed=TRUE), rownames(x@qresult))

##' Convert a list of gene IDs to data.frame object.
##' 
##' @title Convert gene IDs to data.frame object
##' @param inputList A list of gene IDs
##' @return a data.frame object.
##' @noRd
list2df <- function(inputList) {
  # ldf <- lapply(1:length(inputList), function(i) {
  ldf <- lapply(seq_len(length(inputList)), function(i) {
    data.frame(categoryID=rep(names(inputList[i]),
                              length(inputList[[i]])),
               Gene=inputList[[i]])
  })
  
  do.call('rbind', ldf)
}

geneSets <-geneInCategory.enrichResult(res_enrich)
geneSets <- geneSets[order(sapply(lapply(geneSets, length),'[[',1))]
names(geneSets) <- res_enrich@qresult[names(geneSets), "Description"]

labels <- c("Neoplasms", "Neoplasm Metastasis", "Malignant Neoplasms", "Tumor Cell Invasion", "Transitional cell carcinoma of bladder", "Kidney Diseases", "Childhood Renal Cell Carcinoma", "Clear-cell metastatic renal cell carcinoma", "Childhood Kidney Wilms Tumor", "Nephrotic Syndrome")
d <- list2df(geneSets)
d$foldChange <- 0
p_heat  <- ggplot(d, aes_(~Gene, ~categoryID)) + 
  ylab(NULL) +
  scale_y_discrete(limits=rev(labels)) + 
  geom_tile(color = 'white') + 
  geom_tile(aes_(fill = ~foldChange), color = "white") +
  scale_fill_continuous(low="blue", high="red", name = "fold change")  +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))


pdf(file="fig6_DGN.pdf", width = 12, height = 4)
cowplot::plot_grid(
  p_heat, p_enrich,
  ncol = 2,
  labels = "auto",
  rel_widths =c(0.9, 1),
  align = 'h'
)
dev.off()