#  Collinearity Analysis with Correlation

if(!require("corrplot")){install.packages("corrplot")}

sign <-  c("AL353637.1", "AR",  "DPP6", "FOXJ1", "GNB3", "HHLA2", "IL4", "LIMCH1", "LINC01732", "OTX1", "SAA1", "SEMA3G", "ZIC2")

kirc_sign <- data[, c(sign, "age", "status", "ajcc.stage") ]
kirc_sign$ajcc.stage <- as.numeric(sub("T","", kirc_sign$ajcc.stage))

dim(kirc_sign)

kirc_sign %>%
  select_if(is.numeric) %>%
  cor(.) %>%
  as.data.frame() %>%
  mutate(var1 = rownames(.)) %>%
  gather(var2, value, -var1) %>%
  filter(var1 != var2 )  %>%
  arrange(desc(value))

kirc_cor <- cor(kirc_sign)

# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# matrix of the p-value of the correlation
p.mat <- cor.mtest(kirc_sign)
head(p.mat[, 1:5])


col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

pdf("figA05_corrplot.pdf", width = 8, height = 8)
corrplot(kirc_cor, method="color", col=col(200),  
         type="lower", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         tl.cex =0.8,
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)
dev.off()
