
# UALCAN - Expression Survival Curves

#library(gridSVG) #  install.packages("gridSVG")
#library(rsvg) # install.packages("rsvg")

## SURVIVAL ----

library(grImport2) # install.packages("grImport2")
library(ggpubr) #  install.packages("ggpubr")
library(fs)
library(grid)

setwd("~/bio/gene_signature/8_appendix")


# CPTAC ----

ff <- list.files(path = "UALCAN/CPTAC/cancer_stage/", pattern = ".png", include.dirs = FALSE)

img <- list()
i <- 1
for(f in ff){
  p <- png::readPNG(paste0("UALCAN/CPTAC/cancer_stage/",f))
  img[[i]] <- rasterGrob(p, interpolate=TRUE)
  i <- i+1 
  
}

pdf("../figs/figA09_CPTAC.pdf", width = 17, height = 11)
cowplot::plot_grid( img[[1]], img[[2]], NULL, img[[3]], img[[4]], 
                    img[[5]], img[[6]], NULL, img[[7]], img[[8]], 
                    img[[9]], img[[10]], NULL, img[[11]], img[[12]], #align ="l",
                    rel_widths = c(2,1,0.1,2,1),
                    label_size = 13,
                    labels = c("a", "", "", "b", "", 
                               "c", "", "", "d", "", 
                               "e", "", "", "f", ""), ncol=5)
dev.off()



## SURVIVAL ----

ff <- list.files(path = "UALCAN/Survival", pattern = ".svg", include.dirs = FALSE)

img <- list()
i <- 1
for(f in ff){
  p <- grImport2::readPicture(paste0("UALCAN/Survival/",f))
  img[[i]] <- grImport2::pictureGrob(p)
  i <- i+1 
  
}

# f <- function(x){
#   if(is.atomic(x)){
#     list(x)
#   }else{
#     x
#   }
# }
# g <- function(L){
#   out <- unlist(lapply(L, f), recursive=FALSE)
#   while(any(sapply(out, is.list))){
#     out <- g(out)
#   }
#   out
# }


pdf("../figs/fig_supl_surv1.pdf", width = 10, height = 10)

cowplot::plot_grid( img[[1]], img[[2]], 
                    img[[3]], img[[4]], 
                    img[[5]], img[[6]],
                    label_size =12,
                    labels = "auto", ncol=2)
dev.off()



pdf("../figs/fig_supl_surv2.pdf", width = 8, height = 6)
cowplot::plot_grid( img[[7]], img[[8]],
                    img[[9]], img[[10]],
                    label_size = 12,
                    labels = "auto", ncol=2)
dev.off()