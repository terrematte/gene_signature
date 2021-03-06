# "Aalen's additive regression model"
# Installing and Loading Libraries            

# install https://indrajeetpatil.github.io/ggstatsplot/index.html
# remotes::install_github("IndrajeetPatil/ggstatsplot")
# sudo apt-get install libgmp3-dev
# sudo apt-get install libmpfr-dev
if(!require("ggstatsplot")){remotes::install_github("IndrajeetPatil/ggstatsplot")}
if(!require("tidyverse")){install.packages("tidyverse")}
if(!require("ggfortify")){install.packages("ggfortify")}
if(!require("survival")){install.packages("survival")}
if(!require("cowplot")){install.packages("cowplot")}
if(!require("pals")){install.packages("pals")}
if(!require("colorspace")){install.packages("colorspace")}
if(!require("ggsci")){install.packages("ggsci")}


# Load Data ---- 

load("../data/1.data.processed.rda")
sign <- c("age","ZIC2", "SEMA3G", "SAA1", "OTX1", "LINC01732", "LIMCH1", "IL4", "HHLA2", "GNB3", "FOXJ1", "DPP6", "AR", "AL353637.1")
sign <- intersect(sign, colnames(data))
data <- data[!is.na(data$metastasis), ]
data$gender <-ifelse(data$gender == 1, "Male", "Female")

# Aalen's additive regression model for censored data

# The additive model provides an alternative to Cox regression.
# It can provide detailed information about temporal influence of each of the covariates not available in Cox model.

aa_fit <-aareg(as.formula(paste0("Surv(obs.time, status) ~ " , paste(sign, collapse= " + "))),
               data = data)
aa_fit

summary(aa_fit)  # provides a more complete summary of results

labs=c('alphabet','alphabet2', 'glasbey','kelly','polychrome', 'stepped', 'stepped2', 'stepped3', 'tol', 'watlington')
op=par(mar=c(0,5,3,1))
pal.bands(alphabet(), alphabet2(), glasbey(), kelly(),
          polychrome(), stepped(), stepped2(), stepped3(), 
          tol(), watlington(), labels=labs, show.names=FALSE)

  #vars <- c("Intercept", "age", "genderMale",  "metastasisM1", "metastasisMX",  "AL353637.1", "AL138830.2", "DPP6", "FOXJ1", "HHLA2", "MICG", "LIMCH1", "VSX1", "AR", "IL4", "CASR", "CRH", "GNB3", "SAA1")
  
  vars <- c("Intercept", "age", "ZIC2", "SEMA3G", "SAA1", "OTX1", "LINC01732", "LIMCH1", "IL4", "HHLA2", "GNB3", "FOXJ1", "DPP6", "AR", "AL353637.1")
  variables <- rev(factor(vars, levels=vars))
  #-- color decrease of lightness
  cols1 <- as.vector(glasbey(15))
  cols1 <- readhex(file = textConnection(paste(cols1, collapse = "\n")),
                   class = "RGB")
  #transform to hue/lightness/saturation colorspace
  cols1 <- as(cols1, "HLS")
  #additive decrease of lightness
  cols1@coords[, "L"] <- pmax(0, cols1@coords[, "L"] + 0.05)
  cols1 <- as(cols1, "RGB")
  cols1 <- hex(cols1)
  p1 <- ggcoefstats(
    x = aa_fit,
    title = "Aalen's additive regression model",
    subtitle = "(for censored data)",
    only.significant = F,
    #point.args = list(color = "green", shape = 9),
    package = "pals",
    stats.label.args = list(
      max.time = 3,
      direction = "y",
      point.padding = 0.2,
      nudge_x = .15,
      nudge_y = .5,
      segment.curvature = -0.1,
      segment.angle = 10,
      segment.size  = 0.2,
      segment.linetype = 2
      # nudge_x = .15,
      # box.padding = 0.5,
      # nudge_y = 1,
      # segment.curvature = -0.1,
      # segment.ncp = 3,
      # segment.angle = 20
      ),
    palette = "glasbey",
    sort = "none", 
    k = 3
  ) + ggplot2::theme(text = element_text(size=18),
                     axis.text = element_text(size=16)
                     )
  #+  ggplot2::scale_y_discrete(labels = vars) 
  
  p1[["layers"]][[4]][["data"]][["label"]][[1]] <- "list(~widehat(italic(beta))=='7.21 \U00D7 10'^'-3', ~italic(z)=='2.687', ~italic(p)=='0.014')"               
  p1[["layers"]][[4]][["data"]][["label"]][[2]] <- "list(~widehat(italic(beta))=='4.307 \U00D7 10'^'-5', ~italic(z)=='417.530', ~italic(p)=='0.012')" 
  p1[["layers"]][[4]][["data"]][["label"]][[3]] <- "list(~widehat(italic(beta))=='3.991 \U00D7 10'^'-4', ~italic(z)=='32.486', ~italic(p)=='0.054')"  
  p1[["layers"]][[4]][["data"]][["label"]][[4]] <- "list(~widehat(italic(beta))=='-2.954 \U00D7 10'^'-5', ~italic(z)=='-1.515', ~italic(p)=='0.907')" 
  p1[["layers"]][[4]][["data"]][["label"]][[5]] <- "list(~widehat(italic(beta))=='1.395 \U00D7 10'^'-4', ~italic(z)=='61.970', ~italic(p)=='0.068')"  
  p1[["layers"]][[4]][["data"]][["label"]][[6]] <- "list(~widehat(italic(beta))=='5.996 \U00D7 10'^'-4', ~italic(z)=='32.990', ~italic(p)=='0.014')"  
  p1[["layers"]][[4]][["data"]][["label"]][[7]] <- "list(~widehat(italic(beta))=='4.073 \U00D7 10'^'-4', ~italic(z)=='14.360', ~italic(p)=='0.144')"  
  p1[["layers"]][[4]][["data"]][["label"]][[8]] <- "list(~widehat(italic(beta))=='-9.320 \U00D7 10'^'-4', ~italic(z)=='-20.702', ~italic(p)=='0.020')"
  p1[["layers"]][[4]][["data"]][["label"]][[9]] <- "list(~widehat(italic(beta))=='7.920 \U00D7 10'^'-4', ~italic(z)=='27.023', ~italic(p)=='0.006')"  
  p1[["layers"]][[4]][["data"]][["label"]][[10]] <- "list(~widehat(italic(beta))=='-4.194 \U00D7 10'^'-4', ~italic(z)=='-60.500', ~italic(p)=='0.008')"
  p1[["layers"]][[4]][["data"]][["label"]][[11]] <- "list(~widehat(italic(beta))=='2.863 \U00D7 10'^'-4', ~italic(z)=='20.657', ~italic(p)=='0.140')"  
  p1[["layers"]][[4]][["data"]][["label"]][[12]] <- "list(~widehat(italic(beta))=='4.436 \U00D7 10'^'-4', ~italic(z)=='64.615', ~italic(p)=='0.002')"  
  p1[["layers"]][[4]][["data"]][["label"]][[13]] <- "list(~widehat(italic(beta))=='-5.124 \U00D7 10'^'-4', ~italic(z)=='-59.365', ~italic(p)=='0.002')"
  p1[["layers"]][[4]][["data"]][["label"]][[14]] <- "list(~widehat(italic(beta))=='1.693 \U00D7 10'^'-4', ~italic(z)=='6.716', ~italic(p)=='0.554')"   
  p1[["layers"]][[4]][["data"]][["label"]][[15]] <- "list(~widehat(italic(beta))=='-4.680 \U00D7 10'^'-4', ~italic(z)=='-37.743', ~italic(p)=='0.004')"         
  
  
  p2 <- ggplot2::autoplot(aa_fit) +
    theme(legend.position="none", 
          text = element_text(size=16),
          axis.text = element_text(size=10))
  p2$layers[[1]]$data$variable <- factor(p2$layers[[1]]$data$variable,
                                         levels= variables)
  p2$layers[[2]]$data$variable <- factor(p2$layers[[2]]$data$variable,
                                         levels= variables)
  p2$data$variable <- factor(p2$data$variable,
                                         levels= variables)
  p2 <- p2 +
    scale_fill_manual(values=rev(cols1))
  
  pdf("../figs/fig5_mRMR_aareg.pdf", width = 14, height = 8)
    cowplot::plot_grid(p1, p2, labels = "auto", nrow=1, rel_widths = c(0.7,1))
  dev.off()
