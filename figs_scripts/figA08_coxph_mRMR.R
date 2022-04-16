if(!require("tidyverse")){install.packages("tidyverse")}
if(!require("survival")){install.packages("survival")}
if(!require("survminer")){install.packages("survminer")}

explanatory <-c("AR", "AL353637.1", "DPP6", "FOXJ1", "GNB3", "HHLA2", "IL4", "LIMCH1", "LINC01732", "OTX1", "SAA1", "SEMA3G", "ZIC2")

surv_obj <- Surv(time = data$obs.time, event = data$status)

# Fit a Cox proportional hazards model  
fit.coxph <- coxph(as.formula(paste0("surv_obj ~ ", paste(explanatory, collapse= " + "))), data)
fit.coxph

pdf(file = "../figs/figA08_coxph_mRMR", width = 13, height = 6, onefile=FALSE)
ggforest(fit.coxph, data = data[,c("obs.time", "status", explanatory)], fontsize = 1)  
dev.off()