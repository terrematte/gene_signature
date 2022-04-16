# Installing and Loading Libraries            

if(!require('tidyverse')) {install.packages('tidyverse')} 
if(!require('grid')){install.packages('grid')} # 
if(!require('gridExtra')) {install.packages('gridExtra')} #  Provides a number of user-level functions to work with "grid" graphics
if(!require('cowplot')) {install.packages('cowplot')} # provides various features that help with creating publication-quality figures
if(!require('ggpmisc')) {install.packages('ggpmisc')} 

# Load Data

load(file = "../data/2.gsign_kidney.rda")
load(file = "../data/3.data.benchmark_tcgadata_3cv_no_colinear.rda") # job_benchmark_grid.R 3.data.benchmark_tcgadata_3cv.rda

gsign_kidney <- gsign_kidney %>%
  dplyr::filter(!id %in% c("sign.3", "sign.4",  "sign.11", "sign.13", "sign.14", "fs.genetic", "fs.mRMR.status", "fs.mRMR.time")) %>%
  mutate(id =  fct_recode(id,
                          Ha.2015.CaInfo = "sign.1",
                          Zhan.2015.CMMM = "sign.2",
                          Dai.2016.Oncotarget = "sign.5",
                          Wan.2017.IJCancer = "sign.6",
                          Chang.2018.Medicine = "sign.7",
                          YuanLeiChen.2018.JCP = "sign.8",
                          Hu.2019.IJMS = "sign.9",
                          LiangChen.2018.JCP = "sign.10",
                          #Pan.2019.MSM.5genes = "sign.11", # removed
                          Pan.2019.MSM = "sign.12",
                          Wu.2019.FrontiersOnco = "sign.15",
                          Jiang.2020.ACS = "sign.16",
                          Zou.2020.PeerJ = "sign.17",
                          LingChen.2020.Hereditas = "sign.18",
                          DCosta.2020.SciReports = "sign.19",
                          GBM = "filt.gbm",
                          Rpart = "filt.rpart",
                          XGBoost = "filt.xgboost",
                          Boruta = "fs.boruta",
                          RFE = "fs.mlr.rfe",
                          mRMR = "fs.mRMR"
  ))

bmr.measure <- bmr_table_pp %>% 
  dplyr::filter(!task_id %in% c("sign.3", "sign.4",  "sign.11", "sign.13", "sign.14", "fs.genetic", "fs.mRMR.status", "fs.mRMR.time")) %>%
  dplyr::filter(pp %in% "BoxCox") %>%
  dplyr::select(-c("pp")) %>%
  dplyr::mutate(task_id =  fct_recode(task_id,
                                      Ha.2015.CancerInfo = "sign.1",
                                      Zhan.2015.CMMM = "sign.2",
                                      Dai.2016.Oncotarget = "sign.5",
                                      Wan.2017.IJCancer = "sign.6",
                                      Chang.2018.Medicine = "sign.7",
                                      YuanLeiChen.2018.JCP = "sign.8",
                                      Hu.2019.IJMS = "sign.9",
                                      LiangChen.2018.JCP = "sign.10",
                                      #Pan.2019.MSM.5genes = "sign.11",
                                      Pan.2019.MSM = "sign.12",
                                      Wu.2019.FrontiersOnco = "sign.15",
                                      Jiang.2020.ACS = "sign.16",
                                      Zou.2020.PeerJ = "sign.17",
                                      LingChen.2020.Hereditas = "sign.18",
                                      DCosta.2020.SciReports = "sign.19",
                                      GBM = "filt.gbm",
                                      Rpart = "filt.rpart",
                                      XGBoost = "filt.xgboost",
                                      Boruta = "fs.boruta",
                                      RFE = "fs.mlr.rfe",
                                      mRMR = "fs.mRMR"))

rownames(gsign_kidney) <- gsign_kidney$id

gsign_kidney["mRMR", "signature"][[1]] <- list(c("AL353637.1", "AR",  "DPP6", "FOXJ1", "FOXJ1", "GNB3", "HHLA2", "IL4", "LIMCH1", "LINC01732", "OTX1", "SAA1", "SEMA3G", "ZIC2"))

# Heatmap of Signatures metrics for 3cv with 10 repeats
## Arrange signatures means in increasing order

task_id_order <- bmr.measure %>% 
  group_by(task_id) %>% 
  summarise(surv.uno_auc = mean(surv.uno_auc)) %>% 
  arrange(surv.uno_auc)

task_id_order$learner_id <- "mean_signature"
task_id_order$resampling_id <- NA
task_id_order$iters <- NA
task_id_order$surv.harrell_c <- NA

bmr.measure <- rbind(bmr.measure, task_id_order)

bmr.measure$task_id <- factor(bmr.measure$task_id, levels = c("mean_models", as.character(task_id_order$task_id)))

## Arrange learners medians in increasing order

learner_id_order <- bmr.measure %>% 
  group_by(learner_id) %>% 
  summarise(surv.uno_auc = mean(surv.uno_auc)) %>% 
  arrange(surv.uno_auc)

learner_id_order$task_id <- "mean_models"
learner_id_order$resampling_id <- NA
learner_id_order$iters <- NA
learner_id_order$surv.harrell_c <- NA

bmr.measure <- rbind(bmr.measure, learner_id_order)

bmr.measure$learner_id <- factor(bmr.measure$learner_id, levels = c("mean_signature", setdiff(unique(learner_id_order$learner_id), "mean_signature")))

# Plot Heatmap Cross-validation with 100 repetitions 

p1 <- ggplot(bmr.measure, aes(learner_id, task_id)) +
  geom_tile(aes(fill = surv.uno_auc)) + 
  geom_text(aes(label = round(surv.uno_auc, 2)), colour = "white") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  ggtitle("Mean AUC on TCGA-KIRC dataset with 3 CV and 100 repeats")+
  ylab("Gene Signatures & Feature Selection Methods") +
  xlab("Models") +
  theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position="none") +
  scale_x_discrete(expand=c(0,0)) #+ 
# ggsave("../figs/fig4_heatmap_3cv.pdf", width = 14, height = 9, units = c("in"))

rm(list = setdiff(ls(), c("p1","gsign_kidney")))

## Arrange learners means in increasing order

load(file = "../data/3.data.benchmark_train.tcga.test.icgc.rda") # job_benchmark_grid.R 3.data.benchmark_train.tcga.test.icgc.rda

bmr.measure <- bmr_table_pp %>% 
  dplyr::filter(!task_id %in% c("filt.gbm", "filt.rpart", "filt.xgboost", "fs.boruta", "fs.mlr.rfe", "fs.mRMR"))

load(file = "../data/3.data.benchmark_train.tcga.test.icgc_no_colinear.rda") # job_benchmark_grid.R 3.data.benchmark_train.tcga.test.icgc.rda

bmr_table_pp <- rbind(bmr_table_pp,bmr.measure)

bmr.measure <- bmr_table_pp %>% 
  dplyr::filter(!task_id %in% c("sign.3", "sign.4",  "sign.11", "sign.13", "sign.14", "fs.genetic", "fs.mRMR.status", "fs.mRMR.time")) %>%
  dplyr::select(-c("pp")) %>%
  dplyr::mutate(task_id =  fct_recode(task_id,
                                      Ha.2015.CancerInfo = "sign.1",
                                      Zhan.2015.CMMM = "sign.2",
                                      Dai.2016.Oncotarget = "sign.5",
                                      Wan.2017.IJCancer = "sign.6",
                                      Chang.2018.Medicine = "sign.7",
                                      YuanLeiChen.2018.JCP = "sign.8",
                                      Hu.2019.IJMS = "sign.9",
                                      LiangChen.2018.JCP = "sign.10",
                                      #Pan.2019.MSM.5genes = "sign.11",
                                      Pan.2019.MSM = "sign.12",
                                      Wu.2019.FrontiersOnco = "sign.15",
                                      Jiang.2020.ACS = "sign.16",
                                      Zou.2020.PeerJ = "sign.17",
                                      LingChen.2020.Hereditas = "sign.18",
                                      DCosta.2020.SciReports = "sign.19",
                                      GBM = "filt.gbm",
                                      Rpart = "filt.rpart",
                                      XGBoost = "filt.xgboost",
                                      Boruta = "fs.boruta",
                                      RFE = "fs.mlr.rfe",
                                      mRMR = "fs.mRMR"))

task_id_order <- bmr.measure %>% 
  group_by(task_id) %>% 
  summarise(surv.uno_auc = mean(surv.uno_auc)) %>% 
  arrange(surv.uno_auc)

task_id_order$learner_id <- "mean_signature"
task_id_order$resampling_id <- NA
task_id_order$iters <- NA
task_id_order$surv.harrell_c <- NA

bmr.measure <- rbind(bmr.measure, task_id_order)

levs <-  c("mean_models", setdiff(as.character(task_id_order$task_id), c("mRMR", "GBM")), c("mRMR", "GBM"))

bmr.measure$task_id <- factor(bmr.measure$task_id, levels = levs )

learner_id_order <- bmr.measure %>% 
  group_by(learner_id) %>% 
  summarise(surv.uno_auc = mean(surv.uno_auc)) %>% 
  arrange(surv.uno_auc)

learner_id_order$task_id <- "mean_models"
learner_id_order$resampling_id <- NA
learner_id_order$iters <- NA
learner_id_order$surv.harrell_c <- NA

bmr.measure <- rbind(bmr.measure, learner_id_order)

bmr.measure$learner_id <- factor(bmr.measure$learner_id, levels = c("mean_signature", setdiff(unique(learner_id_order$learner_id), "mean_signature")))

# Plot Heatmap Cross-validation with 100 repetitions 

p2 <- ggplot(bmr.measure, aes(learner_id, task_id)) +
  geom_tile(aes(fill = surv.uno_auc)) + 
  geom_text(aes(label = round(surv.uno_auc, 2)), colour = "white") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  ggtitle("Mean AUC on training TCGA-KIRC and testing ICGC-RECA")+
  ylab("Gene Signatures & Feature Selection Methods") +
  xlab("Models") +
  theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position="none")
scale_x_discrete(expand=c(0,0)) #+ 
#ggsave("../figs/fig4_heatmap_trn.tcga.tst.icgc.pdf", width = 14, height = 9, units = c("in"))

rm(learner_id_order, task_id_order)

# Box plot only with significants genes

load("bmk_tcga_3cv_10repeats_no_collinear_significance.rda")

bmk_tcga <- bmk_tcga %>% 
  dplyr::mutate(sign =  fct_recode(sign,
                                   Ha.2015.CancerInfo = "Ha.2014.CancerInfo",
                                   GBM = "filt.gbm",
                                   Rpart = "filt.rpart",
                                   XGBoost = "filt.xgboost",
                                   Boruta = "wrap.boruta",
                                   RFE = "wrap.rfe",
                                   mRMR = "wrap.mRMR"))

bmk_tcga_auc <- bmk_tcga %>% 
  group_by(sign) %>%
  summarise(auc=mean(auc_uno)) %>% 
  arrange(auc) %>%
  as.data.frame(.)

bmk_tcga$sign <- factor(bmk_tcga$sign, levels = c(as.character(bmk_tcga_auc$sign)))

dtF <- rbind(
  data.frame(sign=bmk_tcga$sign, num=log10(bmk_tcga$p.val), cut=0.05, what="p-value of survival risk"),
  data.frame(sign=bmk_tcga$sign, num=bmk_tcga$auc_uno, cut=0.65, what="auc on 10 years"))

dtF$color <- ifelse(dtF$sign %in% bmk_tcga_auc$sign[bmk_tcga_auc$p < 0.05],"#FF1F5B", "gray80")

g1 <- ggplot() +
  geom_boxplot(data=dtF[dtF$what=="p-value of survival risk",], mapping = aes(x = sign, y = num, colour=color), 
               outlier.size=0.2, width=0.2) +
  geom_hline(data=dtF[dtF$what=="p-value of survival risk",], yintercept=log10(0.05), 
             linetype="dashed", color = "red") +
  annotate("text", bmk_tcga_auc$sign[4], log10(0.05), vjust = -0.5, label = "mean of p-value < 0.05", color = "red") +
  scale_x_discrete(NULL, labels = NULL)+
  ggtitle("TCGA-KIRC with CV 3 folds in 100 repeats")+
  theme(plot.title=element_text(), axis.text.x = element_text(),
        axis.ticks.x.bottom = element_line(color = "white"),
        legend.position="none") +
  facet_wrap(~what, strip.position = "right") +
  coord_cartesian(ylim=c(-16,0)) +
  scale_color_manual(values=c("#FF1F5B", "gray60")) +
  xlab("") +
  ylab("log10(p-value)")

g2 <-ggplot() +
  geom_boxplot(data=dtF[dtF$what=="auc on 10 years",], mapping = aes(x = sign, y = num), 
               color = "#009ADE", outlier.size=0.2, width=0.2) +
  geom_text(data=bmk_tcga_auc,
            aes(y=auc, x = sign, label=round(auc,2)), angle=45, vjust = -1.4, size=3,
            color='#0D4A70', position=position_dodge(1)
  ) +
  theme(plot.title=element_text(hjust=0.5, size=12), 
        axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.y = element_text(margin = margin(t = 0, r = 13, b = 0, l = 0))) +
  facet_wrap(~what, strip.position = "right") +
  #scale_y_continuous(breaks=seq(0.4,0.7,0.05)) +
  #coord_cartesian(ylim=c(0.4,0.7))+
  xlab("Gene Signature & Feature Selection Methods") +
  ylab("auc_uno")
#pdf("../figs/KIRC_2_benchmark_fig1_boxplot_tcga.cv3.pdf", width = 10, height = 6)
gE1 <- gridExtra::grid.arrange(g1, g2, ncol=1, heights=c(0.5, 1))
#dev.off()

# Train TCGA Test ICGC - 100 repeats
## Arrange signatures medians in increasing order

# KIRC_2_benchmark_boxplot_job.R 
load("~/gene_signature/4_gene_selection/validation_tcga_icgc_2555.mrmr.rda")

bmk_icgc <- bmk_icgc %>% 
  dplyr::filter(!sign %in% c("sign.3", "sign.4",  "sign.11", "sign.13", "sign.14", "wrap.mRMR")) %>%
  dplyr::mutate(sign =  fct_recode(sign,
                                   Ha.2015.CancerInfo = "Ha.2014.CaInfo",
                                   GBM = "filt.gbm",
                                   Rpart = "filt.rpart",
                                   XGBoost = "filt.xgboost",
                                   Boruta = "wrap.boruta",
                                   RFE = "wrap.rfe",
                                   #mRMR = "wrap.mRMR",
                                   mRMR = "mRMR.search"))

bmk_icgc_auc <- bmk_icgc %>% 
  group_by(sign) %>%
  summarise(auc=median(auc_uno)) %>% 
  arrange(auc) %>%
  as.data.frame(.)

bmk_icgc_pval <- bmk_icgc %>% 
  group_by(sign) %>%
  summarise(p=median(padj)) %>% 
  arrange(desc(p)) %>%
  as.data.frame(.)

bmk_icgc$sign <- factor(bmk_icgc$sign, levels = union(setdiff(as.character(bmk_icgc_auc$sign),"mRMR"), "mRMR"))

dtF <- rbind(
  data.frame(sign=bmk_icgc$sign, num=bmk_icgc$padj, cut=0.05, what="p-adj of survival risk"),
  data.frame(sign=bmk_icgc$sign, num=bmk_icgc$auc_uno, cut=0.65, what="auc on 7 years"))

dtF$color <- ifelse(dtF$sign %in% bmk_icgc_pval$sign[bmk_icgc_pval$p < 0.05], "blue", "gray80") # #FF1F5B

g1 <- ggplot() +
  geom_boxplot(data=dtF[dtF$what=="p-adj of survival risk",], mapping = aes(x = sign, y = num, colour=color), 
               outlier.size=0.2, width=0.2) +
  geom_text(data=bmk_icgc_pval,
            aes(y=p, x = sign, label=ifelse(p<0.05, round(p,3), round(p,2) )),  
            vjust = -1.3, hjust = 0.2, angle=45, size=3,
            color='#0D4A70', position=position_dodge(1)
  ) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "blue") +
  annotate("text", bmk_icgc_auc$sign[4],  0.1, vjust = -0.5, label = "p-adjusted < 0.05", color = "blue") +
  scale_x_discrete(NULL, labels = NULL)+
  ggtitle("Training TCGA-KIRC and Testing ICGC-RECA in 100 repeats")+
  theme(plot.title=element_text(), axis.text.x = element_text(),
        axis.ticks.x.bottom = element_line(color = "white"),
        legend.position="none") +
  facet_wrap(~what, strip.position = "right") +
  scale_y_continuous(breaks=seq(0,1,0.5)) +
  coord_cartesian(ylim=c(0.05,1.1)) +
  scale_color_manual(values=c("#FF1F5B", "gray60")) +
  xlab("") +
  ylab("p-adjusted value")

g2 <-ggplot() +
  geom_boxplot(data=dtF[dtF$what=="auc on 7 years",], mapping = aes(x = sign, y = num), 
               color = "#009ADE", outlier.size=0.2, width=0.2) +
  geom_text(data=bmk_icgc_auc,
            aes(y=auc, x = sign ,label=round(auc,2)),  vjust = -1, angle=45, size=3,
            color='#0D4A70', position=position_dodge(1)
  ) +
  theme(plot.title=element_text(hjust=0.5, size=12), 
        axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_wrap(~what, strip.position = "right") +
  scale_y_continuous(breaks=seq(0.4,0.8,0.05)) +
  coord_cartesian(ylim=c(0.4,0.8))+
  xlab("Survival Predition of Gene Signatures & Feature Selection Methods") +
  ylab("AUC-Uno")
gE2 <- gridExtra::grid.arrange(g1, g2, ncol=1, heights=c(0.5, 1))

pdf("../figs/fig3_heatmap_boxplots.pdf", width = 9, height = 12)
  cowplot::plot_grid(p1, gE2, labels = c("a","b"), nrow=2)
dev.off()

# pdf("../figs/fig_supl_heatmap_ICGC.pdf", width = 8, height = 6)
# cowplot::plot_grid(p2, labels = "", nrow=1)
# dev.off()
# 
# pdf("../figs/fig_supl_boxplots_TCGA.pdf", width = 8, height = 6)
# cowplot::plot_grid(gE1, labels = "", nrow=1)
# dev.off()
