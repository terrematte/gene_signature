# Time Run : 18:33:01
# Installing and Loading Libraries      -----      

setwd("~/gene_signature/4_gene_selection")

packages_cran = c("tidyverse", "caret","mlr3", "mlr3verse", "mlr3proba", "survAUC")

package.check <- lapply(packages_cran, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

rm(packages_cran, package.check)

# Load Learners ----

load(file = "../data/1.data.processed.rda") 
load(file = "../data/2.gsign_kidney.rda")

#load("../data/2.data.full.count.rda")
# 
# 
# data <- full.data
# data$pp <- "raw"
# col_clini <- c("obs.time", "status", "age", "gender", "neoplasm", "metastasis", "ajcc.stage", "dataset", "pp")
# col_genes <- setdiff(colnames(data), col_clini)
# # To avoid inf exp on log-normalization
# data[, col_genes] <- data[, col_genes] + 1
# 
# # DESeq+vst ===========================================================================
# 
# vst1 <- function(countdata){
#   library(DESeq)
#   condition <- factor(rep("Tumour", nrow(countdata)))
#   countdata <- DESeq::newCountDataSet(as.matrix(t(countdata)),condition )
#   countdata <- DESeq::estimateSizeFactors( countdata )
#   cdsBlind <- DESeq::estimateDispersions( countdata, method="blind")
#   vstdata <- DESeq::varianceStabilizingTransformation( cdsBlind )
#   return(t(exprs(vstdata)))
# }
# 
# data.vst1 <- vst1(data[, col_genes])
# data.vst1 <- cbind(data[, col_clini], data.vst1)
# data.vst1$pp <- "DESeq1+vst"
# 
# # DESeq2+vst ===========================================================================
# 
# vst2 <- function(countdata, condition){
#   library(DESeq2)
#   colData <- data.frame(matrix( nrow=nrow(countdata) ))
#   colData$condition <- as.factor(condition)
#   rownames(colData) <- rownames(countdata)
#   countdata <-  DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(t(countdata)), 
#                                                colData = colData, design = ~ condition)
#   countdata <- DESeq2::estimateSizeFactors( countdata )
#   cdsBlind <- DESeq2::estimateDispersions( countdata, fitType="mean")
#   vstdata <- DESeq2::varianceStabilizingTransformation( cdsBlind , fitType="mean")
#   return( t(assay(vstdata)))
# }
# 
# data.vst2 <- vst2(data[, col_genes], data$dataset)
# 
# data.vst2 <- cbind(data[, col_clini], data.vst2)
# data.vst2$pp <- "DESeq2+vst"
# 
# # caret+rang ===========================================================================
# data.range <- data
# 
# # calculate the pre-process parameters from the dataset
# #preprocessParams <- preProcess(data[data$dataset == "TCGA",col_genes], method=c("BoxCox"), fudge=1)
# preprocessParams <- preProcess(data[data$dataset == "TCGA",col_genes], method=c("range"))
# 
# # transform the dataset using the parameters
# data.range[data.range$dataset == "TCGA",col_genes] <- predict(preprocessParams, data[data$dataset == "TCGA",col_genes])
# data.range[data.range$dataset == "ICGC",col_genes] <- predict(preprocessParams, data[data$dataset == "ICGC", col_genes])
# data.range$pp <- "caret+range"
# 
# # center+scale+nzv ===========================================================================
# data.csn <- data
# 
# # calculate the pre-process parameters from the dataset
# preprocessParams <- preProcess(data[data$dataset == "TCGA",col_genes], method=c("center","scale","nzv"))
# 
# # transform the dataset using the parameters
# data.csn[data.csn$dataset == "TCGA",col_genes] <- predict(preprocessParams, data[data$dataset == "TCGA",col_genes])
# data.csn[data.csn$dataset == "ICGC",col_genes] <- predict(preprocessParams, data[data$dataset == "ICGC", col_genes])
# data.csn$pp  <- "center+scale+nzv"
# 
# # BoxCox ===========================================================================
# data.bc <- data
# 
# # calculate the pre-process parameters from the dataset
# preprocessParams <- preProcess(data[data$dataset == "TCGA",col_genes], method=c("BoxCox"), fudge=1)
# 
# # transform the dataset using the parameters
# data.bc[data.bc$dataset == "TCGA",col_genes] <- predict(preprocessParams, data[data$dataset == "TCGA",col_genes])
# data.bc[data.bc$dataset == "ICGC",col_genes] <- predict(preprocessParams, data[data$dataset == "ICGC", col_genes])
# data.bc$pp  <- "BoxCox"

#save(data, data.bc, data.csn, data.range, data.vst1, data.vst2, file = "../data/1.norms.processed.rda", compress = T)

load("../data/1.norms.processed.rda")

rm(data)

  data0 <- rbind(data.bc, data.csn, data.range, data.vst1, data.vst2)

  data.m20 <- data0[data0$obs.time > 20, ]
  data.m20$pp <- paste0(data.m20$pp, "+>20")
  data.m20l10y <- data0[data0$obs.time > 20 & data0$obs.time < 3650, ]
  data.m20l10y$pp <- paste0(data.m20l10y$pp, "+>20d&<10y")
  data.m20l9y <- data0[data0$obs.time > 20 & data0$obs.time < 3285, ]
  data.m20l9y$pp <- paste0(data.m20l9y$pp, "+>20d&<9y")
  data.m20l8y <- data0[data0$obs.time > 20 & data0$obs.time < 2920, ]
  data.m20l8y$pp <- paste0(data.m20l8y$pp, "+>20d&<8y")
  data.m20l7y <- data0[data0$obs.time > 20 & data0$obs.time < 2555, ]
  data.m20l7y$pp <- paste0(data.m20l7y$pp, "+>20d&<7y")
  # data.m20l6y <- data0[data0$obs.time > 20 & data0$obs.time < 2190, ]
  # data.m20l6y$pp <- paste0(data.m20l6y$pp, "+>20d&<6y")
  # data.m20l5y <- data0[data0$obs.time > 20 & data0$obs.time < 1825, ]
  # data.m20l5y$pp <- paste0(data.m20l5y$pp, "+>20d&<6y")

  data.m30 <- data0[data0$obs.time > 30, ]
  data.m30$pp <- paste0(data.m30$pp, "+>30")
  data.m30l10y <- data0[data0$obs.time > 30 & data0$obs.time < 3650, ]
  data.m30l10y$pp <- paste0(data.m30l10y$pp, "+>30d&<10y")
  data.m30l9y <- data0[data0$obs.time > 30 & data0$obs.time < 3285, ]
  data.m30l9y$pp <- paste0(data.m30l9y$pp, "+>30d&<9y")
  data.m30l8y <- data0[data0$obs.time > 30 & data0$obs.time < 2930, ]
  data.m30l8y$pp <- paste0(data.m30l8y$pp, "+>30d&<8y")
  data.m30l7y <- data0[data0$obs.time > 30 & data0$obs.time < 2555, ]
  data.m30l7y$pp <- paste0(data.m30l7y$pp, "+>30d&<7y")
  # data.m30l6y <- data0[data0$obs.time > 30 & data0$obs.time < 2190, ]
  # data.m30l6y$pp <- paste0(data.m30l6y$pp, "+>30d&<6y")
  # data.m30l5y <- data0[data0$obs.time > 30 & data0$obs.time < 1825, ]
  # data.m30l5y$pp <- paste0(data.m30l5y$pp, "+>30d&<6y")

  data.m60 <- data0[data0$obs.time > 60, ]
  data.m60$pp <- paste0(data.m60$pp, "+>60")
  data.m60l10y <- data0[data0$obs.time > 60 & data0$obs.time < 3650, ]
  data.m60l10y$pp <- paste0(data.m60l10y$pp, "+>60d&<10y")
  data.m60l9y <- data0[data0$obs.time > 60 & data0$obs.time < 3285, ]
  data.m60l9y$pp <- paste0(data.m60l9y$pp, "+>60d&<9y")
  data.m60l8y <- data0[data0$obs.time > 60 & data0$obs.time < 2930, ]
  data.m60l8y$pp <- paste0(data.m60l8y$pp, "+>60d&<8y")
  data.m60l7y <- data0[data0$obs.time > 60 & data0$obs.time < 2555, ]
  data.m60l7y$pp <- paste0(data.m60l7y$pp, "+>60d&<7y")
  # data.m60l6y <- data0[data0$obs.time > 60 & data0$obs.time < 2190, ]
  # data.m60l6y$pp <- paste0(data.m60l6y$pp, "+>60d&<6y")
  # data.m60l5y <- data0[data0$obs.time > 60 & data0$obs.time < 1825, ]
  # data.m60l5y$pp <- paste0(data.m60l5y$pp, "+>60d&<6y")

  data.m120 <- data0[data0$obs.time > 120, ]
  data.m120$pp <- paste0(data.m120$pp, "+>120")
  data.m120l10y <- data0[data0$obs.time > 120 & data0$obs.time < 3650, ]
  data.m120l10y$pp <- paste0(data.m120l10y$pp, "+>120d&<10y")
  data.m120l9y <- data0[data0$obs.time > 120 & data0$obs.time < 3285, ]
  data.m120l9y$pp <- paste0(data.m120l9y$pp, "+>120d&<9y")
  data.m120l8y <- data0[data0$obs.time > 120 & data0$obs.time < 2930, ]
  data.m120l8y$pp <- paste0(data.m120l8y$pp, "+>120d&<8y")
  data.m120l7y <- data0[data0$obs.time > 120 & data0$obs.time < 2555, ]
  data.m120l7y$pp <- paste0(data.m120l7y$pp, "+>120d&<7y")
  # data.m120l6y <- data0[data0$obs.time > 120 & data0$obs.time < 2190, ]
  # data.m120l6y$pp <- paste0(data.m120l6y$pp, "+>120d&<6y")
  # data.m120l5y <- data0[data0$obs.time > 120 & data0$obs.time < 1825, ]
  # data.m120l5y$pp <- paste0(data.m120l5y$pp, "+>120d&<6y")

  rm(data0)
  #ls(pattern = "(^data.)"))
  data <- rbind(data.bc, data.csn, data.m120, data.m120l10y, data.m120l7y, data.m120l8y, data.m120l9y, data.m20, data.m20l10y, data.m20l7y, data.m20l8y, data.m20l9y, data.m30, data.m30l10y, data.m30l7y, data.m30l8y, data.m30l9y, data.m60, data.m60l10y, data.m60l7y, data.m60l8y, data.m60l9y, data.range, data.vst1, data.vst2)
  
print(table(data$pp))

rm(list=setdiff(ls(), c("data", "gsign_kidney")))

bmr_table_pp <- data.frame(matrix(0, nrow=0, ncol=7))

colnames(bmr_table_pp) <- c("task_id", "learner_id", "resampling_id", "iters", "surv.harrell_c", "surv.uno_auc", "pp")

for(pp in unique(data$pp)){ 
  
print(pp)
  
data_Task <- data[data$pp %in% pp, ]

# Create TaskSurv  ===========================================================================
for(s in c(1:nrow(gsign_kidney))){ 
    sign <- gsign_kidney$id[s]
    cols <- colnames(data) %in% unlist(c("obs.time", "status", intersect(colnames(data), unlist(gsign_kidney$signature[s]))))
    assign(sign,
           TaskSurv$new(id = sign, backend = data_Task[, cols],
                        time = "obs.time", event = "status", type = "right")
    ) 
  }
  rm(s, sign, cols)
  
  list.surv.task <- list(filt.gbm, fs.mRMR.status, sign.10)
  
  # list.surv.task <- list(filt.gbm, filt.rpart, filt.xgboost, fs.boruta,
  #                        fs.genetic, fs.mlr.rfe, fs.mRMR, fs.mRMR.status, fs.mRMR.time,
  #                        sign.1,  sign.2,   sign.3,   sign.4,
  #                        sign.5,   sign.6,   sign.7,   sign.8,   sign.9,
  #                        sign.10,  sign.11,  sign.12,  sign.13,  sign.14,
  #                        sign.15,  sign.16,  sign.17,  sign.18,  sign.19)
  
  rm(list=ls(pattern = "(^filt.)|(^fs.)|(^sign.)"))
  
  
  # Load metrics ----
  
  bmr_table <- data.frame(matrix(0, nrow=0, ncol=6))
  colnames(bmr_table) <- c("task_id", "learner_id", "resampling_id", "iters", "surv.harrell_c", "surv.uno_auc")
  
  measures = list(
    msr("surv.cindex", id = "surv.harrell_c", predict_sets =  "test"),
    msr("surv.uno_auc", id = "surv.uno_auc", predict_sets =  "test")
  )
  
  measure = msr("surv.cindex")
  
  # Execute benchmark ----
  for(i in c(1:length(list.surv.task))){
    list.surv.task[[i]]$print() 
    
    #  Load Tunning ---
    tuner = tnr("random_search")
    terminator = trm("evals", n_evals = 20)
    
    # glmnet ----
    search_space = ps(
      s = p_dbl(lower = 0.001, upper = 0.1),
      alpha = p_dbl(lower = 0, upper = 1)
    )
    lrn.glmnet <- AutoTuner$new(
      learner = lrn("surv.glmnet", id="glmnet"),
      resampling = rsmp("holdout"),
      measure = measure,
      search_space = search_space,
      terminator = terminator,
      tuner = tuner
    )
    
    # xgb ----
    search_space = ps(
      alpha = p_dbl(lower = 0, upper = 1),
      eta = p_dbl(lower = 0.001, upper = 0.1),
      nrounds = p_int(lower = 1, upper = 1000)#,
      #max_depth = p_int(lower = 3, upper = 10), # -gbtree
      #subsample = p_dbl(lower = 0.3, upper = 0.8), # -gbtree
      #gamma = p_int(lower = 0, upper = 5) # -gbtree
    )
    lrn.xgb <- AutoTuner$new(
      learner = lrn("surv.xgboost", id="xgb", booster="gblinear"),
      resampling = rsmp("holdout"),
      measure = measure,
      search_space = search_space,
      terminator = terminator,
      tuner = tuner
    )
    
    # glmboost ----
    # lrn("surv.glmboost")$param_set
    search_space = ps(
      nu = p_dbl(lower = 0, upper = 1),
      sigma = p_dbl(lower = 0, upper = 1)
    )
    lrn.glmboost <- AutoTuner$new(
      learner = lrn("surv.glmboost", id="glmboost", family="cindex"),
      resampling = rsmp("holdout"),
      measure = measure,
      search_space = search_space,
      terminator = terminator,
      tuner = tuner
    )
    
    # surv.gbm ---- 
    # lrn("surv.gbm")$param_set
    search_space = ps(
      n.trees = p_int(lower = 300, upper = 1500),
      interaction.depth = p_int(lower = 3, upper = 8),
      shrinkage = p_dbl(lower = 0.0001, upper = 0.1)
    )
    lrn.gbm <- AutoTuner$new(
      learner = lrn("surv.gbm", id="gbm"),
      resampling = rsmp("holdout"),
      measure = measure,
      search_space = search_space,
      terminator = terminator,
      tuner = tuner
    )

    lrn.glmnet.lasso <- lrn("surv.glmnet", id="glmnet.lasso", alpha = 1)
    lrn.glmnet.ridge <- lrn("surv.glmnet", id="glmnet.ridge", alpha = 0)
    lrn.glmnet.en <- lrn("surv.glmnet", id="glmnet.en", alpha = 0.5)
    #learners = c("surv.parametric", "surv.coxboost", "surv.xgboost", "surv.coxph", "surv.blackboost", "surv.flexible", "surv.gbm", "surv.mboost", "surv.glmboost","surv.gamboost")
    #learners = lrns(c("surv.coxph", "surv.blackboost", "surv.gbm", "surv.ranger"))
    learners = c("surv.coxph", "surv.blackboost")
    learners = lapply(learners, lrn, predict_type = "lp", predict_sets = c("test"))
    
    resampling = rsmp("custom")
    resampling$instantiate(list.surv.task[[i]],
                           train = list(c(1:nrow(data_Task))[data_Task$dataset == "TCGA"]),
                           test = list(c(1:nrow(data_Task))[data_Task$dataset == "ICGC"])
    )
    
    # Set multicore or multisession  ----
    if (future::supportsMulticore()) {
      print("supportsMulticore")
      print(future::availableCores())
      future::plan(list(
        future::tweak(future::multicore, workers = 2),
        future::tweak(future::multisession, workers = 4)
      ))
    } else {
      print("notSupportsMulticore")
      print(future::availableCores())
      # 9 cores in outer and 4 inner
      future::plan(list(
        future::tweak(future::multisession, workers = 2),
        future::tweak(future::multisession, workers = 4)
      ))
    }
    
    # create a benchmark design object
    set.seed(1)
    design = benchmark_grid(list.surv.task[i], 
                            c(lrn.xgb, lrn.glmnet.lasso, lrn.glmnet.ridge, lrn.glmnet.en, lrn.glmboost, lrn.gbm, learners),
                            resampling = resampling
    )
    
    try({
      bmr = benchmark(design) 
      bmr <- bmr$aggregate(measures) %>%
        dplyr::select(colnames(bmr_table)) %>%
        as.data.frame()
      
      print(bmr)
      bmr_table <- rbind(bmr_table, bmr)

    }, silent = T)
    bmr <- NULL
  }

  bmr_table$pp <- pp
  bmr_table_pp <- rbind(bmr_table_pp, bmr_table[,colnames(bmr_table_pp)])
  
  #print(bmr_table)
}

library(knitr)
library(rmarkdown)

bmr_table_pp <- bmr_table_pp %>%
  dplyr::arrange(desc(surv.harrell_c))

table <- kable(bmr_table_pp, format="markdown")
save(bmr_table_pp, file = "../data/job_benchmark_normalizations.rda")

# 
# title <- "#### Supplementary Table 2. Survival Analysis Validation KIRC-TCGA and ICGC"
# 
# text <- "Train dataset KIRC-TCGA and test with ICGC.\n \n"
# 
# cat(title, text, table, sep="\n", file="~/public_html/results/train.tcga.test.icgc.Rmd")
# render("~/public_html/results/train.tcga.test.icgc.Rmd", output_format = "html_document")




# bmr_table_pp %>% group_by(pp) %>% summarise(ci=mean(surv.harrell_c)) %>% arrange(desc( ci))
# # A tibble: 125 x 2
# pp                            ci
# <chr>                      <dbl>
#   1 DESeq1+vst+>60d&<7y        0.619
# 2 DESeq2+vst+>60d&<7y        0.616
# 3 center+scale+nzv+>60d&<6y  0.615
# 4 BoxCox+>60d&<7y            0.614
# 5 center+scale+nzv+>120d&<6y 0.614
# 6 caret+range+>60d&<7y       0.611
# 7 center+scale+nzv+>60d&<10y 0.610
# 8 center+scale+nzv+>60       0.606
# 9 center+scale+nzv+>20d&<6y  0.606
# 10 center+scale+nzv+>20d&<7y  0.606


# > bmr_table_pp %>% group_by(pp) %>% summarise(uno_auc=mean(surv.uno_auc)) %>% arrange(desc( uno_auc))
# # A tibble: 125 x 2
# pp                  uno_auc
# <chr>                 <dbl>
#   1 DESeq1+vst+>60d&<7y   0.650
# 2 BoxCox+>60d&<7y       0.648
# 3 DESeq2+vst+>60d&<7y   0.648
# 4 BoxCox                0.633
# 5 DESeq1+vst+>30d&<7y   0.632
# 6 DESeq2+vst+>30d&<7y   0.631
# 7 DESeq1+vst+>60        0.630
# 8 BoxCox+>30d&<7y       0.629
# 9 BoxCox+>20d&<7y       0.629
# 10 DESeq1+vst+>20d&<7y   0.628


