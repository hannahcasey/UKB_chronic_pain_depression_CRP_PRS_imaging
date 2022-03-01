## Load data

BEST_SBayesR_CRP <- read.table("~/Desktop/PhD/projects/UKBCRPImagingPRS/output/PRS_analysis/CRP_PRS_SBayesR_noMHC.best", header = T)
BEST_SBayesR_CRP_MHC <- read.table("~/Desktop/PhD/projects/UKBCRPImagingPRS/output/PRS_analysis/CRP_PRS_SBayesR.best", header = T)
BEST_SBayesR_CRP_MHC <- read.table("~/Desktop/PhD/projects/UKBCRPImagingPRS/output/PRS_analysis/CRP_PRS", header = T)


UKB_assay_results <- readRDS("/Volumes/GenScotDepression/data/ukb/phenotypes/fields/2020-09-phenotypes-ukb43743/AssayResults.rds")
## f.30710 = C-reactive protein
## Extract ID column and columns pretaining to CRP measurment
UKB_assay_results_filtered <- cbind(UKB_assay_results[,1], UKB_assay_results[grepl("f.30710.", names(UKB_assay_results))])

## rename FID and CRP columns
UKB_assay_results_filtered <- UKB_assay_results_filtered %>%
  rename(FID = "UKB_assay_results[, 1]",
         CRP = f.30710.0.0)

## Join dataframes
CRP_serum_PRS <- inner_join(BEST_SBayesR_CRP, UKB_assay_results_filtered, by = "FID")

## Get correlation
cor(CRP_serum_PRS$CRP, CRP_serum_PRS$PRS, method = "spearman", use = "complete.obs")

## Log transform CRP
CRP_serum_PRS$`Log(CRP)` <- log(CRP_serum_PRS$CRP)

## Plot correlation
library("ggpubr")
ggscatter(CRP_serum_PRS, x = "Log(CRP)", y = "PRS", 
          add = "reg.line", conf.int = TRUE,  alpha = 0.03,
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Log(CRP mg/L)", ylab = "CRP PRS", title = "Correlation Between Serum CRP and CRP PRS")
