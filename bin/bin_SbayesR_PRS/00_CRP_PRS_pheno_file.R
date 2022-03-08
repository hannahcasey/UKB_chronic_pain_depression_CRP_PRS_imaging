## Load packages
library(dplyr)

## Load in assay data
UKBAssayResults <- readRDS("/Volumes/GenScotDepression/data/ukb/phenotypes/fields/2020-09-phenotypes-ukb43743/AssayResults.rds")
## f.30710 = C-reactive protein
## Extract ID column and columns pretaining to CRP measurment
UKBAssayResultsCRP <- cbind(UKBAssayResults[,1], UKBAssayResults[grepl("f.30710.", names(UKBAssayResults))])
## Remove redundant dataframe
rm(UKBAssayResults)

## rename FID and CRP columns
UKBAssayResultsCRP <- UKBAssayResultsCRP %>%
  rename(IID = "UKBAssayResults[, 1]",
         CRP = f.30710.0.0)


## load in final QC fam file 
FAM <- read.table("~/Desktop/PhD/projects/UKB_inflammation_imaging/resources//ukb_imp_v3.qc.fam", header = F)

## rename FID and IID columns
FAM <- FAM %>%
  rename(FID = V1,
         IID = V2)

## limit covar file to QCd individuals with imaging data
UKB_CRP_imaging.pheno <- UKBAssayResultsCRP %>%
  filter(IID %in% FAM$IID)

## limit covar file to PCs and array and add FID column
UKB_CRP_imaging.pheno$FID <- UKB_CRP_imaging.pheno$IID
UKB_CRP_imaging.pheno <- UKB_CRP_imaging.pheno[c("FID", "IID", "CRP")]

## remove participants w/o CRP measurment
UKB_CRP_imaging.pheno <- na.omit(UKB_CRP_imaging.pheno)

## remove participants w/o CRP measurment
UKB_CRP_log_imaging.pheno <- UKB_CRP_imaging.pheno %>%
  mutate(CRP_log = log(CRP)) %>%
  select(FID, IID, CRP_log)

## write to resource directory
write.table(UKB_CRP_imaging.pheno, "~/Desktop/PhD/projects/UKB_inflammation_imaging/resources/UKB_CRP_imaging.pheno", row.names = FALSE, quote = FALSE)
write.table(UKB_CRP_log_imaging.pheno, "~/Desktop/PhD/projects/UKB_inflammation_imaging/resources/UKB_CRP_log_imaging.pheno", row.names = FALSE, quote = FALSE)
