## Load packages
library(dplyr)

## Load in genotype data
UKB_genotype <- readRDS("/Volumes/GenScotDepression/data/ukb/phenotypes/fields/2021-04-phenotypes-ukb44797/Genotypes.rds")

## load in final QC fam file 

FAM <- read.table("~/Desktop/PhD/projects/UKBCRPImagingPRS/resources/ukb_imp_autosome_v3_imaging_rmSNP.QC.fam", header = F)

## rename FID and IID columns
FAM <- FAM %>%
  rename(FID = V1,
         IID = V2)

## Retain genotype measurement batch and principal component
## f.22000 = Genotype measurement batch
## f.22009 = Genetic principal components
UKB_genotype_batch_PC <- cbind(UKB_genotype$f.eid,
                                   UKB_genotype[grepl("f.22000.", names(UKB_genotype))],
                                   UKB_genotype[grepl("f.22009.", names(UKB_genotype))])


## Identify genotype array based on genotype batch
## Identify participants genotyped using the UK BiLEVE array, coded as -11....-1 (n = ~50,000)
## BiLEVE array = 1
## Create new column to store array data
UKB_genotype_batch_PC$array <- NA

nrow(UKB_genotype_batch_PC[(UKB_genotype_batch_PC$f.22000.0.0 < 0 & UKB_genotype_batch_PC$f.22000.0.0 > -11),])

UKB_genotype_batch_PC$array[(UKB_genotype_batch_PC$f.22000.0.0 < 0 & UKB_genotype_batch_PC$f.22000.0.0 > -11)] <- 1

## Identify participants genotyped using the UKB Axiom array, coded as 1....95 (n = ~450,000)
## Axiom array = 2

nrow(UKB_genotype_batch_PC[(UKB_genotype_batch_PC$f.22000.0.0 > 0 & UKB_genotype_batch_PC$f.22000.0.0 < 95),])

UKB_genotype_batch_PC$array[(UKB_genotype_batch_PC$f.22000.0.0 > 0 & UKB_genotype_batch_PC$f.22000.0.0 < 95)] <- 2


## change column names
names(UKB_genotype_batch_PC)

UKB_genotype_batch_PC <- UKB_genotype_batch_PC %>%
  rename(IID = "UKB_genotype$f.eid",
         PC1 = f.22009.0.1,
         PC2 = f.22009.0.2,
         PC3 = f.22009.0.3,
         PC4 = f.22009.0.4,
         PC5 = f.22009.0.5,
         PC6 = f.22009.0.6,
         PC7 = f.22009.0.7,
         PC8 = f.22009.0.8,
         PC9 = f.22009.0.9,
         PC10 = f.22009.0.10)

## limit covar file to QCd individuals with imaging data
UKB_CRP_imaging.cov <- UKB_genotype_batch_PC %>%
  filter(IID %in% FAM$IID)

## limit covar file to PCs and array and add FID column
UKB_CRP_imaging.cov$FID <- UKB_CRP_imaging.cov$IID
UKB_CRP_imaging.cov <- UKB_CRP_imaging.cov[c("FID", "IID", paste0("PC",1:10), "array")]

## write to resource directory
write.table(UKB_CRP_imaging.cov, "~/Desktop/PhD/projects/UKBCRPImagingPRS/resources/UKB_CRP_imaging.cov", row.names = FALSE, quote = FALSE)

