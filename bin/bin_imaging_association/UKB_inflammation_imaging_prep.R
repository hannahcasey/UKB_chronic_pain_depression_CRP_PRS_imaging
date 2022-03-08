## Load in Packages
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)
library(rebus) #to build regular expressions


## Load in data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Load in CRP polygenic risk scores for UKB participants
## More information on PRS calculation: <github link>
UKB_CRP_PRS <- read.table("~/Desktop/PhD/projects/UKB_inflammation_imaging/output/SBayesR_CRP/CRP_log_PRS_noMHC_SBayesR.best", header = T)

## Rename ID column
UKB_CRP_PRS <- UKB_CRP_PRS %>%
  rename(f.eid = FID)

## Load in serum CRP measures
## Load in assay data
UKB_assay_results <- readRDS("/Volumes/GenScotDepression/data/ukb/phenotypes/fields/2020-09-phenotypes-ukb43743/AssayResults.rds")
## f.30710 = C-reactive protein
## Extract ID column and columns pretaining to CRP measurment
UKB_CRP_serum <- cbind(UKB_assay_results[,1], UKB_assay_results[grepl("f.30710.", names(UKB_assay_results))])

## rename FID and CRP columns
UKB_CRP_serum <- UKB_CRP_serum %>%
  rename(f.eid = "UKB_assay_results[, 1]",
         CRP = f.30710.0.0)

## Remove redundant dataframes 
rm(UKB_assay_results)


## Load in NMR metabolomics data
UKB_NMR <- read.table("/Volumes/GenScotDepression/data/ukb/phenotypes/fields/2021-10-nmr-metabolomics-ukb48936/NMRMetabolomics.tsv.gz", header = TRUE)

## Extract ID column and columns pretaining to GlycA measurment
UKB_glycA <- cbind(UKB_NMR[,1], UKB_NMR[grepl("f.23480.", names(UKB_NMR))])

## Rename ID column
UKB_glycA <- UKB_glycA %>%
  rename(f.eid = "UKB_NMR[, 1]",
         glycA_1 = f.23480.0.0,
         glycA_2 = f.23480.1.0)

## Remove redundant dataframes 
rm(UKB_NMR)

## Load in covariate information:
## f.21022 = Age at recruitment -
## f.31 = Sex -
## f.21001 = BMI -
## f.54 = assessment centre -
## Volume of EstimatedTotalIntraCranial (whole brain)

## Load in baseline characteristics - sex
UKB_baseline_characteristics <- readRDS("/Volumes/GenScotDepression/data/ukb/phenotypes/fields/2020-09-phenotypes-ukb43743/BaselineCharacteristics.rds")

## Extract sex and sex information
UKB_sex_age <- cbind(UKB_baseline_characteristics[,1],
                                               UKB_baseline_characteristics[grepl("f.31.0.0", names(UKB_baseline_characteristics))],
                                               UKB_baseline_characteristics[grepl("f.21022.0.0", names(UKB_baseline_characteristics))])

## Rename ID and sex column in new filtered dataframe
UKB_sex_age <- UKB_sex_age %>%
  rename(f.eid = "UKB_baseline_characteristics[, 1]",
         sex = f.31.0.0,
         age = f.21022.0.0)

## Create new variable, age^2
UKB_sex_age$age_squared <- (UKB_sex_age$age^2)

## Remove redundant dataframes 
rm(UKB_baseline_characteristics)

## Load in physical measures data - BMI
UKB_physical_measures <- readRDS("/Volumes/GenScotDepression/data/ukb/phenotypes/fields/2021-04-phenotypes-ukb44797/PhysicalMeasures.rds")

## Extract BMI data
UKB_BMI <- cbind(UKB_physical_measures[,1], UKB_physical_measures[grepl("f.21001.0.0", names(UKB_physical_measures))])

## Rename ID and BMI column in new filtered dataframe
UKB_BMI <- UKB_BMI %>%
  rename(f.eid = "UKB_physical_measures[, 1]",
         BMI = f.21001.0.0)

## Remove redundant dataframes 
rm(UKB_physical_measures)

## Load in recruitment data - assessment centre
UKB_recruitment <- readRDS("/Volumes/GenScotDepression/data/ukb/phenotypes/fields/2021-04-phenotypes-ukb44797/Recruitment.rds")
UKB_assessment_centre_first_imaging <- cbind(UKB_recruitment[,1],
                                  UKB_recruitment[grepl("f.54.2.0", names(UKB_recruitment))])

## Rename ID and 2nd instance assessment centre column in new filtered dataframe
UKB_assessment_centre_first_imaging <- UKB_assessment_centre_first_imaging %>%
  rename(f.eid = "UKB_recruitment[, 1]",
         assessment_centre_first_imaging = f.54.2.0)

## Remove redundant dataframes 
rm(UKB_recruitment)

## Load in structural imaging data
UKB_imaging <- readRDS("/Volumes/GenScotDepression/data/ukb/phenotypes/fields/2020-02-imaging-ukb40531/Imaging.rds")

## Extract total ICV data - sum of volume of white matter(f.25008), volume of grey matter (f.25006) and volume of ventricular cerebrospinal fluid (25004)
UKB_imaging_ICV <- cbind(UKB_imaging[,1],
                                     UKB_imaging[grepl("f.25008.", names(UKB_imaging))],
                                     UKB_imaging[grepl("f.25006.", names(UKB_imaging))],
                                     UKB_imaging[grepl("f.25004.", names(UKB_imaging))])
## Add white matter, grey matter and cerebrospinal fluid volume to calcuate total ICV
UKB_imaging_ICV$ICV <- UKB_imaging_ICV$f.25008.2.0 + UKB_imaging_ICV$f.25006.2.0 + UKB_imaging_ICV$f.25004.2.0


## Rename ID in new filtered dataframe
UKB_imaging_ICV <- UKB_imaging_ICV %>%
  rename(f.eid = "UKB_imaging[, 1]")

## Create regular expression matchinng number ranges identifying field IDs containing regional cortical volumes, regional FA measures, regional MD measures and subcortical volumes
rx_cort <- number_range(26721, 26922)
rx_FA <- number_range(25488, 25514)
rx_MA <- number_range(25515, 25541)
rx_subcort <- number_range(25011, 25024)

UKB_imaging_small <- cbind(UKB_imaging[,1],
                               UKB_imaging[grepl(rx_cort, names(UKB_imaging))],
                               UKB_imaging[grepl(rx_FA, names(UKB_imaging))],
                               UKB_imaging[grepl(rx_MA, names(UKB_imaging))],
                               UKB_imaging[grepl(rx_subcort, names(UKB_imaging))])

## Rename ID column in new filtered dataframe
UKB_imaging_small <- UKB_imaging_small %>%
  rename(f.eid = "UKB_imaging[, 1]")

## Remove redundant dataframes 
rm(UKB_imaging)


## Combine ID, covariates, CRP PRS and DKT structural measures into single dataframe
UKB_inflammation_imaging_covariates <- left_join(UKB_CRP_serum, UKB_CRP_PRS, by = "f.eid")
UKB_inflammation_imaging_covariates <- left_join(UKB_inflammation_imaging_covariates, UKB_glycA, by = "f.eid")
UKB_inflammation_imaging_covariates <- left_join(UKB_inflammation_imaging_covariates, UKB_sex_age, by = "f.eid")
UKB_inflammation_imaging_covariates <- left_join(UKB_inflammation_imaging_covariates, UKB_BMI, by = "f.eid")
UKB_inflammation_imaging_covariates <- left_join(UKB_inflammation_imaging_covariates, UKB_assessment_centre_first_imaging, by = "f.eid")
UKB_inflammation_imaging_covariates <- left_join(UKB_inflammation_imaging_covariates, UKB_imaging_ICV, by = "f.eid")
UKB_inflammation_imaging_covariates <- left_join(UKB_inflammation_imaging_covariates, UKB_imaging_small, by = "f.eid")

## Calculate Global and lobar meaurments for imaging features
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## General Cortical Volume
## General Cortical Area
## General Cortical Thickness


## Lobar Cortical Volumes

## Frontal lobe FIDs:
## Superior frontal gyrus = 26815 + 26916
## rostralmiddlefrontal = 26814 + 26915
## Caudalmiddlefrontal = 26791 + 26892
## Pars orbitalis = 26806 + 26907
## pars triangularis = 26807 + 26908
## pars opercularis = 26805 + 26906
## Frontal pole = 26819 + 26920
## Lateral orbitofrontal = 26799 + 26900
## medial orbitofrontal = 26801 + 26902
## Precentral gyrus = 26811 + 26912
## paracentral cortex = 26804 + 26905

## Sum volume of individual structures in frontal lobe
UKB_inflammation_imaging_covariates$frontal_lobe_volume <- UKB_inflammation_imaging_covariates$f.26815.2.0 + UKB_inflammation_imaging_covariates$f.26916.2.0 + UKB_inflammation_imaging_covariates$f.26814.2.0 +
  UKB_inflammation_imaging_covariates$f.26915.2.0 + UKB_inflammation_imaging_covariates$f.26791.2.0 + UKB_inflammation_imaging_covariates$f.26892.2.0 + UKB_inflammation_imaging_covariates$f.26806.2.0 +
  UKB_inflammation_imaging_covariates$f.26907.2.0 + UKB_inflammation_imaging_covariates$f.26807.2.0 + UKB_inflammation_imaging_covariates$f.26908.2.0 + UKB_inflammation_imaging_covariates$f.26805.2.0 + 
  UKB_inflammation_imaging_covariates$f.26906.2.0 + UKB_inflammation_imaging_covariates$f.26819.2.0 + UKB_inflammation_imaging_covariates$f.26920.2.0 + UKB_inflammation_imaging_covariates$f.26799.2.0 +
  UKB_inflammation_imaging_covariates$f.26900.2.0 + UKB_inflammation_imaging_covariates$f.26801.2.0 + UKB_inflammation_imaging_covariates$f.26902.2.0 + UKB_inflammation_imaging_covariates$f.26811.2.0 + 
  UKB_inflammation_imaging_covariates$f.26912.2.0 + UKB_inflammation_imaging_covariates$f.26804.2.0 + UKB_inflammation_imaging_covariates$f.26905.2.0

## Temporal lobe FIDs:
## Insula= 26821 + 26922
## Superior temporal = 26817 + 26918
## transverse temporal = 26820 + 26921
## banks STS = 	26789 + 26890
## Middle temporal gyrus = 26802 + 26903
## Inferior temporal gyrus = 26796 + 26897
## Fusiform = 26794 + 26895
## parahippocampal = 26803 + 26904
## entorhinal = 26793 + 26894

## Sum volume of individual structures in temporal lobe
UKB_inflammation_imaging_covariates$temporal_lobe_volume <- UKB_inflammation_imaging_covariates$f.26821.2.0 + UKB_inflammation_imaging_covariates$f.26922.2.0 + UKB_inflammation_imaging_covariates$f.26817.2.0 +
  UKB_inflammation_imaging_covariates$f.26918.2.0 + UKB_inflammation_imaging_covariates$f.26820.2.0 + UKB_inflammation_imaging_covariates$f.26921.2.0 + UKB_inflammation_imaging_covariates$f.26789.2.0 +
  UKB_inflammation_imaging_covariates$f.26890.2.0 + UKB_inflammation_imaging_covariates$f.26802.2.0 + UKB_inflammation_imaging_covariates$f.26903.2.0 + UKB_inflammation_imaging_covariates$f.26796.2.0 +
  UKB_inflammation_imaging_covariates$f.26897.2.0 + UKB_inflammation_imaging_covariates$f.26794.2.0 + UKB_inflammation_imaging_covariates$f.26895.2.0 + UKB_inflammation_imaging_covariates$f.26803.2.0 +
  UKB_inflammation_imaging_covariates$f.26904.2.0 + UKB_inflammation_imaging_covariates$f.26793.2.0 + UKB_inflammation_imaging_covariates$f.26894.2.0

## Parietal lobe:
## Postcentral gyrus = 26809 + 26910
## paracentral cortex = 26804 + 26905
## Superior parietal cortex = 26816 + 26917
## Inferior parietal cortex = 26795 + 26896
## Supramarginal gyrus = 26818 + 26919
## Precuneus = 26812 + 26913

## Sum volume of individual structures in parietal lobe
UKB_inflammation_imaging_covariates$parietal_lobe_volume <- UKB_inflammation_imaging_covariates$f.26809.2.0 + UKB_inflammation_imaging_covariates$f.26910.2.0 + UKB_inflammation_imaging_covariates$f.26804.2.0 +
  UKB_inflammation_imaging_covariates$f.26905.2.0 + UKB_inflammation_imaging_covariates$f.26816.2.0 + UKB_inflammation_imaging_covariates$f.26917.2.0 + UKB_inflammation_imaging_covariates$f.26795.2.0 +
  UKB_inflammation_imaging_covariates$f.26896.2.0 + UKB_inflammation_imaging_covariates$f.26818.2.0 + UKB_inflammation_imaging_covariates$f.26919.2.0 + UKB_inflammation_imaging_covariates$f.26812.2.0 +
  UKB_inflammation_imaging_covariates$f.26913.2.0

## Occipital lobe:
## Lateral occipital cortex = 26798, 26899
## Cuneus = 26792, 26893
## Pericalcarine cortex = 26808, 26909
## Lingual gyrus = 26800, 26901

UKB_inflammation_imaging_covariates$occipital_lobe_volume <- UKB_inflammation_imaging_covariates$f.26798.2.0 + UKB_inflammation_imaging_covariates$f.26899.2.0 + UKB_inflammation_imaging_covariates$f.26792.2.0 +
  UKB_inflammation_imaging_covariates$f.26893.2.0 + UKB_inflammation_imaging_covariates$f.26808.2.0 + UKB_inflammation_imaging_covariates$f.26909.2.0 + UKB_inflammation_imaging_covariates$f.26800.2.0 +
  UKB_inflammation_imaging_covariates$f.26901.2.0

## Cingulate lobe:
## Rostral ACC = 26813, 26914
## Caudal ACC = 26790, 26891
## Posterior cingulate cortex = 26810 + 26911
## Cingulate isthmus = 26797 + 26898

UKB_inflammation_imaging_covariates$cingulate_lobe_volume <- UKB_inflammation_imaging_covariates$f.26813.2.0 + UKB_inflammation_imaging_covariates$f.26914.2.0 + UKB_inflammation_imaging_covariates$f.26790.2.0 +
  UKB_inflammation_imaging_covariates$f.26891.2.0 + UKB_inflammation_imaging_covariates$f.26810.2.0 + UKB_inflammation_imaging_covariates$f.26911.2.0 + UKB_inflammation_imaging_covariates$f.26797.2.0 +
  UKB_inflammation_imaging_covariates$f.26898.2.0

## Sum lobar volumes to get global cortical volume
UKB_inflammation_imaging_covariates$global_cortical_volume <- UKB_inflammation_imaging_covariates$frontal_lobe_volume + UKB_inflammation_imaging_covariates$temporal_lobe_volume + UKB_inflammation_imaging_covariates$parietal_lobe_volume +
  UKB_inflammation_imaging_covariates$occipital_lobe_volume +UKB_inflammation_imaging_covariates$cingulate_lobe_volume 



## Calculate cortical lobar area measures
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Frontal lobe area FIDs:
## Superior frontal gyrus = 26748 + 26849
## rostralmiddlefrontal = 26747 + 26848
## Caudalmiddlefrontal = 26724 + 26825
## Pars orbitalis = 26739 + 26840
## pars triangularis = 26740 + 26841
## pars opercularis = 26805 + 26839
## Frontal pole = 26752 + 26853
## Lateral orbitofrontal = 26732 + 26833
## medial orbitofrontal = 26734 + 26835
## Precentral gyrus = 26744 + 26845
## paracentral cortex = 26737 + 26838

## Sum area of individual structures in frontal lobe
UKB_inflammation_imaging_covariates$frontal_lobe_area <- UKB_inflammation_imaging_covariates$f.26748.2.0 + UKB_inflammation_imaging_covariates$f.26849.2.0 + UKB_inflammation_imaging_covariates$f.26747.2.0 +
  UKB_inflammation_imaging_covariates$f.26848.2.0 + UKB_inflammation_imaging_covariates$f.26724.2.0 + UKB_inflammation_imaging_covariates$f.26825.2.0 + UKB_inflammation_imaging_covariates$f.26739.2.0 +
  UKB_inflammation_imaging_covariates$f.26840.2.0 + UKB_inflammation_imaging_covariates$f.26740.2.0 + UKB_inflammation_imaging_covariates$f.26841.2.0 + UKB_inflammation_imaging_covariates$f.26805.2.0 + 
  UKB_inflammation_imaging_covariates$f.26839.2.0 + UKB_inflammation_imaging_covariates$f.26752.2.0 + UKB_inflammation_imaging_covariates$f.26853.2.0 + UKB_inflammation_imaging_covariates$f.26732.2.0 +
  UKB_inflammation_imaging_covariates$f.26833.2.0 + UKB_inflammation_imaging_covariates$f.26734.2.0 + UKB_inflammation_imaging_covariates$f.26835.2.0 + UKB_inflammation_imaging_covariates$f.26744.2.0 + 
  UKB_inflammation_imaging_covariates$f.26845.2.0 + UKB_inflammation_imaging_covariates$f.26737.2.0 + UKB_inflammation_imaging_covariates$f.26838.2.0


## Temporal lobe area FIDs:
## Insula= 26754 + 26855
## Superior temporal = 26750 + 26851
## transverse temporal = 26753 + 26854
## banks STS = 	26722 + 26823
## Middle temporal gyrus = 26735 + 26836
## Inferior temporal gyrus = 26796 + 26897
## Fusiform = 26727 + 26828
## parahippocampal = 26736 + 26837
## entorhinal = 26726 + 26827

## Sum area of individual structures in temporal lobe
UKB_inflammation_imaging_covariates$temporal_lobe_area <- UKB_inflammation_imaging_covariates$f.26754.2.0 + UKB_inflammation_imaging_covariates$f.26855.2.0 + UKB_inflammation_imaging_covariates$f.26750.2.0 +
  UKB_inflammation_imaging_covariates$f.26851.2.0 + UKB_inflammation_imaging_covariates$f.26753.2.0 + UKB_inflammation_imaging_covariates$f.26854.2.0 + UKB_inflammation_imaging_covariates$f.26722.2.0 +
  UKB_inflammation_imaging_covariates$f.26823.2.0 + UKB_inflammation_imaging_covariates$f.26735.2.0 + UKB_inflammation_imaging_covariates$f.26836.2.0 + UKB_inflammation_imaging_covariates$f.26796.2.0 +
  UKB_inflammation_imaging_covariates$f.26897.2.0 + UKB_inflammation_imaging_covariates$f.26727.2.0 + UKB_inflammation_imaging_covariates$f.26828.2.0 + UKB_inflammation_imaging_covariates$f.26736.2.0 +
  UKB_inflammation_imaging_covariates$f.26837.2.0 + UKB_inflammation_imaging_covariates$f.26726.2.0 + UKB_inflammation_imaging_covariates$f.26827.2.0

## Parietal area lobe:
## Postcentral gyrus = 26742 + 26843
## paracentral cortex = 26737 + 26838
## Superior parietal cortex = 26749 + 26850
## Inferior parietal cortex = 26728 + 26829
## Supramarginal gyrus = 26751 + 26852
## Precuneus = 26745 + 26846

## Sum area of individual structures in parietal lobe
UKB_inflammation_imaging_covariates$parietal_lobe_area <- UKB_inflammation_imaging_covariates$f.26742.2.0 + UKB_inflammation_imaging_covariates$f.26843.2.0 + UKB_inflammation_imaging_covariates$f.26737.2.0 +
  UKB_inflammation_imaging_covariates$f.26838.2.0 + UKB_inflammation_imaging_covariates$f.26749.2.0 + UKB_inflammation_imaging_covariates$f.26850.2.0 + UKB_inflammation_imaging_covariates$f.26728.2.0 +
  UKB_inflammation_imaging_covariates$f.26829.2.0 + UKB_inflammation_imaging_covariates$f.26751.2.0 + UKB_inflammation_imaging_covariates$f.26852.2.0 + UKB_inflammation_imaging_covariates$f.26745.2.0 +
  UKB_inflammation_imaging_covariates$f.26846.2.0

## Occipital area lobe:
## Lateral occipital cortex = 26731, 26832
## Cuneus = 26725, 26826
## Pericalcarine cortex = 26741, 26842
## Lingual gyrus = 26733, 26834

UKB_inflammation_imaging_covariates$occipital_lobe_area <- UKB_inflammation_imaging_covariates$f.26731.2.0 + UKB_inflammation_imaging_covariates$f.26832.2.0 + UKB_inflammation_imaging_covariates$f.26725.2.0 +
  UKB_inflammation_imaging_covariates$f.26826.2.0 + UKB_inflammation_imaging_covariates$f.26741.2.0 + UKB_inflammation_imaging_covariates$f.26842.2.0 + UKB_inflammation_imaging_covariates$f.26733.2.0 +
  UKB_inflammation_imaging_covariates$f.26834.2.0

## Cingulate area lobe:
## Rostral ACC = 26746, 26847
## Caudal ACC = 26723, 26824
## Posterior cingulate cortex = 26743 + 26844
## Cingulate isthmus = 26730 + 26831

UKB_inflammation_imaging_covariates$cingulate_lobe_area <- UKB_inflammation_imaging_covariates$f.26746.2.0 + UKB_inflammation_imaging_covariates$f.26847.2.0 + UKB_inflammation_imaging_covariates$f.26723.2.0 +
  UKB_inflammation_imaging_covariates$f.26824.2.0 + UKB_inflammation_imaging_covariates$f.26743.2.0 + UKB_inflammation_imaging_covariates$f.26844.2.0 + UKB_inflammation_imaging_covariates$f.26730.2.0 +
  UKB_inflammation_imaging_covariates$f.26831.2.0

## Sum lobar areas to get global cortical area
UKB_inflammation_imaging_covariates$global_cortical_area <- UKB_inflammation_imaging_covariates$frontal_lobe_area + UKB_inflammation_imaging_covariates$temporal_lobe_area + UKB_inflammation_imaging_covariates$parietal_lobe_area +
  UKB_inflammation_imaging_covariates$occipital_lobe_area +UKB_inflammation_imaging_covariates$cingulate_lobe_area

## Calculate cortical lobar thickness measures
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Frontal lobe thickness FIDs:
## Superior frontal gyrus = 26782 + 26883
## rostralmiddlefrontal = 26781 + 26882
## Caudalmiddlefrontal = 26758 + 26859
## Pars orbitalis = 26773 + 26874
## pars triangularis = 26774 + 26875
## pars opercularis = 26772 + 26873
## Frontal pole = 26786 + 26887
## Lateral orbitofrontal = 26766 + 26867
## medial orbitofrontal = 26768 + 26869
## Precentral gyrus = 26778 + 26879
## paracentral cortex = 26771 + 26872

## Sum thickness of individual structures in frontal lobe
UKB_inflammation_imaging_covariates$frontal_lobe_thickness <- UKB_inflammation_imaging_covariates$f.26782.2.0 + UKB_inflammation_imaging_covariates$f.26883.2.0 + UKB_inflammation_imaging_covariates$f.26781.2.0 +
  UKB_inflammation_imaging_covariates$f.26882.2.0 + UKB_inflammation_imaging_covariates$f.26758.2.0 + UKB_inflammation_imaging_covariates$f.26859.2.0 + UKB_inflammation_imaging_covariates$f.26773.2.0 +
  UKB_inflammation_imaging_covariates$f.26874.2.0 + UKB_inflammation_imaging_covariates$f.26774.2.0 + UKB_inflammation_imaging_covariates$f.26875.2.0 + UKB_inflammation_imaging_covariates$f.26772.2.0 + 
  UKB_inflammation_imaging_covariates$f.26873.2.0 + UKB_inflammation_imaging_covariates$f.26786.2.0 + UKB_inflammation_imaging_covariates$f.26887.2.0 + UKB_inflammation_imaging_covariates$f.26766.2.0 +
  UKB_inflammation_imaging_covariates$f.26867.2.0 + UKB_inflammation_imaging_covariates$f.26768.2.0 + UKB_inflammation_imaging_covariates$f.26869.2.0 + UKB_inflammation_imaging_covariates$f.26778.2.0 + 
  UKB_inflammation_imaging_covariates$f.26879.2.0 + UKB_inflammation_imaging_covariates$f.26771.2.0 + UKB_inflammation_imaging_covariates$f.26872.2.0

## Temporal lobe thickness FIDs:
## Insula= 26788 + 26889
## Superior temporal = 26784 + 26885
## transverse temporal = 26787 + 26888
## banks STS = 	26756 + 26857
## Middle temporal gyrus = 26769 + 26870
## Inferior temporal gyrus = 26796 + 26897
## Fusiform = 26761 + 26862
## parahippocampal = 26770 + 26871
## entorhinal = 26760 + 26760

## Sum thickness of individual structures in temporal lobe
UKB_inflammation_imaging_covariates$temporal_lobe_thickness <- UKB_inflammation_imaging_covariates$f.26788.2.0 + UKB_inflammation_imaging_covariates$f.26889.2.0 + UKB_inflammation_imaging_covariates$f.26784.2.0 +
  UKB_inflammation_imaging_covariates$f.26885.2.0 + UKB_inflammation_imaging_covariates$f.26787.2.0 + UKB_inflammation_imaging_covariates$f.26888.2.0 + UKB_inflammation_imaging_covariates$f.26756.2.0 +
  UKB_inflammation_imaging_covariates$f.26857.2.0 + UKB_inflammation_imaging_covariates$f.26769.2.0 + UKB_inflammation_imaging_covariates$f.26870.2.0 + UKB_inflammation_imaging_covariates$f.26796.2.0 +
  UKB_inflammation_imaging_covariates$f.26897.2.0 + UKB_inflammation_imaging_covariates$f.26761.2.0 + UKB_inflammation_imaging_covariates$f.26862.2.0 + UKB_inflammation_imaging_covariates$f.26770.2.0 +
  UKB_inflammation_imaging_covariates$f.26871.2.0 + UKB_inflammation_imaging_covariates$f.26760.2.0 + UKB_inflammation_imaging_covariates$f.26760.2.0

## Parietal thickness lobe:
## Postcentral gyrus = 26776 + 26877
## paracentral cortex = 26771 + 26872
## Superior parietal cortex = 26783 + 26884
## Inferior parietal cortex = 26762 + 26863
## Supramarginal gyrus = 26785 + 26886
## Precuneus = 26779 + 26880

## Sum thickness of individual structures in parietal lobe
UKB_inflammation_imaging_covariates$parietal_lobe_thickness <- UKB_inflammation_imaging_covariates$f.26776.2.0 + UKB_inflammation_imaging_covariates$f.26877.2.0 + UKB_inflammation_imaging_covariates$f.26771.2.0 +
  UKB_inflammation_imaging_covariates$f.26872.2.0 + UKB_inflammation_imaging_covariates$f.26783.2.0 + UKB_inflammation_imaging_covariates$f.26884.2.0 + UKB_inflammation_imaging_covariates$f.26762.2.0 +
  UKB_inflammation_imaging_covariates$f.26863.2.0 + UKB_inflammation_imaging_covariates$f.26785.2.0 + UKB_inflammation_imaging_covariates$f.26886.2.0 + UKB_inflammation_imaging_covariates$f.26779.2.0 +
  UKB_inflammation_imaging_covariates$f.26880.2.0

## Occipital thickness lobe:
## Lateral occipital cortex = 26765, 26866
## Cuneus = 26759, 26860
## Pericalcarine cortex = 26775, 26876
## Lingual gyrus = 26767, 26868

UKB_inflammation_imaging_covariates$occipital_lobe_thickness <- UKB_inflammation_imaging_covariates$f.26765.2.0 + UKB_inflammation_imaging_covariates$f.26866.2.0 + UKB_inflammation_imaging_covariates$f.26759.2.0 +
  UKB_inflammation_imaging_covariates$f.26860.2.0 + UKB_inflammation_imaging_covariates$f.26775.2.0 + UKB_inflammation_imaging_covariates$f.26876.2.0 + UKB_inflammation_imaging_covariates$f.26767.2.0 +
  UKB_inflammation_imaging_covariates$f.26868.2.0

## Cingulate thickness lobe:
## Rostral ACC = 26780, 26881
## Caudal ACC = 26757, 26858
## Posterior cingulate cortex = 26777 + 26878
## Cingulate isthmus = 26764 + 26865

UKB_inflammation_imaging_covariates$cingulate_lobe_thickness <- UKB_inflammation_imaging_covariates$f.26780.2.0 + UKB_inflammation_imaging_covariates$f.26881.2.0 + UKB_inflammation_imaging_covariates$f.26757.2.0 +
  UKB_inflammation_imaging_covariates$f.26858.2.0 + UKB_inflammation_imaging_covariates$f.26777.2.0 + UKB_inflammation_imaging_covariates$f.26878.2.0 + UKB_inflammation_imaging_covariates$f.26764.2.0 +
  UKB_inflammation_imaging_covariates$f.26865.2.0

## Sum lobar thicknesss to get global cortical thickness
UKB_inflammation_imaging_covariates$global_cortical_thickness <- UKB_inflammation_imaging_covariates$frontal_lobe_volume + UKB_inflammation_imaging_covariates$temporal_lobe_volume + UKB_inflammation_imaging_covariates$parietal_lobe_volume +
  UKB_inflammation_imaging_covariates$occipital_lobe_volume +UKB_inflammation_imaging_covariates$cingulate_lobe_volume 


## Standardise all cortical global and lobar volumes, area and thickness
UKB_inflammation_imaging_covariates <- UKB_inflammation_imaging_covariates %>%
  mutate_at(vars(contains('_lobe_')), scale) %>%
  mutate_at(vars(contains('global_cortical')), scale)

## Calculate general MD measures
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Association/Commissural fibres

# Cingulate gyrus part of cingulum (left) - f.25519
# Cingulate gyrus part of cingulum (right) - f.25520
# Inferior fronto-occipital fasciculus (left) - f.25527
# Inferior fronto-occipital fasciculus (right) - f.25528
# Inferior longitudinal fasciculus (left) - f.25529
# Inferior longitudinal fasciculus (right) - f.25530
# Parahippocampal part of cingulum (left) - f.25521
# Parahippocampal part of cingulum (right) - f.25522
# Superior longitudinal fasciculus (left) - f.25536
# Superior longitudinal fasciculus (right) - f.25537
# Uncinate fasciculus (left) - f.25540
# Uncinate fasciculus (right) - f.25541
# Forceps major - f.25525
# Forceps minor - f.25526


association_fibres_MD <- c("f.25519.2.0", "f.25520.2.0", "f.25527.2.0", "f.25528.2.0", "f.25529.2.0", "f.25530.2.0",
                           "f.25521.2.0", "f.25522.2.0", "f.25536.2.0", "f.25537.2.0", "f.25540.2.0", "f.25541.2.0",
                           "f.25525.2.0", "f.25526.2.0")


## Thalamic radiations

# Anterior thalamic radiation (left)  - f.25517
# Anterior thalamic radiation (right) - f.25518
# Posterior thalamic radiation (left) - f.25534
# Posterior thalamic radiation (right) - f.25535
# Superior thalamic radiation (left) - f.25538
# Superior thalamic radiation (right) - f.25539

thalamic_radiations_MD <- c("f.25517.2.0", "f.25518.2.0", "f.25534.2.0", "f.25535.2.0", "f.25538.2.0", "f.25539.2.0")


## Projection fibres

# AcouDTIc radiation (left) - f.25516
# AcouDTIc radiation (right) - f.25517
# Corticospinal tract (left) - f.25523
# Corticospinal tract (right) - f.25524
# Medial lemniscus (left) - f.25532
# Medial lemniscus (right) - f.25533
# Middle cerebellar peduncle - f.25531


projection_fibres_MD <- c("f.25516.2.0", "f.25517.2.0", "f.25523.2.0", "f.25524.2.0", "f.25532.2.0", "f.25533.2.0", "f.25531.2.0")

## Global_MD
global_MD <- c(association_fibres_MD, thalamic_radiations_MD, projection_fibres_MD)


## Remove participants without DTI information and store in temporary file (PCA does not allow NAs)
UKB_inflammation_imaging_covariates_temp <- UKB_inflammation_imaging_covariates[!is.na(UKB_inflammation_imaging_covariates$f.25516.2.0),]

## Run PCA in each general MD measure
PC_res_MD_association_fibres <- prcomp(UKB_inflammation_imaging_covariates_temp[,association_fibres_MD],scale = TRUE)
PC_res_MD_thalamic_radiations <- prcomp(UKB_inflammation_imaging_covariates_temp[,thalamic_radiations_MD],scale = TRUE)
PC_res_MD_projection_fibres <- prcomp(UKB_inflammation_imaging_covariates_temp[,projection_fibres_MD],scale = TRUE)
PC_res_MD_global <- prcomp(UKB_inflammation_imaging_covariates_temp[,global_MD],scale = TRUE)

## Add 1st PC to temporary dataframe
UKB_inflammation_imaging_covariates_temp <- cbind(UKB_inflammation_imaging_covariates_temp, PC_res_MD_association_fibres$x[,1], PC_res_MD_thalamic_radiations$x[,1], PC_res_MD_projection_fibres$x[,1],
                                                             PC_res_MD_global$x[,1])

## Rename column containing 1st PC of each general MD meaure and keep only ID and general MD meaures
UKB_inflammation_imaging_covariates_temp <- UKB_inflammation_imaging_covariates_temp %>%
  rename(global_MD= "PC_res_MD_global$x[, 1]",
         association_fibres_MD= "PC_res_MD_association_fibres$x[, 1]",
         thalamic_radiations_MD= "PC_res_MD_thalamic_radiations$x[, 1]",
         projection_fibres_MD= "PC_res_MD_projection_fibres$x[, 1]") %>%
  select(f.eid, global_MD, association_fibres_MD, thalamic_radiations_MD, projection_fibres_MD)


## Merge with full dataframe 
UKB_inflammation_imaging_covariates <- left_join(UKB_inflammation_imaging_covariates, UKB_inflammation_imaging_covariates_temp, by = "f.eid")



## Calculate general FA measures
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Association/Commissural fibres


## Association/Commissural fibres

# Cingulate gyrus part of cingulum (left) - f.25492
# Cingulate gyrus part of cingulum (right) - f.25493
# Inferior fronto-occipital fasciculus (left) - f.25500
# Inferior fronto-occipital fasciculus (right) - f.25501
# Inferior longitudinal fasciculus (left) - f.25502
# Inferior longitudinal fasciculus (right) - f.25503
# Parahippocampal part of cingulum (left) - f.25494
# Parahippocampal part of cingulum (right) - f.25495
# Superior longitudinal fasciculus (left) - f.25509
# Superior longitudinal fasciculus (right) - f.25510
# Uncinate fasciculus (left) - f.25513
# Uncinate fasciculus (right) - f.25514
# Forceps major - f.25498
# Forceps minor - f.25499


association_fibres_FA <- c("f.25492.2.0", "f.25493.2.0", "f.25500.2.0", "f.25501.2.0", "f.25502.2.0", "f.25503.2.0",
                           "f.25494.2.0", "f.25495.2.0", "f.25509.2.0", "f.25510.2.0", "f.25513.2.0", "f.25514.2.0",
                           "f.25498.2.0", "f.25499.2.0")

## Thalamic radiations

# Anterior thalamic radiation (left)  - f.25490
# Anterior thalamic radiation (right) - f.25491
# Posterior thalamic radiation (left) - f.25507
# Posterior thalamic radiation (right) - f.25508
# Superior thalamic radiation (left) - f.25511
# Superior thalamic radiation (right) - f.25512

thalamic_radiations_FA <- c("f.25490.2.0", "f.25491.2.0", "f.25507.2.0", "f.25508.2.0", "f.25511.2.0", "f.25512.2.0")


## Projection fibres

# AcouDTIc radiation (left) - f.25488
# AcouDTIc radiation (right) - f.25489
# Corticospinal tract (left) - f.25496
# Corticospinal tract (right) - f.25497
# Medial lemniscus (left) - f.25505
# Medial lemniscus (right) - f.25506
# Middle cerebellar peduncle - f.25504



projection_fibres_FA <- c("f.25488.2.0", "f.25489.2.0", "f.25496.2.0", "f.25497.2.0", "f.25505.2.0", "f.25506.2.0", "f.25504.2.0")

## Global_FA
global_FA <- c(association_fibres_FA, thalamic_radiations_FA, projection_fibres_FA)

## Remove participants without DTI information and store in temporary file (PCA does not allow NAs)
UKB_inflammation_imaging_covariates_temp <- UKB_inflammation_imaging_covariates[!is.na(UKB_inflammation_imaging_covariates$f.25516.2.0),]

## Run PCA in each general FA measure
PC_res_FA_association_fibres <- prcomp(UKB_inflammation_imaging_covariates_temp[,association_fibres_FA],scale = TRUE)
PC_res_FA_thalamic_radiations <- prcomp(UKB_inflammation_imaging_covariates_temp[,thalamic_radiations_FA],scale = TRUE)
PC_res_FA_projection_fibres <- prcomp(UKB_inflammation_imaging_covariates_temp[,projection_fibres_FA],scale = TRUE)
PC_res_FA_global <- prcomp(UKB_inflammation_imaging_covariates_temp[,global_FA],scale = TRUE)

## Add 1st PC to temporary dataframe
UKB_inflammation_imaging_covariates_temp <- cbind(UKB_inflammation_imaging_covariates_temp, PC_res_FA_association_fibres$x[,1], PC_res_FA_thalamic_radiations$x[,1], PC_res_FA_projection_fibres$x[,1],
                                                  PC_res_FA_global$x[,1])

## Rename column containing 1st PC of each general FA meaure and keep only ID and general FA meaures
UKB_inflammation_imaging_covariates_temp <- UKB_inflammation_imaging_covariates_temp %>%
  rename(global_FA= "PC_res_FA_global$x[, 1]",
         association_fibres_FA= "PC_res_FA_association_fibres$x[, 1]",
         thalamic_radiations_FA= "PC_res_FA_thalamic_radiations$x[, 1]",
         projection_fibres_FA= "PC_res_FA_projection_fibres$x[, 1]") %>%
  select(f.eid, global_FA, association_fibres_FA, thalamic_radiations_FA, projection_fibres_FA)


## Merge with full dataframe 
UKB_inflammation_imaging_covariates <- left_join(UKB_inflammation_imaging_covariates, UKB_inflammation_imaging_covariates_temp, by = "f.eid")


## Save full dataframe to output directory
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write.csv(UKB_inflammation_imaging_covariates, "~/Desktop/PhD/projects/UKB_inflammation_imaging/resources//UKB_inflammation_imaging_covariates.csv", quote = F)

