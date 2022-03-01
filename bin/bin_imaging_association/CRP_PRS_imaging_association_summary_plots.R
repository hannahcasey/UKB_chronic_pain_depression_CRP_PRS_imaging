## Load Packages
#install.packages("remotes")
#remotes::install_github("ggseg/ggseg")
library(ggplot2)
#install.packages("RColorBrewer")
library(RColorBrewer)
library(dplyr)
#install.packages("gplots")
library(gplots)
library(patchwork)


## load in all data
file_paths = list.files(path = "~/Desktop/PhD/projects/UKB_inflammation_imaging/output/PRS_imaging_association/", pattern ="*.csv", full.names = T)
file_paths_short = list.files(path = "~/Desktop/PhD/projects/UKB_inflammation_imaging/output/PRS_imaging_association/", pattern ="*.csv", full.names = F)
myfiles = lapply(file_paths, read.csv)
names(myfiles) <- file_paths_short
list2env(myfiles,envir=.GlobalEnv)

file_paths = list.files(path = "~/Desktop/PhD/projects/UKB_inflammation_imaging/output/CRP_pheno_imaging_association/", pattern ="*.csv", full.names = T)
file_paths_short = list.files(path = "~/Desktop/PhD/projects/UKB_inflammation_imaging/output/CRP_pheno_imaging_association/", pattern ="*.csv", full.names = F)
myfiles = lapply(file_paths, read.csv)
names(myfiles) <- file_paths_short
list2env(myfiles,envir=.GlobalEnv)

file_paths = list.files(path = "~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/", pattern ="*.csv", full.names = T)
file_paths_short = list.files(path = "~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/", pattern ="*.csv", full.names = F)
myfiles = lapply(file_paths, read.csv)
names(myfiles) <- file_paths_short
list2env(myfiles,envir=.GlobalEnv)


## Plot all CRP PRS results
CRP_PRS_general_cortical_volumes_small <- CRP_PRS_cortical_general.csv %>%
  rename(feature = cortical_volume, p.value.adjust = p.adjust) %>%
  select(feature, p.value, p.value.adjust)
CRP_PRS_general_cortical_volumes_small$Group <- "General Cortical Volume"

CRP_PRS_general_cortical_area_small <- CRP_PRS_cortical_area_general.csv %>%
  rename(feature = cortical_area, p.value.adjust = p.adjust) %>%
  select(feature, p.value, p.value.adjust)
CRP_PRS_general_cortical_area_small$Group <- "General Cortical Area"

CRP_PRS_general_cortical_thickness_small <- CRP_PRS_cortical_thickness_general.csv %>%
  rename(feature = cortical_thickness, p.value.adjust = p.adjust) %>%
  select(feature, p.value, p.value.adjust)
CRP_PRS_general_cortical_thickness_small$Group <- "General Cortical Thickness"

CRP_PRS_regional_cortical_volumes_small <- CRP_regional_cortical_volumes.csv %>%
  rename(feature = cortical_volume) %>%
  select(feature, p.value, p.value.adjust)
CRP_PRS_regional_cortical_volumes_small$Group <- "Regional Cortical Volume"

CRP_PRS_regional_cortical_area_small <- lme_cortical_area_CRP_PRS.csv %>%
  rename(feature = cortical_area) %>%
  select(feature, p.value, p.value.adjust)
CRP_PRS_regional_cortical_area_small$Group <- "Regional Cortical Area"

CRP_PRS_regional_cortical_thickness_small <- lme_cortical_thickness_CRP_PRS.csv %>%
  rename(feature = cortical_thickness) %>%
  select(feature, p.value, p.value.adjust)
CRP_PRS_regional_cortical_thickness_small$Group <- "Regional Cortical Thickness"

CRP_PRS_subcortical_volumes_small <- CRP_regional_subcortical_volumes.csv %>%
  rename(feature = subcortical_volume) %>%
  select(feature, p.value, p.value.adjust)
CRP_PRS_subcortical_volumes_small$Group <- "Subcortical Volume"

CRP_bilateral_FA_tracts_small <- CRP_bilateral_FA_tracts.csv %>%
  rename(feature = FA_tract) %>%
  select(feature, p.value, p.value.adjust)
CRP_bilateral_FA_tracts_small$Group <- "FA Tracts"

CRP_bilateral_MD_tracts_small <- CRP_bilateral_MD_tracts.csv %>%
  rename(feature = MD_tract) %>%
  select(feature, p.value, p.value.adjust)
CRP_bilateral_MD_tracts_small$Group <- "MD Tracts"

CRP_unilateral_FA_tracts_small <- CRP_unilateral_FA_tracts.csv %>%
  rename(feature = FA_tract) %>%
  select(feature, p.value, p.value.adjust)
CRP_unilateral_FA_tracts_small$Group <- "FA Tracts"

CRP_unilateral_MD_tracts_small <- CRP_unilateral_MD_tracts.csv %>%
  rename(feature = MD_tract) %>%
  select(feature, p.value, p.value.adjust)
CRP_unilateral_MD_tracts_small$Group <- "MD Tracts"

CRP_PRS_FA_general_small <- CRP_PRS_FA_general.csv %>%
  rename(feature = general_DTI_FA_measure,  p.value.adjust = p.adjust) %>%
  select(feature, p.value, p.value.adjust)
CRP_PRS_FA_general_small$Group <- "General FA Measures"

CRP_PRS_MD_general_small <- CRP_PRS_MD_general.csv %>%
  rename(feature = general_DTI_MD_measure, p.value.adjust = p.adjust) %>%
  select(feature, p.value, p.value.adjust)
CRP_PRS_MD_general_small$Group <- "General MD Measures"
## Combine all dataframes 
CRP_PRS_all <- rbind(CRP_PRS_general_cortical_volumes_small, CRP_PRS_general_cortical_area_small, CRP_PRS_general_cortical_thickness_small, CRP_PRS_regional_cortical_volumes_small, 
                     CRP_PRS_regional_cortical_area_small, CRP_PRS_regional_cortical_thickness_small, CRP_PRS_subcortical_volumes_small, CRP_bilateral_FA_tracts_small, 
                     CRP_bilateral_MD_tracts_small, CRP_unilateral_FA_tracts_small, CRP_unilateral_MD_tracts_small, CRP_PRS_FA_general_small, CRP_PRS_MD_general_small)

## If p.adjust = NA (no FDR correction) replace with pvalue
CRP_PRS_all$p.value.adjust[is.na(CRP_PRS_all$p.value.adjust)] <- CRP_PRS_all$p.value[is.na(CRP_PRS_all$p.value.adjust)]

## Log tramsform p-value 
CRP_PRS_all$p_minus_log_10 <- -log10(CRP_PRS_all$p.value)
## Plot results

## Reorder feature type categories
CRP_PRS_all$Group <- factor(CRP_PRS_all$Group, levels = c('General Cortical Volume', 'General Cortical Area', 'General Cortical Thickness', 'Regional Cortical Volume', 'Regional Cortical Area',
                                                      'Regional Cortical Thickness', 'Subcortical Volume', 'General FA Measures', 'General MD Measures', 'FA Tracts', 'MD Tracts'))


jpeg("~/Desktop/PhD/projects/UKB_inflammation_imaging/output/summary_plots/all_CRP_PRS.jpeg", width = 1000, height = 1000)

if (nrow(CRP_PRS_all[which(CRP_PRS_all$p.value.adjust<0.05),])>0){
ggplot(CRP_PRS_all, aes(x = Group, y = p_minus_log_10)) +
  geom_point(data = CRP_PRS_all[which(CRP_PRS_all$p.value.adjust>0.05),], 
             aes(colour =  Group), size = 5, alpha = 0.8, shape = 16,
             position="jitter") +
    geom_point(data = CRP_PRS_all[which(CRP_PRS_all$p.value.adjust<0.05),], 
                aes(colour =  Group), size = 5, alpha = 0.8, shape = 18,
                position="jitter") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", 
             color = "red", size=1, alpha = 0.4) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y=element_blank()) +
  xlab("") + ylab("-log10(p)")
} else{
  ggplot(CRP_PRS_all, aes(x = Group, y = p_minus_log_10)) +
    geom_point(data = CRP_PRS_all[which(CRP_PRS_all$p.value.adjust>0.05),], 
               aes(colour =  Group), size = 5, alpha = 0.8, shape = 16,
               position="jitter") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", 
             color = "red", size=1, alpha = 0.4) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y=element_blank(), text = element_text(size = 20)) +
  xlab("") + ylab("-log10(p)")
}
dev.off()

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
## Plot all Serum CRP (wothin week) results
##

CRP_pheno_within_week_regional_cortical_volumes_small <- CRP_pheno_within_week_cortical_regions.csv %>%
  rename(feature = cortical_volume) %>%
  select(feature, p.value, p.value.adjust)
CRP_pheno_within_week_regional_cortical_volumes_small$Group <- "Regional Cortical Volume"

CRP_pheno_within_week_general_cortical_volumes_small <- CRP_pheno_within_week_cortical_general.csv %>%
  rename(feature = cortical_volume, p.value.adjust = p.adjust) %>%
  select(feature, p.value, p.value.adjust)
CRP_pheno_within_week_general_cortical_volumes_small$Group <- "General Cortical Volume"

CRP_pheno_within_week_subcortical_volumes_small <- CRP_pheno_within_week_subcortical.csv %>%
  rename(feature = subcortical_volume) %>%
  select(feature, p.value, p.value.adjust)
CRP_pheno_within_week_subcortical_volumes_small$Group <- "Subortical Volume"

CRP_pheno_within_week_bilateral_FA_tracts_small <- CRP_pheno_within_week_FA_bilateral_tracts.csv %>%
  rename(feature = FA_tract) %>%
  select(feature, p.value, p.value.adjust)
CRP_pheno_within_week_bilateral_FA_tracts_small$Group <- "FA Tracts"

CRP_pheno_within_week_bilateral_MD_tracts_small <- CRP_pheno_within_week_MD_bilateral_tracts.csv %>%
  rename(feature = MD_tract) %>%
  select(feature, p.value, p.value.adjust)
CRP_pheno_within_week_bilateral_MD_tracts_small$Group <- "MD Tracts"

CRP_pheno_within_week_unilateral_FA_tracts_small <- CRP_pheno_within_week_FA_unilateral_tracts.csv %>%
  rename(feature = FA_tract) %>%
  select(feature, p.value, p.value.adjust)
CRP_pheno_within_week_unilateral_FA_tracts_small$Group <- "FA Tracts"

CRP_pheno_within_week_unilateral_MD_tracts_small <- CRP_pheno_within_week_MD_unilateral_tracts.csv %>%
  rename(feature = MD_tract) %>%
  select(feature, p.value, p.value.adjust)
CRP_pheno_within_week_unilateral_MD_tracts_small$Group <- "MD Tracts"

CRP_pheno_within_week_FA_general_small <- CRP_pheno_within_week_FA_general.csv %>%
  rename(feature = general_DTI_FA_measure,  p.value.adjust = p.adjust) %>%
  select(feature, p.value, p.value.adjust)
CRP_pheno_within_week_FA_general_small$Group <- "General FA Measures"

CRP_pheno_within_week_MD_general_small <- CRP_pheno_within_week_MD_general.csv %>%
  rename(feature = general_DTI_MD_measure, p.value.adjust = p.adjust) %>%
  select(feature, p.value, p.value.adjust)
CRP_pheno_within_week_MD_general_small$Group <- "General MD Measures"

## Combine all dataframes 
CRP_pheno_within_week_all <- rbind(CRP_pheno_within_week_regional_cortical_volumes_small, CRP_pheno_within_week_general_cortical_volumes_small, CRP_pheno_within_week_subcortical_volumes_small, CRP_pheno_within_week_bilateral_FA_tracts_small, 
              CRP_pheno_within_week_bilateral_MD_tracts_small, CRP_pheno_within_week_unilateral_FA_tracts_small, CRP_pheno_within_week_unilateral_MD_tracts_small, CRP_pheno_within_week_FA_general_small, CRP_pheno_within_week_MD_general_small)

## If p.adjust = NA (no FDR correction) replace with pvalue
CRP_pheno_within_week_all$p.value.adjust[is.na(CRP_pheno_within_week_all$p.value.adjust)] <- CRP_pheno_within_week_all$p.value[is.na(CRP_pheno_within_week_all$p.value.adjust)]

## Log tramsform p-value 
CRP_pheno_within_week_all$p_minus_log_10 <- -log10(CRP_pheno_within_week_all$p.value)
## Plot results

jpeg("~/Desktop/PhD/projects/UKB_inflammation_imaging/output/summary_plots/all_CRP_pheno_within_week.jpeg", width = 1000, height = 1000)

if (nrow(CRP_pheno_within_week_all[which(CRP_pheno_within_week_all$p.value.adjust<0.05),])>0){
  ggplot(CRP_pheno_within_week_all, aes(x = Group, y = p_minus_log_10)) +
    geom_point(data = CRP_pheno_within_week_all[which(CRP_pheno_within_week_all$p.value.adjust>0.05),], 
               aes(colour =  Group), size = 5, alpha = 0.8, shape = 16,
               position="jitter") +
    geom_point(data = CRP_pheno_within_week_all[which(CRP_pheno_within_week_all$p.value.adjust<0.05),], 
               aes(colour =  Group), size = 5, alpha = 0.8, shape = 18,
               position="jitter") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", 
               color = "red", size=1, alpha = 0.4) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y=element_blank()) +
    xlab("") + ylab("-log10(p)")
} else{
  ggplot(CRP_pheno_within_week_all, aes(x = Group, y = p_minus_log_10)) +
    geom_point(data = CRP_pheno_within_week_all[which(CRP_pheno_within_week_all$p.value.adjust>0.05),], 
               aes(colour =  Group), size = 5, alpha = 0.8, shape = 16,
               position="jitter") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", 
               color = "red", size=1, alpha = 0.4) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y=element_blank(), text = element_text(size = 20)) +
    xlab("") + ylab("-log10(p)")
}
dev.off()

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
## Plot all Serum CRP (wothin month) results
##

CRP_pheno_within_month_regional_cortical_volumes_small <- CRP_pheno_within_month_cortical_regions.csv %>%
  rename(feature = cortical_volume) %>%
  select(feature, p.value, p.value.adjust)
CRP_pheno_within_month_regional_cortical_volumes_small$Group <- "Regional Cortical Volumes"

CRP_pheno_within_month_general_cortical_volumes_small <- CRP_pheno_within_month_cortical_general.csv %>%
  rename(feature = cortical_volume, p.value.adjust = p.adjust) %>%
  select(feature, p.value, p.value.adjust)
CRP_pheno_within_month_general_cortical_volumes_small$Group <- "General Cortical Volumes"

CRP_pheno_within_month_subcortical_volumes_small <- CRP_pheno_within_month_subcortical.csv %>%
  rename(feature = subcortical_volume) %>%
  select(feature, p.value, p.value.adjust)
CRP_pheno_within_month_subcortical_volumes_small$Group <- "Subortical Volumes"

CRP_pheno_within_month_bilateral_FA_tracts_small <- CRP_pheno_within_month_FA_bilateral_tracts.csv %>%
  rename(feature = FA_tract) %>%
  select(feature, p.value, p.value.adjust)
CRP_pheno_within_month_bilateral_FA_tracts_small$Group <- "FA Tracts"

CRP_pheno_within_month_bilateral_MD_tracts_small <- CRP_pheno_within_month_MD_bilateral_tracts.csv %>%
  rename(feature = MD_tract) %>%
  select(feature, p.value, p.value.adjust)
CRP_pheno_within_month_bilateral_MD_tracts_small$Group <- "MD Tracts"

CRP_pheno_within_month_unilateral_FA_tracts_small <- CRP_pheno_within_month_FA_unilateral_tracts.csv %>%
  rename(feature = FA_tract) %>%
  select(feature, p.value, p.value.adjust)
CRP_pheno_within_month_unilateral_FA_tracts_small$Group <- "FA Tracts"

CRP_pheno_within_month_unilateral_MD_tracts_small <- CRP_pheno_within_month_MD_unilateral_tracts.csv %>%
  rename(feature = MD_tract) %>%
  select(feature, p.value, p.value.adjust)
CRP_pheno_within_month_unilateral_MD_tracts_small$Group <- "MD Tracts"

CRP_pheno_within_month_FA_general_small <- CRP_pheno_within_month_FA_general.csv %>%
  rename(feature = general_DTI_FA_measure,  p.value.adjust = p.adjust) %>%
  select(feature, p.value, p.value.adjust)
CRP_pheno_within_month_FA_general_small$Group <- "General FA Measures"

CRP_pheno_within_month_MD_general_small <- CRP_pheno_within_month_MD_general.csv %>%
  rename(feature = general_DTI_MD_measure, p.value.adjust = p.adjust) %>%
  select(feature, p.value, p.value.adjust)
CRP_pheno_within_month_MD_general_small$Group <- "General MD Measures"
## Combine all dataframes 
CRP_pheno_within_month_all <- rbind(CRP_pheno_within_month_regional_cortical_volumes_small, CRP_pheno_within_month_general_cortical_volumes_small, CRP_pheno_within_month_subcortical_volumes_small, CRP_pheno_within_month_bilateral_FA_tracts_small, 
                                   CRP_pheno_within_month_bilateral_MD_tracts_small, CRP_pheno_within_month_unilateral_FA_tracts_small, CRP_pheno_within_month_unilateral_MD_tracts_small, CRP_pheno_within_month_FA_general_small, CRP_pheno_within_month_MD_general_small)

## If p.adjust = NA (no FDR correction) replace with pvalue
CRP_pheno_within_month_all$p.value.adjust[is.na(CRP_pheno_within_month_all$p.value.adjust)] <- CRP_pheno_within_month_all$p.value[is.na(CRP_pheno_within_month_all$p.value.adjust)]

## Log tramsform p-value 
CRP_pheno_within_month_all$p_minus_log_10 <- -log10(CRP_pheno_within_month_all$p.value)
## Plot results

jpeg("~/Desktop/PhD/projects/UKB_inflammation_imaging/output/summary_plots/all_CRP_pheno_within_month.jpeg", width = 1000, height = 1000)

if (nrow(CRP_pheno_within_month_all[which(CRP_pheno_within_month_all$p.value.adjust<0.05),])>0){
  ggplot(CRP_pheno_within_month_all, aes(x = Group, y = p_minus_log_10)) +
    geom_jitter(data = CRP_pheno_within_month_all[which(CRP_pheno_within_month_all$p.value.adjust>0.05),], 
               aes(colour =  Group), size = 5, alpha = 0.8, shape = 16,
               position="jitter") +
    geom_jitter(data = CRP_pheno_within_month_all[which(CRP_pheno_within_month_all$p.value.adjust<0.05),], 
               aes(colour =  Group), size = 5, alpha = 0.8, shape = 15,
               height= 0) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", 
               color = "red", size=1, alpha = 0.4) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y=element_blank()) +
    xlab("") + ylab("-log10(p)")
} else{
  ggplot(CRP_pheno_within_month_all, aes(x = Group, y = p_minus_log_10)) +
    geom_jitter(data = CRP_pheno_within_month_all[which(CRP_pheno_within_month_all$p.value.adjust>0.05),], 
               aes(colour =  Group), size = 5, alpha = 0.8, shape = 16,
               position="jitter") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", 
               color = "red", size=1, alpha = 0.4) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y=element_blank(), text = element_text(size = 20)) +
    xlab("") + ylab("-log10(p)")
}
dev.off()


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
## Plot all Serum CRP (wothin 3 months) results
##

CRP_pheno_within_3_months_regional_cortical_volumes_small <- CRP_pheno_within_3_months_cortical_regions.csv %>%
  rename(feature = cortical_volume) %>%
  select(feature, p.value, p.value.adjust)
CRP_pheno_within_3_months_regional_cortical_volumes_small$Group <- "Regional Cortical Volumes"

CRP_pheno_within_3_months_general_cortical_volumes_small <- CRP_pheno_within_3_months_cortical_general.csv %>%
  rename(feature = cortical_volume, p.value.adjust = p.adjust) %>%
  select(feature, p.value, p.value.adjust)
CRP_pheno_within_3_months_general_cortical_volumes_small$Group <- "General Cortical Volumes"

CRP_pheno_within_3_months_subcortical_volumes_small <- CRP_pheno_within_3_months_subcortical.csv %>%
  rename(feature = subcortical_volume) %>%
  select(feature, p.value, p.value.adjust)
CRP_pheno_within_3_months_subcortical_volumes_small$Group <- "Subortical Volumes"

CRP_pheno_within_3_months_bilateral_FA_tracts_small <- CRP_pheno_within_3_months_FA_bilateral_tracts.csv %>%
  rename(feature = FA_tract) %>%
  select(feature, p.value, p.value.adjust)
CRP_pheno_within_3_months_bilateral_FA_tracts_small$Group <- "FA Tracts"

CRP_pheno_within_3_months_bilateral_MD_tracts_small <- CRP_pheno_within_3_months_MD_bilateral_tracts.csv %>%
  rename(feature = MD_tract) %>%
  select(feature, p.value, p.value.adjust)
CRP_pheno_within_3_months_bilateral_MD_tracts_small$Group <- "MD Tracts"

CRP_pheno_within_3_months_unilateral_FA_tracts_small <- CRP_pheno_within_3_months_FA_unilateral_tracts.csv %>%
  rename(feature = FA_tract) %>%
  select(feature, p.value, p.value.adjust)
CRP_pheno_within_3_months_unilateral_FA_tracts_small$Group <- "FA Tracts"

CRP_pheno_within_3_months_unilateral_MD_tracts_small <- CRP_pheno_within_3_months_MD_unilateral_tracts.csv %>%
  rename(feature = MD_tract) %>%
  select(feature, p.value, p.value.adjust)
CRP_pheno_within_3_months_unilateral_MD_tracts_small$Group <- "MD Tracts"

CRP_pheno_within_3_months_FA_general_small <- CRP_pheno_within_3_months_FA_general.csv %>%
  rename(feature = general_DTI_FA_measure,  p.value.adjust = p.adjust) %>%
  select(feature, p.value, p.value.adjust)
CRP_pheno_within_3_months_FA_general_small$Group <- "General FA Measures"

CRP_pheno_within_3_months_MD_general_small <- CRP_pheno_within_3_months_MD_general.csv %>%
  rename(feature = general_DTI_MD_measure, p.value.adjust = p.adjust) %>%
  select(feature, p.value, p.value.adjust)
CRP_pheno_within_3_months_MD_general_small$Group <- "General MD Measures"
## Combine all dataframes 
CRP_pheno_within_3_months_all <- rbind(CRP_pheno_within_3_months_regional_cortical_volumes_small, CRP_pheno_within_3_months_general_cortical_volumes_small, CRP_pheno_within_3_months_subcortical_volumes_small, CRP_pheno_within_3_months_bilateral_FA_tracts_small, 
                                    CRP_pheno_within_3_months_bilateral_MD_tracts_small, CRP_pheno_within_3_months_unilateral_FA_tracts_small, CRP_pheno_within_3_months_unilateral_MD_tracts_small, CRP_pheno_within_3_months_FA_general_small, CRP_pheno_within_3_months_MD_general_small)

## If p.adjust = NA (no FDR correction) replace with pvalue
CRP_pheno_within_3_months_all$p.value.adjust[is.na(CRP_pheno_within_3_months_all$p.value.adjust)] <- CRP_pheno_within_3_months_all$p.value[is.na(CRP_pheno_within_3_months_all$p.value.adjust)]

## Log tramsform p-value 
CRP_pheno_within_3_months_all$p_minus_log_10 <- -log10(CRP_pheno_within_3_months_all$p.value)
## Plot results

jpeg("~/Desktop/PhD/projects/UKB_inflammation_imaging/output/summary_plots/all_CRP_pheno_within_3_months.jpeg", width = 1000, height = 1000)

if (nrow(CRP_pheno_within_3_months_all[which(CRP_pheno_within_3_months_all$p.value.adjust<0.05),])>0){
  ggplot(CRP_pheno_within_3_months_all, aes(x = Group, y = p_minus_log_10)) +
    geom_jitter(data = CRP_pheno_within_3_months_all[which(CRP_pheno_within_3_months_all$p.value.adjust>0.05),], 
                aes(colour =  Group), size = 5, alpha = 0.8, shape = 16,
                position="jitter") +
    geom_jitter(data = CRP_pheno_within_3_months_all[which(CRP_pheno_within_3_months_all$p.value.adjust<0.05),], 
                aes(colour =  Group), size = 5, alpha = 0.8, shape = 15,
                height= 0) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", 
               color = "red", size=1, alpha = 0.4) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y=element_blank()) +
    xlab("") + ylab("-log10(p)")
} else{
  ggplot(CRP_pheno_within_3_months_all, aes(x = Group, y = p_minus_log_10)) +
    geom_jitter(data = CRP_pheno_within_3_months_all[which(CRP_pheno_within_3_months_all$p.value.adjust>0.05),], 
                aes(colour =  Group), size = 5, alpha = 0.8, shape = 16,
                position="jitter") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", 
               color = "red", size=1, alpha = 0.4) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y=element_blank(), text = element_text(size = 20)) +
    xlab("") + ylab("-log10(p)")
}
dev.off()



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
## Plot all Serum CRP (wothin 3 months) results
##

glycA_general_cortical_volumes_small <- glycA_cortical_volume_general.csv %>%
  rename(feature = cortical_volume, p.value.adjust = p.adjust) %>%
  select(feature, p.value, p.value.adjust)
glycA_general_cortical_volumes_small$Group <- "General Cortical Volume"

glycA_general_cortical_thickness_small <- glycA_cortical_thickness_general.csv%>%
  rename(feature = cortical_thickness, p.value.adjust = p.adjust) %>%
  select(feature, p.value, p.value.adjust)
glycA_general_cortical_thickness_small$Group <- "General Cortical Thickness"

glycA_general_cortical_area_small <- glycA_cortical_area_general.csv %>%
  rename(feature = cortical_area, p.value.adjust = p.adjust) %>%
  select(feature, p.value, p.value.adjust)
glycA_general_cortical_area_small$Group <- "General Cortical Area"

glycA_regional_cortical_volume_small <- lme_cortical_volume.csv %>%
  rename(feature = cortical_volume) %>%
  select(feature, p.value, p.value.adjust)
glycA_regional_cortical_volume_small$Group <- "Regional Cortical Volume"

glycA_regional_cortical_area_small <- lme_cortical_area_glycA.csv %>%
  rename(feature = cortical_area) %>%
  select(feature, p.value, p.value.adjust)
glycA_regional_cortical_area_small$Group <- "Regional Cortical Area"

glycA_regional_cortical_thickness_small <- lme_cortical_thickness_glycA.csv %>%
  rename(feature = cortical_thickness) %>%
  select(feature, p.value, p.value.adjust)
glycA_regional_cortical_thickness_small$Group <- "Regional Cortical Thickness"

glycA_subcortical_volumes_small <- glycA_subcortical.csv %>%
  rename(feature = subcortical_volume) %>%
  select(feature, p.value, p.value.adjust)
glycA_subcortical_volumes_small$Group <- "Subortical Volumes"

glycA_FA_bilateral_tracts_small <- glycA_FA_bilateral_tracts.csv %>%
  rename(feature = FA_tract) %>%
  select(feature, p.value, p.value.adjust)
glycA_FA_bilateral_tracts_small$Group <- "FA Tracts"

glycA_MD_bilateral_tracts_small <- glycA_MD_bilateral_tracts.csv %>%
  rename(feature = MD_tract) %>%
  select(feature, p.value, p.value.adjust)
glycA_MD_bilateral_tracts_small$Group <- "MD Tracts"

glycA_FA_unilateral_tracts_small <- glycA_FA_unilateral_tracts.csv %>%
  rename(feature = FA_tract) %>%
  select(feature, p.value, p.value.adjust)
glycA_FA_unilateral_tracts_small$Group <- "FA Tracts"

glycA_FA_bilateral_tracts_hemisphere_interaction_small <- glycA_FA_bilateral_tracts_hemisphere_interaction.csv %>%
  rename(feature = FA_tract) %>%
  select(feature, p.value, p.value.adjust)
glycA_FA_bilateral_tracts_hemisphere_interaction_small$Group <- "FA Tracts"

glycA_MD_unilateral_tracts_small <- glycA_MD_unilateral_tracts.csv %>%
  rename(feature = MD_tract) %>%
  select(feature, p.value, p.value.adjust)
glycA_MD_unilateral_tracts_small$Group <- "MD Tracts"

glycA_FA_general_small <- glycA_FA_general.csv %>%
  rename(feature = general_DTI_FA_measure,  p.value.adjust = p.adjust) %>%
  select(feature, p.value, p.value.adjust)
glycA_FA_general_small$Group <- "General FA Measures"

glycA_MD_general_small <- glycA_MD_general.csv %>%
  rename(feature = general_DTI_MD_measure, p.value.adjust = p.adjust) %>%
  select(feature, p.value, p.value.adjust)
glycA_MD_general_small$Group <- "General MD Measures"
## Combine all dataframes 

glycA_all <- rbind(glycA_general_cortical_volumes_small, glycA_general_cortical_thickness_small, glycA_general_cortical_area_small, glycA_regional_cortical_volume_small, 
                   glycA_regional_cortical_area_small, glycA_regional_cortical_thickness_small, glycA_subcortical_volumes_small, glycA_FA_bilateral_tracts_small, 
                   glycA_MD_bilateral_tracts_small, glycA_FA_bilateral_tracts_hemisphere_interaction_small, glycA_FA_unilateral_tracts_small, glycA_MD_unilateral_tracts_small,
                   glycA_FA_general_small, glycA_MD_general_small)

## If p.adjust = NA (no FDR correction) replace with pvalue
glycA_all$p.value.adjust[is.na(glycA_all$p.value.adjust)] <- glycA_all$p.value[is.na(glycA_all$p.value.adjust)]

## Log tramsform p-value 
glycA_all$p_minus_log_10 <- -log10(glycA_all$p.value)
## Plot results

## Reorder feature type categories
glycA_all$Group <- factor(glycA_all$Group, levels = c('General Cortical Volume', 'General Cortical Area', 'General Cortical Thickness', 'Regional Cortical Volume', 'Regional Cortical Area',
                                                      'Regional Cortical Thickness', 'Subortical Volumes', 'General FA Measures', 'General MD Measures', 'FA Tracts', 'MD Tracts'))

jpeg("~/Desktop/PhD/projects/UKB_inflammation_imaging/output/summary_plots/all_glycA.jpeg", width = 1000, height = 1000)

if (nrow(glycA_all[which(glycA_all$p.value.adjust<0.05),])>0){
  ggplot(glycA_all, aes(x = Group, y = p_minus_log_10)) +
    geom_jitter(data = glycA_all[which(glycA_all$p.value.adjust>0.05),], 
                aes(colour =  Group), size = 5, alpha = 0.8, shape = 16,
                position="jitter") +
    geom_jitter(data = glycA_all[which(glycA_all$p.value.adjust<0.05),], 
                aes(colour =  Group), size = 5, alpha = 0.8, shape = 15,
                height= 0) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", 
               color = "red", size=1, alpha = 0.4) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y=element_blank()) +
    xlab("") + ylab("-log10(p)")
} else{
  ggplot(glycA_all, aes(x = Group, y = p_minus_log_10)) +
    geom_jitter(data = glycA_all[which(glycA_all$p.value.adjust>0.05),], 
                aes(colour =  Group), size = 5, alpha = 0.8, shape = 16,
                position="jitter") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", 
               color = "red", size=1, alpha = 0.4) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y=element_blank(), text = element_text(size = 20)) +
    xlab("") + ylab("-log10(p)")
}
dev.off()



