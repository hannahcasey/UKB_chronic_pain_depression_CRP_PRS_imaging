# UKB_chronic_pain_depression_CRP_PRS_imaging
Association analysis of CRP PRS and imaging features in UKB.

## bin direcctory:

### bin_SbayesR_PRS
Scripts used to calculate PRS of CRP in UKB using C + T and SBayesR methods (w/ and w/o MHC)

### bin_imaging_association
Scripts used to prepare UKB data and perform analysis of association analysis between inflammation markers and imaging features in UKB.

#### UKB_inflammation_imaging_prep.R
Prepare CRP PRS, glycA, covariates and imaging measures into single dataframe

#### CRP_PRS_imaging_assoc.R
Perform association analysis between inflammation marker and imaging features.

### CRP_PRS_imaging_association_summary_plots.R
Plot summary of imaging association analysis using manhattan plots

## func directory:
Contains function used in scripts

### lme_Shen.R
Function written by Shen used to perfrom LME modelling


