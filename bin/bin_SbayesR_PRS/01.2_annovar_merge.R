# Merge output from ANNOVAR with formatted CRP sum stats Hannah made
# ensure col names are in the correct order and named correctly for SBayesR
# This might also be a good time to do the QC whilst the sum stats are loaded in R.


setwd("/exports/eddie/scratch/s2093301/CRP_SBayesR_PRS/")

library(dplyr)

# Read in output from ANNOVAR:
# 2nd column = allele frequency
sumstats <- read.table("resources/annovar/ANNOVAR/CRP_formatted_MAF.hg19_EUR.sites.2015_08_dropped", header = F)

# Read in col names
colnames <- read.table( "resources/annovar/ANNOVAR/CRP_formatted.txt.avinput.header", header = T)
colnames(sumstats) <- c("Build", "freq", colnames(colnames))

# Check for duplicate SNPs
length( unique(sumstats$rsID) ) == nrow(sumstats)

# Add col for sample size "N"
sumstats$N <- 148164 # obtained from methods of Ligthart et al. 2018: "1KG (148,164 individuals from 49 studies) imputed-genotype GWASs"

# Check format and col names:
# SNP A1 A2 freq b se p N
sumstats <- sumstats %>%
			select(SNP = rsID, A1 = Allele1, A2 = Allele2, freq, b = Effect, se = StdErr, p = P, N = N,
				CHR = ChromNum, BP = Start) %>%
			filter(freq > 0.01) # Remove SNPs with MAF < 0.01



head(sumstats)
nrow(sumstats) # 3,722,473 SNPs

write.table(sumstats, "resources/CRP_formatted_MAF_QC.txt" , quote = F, row.names = F, col.names = T)

