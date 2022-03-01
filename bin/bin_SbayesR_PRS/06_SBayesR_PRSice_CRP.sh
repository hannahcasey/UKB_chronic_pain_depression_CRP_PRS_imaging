#!/bin/sh
#$ -e /exports/eddie/scratch/s2093301/CRP_SBayesR_PRS/err/
#$ -o /exports/eddie/scratch/s2093301/CRP_SBayesR_PRS/out/
#$ -l h_vmem=32G
#$ -pe sharedmem 8

PACKAGES=/exports/igmm/eddie/GenScotDepression/users/hcasey/packages
DIR=/exports/eddie/scratch/s2093301/CRP_SBayesR_PRS

    
# -------------------------------------------------
# Create PRSs
# -------------------------------------------------
# From SBayesR output:
# Insert name of beta, A1, A2, bp, chr, SNP, BETA cols

# CRP with MHC
Rscript /exports/igmm/eddie/GenScotDepression/users/hcasey/packages/PRSice.v2.3.5/PRSice.R \
    --prsice /exports/igmm/eddie/GenScotDepression/users/hcasey/packages/PRSice.v2.3.5/PRSice_linux \
    --dir . \
    --base $DIR/output/CRP_summstats_noMHC.SBayesR.snpRes \
    --target $DIR/resources/ukb_imp_v3_imaging.qc \
    --beta \
    --A1 A1 \
    --A2 A2 \
    --bp Position \
    --chr Chrom \
    --pvalue PIP \
    --snp Name \
    --stat A1Effect \
    --no-default \
    --binary-target F\
    --pheno-file resources/UKB_CRP_imaging.pheno \
    --pheno-col CRP \
    --cov resources/UKB_CRP_imaging.cov \
    --fastscore \
    --bar-levels 1 \
    --all-score  \
    --no-clump \
    --thread 8 \
    --out output/SBayesR_CRP/CRP_PRS_SBayesR


# CRP without MHC
Rscript /exports/igmm/eddie/GenScotDepression/users/hcasey/packages/PRSice.v2.3.5/PRSice.R \
    --prsice /exports/igmm/eddie/GenScotDepression/users/hcasey/packages/PRSice.v2.3.5/PRSice_linux \
    --dir . \
    --base $DIR/output/CRP_summstats_noMHC.SBayesR.snpRes \
    --target $DIR/resources/ukb_imp_v3_imaging.qc \
    --beta \
    --A1 A1 \
    --A2 A2 \
    --bp Position \
    --chr Chrom \
    --pvalue PIP \
    --snp Name \
    --stat A1Effect \
    --no-default \
    --binary-target F\
    --pheno-file resources/UKB_CRP_imaging.pheno \
    --pheno-col CRP \
    --cov resources/UKB_CRP_imaging.cov \
    --fastscore \
    --bar-levels 1 \
    --all-score  \
    --no-clump \
    --thread 8 \
    --out output/SBayesR_CRP/CRP_PRS_SBayesR
