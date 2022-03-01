#!/bin/sh
#$ -N keep_imaging
#$ -o /exports/eddie/scratch/s2093301/CRP_PRS/out/
#$ -e /exports/eddie/scratch/s2093301/CRP_PRS/err/
#$ -l h_vmem=16G


DIR=/exports/eddie/scratch/s2093301/CRP_SBayesR_PRS/resources


# Keep participants who have only imaging available

/exports/igmm/eddie/GenScotDepression/local/bin/plink2 \
        --bfile $DIR/ukb_imp_v3.qc \
        --keep $DIR/UKB_imaging.keep \
        --make-bed \
        --out $DIR/ukb_imp_v3_imaging.qc