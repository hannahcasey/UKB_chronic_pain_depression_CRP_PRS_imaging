#$ -N target_MHC_get_SNP_list
#$ -l h_vmem=4G
# -pe sharedmem 4
#$ -l h_rt=4:00:00
#$ -e logs
#$ -o logs
#$ -cwd

DIR=/exports/eddie/scratch/s2093301/CRP_SBayesR_PRS/resources

## Get list of SNPs in MHC region

/exports/igmm/eddie/GenScotDepression/local/bin/plink2 \
    --bfile $DIR/ukb_imp_v3.qc\
    --chr 6\
    --from-kb 28000\
    --to-kb 34000\
    --write-snplist\
    --out $DIR/ukb_imp_v3_MHC.qc
