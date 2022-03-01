#$ -N target_qc
#$ -l h_vmem=4G
# -pe sharedmem 4
#$ -l h_rt=4:00:00
#$ -e logs
#$ -o logs
#$ -cwd

# Example job to QC UKB imputed PGEN data for PRS target 

IMPV3=/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/impv3
DIR=/exports/eddie/scratch/s2093301/CRP_SBayesR_PRS/resources

## Extract list of SNPs in base data GWAS
cat $DIR/CRP_formatted_MAF_QC_unambig.txt | tail -n +2 | awk '{print $1}' > $DIR/CRP.snps


# European ancestries analysis

/exports/igmm/eddie/GenScotDepression/local/bin/plink2 \
--pfile $IMPV3/pgen/ukb_imp_v3.autosome.qc.mac100.info9.hwe10 \
--keep /exports/igmm/eddie/GenScotDepression/data/ukb/genetics/input_filters/v2/ukb_sqc_qc_WhiteBritishPCs_addPrunedRels_noPGC_noGenScot_v2.id \
--extract $DIR/CRP.snps \
--maf 0.01 \
--mac 100 \
--geno 0.02 \
--mind 0.02 \
--make-bed \
--out $DIR/ukb_imp_v3.qc \
--memory 15360 \
--threads 4