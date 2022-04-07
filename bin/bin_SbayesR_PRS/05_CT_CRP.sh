#$ -N CT_prsice_CRP
#$ -l h_vmem=2G
#$ -pe sharedmem 8
#$ -l h_rt=24:00:00
#$ -e err
#$ -o out
#$ -cwd

# initialize module
. /etc/profile.d/modules.sh

module add igmm/apps/R/4.0.3
module add igmm/apps/PRSice/2.1.11

DIR=/exports/eddie/scratch/s2093301/CRP_SBayesR_PRS

## Run PRSice with MHC 

Rscript /exports/igmm/eddie/GenScotDepression/users/hcasey/packages/PRSice.v2.3.5/PRSice.R \
    --prsice /exports/igmm/eddie/GenScotDepression/users/hcasey/packages/PRSice.v2.3.5/PRSice_linux \
    --base $DIR/resources/CRP_formatted_MAF_QC_unambig.txt \
    --target $DIR/resources/ukb_imp_v3.qc \
    --beta \
    --A1 A1 \
    --A2 A2 \
    --bp BP \
    --chr CHR \
    --pvalue p \
    --snp SNP \
    --stat b \
    --binary-target F\
    --pheno-file $DIR/resources/UKB_CRP_log_imaging.pheno \
    --cov-file $DIR/resources/UKB_CRP_imaging.cov \
    --fastscore \
    --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \
    --all-score  \
    --clump-r2 0.25 \
    --clump-kb 500 \
    --thread 8 \
    --print-snp \
    --out $DIR/output/CRP_PRSice_MHC


## Run PRSice witout MHC 

Rscript /exports/igmm/eddie/GenScotDepression/users/hcasey/packages/PRSice.v2.3.5/PRSice.R \
    --prsice /exports/igmm/eddie/GenScotDepression/users/hcasey/packages/PRSice.v2.3.5/PRSice_linux \
    --base $DIR/resources/CRP_formatted_MAF_QC_unambig.txt \
    --target $DIR/resources/ukb_imp_v3.qc \
    --beta \
    --A1 A1 \
    --A2 A2 \
    --bp BP \
    --chr CHR \
    --pvalue p \
    --snp SNP \
    --stat b \
    --binary-target F\
    --pheno-file $DIR/resources/UKB_CRP_log_imaging.pheno \
    --cov-file $DIR/resources/UKB_CRP_imaging.cov \
    --fastscore \
    --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \
    --all-score  \
    --clump-r2 0.25 \
    --clump-kb 500 \
    --thread 8 \
    --print-snp \
    --exclude $DIR/resources/ukb_imp_v3_MHC.qc.snplist \
    --out $DIR/output/CRP_PRSice_noMHC

