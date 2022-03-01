PACKAGES=/exports/igmm/eddie/GenScotDepression/users/hcasey/packages
DIR=/exports/eddie/scratch/s2093301/CRP_SBayesR_PRS

## Include MHC

$PACKAGES/GCTB/gctb_2.02_Linux/gctb --sbayes R \
     --mldm $DIR/resources/LD_matrices/ukb_50k_bigset_2.8M/ukb50k_2.8M_shrunk_sparse.new.mldmlist \
     --pi 0.95,0.02,0.02,0.01 \
     --gamma 0.0,0.01,0.1,1 \
     --ambiguous-snp \
     --impute-n \
     --gwas-summary $DIR/resources/CRP_formatted_MAF_QC_unambig.txt \
     --chain-length 10000 \
     --burn-in 2000 \
     --out-freq 10 \
     --out $DIR/output/CRP_summstats.SBayesR

## Exclude MHC

$PACKAGES/GCTB/gctb_2.02_Linux/gctb --sbayes R \
     --mldm $DIR/resources/LD_matrices/ukb_50k_bigset_2.8M/ukb50k_2.8M_shrunk_sparse.new.mldmlist \
     --pi 0.95,0.02,0.02,0.01 \
     --gamma 0.0,0.01,0.1,1 \
     --ambiguous-snp \
     --impute-n \
     --gwas-summary $DIR/resources/CRP_formatted_MAF_QC_unambig.txt \
     --chain-length 10000 \
     --burn-in 2000 \
     --out-freq 10 \
     --out $DIR/output/CRP_summstats_noMHC.SBayesR