
## Annotate Sumstats with MAF

. /etc/profile.d/modules.sh
module load igmm/apps/R/3.6.1

# 1. Download ANNOVAR from:
# https://annovar.openbioinformatics.org/en/latest/user-guide/download/
# Wang K, Li M, Hakonarson H. ANNOVAR: Functional annotation of genetic variants from next-generation sequencing data, Nucleic Acids Research, 38:e164, 2010
# cd /exports/igmm/eddie/GenScotDepression/amelia/packages
# wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz
# tar xvfz annovar.latest.tar.gz
PATH=$PATH:/exports/eddie/scratch/s2093301/CRP_SBayesR_PRS/resources/annovar
cd -

# 2. use ANNOVAR's command below to download relevant reference files from 1000 genomes project (August 2015 European ancestry used here)
# mkdir ANNOVAR
cd ANNOVAR
annotate_variation.pl -downdb 1000g2015aug humandb -buildver hg19
cd -

# 3. Format GWAS sum stats for ANNOVAR to deal with
# http://annovar.openbioinformatics.org/en/latest/user-guide/input/
#"The first five space or tab delimited fields are Chromosome ("chr" prefix is optional), Start, End, Reference Allele, Alternative Allele. The rest of the columns are completely optional"
# text format with no quotes, no row IDs, no column headings
# first line of example input file is: 1        948921  948921  T       C
## Copy over Amelia's annovar.R script 
#cp annovar.R /exports/eddie/scratch/s2093301/CRP_SBayesR_PRS/resources/annovar/
Rscript annovar.R

# 4. run the annovar_variation.pl command with -filter option, which annotates the provided SNP with allele frequencies taken from the 1000 genome datatbase (August 2015 release, European ancestry, hg19 genome build)

# it outputs 3 files all starting with "CRP_formatted_MAF.txt" to the working directory
# ..._dropped: SNP annotated with allele frequency (1st column = ref genome, 2nd column = allele frequency)
# ...filtered: SNP which cannot be annotated (presumably because not in database, but not sure)
# .invalid_input: SNP which cannot be annotated because they are insertions or deletions
# for more info: https://doc-openbio.readthedocs.io/projects/annovar/en/latest/user-guide/filter/

annotate_variation.pl -filter \
 -dbtype 1000g2015aug_eur \
 -buildver hg19 \
 -out ANNOVAR/CRP_formatted_MAF ANNOVAR/CRP_formatted.txt.avinput ANNOVAR/humandb/

# 4. run the annovar_variation.pl command with -filter option, which annotates the provided SNP with allele frequencies taken from the 1000 genome datatbase (August 2015 release, European ancestry, hg19 genome build)

# it outputs 3 files all starting with "CRP_formatted_MAF.txt" to the working directory
# ..._dropped: SNP annotated with allele frequency (1st column = ref genome, 2nd column = allele frequency)
# ...filtered: SNP which cannot be annotated (presumably because not in database, but not sure)
# .invalid_input: SNP which cannot be annotated because they are insertions or deletions
# for more info: https://doc-openbio.readthedocs.io/projects/annovar/en/latest/user-guide/filter/

annotate_variation.pl -filter \
 -dbtype 1000g2015aug_eur \
 -buildver hg19 \
 -out ANNOVAR/CRP_formatted_MAF ANNOVAR/CRP_formatted.txt.avinput ANNOVAR/humandb/

head ANNOVAR/CRP_formatted_MAF.hg19_EUR.sites.2015_08_dropped

# 5. Merge output from ANNOVAR with formatted CRP sum stats Hannah made
# ensure col names are in the correct order and named correctly for SBayesR
# This might also be a good time to do the QC whilst the sum stats are loaded in R.
## copy over Amelia's annovar_merge.R script
Rscript annovar_merge.R

 Duplicate SNPs

cat /exports/eddie/scratch/s2093301/CRP_SBayesR_PRS/resources/CRP_sumstats.txt |\
awk '{seen[$6]++; if(seen[$6]==1){ print}}' > /exports/eddie/scratch/s2093301/CRP_SBayesR_PRS/resources/CRP_sumstats.nodup.txt



## remove ambigious SNPs
cat /exports/eddie/scratch/s2093301/CRP_SBayesR_PRS/resources/CRP_sumstats.nodup.txt |\
awk '!( ($1=="A" && $2=="T") || \
        ($1=="T" && $2=="A") || \
        ($1=="G" && $2=="C") || \
        ($1=="C" && $2=="G")) {print}' > /exports/eddie/scratch/s2093301/CRP_SBayesR_PRS/resources/CRP_sumstats.QC.txt
