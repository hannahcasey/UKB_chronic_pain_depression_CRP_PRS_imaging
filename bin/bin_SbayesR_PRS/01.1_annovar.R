# Format GWAS sum stats for ANNOVAR to deal with
# http://annovar.openbioinformatics.org/en/latest/user-guide/input/
#"The first five space or tab delimited fields are:
# Chromosome ("chr" prefix is optional), Start, End, Reference Allele, Alternative Allele. The rest of the columns are completely optional"
# text format with no quotes, no row IDs, no column headings
# first line of example input file is: 1	948921	948921	T	C

setwd("/exports/eddie/scratch/s2093301/CRP_SBayesR_PRS/")

library(dplyr)

CRP <- read.table("resources/CRP_formatted.txt", header = T)

CRP_formatted <- CRP %>%
	mutate(Start = BP) %>%
	rename(End = BP) %>%
	select(c(ChromNum, Start, End, Allele1, Allele2, Effect, StdErr, P, rsID)) 

head(CRP_formatted)

write.table(CRP_formatted, "resources/annovar/ANNOVAR/CRP_formatted.txt.avinput", row.names = F, col.names = F, quote = F)

write.table(slice(CRP_formatted, 1), "resources/annovar/ANNOVAR/CRP_formatted.txt.avinput.header", row.names = F, col.names = T, quote = F)

