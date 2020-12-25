#### Create Genotype Summary Table ####

library(gdsfmt)
library(SNPRelate)
library(plyr)
library(tidyverse)

options(stringsAsFactors=F)

# Set VCF name:
vcf.fn = "./03-callvariants/raw_variants.vcf"

# Convert ALL genotypes to GDS:
snpgdsVCF2GDS(vcf.fn, "./03-callvariants/raw_variants.gds", method="copy.num.of.ref")

## Open file for access:
genofile = snpgdsOpen("./03-callvariants/raw_variants.gds")

chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
pos <- read.gdsn(index.gdsn(genofile, "snp.position"))
snpid <- read.gdsn(index.gdsn(genofile, "snp.id"))

chr2 = ordered(chr, levels=c(1:19, 'X', 'Y'))
col.idx = order(chr2, pos)
chr.ord = chr[col.idx]
pos.ord = pos[col.idx]

# Get the genotype matrix
gmat = snpgdsGetGeno(genofile, sample.id=NULL, snp.id=NULL, snpfirstdim=NA, .snpread=NA, with.id=T, verbose=TRUE)

genotype = gmat$genotype[, col.idx]
#genotype.df <- as.data.frame(genotype, row.names = gmat$sample.id, col.names = paste(chr.ord, pos.ord, sep='.'))

samplenames = gmat$sample.id
locusnames = paste(chr.ord, pos.ord, sep='.')

genotype.df <- as.data.frame(t(genotype), col.names = samplenames, row.names = locusnames)
genotype.df[is.na(genotype.df)] <- 3
colnames(genotype.df) <- samplenames

#Isolate autosomes by removing the X & Y markers
Yloci <- locusnames %>% .[matches("[Y^]", vars=.)]
Xloci <- locusnames %>% .[matches("[X^]", vars=.)]

genotype.noY = genotype.df[ !(rownames(genotype.df) %in% Yloci), ]
#genotype.autosomes  = genotype.noY[ !(rownames(genotype.noY) %in% Xloci), ]

genotype_sum0 <- ldply(genotype.noY, function(r) sum(r=="0"))
genotype_sum1 <- ldply(genotype.noY, function(r) sum(r=="1"))
genotype_sum2 <- ldply(genotype.noY, function(r) sum(r=="2"))
genotype_sum3 <- ldply(genotype.noY, function(r) sum(r=="3"))
genotype_sumAll <- join_all(list(genotype_sum0,genotype_sum1,genotype_sum2, genotype_sum3), by='.id', type='left')
#colnames(genotype_sumAll)
names(genotype_sumAll)[1] <- "Sample"
names(genotype_sumAll)[2] <- "BB"
names(genotype_sumAll)[3] <- "AB"
names(genotype_sumAll)[4] <- "AA"
names(genotype_sumAll)[5] <- "Failed"

genotype_sumAll <- genotype_sumAll  %>% 
  mutate(Percent_Failed = 100*(1-((BB+AB+AA)/1640)), Percent_Failed = round(Percent_Failed, 2)) %>% 
  mutate(Percent_C57BL6J = 100*((2*AA+AB)/(2*(BB+AB+AA))), Percent_C57BL6J = round(Percent_C57BL6J, 2))
names(genotype_sumAll)[6] <- "Failed (%)"
names(genotype_sumAll)[7] <- "C57BL6J alleles (%)"

write_tsv(genotype_sumAll, "./04-deliverables/Genotype_Summary_Table.tsv", quote_escape = "none")
