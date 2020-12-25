#### Output genotypes and plots of per-sample depth of coverage and genotypes ####

library("gdsfmt")
library("SNPRelate")
library("pheatmap")
library("RColorBrewer")

options(stringsAsFactors=F)

# Set VCF name:
vcf.fn = "./03-callvariants/raw_variants.vcf"

# Convert ALL genotypes to GDS:
snpgdsVCF2GDS(vcf.fn, "./03-callvariants/raw_variants.gds", method="copy.num.of.ref")

## Open file for access:
genofile = snpgdsOpen("./03-callvariants/raw_variants.gds")

# Extract the genotype information and reorder by genome position:
# There are four possible values stored in the variable genotype: 0, 1, 2 and 3. 
# For bi-allelic SNP sites, “0” indicates two B alleles, 
#                           “1” indicates one A allele and one B allele, 
#                           “2” indicates two A alleles, and 
#                           “3” is a missing genotype. 
#   For multi-allelic sites, it is a count of the reference allele (3 meaning no call). 
#   “Bit2” indicates that each byte encodes up to four SNP genotypes since one byte consists of eight bits.
## This code helps with ordering SNPs by genome location for plotting:
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
rownames(genotype) = gmat$sample.id
colnames(genotype) = paste(chr.ord, pos.ord, sep='.')


# Read in targets and coverage data and see how many were called:
targets = read.table("ST2181G_1_target_ST2181G_1_2053-noextend.bed", as.is=T, sep='\t')
coverage = read.table("target_depth.tsv", header=F, as.is=T, sep='\t')

bams = read.table("samples_bam.txt", as.is=T, header=F)
colnames(coverage) = c("chr", "pos", gsub(".bam", '', gsub('./02-mapped/', '', bams$V9)))


## Create identifiers for matching up targets, SNPs, and coverage:
targetpos = paste(targets$V1, targets$V3, sep=':')
snppos = paste(paste("chr", chr.ord, sep=''), pos.ord, sep=":")
coveragepos = paste(coverage$chr, coverage$pos, sep=":")

### How many of the snps are in the targets list?
table(snppos %in% targetpos)

### Get coverage values for ONLY the targeted loci:
coverage.targeted = coverage[coveragepos %in% targetpos, ]
rownames(coverage.targeted) = coveragepos[coveragepos %in% targetpos]
#colnames(coverage.targeted) = c("chr","pos",sort(gmat$sample.id))

### Make a coverage assessment plot: 
pdf(file="./04-deliverables/coverage_assessment.pdf", w=12, h=8)
   par(las=2)
   par(mar=c(6,5,4,2))
   boxplot(coverage.targeted[,3:50]+1, ylab='Mapping coverage')
   abline(h=10, col='red')

   boxplot(log2(coverage.targeted[,3:50]+1), ylab='log2 mapping coverage')
   abline(h=log2(10), col='red')

   #breaksList = c(0, 20, 40, 60, max(coverage.targeted[,3:50]))
   pheatmap(coverage.targeted[,3:50], 
      cluster_rows=F, cluster_cols=F, show_rownames=F)
dev.off()


## Write a table of poor targets:
# Coverage < 5 in more than 30 samples
write.table(file="./04-deliverables/poor_targets.tsv", row.names=F, col.names=F, sep='\t', 
   coverage.targeted[rowSums(coverage.targeted[,3:ncol(coverage.targeted)] < 5) > 30, ])


# Reorder genotype matrix, rows by strain, cols by chromosome and loc:

# Create groups for colums and rows
ordered(chr.ord, levels=c(1:19, 'X', 'Y'))

annotation_col = data.frame(chr=ordered(chr.ord, levels=c(1:19, 'X', 'Y')))
rownames(annotation_col) = colnames(genotype)

png(file="./04-deliverables/genotype_plots-ordered.png", w=32000, h=1200)
pheatmap(genotype,  cluster_rows=F, cluster_cols=F,
   na_col="black",
   show_rownames=T, display_numbers=T, number_format="%i", number_color="black", fontsize_number=12, 
   fontsize_row=14, fontsize_col=14)
dev.off()

## Write the calls out to a spreadsheet##
write.table(file="./04-deliverables/genotype.tsv", data.frame(chr.pos = colnames(genotype), t(genotype)), row.names=F, sep='\t', quote=F)

