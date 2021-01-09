# Pipeline for mouse genotyping assay

### Analyses implemented in the following study:

Andrews KR, Hunter SS, Torrevillas BK, Cespedes N, New DD, Fagnan MW, Luckhart S. A new SNP genotyping assay for speed congenics: combining flexibility, affordability, and power. In review.

## Requirements
### Programs
* Python v.2.7
* Fastqc
* MultiQC
* HTStream v1.1.0
* BWA v.0.7.17
* Samtools v.1.5 
* GATK v.4.1.3.0 
* R v.3.6.0
* R packages:
  * gdsfmt
  * SNPRelate
  * pheatmap
  * RColorBrewer
  * plyr
  * tidyverse

### Scripts
* 01-cleaning.py  
* 02-map_reads.py  
* 03.0-gatk.py  
* 03.1-combine-gatk.py  
* 03.2-call-variants.py  
* 04-analyze_VCF.R  
* 05-GenotypeSummary.R  

### Other Files
* Raw fastq.gz files, demultiplexed by sample and stored in ./00-RawData
* Reference genome (indexed with BWA), stored as ./ref/reference.fasta
* Bed file of SNP targets: ST2181G_1_target_ST2181G_1_2053-noextend.bed  
* Illumina adapter sequences: adapters.fa  

## Genotyping Pipeline

### Quality assessment  

```
mkdir -p 00-fastqc
fastqc -o ./00-fastqc ./00-RawData/*
multiqc -i fastqc ./00-fastqc/
```
### Cleaning
```
python 01-cleaning.py
bash 01-cleaning_commands.sh
```
### Mapping
```
python 02-map_reads.py
bash 02-mapping_commands.sh
```
### Quality assessment of mapping  

```
multiqc -d ./02-mapped/ -i Mapping -o ./02-mapped/
```
### Build samtools stats for all bam files
```
for f in ./02-mapped/*.bam
do
samtools stats $f > $f.stats &
done
```
### Index all bam files
```
cd ./02-mapped
for f in *.bam
do
samtools index $f &
done
cd ../
```
### Genotype with GATK ##
```
python 03.0-gatk.py
bash 03.0-gatk_command.sh

python 03.1-combine-gatk.py
bash 03.1-gatk_combine.sh

python 03.2-call-variants.py
bash 03.2-gatk_call.sh
```
### Calculate depth by SNP
```
samtools depth -a -b ST2181G_1_target_ST2181G_1_2053-noextend.bed ./02-mapped/*.bam > target_depth.tsv
ls -la ./02-mapped/*.bam > samples_bam.txt
```
### Analyze VCF
```
mkdir -p 04-deliverables
Rscript 04-analyze_VCF.R
```
Outputs:   
* Genotypes for each sample and SNP: genotype.tsv  
* Plots of coverage depth per SNP per sample: coverage_assessment.pdf  
* Plot of genotypes for each sample and SNP: genotype_plots-ordered.png  
* List of poorly performing SNPs (SNPs with coverage < 5 in more than 30 samples): poor_targets.tsv  

```
Rscript 05-GenotypeSummary.R
```
Outputs:  
*  Ancestry proportions for each sample: Genotype_Summary_Table.tsv



