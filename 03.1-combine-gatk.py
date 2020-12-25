######### Combine VCF using GenomicsDBImport ####
## See https://software.broadinstitute.org/gatk/documentation/article?id=11813
from glob import glob
import os
from os.path import join as jp

ref = './ref/MM10_combined.fasta'
gatk = 'gatk --java-options "-Xmx4g"'
bed = 'ST2181G_1_target_ST2181G_1_2053-noextend.bed'

outf = open("03.1-gatk_combine.sh", 'w')
cmd = gatk + ' GenomicsDBImport --genomicsdb-workspace-path ./03-callvariants/combined.db -L ' + bed
for gvcf in glob("./03-callvariants/*.g.vcf"):
    cmd += " -V " + gvcf
outf.write(cmd + '\n')
outf.close()



