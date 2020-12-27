from glob import glob
import os
from os.path import join as jp

ref = './ref/reference.fasta'
gatk = 'gatk --java-options "-Xmx4g"'
bed = 'ST2181G_1_target_ST2181G_1_2053-noextend.bed'

# gatk GenomicsDBImport \
#     -V data/gvcfs/mother.g.vcf \
#     -V data/gvcfs/father.g.vcf \
#     -V data/gvcfs/son.g.vcf \
#     --genomicsdb-workspace-path my_database \
#     --intervals chr20,chr21


######### Generate multi-sample VCF ######
outf = open("03.2-gatk_call.sh", 'w')
cmd = gatk + ' GenotypeGVCFs --all-sites true -R ' + ref + ' -V gendb://./03-callvariants/combined.db -L ' + bed 
cmd += " -O ./03-callvariants/raw_variants.vcf"
outf.write(cmd+'\n')
outf.close()
