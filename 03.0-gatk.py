####### Revised for single-bam GVCF mode to see if we can get it to emit all sites #####
from glob import glob
import os
from os.path import join as jp

os.system("mkdir -p 03-callvariants")

ref = './ref/MM10_combined.fasta'
gatk = 'gatk --java-options "-Xmx4g"'
bed = 'ST2181G_1_target_ST2181G_1_2053-noextend.bed'

outf = open("03.0-gatk_command.sh", 'w')

for bam in glob("./02-mapped/*.bam"):
    s = bam.split('/')[-1].replace('.bam','')
    out_gvcf = './03-callvariants/' + s + '.g.vcf'
    cmd = gatk + ' HaplotypeCaller --output-mode EMIT_ALL_SITES --all-site-pls true --ERC GVCF '
    cmd += ' -R ' + ref + ' -O ' + out_gvcf + " --max-assembly-region-size 500 -L " + bed
    cmd += " -I " + bam + " &> ./03-callvariants/" + s + "_gvcf.log"
    outf.write(cmd+'\n')

outf.close()


