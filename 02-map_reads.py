from glob import glob
import os
from os.path import join as jp

mapping = open("02-mapping_commands.sh", 'w')
ref = './ref/MM10_combined.fasta'
out_prefix = './02-mapped/'
os.system('mkdir -p 02-mapped')

for r1 in glob("./01-cleaned/*_R1.fastq.gz"):
    r2 = r1.replace("_R1", "_R2")
    s = r1.split('/')[-1].replace("_R1.fastq.gz", '')
    cmd = "bwa mem -t 50 -R '@RG\\tID:bwa\\tSM:" + s + "\\tPL:ILLUMINA' " + ref + ' ' + r1 + ' ' + r2
    cmd += ' 2> ' + jp(out_prefix, s) + '.log' + ' | samtools view -b - | '
    cmd += 'samtools sort - > ' + jp(out_prefix, s) + '.bam'
    mapping.write(cmd+'\n')
mapping.close()

