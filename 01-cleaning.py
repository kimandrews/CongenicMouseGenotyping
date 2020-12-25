from glob import glob
import os

cleaning = open("01-cleaning_commands.sh", 'w')
os.system("mkdir -p 01-cleaned")

for r1 in glob("./00-RawData/*_R1_*.fastq.gz"):
    r2 = r1.replace("_R1_", "_R2_")
    s = r1.split('/')[-1].replace("_L001_R1_001.fastq.gz", '')
    log = "./01-cleaned/" + s + "_stats.log"
    cmd =  "hts_SuperDeduper -e 10000 -L " + log + " -1 " + r1 + " -2 " + r2 + " | "
    cmd += "hts_AdapterTrimmer -m 90 -AL " + log + " | "
    cmd += "hts_CutTrim -a 40 -m 90 -AL " + log + " | "  # Cut off probe
    #cmd = "hts_NTrimmer -n -m 100 -L " + log + " -1 " + r1 + " -2 " + r2 + " -O |"
    #cmd += "hts_Overlapper -m 300 -n -S -o 15 -e 0.2 -A -L " + log + " -O |"
    cmd += "hts_SeqScreener -k 12 -AL " + log + " | "
    cmd += "hts_SeqScreener -AL " + log + " -k 15 -x .01 --seq adapters.fa -fgp ./01-cleaned/" + s
    cleaning.write(cmd+'\n')

cleaning.close()
