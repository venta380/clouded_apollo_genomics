import sys
import os
import string
import pandas
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
import math
import time
import sys
import subprocess



def run_SRA_BWA_GATK_lib_1(Sample_ID, genome, out_dir, fastq_1, fastq_2, outFile):
	outFileName = outFile
	outFile = open(outFileName, "w")
	outFile.write('#!/bin/sh'+'\n')
	outFile.write('#PBS -S /bin/bash'+'\n')
	outFile.write('#PBS -q batch'+'\n')
	outFile.write('#PBS -l nodes=1:ppn=20'+'\n')
	outFile.write('#PBS -N mapping-'+Sample_ID+'\n')
	outFile.write('#PBS -l walltime=80:00:00'+'\n')
	outFile.write('#PBS -o mapping-'+Sample_ID+'.out'+'\n')
	outFile.write('#PBS -e mapping-'+Sample_ID+'.err'+'\n')
	outFile.write('#PBS -l mem=100gb'+'\n')
	outFile.write('\n')
	outFile.write("module load SAMtools"+'\n')
	outFile.write("module load BWA"+'\n')
	outFile.write("module load Stampy"+'\n')
	outFile.write('wait'+'\n')
	outFile.write('cd $TMPDIR'+'\n')
	outFile.write('\n')
	outFile.write("zcat "+fastq_1+" > reads.fastq"+'\n')
	outFile.write("zcat "+fastq_2+" >> reads.fastq"+'\n')
	outFile.write("stampy.py -g "+genome+" -h "+genome+" -t 20 --fast --overwrite --readgroup=ID:group1,SM:"+Sample_ID+",PL:illumina,LB:lib1,PU:unit1  --substitutionrate=0.025 -o "+Sample_ID+".sam -M reads.fastq"+'\n')
	outFile.write("samtools sort $TMPDIR/"+Sample_ID+".sam -o "+out_dir+Sample_ID+".bam"+'\n') 
	outFile.write("samtools index "+out_dir+Sample_ID+".bam"+'\n') 
	outFile.write('\n')





samples_1_libraries=set(['S121','S122','S123','S124','S125','S126','S127','S128','S129','S130','S131','S132','S133','S134','S135','S136','S137','S138','S139','S140','S141','S142','S143','S144','S145','S146','S147','S148','S149','S150','S151','S152','S153','S154','S155','S156','S157','S158','S159','S160'])

files={}
files['S121']=['/scratch/vt20265/clouded_apollo/P21213_101_S121_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_101_S121_L003_R2_001.fastq.gz']
files['S122']=['/scratch/vt20265/clouded_apollo/P21213_102_S122_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_102_S122_L003_R2_001.fastq.gz']
files['S123']=['/scratch/vt20265/clouded_apollo/P21213_103_S123_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_103_S123_L003_R2_001.fastq.gz']
files['S124']=['/scratch/vt20265/clouded_apollo/P21213_104_S124_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_104_S124_L003_R2_001.fastq.gz']
files['S125']=['/scratch/vt20265/clouded_apollo/P21213_105_S125_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_105_S125_L003_R2_001.fastq.gz']
files['S126']=['/scratch/vt20265/clouded_apollo/P21213_106_S126_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_106_S126_L003_R2_001.fastq.gz']
files['S127']=['/scratch/vt20265/clouded_apollo/P21213_107_S127_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_107_S127_L003_R2_001.fastq.gz']
files['S128']=['/scratch/vt20265/clouded_apollo/P21213_108_S128_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_108_S128_L003_R2_001.fastq.gz']
files['S129']=['/scratch/vt20265/clouded_apollo/P21213_109_S129_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_109_S129_L003_R2_001.fastq.gz']
files['S130']=['/scratch/vt20265/clouded_apollo/P21213_110_S130_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_110_S130_L003_R2_001.fastq.gz']
files['S131']=['/scratch/vt20265/clouded_apollo/P21213_111_S131_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_111_S131_L003_R2_001.fastq.gz']
files['S132']=['/scratch/vt20265/clouded_apollo/P21213_112_S132_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_112_S132_L003_R2_001.fastq.gz']
files['S133']=['/scratch/vt20265/clouded_apollo/P21213_113_S133_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_113_S133_L003_R2_001.fastq.gz']
files['S134']=['/scratch/vt20265/clouded_apollo/P21213_114_S134_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_114_S134_L003_R2_001.fastq.gz']
files['S135']=['/scratch/vt20265/clouded_apollo/P21213_115_S135_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_115_S135_L003_R2_001.fastq.gz']
files['S136']=['/scratch/vt20265/clouded_apollo/P21213_116_S136_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_116_S136_L003_R2_001.fastq.gz']
files['S137']=['/scratch/vt20265/clouded_apollo/P21213_117_S137_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_117_S137_L003_R2_001.fastq.gz']
files['S138']=['/scratch/vt20265/clouded_apollo/P21213_118_S138_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_118_S138_L003_R2_001.fastq.gz']
files['S139']=['/scratch/vt20265/clouded_apollo/P21213_119_S139_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_119_S139_L003_R2_001.fastq.gz']
files['S140']=['/scratch/vt20265/clouded_apollo/P21213_120_S140_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_120_S140_L003_R2_001.fastq.gz']
files['S141']=['/scratch/vt20265/clouded_apollo/P21213_121_S141_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_121_S141_L003_R2_001.fastq.gz']
files['S142']=['/scratch/vt20265/clouded_apollo/P21213_122_S142_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_122_S142_L003_R2_001.fastq.gz']
files['S143']=['/scratch/vt20265/clouded_apollo/P21213_123_S143_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_123_S143_L003_R2_001.fastq.gz']
files['S144']=['/scratch/vt20265/clouded_apollo/P21213_124_S144_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_124_S144_L003_R2_001.fastq.gz']
files['S145']=['/scratch/vt20265/clouded_apollo/P21213_125_S145_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_125_S145_L003_R2_001.fastq.gz']
files['S146']=['/scratch/vt20265/clouded_apollo/P21213_126_S146_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_126_S146_L003_R2_001.fastq.gz']
files['S147']=['/scratch/vt20265/clouded_apollo/P21213_127_S147_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_127_S147_L003_R2_001.fastq.gz']
files['S148']=['/scratch/vt20265/clouded_apollo/P21213_128_S148_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_128_S148_L003_R2_001.fastq.gz']
files['S149']=['/scratch/vt20265/clouded_apollo/P21213_129_S149_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_129_S149_L003_R2_001.fastq.gz']
files['S150']=['/scratch/vt20265/clouded_apollo/P21213_130_S150_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_130_S150_L003_R2_001.fastq.gz']
files['S151']=['/scratch/vt20265/clouded_apollo/P21213_131_S151_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_131_S151_L003_R2_001.fastq.gz']
files['S152']=['/scratch/vt20265/clouded_apollo/P21213_132_S152_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_132_S152_L003_R2_001.fastq.gz']
files['S153']=['/scratch/vt20265/clouded_apollo/P21213_133_S153_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_133_S153_L003_R2_001.fastq.gz']
files['S154']=['/scratch/vt20265/clouded_apollo/P21213_134_S154_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_134_S154_L003_R2_001.fastq.gz']
files['S155']=['/scratch/vt20265/clouded_apollo/P21213_135_S155_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_135_S155_L003_R2_001.fastq.gz']
files['S156']=['/scratch/vt20265/clouded_apollo/P21213_136_S156_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_136_S156_L003_R2_001.fastq.gz']
files['S157']=['/scratch/vt20265/clouded_apollo/P21213_137_S157_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_137_S157_L003_R2_001.fastq.gz']
files['S158']=['/scratch/vt20265/clouded_apollo/P21213_138_S158_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_138_S158_L003_R2_001.fastq.gz']
files['S159']=['/scratch/vt20265/clouded_apollo/P21213_139_S159_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_139_S159_L003_R2_001.fastq.gz']
files['S160']=['/scratch/vt20265/clouded_apollo/P21213_140_S160_L003_R1_001.fastq.gz','/scratch/vt20265/clouded_apollo/P21213_140_S160_L003_R2_001.fastq.gz']


for i in samples_1_libraries:
	fastq_1 = files[i][0]
	fastq_2 = files[i][1]
	Sample_ID=str(i)
	out_dir='/scratch/vt20265/clouded_apollo/mapping/'
	genome='/scratch/vt20265/clouded_apollo/Genome/GCA_907164705.1_Parnassius_apollo_genomic.fasta'
	outFile='/home/vt20265/Scripts/'+Sample_ID+'_Appolo_BWA_GATK'+'.sh'
	run_SRA_BWA_GATK_lib_1(Sample_ID, genome, out_dir, fastq_1, fastq_2, outFile)



