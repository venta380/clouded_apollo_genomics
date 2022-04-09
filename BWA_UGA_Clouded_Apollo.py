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
	outFile.write('#PBS -l nodes=1:ppn=10'+'\n')
	outFile.write('#PBS -N mapping-'+Sample_ID+'\n')
	outFile.write('#PBS -l walltime=60:00:00'+'\n')
	outFile.write('#PBS -o mapping-'+Sample_ID+'.out'+'\n')
	outFile.write('#PBS -e mapping-'+Sample_ID+'.err'+'\n')
	outFile.write('#PBS -l mem=100gb'+'\n')
	outFile.write('\n')
	outFile.write("module load cutadapt"+'\n')
	outFile.write("module load pigz"+'\n')
	outFile.write("module load SRA-Toolkit"+'\n')
	outFile.write('wait'+'\n')
	outFile.write('cd $TMPDIR'+'\n')
	outFile.write('wait'+'\n')
	outFile.write('\n')
	outFile.write(""+'\n')
	outFile.write("ulimit -c unlimited"+'\n')
	outFile.write('\n')
	outFile.write("cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -j 10 -q 30,30 -o "+fastq_1.split('.')[0]+'_trimmed.fastq.gz -p '+fastq_2.split('.')[0]+'_trimmed.fastq.gz '+' /scratch/vt20265/bam_files/DNA_seq_temp/'+fastq_1+' /scratch/vt20265/bam_files/DNA_seq_temp/'+fastq_2+' '+'\n')
	outFile.write('\n')
	outFile.write("module load SAMtools"+'\n')
	outFile.write("module load BWA"+'\n')
	outFile.write("module load BCFtools"+'\n')
	outFile.write("module load VCFtools"+'\n')
	outFile.write('module load Java/1.8.0_92'+'\n')
	outFile.write('\n')
	outFile.write('bwa mem -t 10 -M '+genome+' -R '+' '+'"@RG\\tLB:Lib1\\tID:1\\tSM:'+Sample_ID+'\\tPL:ILLUMINA"'+' $TMPDIR/'+fastq_1.split('.')[0]+'_trimmed.fastq.gz '+' $TMPDIR/'+fastq_2.split('.')[0]+'_trimmed.fastq.gz '+' > '+'$TMPDIR/'+Sample_ID+'.sam'+'\n')
	outFile.write('samtools sort '+'$TMPDIR/'+Sample_ID+'.sam -o '+'$TMPDIR/'+Sample_ID+'.bam'+'\n')
	outFile.write('samtools index '+'$TMPDIR/'+Sample_ID+'.bam'+'\n')
	outFile.write('module load picard/2.21.6-Java-11'+'\n')
	outFile.write('java -Xmx32g -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT='+'$TMPDIR/'+Sample_ID+'.bam '+'OUTPUT='+'$TMPDIR/'+Sample_ID+'.bam.dedup.bam METRICS_FILE='+'$TMPDIR/'+Sample_ID+'.bam.metrics REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT'+'\n')
	outFile.write('samtools index '+'$TMPDIR/'+Sample_ID+'.bam.dedup.bam'+'\n')
	outFile.write('module load GATK/3.8-1-Java-1.8.0_144'+'\n')
	outFile.write('java -Xmx32g -jar $EBROOTGATK/GenomeAnalysisTK.jar -I ' + '$TMPDIR/'+Sample_ID+'.bam.dedup.bam -R '+genome+' -T RealignerTargetCreator -o '+'$TMPDIR/'+Sample_ID+'.intervals'+'\n')
	outFile.write('java -Xmx32g -jar $EBROOTGATK/GenomeAnalysisTK.jar -I ' + '$TMPDIR/'+Sample_ID+'.bam.dedup.bam -R '+genome+' -T IndelRealigner --filter_bases_not_stored -o '+out_dir+Sample_ID+'.realign.bam -targetIntervals '+'$TMPDIR/'+Sample_ID+'.intervals'+'\n')
	outFile.write("#java -Xmx30g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T BaseRecalibrator -I "+out_dir+Sample_ID+'.realign.bam -R '+genome+' -o '+out_dir+Sample_ID+'.bam.dedup.realign.calibration.csv -knownSites '+out_dir+'known.vcf.gz'+'\n')
	outFile.write("#java -Xmx30g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T PrintReads -I "+out_dir+Sample_ID+'.realign.bam -R '+genome+' -BQSR '+out_dir+Sample_ID+'.bam.dedup.realign.calibration.csv -o '+out_dir+Sample_ID+'.final.bam'+'\n')
	outFile.write("#samtools flagstat "+out_dir+Sample_ID+'.final.bam > '+out_dir+Sample_ID+'.flagstat'+'\n')
	outFile.write("#java -Xmx30g -jar  $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R "+genome+' -I '+out_dir+Sample_ID+'.final.bam --emitRefConfidence GVCF  -o '+out_dir+Sample_ID+'.final.bam.g.vcf'+'\n'+'\n')	
	outFile.write('\n')
	outFile.write('\n')





Knownsites="known.vcf.gz"
srafiles={"15_77":["15_77_S3_L001_R2_001.fastq.gz", "15_77_S3_L001_R1_001.fastq.gz"],"11_77":["11_77_S8_L001_R2_001.fastq.gz", "11_77_S8_L001_R1_001.fastq.gz"],"09_77":["09_77_S8_L001_R2_001.fastq.gz", "09_77_S8_L001_R1_001.fastq.gz"],"16_77":["16_77_S8_L001_R2_001.fastq.gz", "16_77_S8_L001_R1_001.fastq.gz"],"01_77":["01_77_S8_L001_R2_001.fastq.gz","01_77_S8_L001_R1_001.fastq.gz"]}
srafiles={"ST_103": ["ST_103_S9_L001_R1_001.fastq.gz", "ST_103_S9_L001_R2_001.fastq.gz"],"ST_37": ["ST_37_S9_L001_R1_001.fastq.gz", "ST_37_S9_L001_R2_001.fastq.gz"],"ST_34": ["ST_34_S4_L001_R1_001.fastq.gz", "ST_34_S4_L001_R2_001.fastq.gz"],"ST_54": ["ST_54_S4_L001_R1_001.fastq.gz", "ST_54_S4_L001_R2_001.fastq.gz"],"ST_29": ["ST_29_S4_L001_R1_001.fastq.gz", "ST_29_S4_L001_R2_001.fastq.gz"],"ST_113": ["ST_113_S4_L001_R1_001.fastq.gz", "ST_113_S4_L001_R2_001.fastq.gz"],"ST_111": ["ST_111_S4_L001_R2_001.fastq.gz", "ST_111_S4_L001_R1_001.fastq.gz"],"ST_38": ["ST_38_S4_L001_R1_001.fastq.gz", "ST_38_S4_L001_R2_001.fastq.gz"],"ST_105": ["ST_105_S8_L001_R2_001.fastq.gz", "ST_105_S8_L001_R1_001.fastq.gz"],"ST_51": ["ST_51_S8_L001_R1_001.fastq.gz", "ST_51_S8_L001_R2_001.fastq.gz"]}


srafiles={"AZTC03":["AZTC03_S1_L001_R1_001.fastq.gz", "AZTC03_S1_L001_R2_001.fastq.gz"],"AZTC05":["AZTC05_S1_L001_R1_001.fastq.gz", "AZTC05_S1_L001_R2_001.fastq.gz"],"AZTC07":["AZTC07_S1_L001_R1_001.fastq.gz", "AZTC07_S1_L001_R2_001.fastq.gz"],"AZTC10":["AZTC10_S1_L001_R1_001.fastq.gz", "AZTC10_S1_L001_R2_001.fastq.gz"],"AZTC11":["AZTC11_S1_L001_R1_001.fastq.gz", "AZTC11_S1_L001_R2_001.fastq.gz"],"MAD03":["MAD03_S1_L001_R1_001.fastq.gz", "MAD03_S1_L001_R2_001.fastq.gz"],"MAD08":["MAD08_S6_L001_R1_001.fastq.gz", "MAD08_S6_L001_R2_001.fastq.gz"],"MAD13":["MAD13_S6_L001_R1_001.fastq.gz", "MAD13_S6_L001_R2_001.fastq.gz"],"MAD18":["MAD18_S6_L001_R1_001.fastq.gz", "MAD18_S6_L001_R2_001.fastq.gz"],"MAD22":["MAD22_S6_L001_R1_001.fastq.gz", "MAD22_S6_L001_R2_001.fastq.gz"],"MAD25":["MAD25_S6_L001_R1_001.fastq.gz", "MAD25_S6_L001_R2_001.fastq.gz"],"MAD28":["MAD28_S6_L001_R1_001.fastq.gz", "MAD28_S6_L001_R2_001.fastq.gz"],"CNTN01":["CNTN01_S2_L001_R1_001.fastq.gz", "CNTN01_S2_L001_R2_001.fastq.gz"],"CNTN02":["CNTN02_S2_L001_R1_001.fastq.gz", "CNTN02_S2_L001_R2_001.fastq.gz"],"CNTN03":["CNTN03_S2_R1_001.fastq.gz", "CNTN03_S2_R2_001.fastq.gz"],"CNTN04":["CNTN04_S2_L001_R1_001.fastq.gz", "CNTN04_S2_L001_R2_001.fastq.gz"],"CNTN05":["CNTN05_S2_L001_R1_001.fastq.gz", "CNTN05_S2_L001_R2_001.fastq.gz"],"CNTN06":["CNTN06_S2_L001_R1_001.fastq.gz", "CNTN06_S2_L001_R2_001.fastq.gz"],"CNTN07":["CNTN07_S7_L001_R1_001.fastq.gz", "CNTN07_S7_L001_R2_001.fastq.gz"]}

srafiles={'S141':['P21213_121_S141_L003_R1_001.fastq.gz','P21213_121_S141_L003_R2_001.fastq.gz']}

for i in srafiles:
	fastq_1 = srafiles[i][0]
	fastq_2 = srafiles[i][1]
	Sample_ID=str(i)
	out_dir='/scratch/vt20265/clouded_apollo/mapping/'
	genome='/scratch/vt20265/clouded_apollo/Genome/pt_114_001.Clouded_apollo.ipa1.5.0.purged.primary.fasta'
	outFile='/home/vt20265/Scripts/'+Sample_ID+'_BWA_GATK'+'.sh'
	run_SRA_BWA_GATK_lib_1(Sample_ID, genome, out_dir, fastq_1, fastq_2, outFile)






