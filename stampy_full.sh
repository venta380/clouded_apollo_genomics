#!/bin/bash -l
#SBATCH -A snic2021-5-20
#SBATCH -p core -n 20
#SBATCH -J Test_mapping
#SBATCH -t 200:00:00

module load bioinfo-tools
module load samtools/1.2
module load bwa/0.7.12
module load bcftools
module load vcftools
module load Stampy
module load cutadapt

cd $TMPDIR

PE1='/proj/uppstore2017185/b2014034_nobackup/Veronika/Clouded_Apollo_Project/data/Raw_Data/INBOX/P21213/P21213_121/02-FASTQ/210831_A00689_0336_AHGNCJDSX2/P21213_121_S141_L003_R1_001.fastq.gz'
PE2='/proj/uppstore2017185/b2014034_nobackup/Veronika/Clouded_Apollo_Project/data/Raw_Data/INBOX/P21213/P21213_121/02-FASTQ/210831_A00689_0336_AHGNCJDSX2/P21213_121_S141_L003_R2_001.fastq.gz'
outdir='/proj/uppstore2017185/b2014034_nobackup/Venkat/Clouded_Apollo_Project/mapping/'
ref='/proj/uppstore2017185/b2014034_nobackup/Venkat/Clouded_Apollo_Project/Genome/GCA_907164705.1_Parnassius_apollo_genomic.fasta'

#cutadapt -j 10 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -q 30,30 -o P21213_121_S141_L003_R1_001.fastq_trimmed.fastq.gz -p P21213_121_S141_L003_R2_001.fastq_trimmed.fastq.gz  $PE1 $PE2 &

zcat $PE1 > reads.fastq
zcat $PE2 >> reads.fastq
#bwa index $ref
#
#stampy.py -G $ref $ref
#stampy.py -g $ref -H $ref
stampy.py -g $ref -h $ref -t 5 --fast --overwrite --readgroup=ID:group1,SM:P21213_121_S141,PL:illumina,LB:lib1,PU:unit1  --substitutionrate=0.025 -o $outdir'/P21213_121_S141_0.025.sam' -M reads.fastq &
stampy.py -g $ref -h $ref -t 5 --fast --overwrite --readgroup=ID:group1,SM:P21213_121_S141,PL:illumina,LB:lib1,PU:unit1  --substitutionrate=0.06 -o $outdir'/P21213_121_S141_0.06.sam' -M reads.fastq &
stampy.py -g $ref -h $ref -t 5 --fast --overwrite --readgroup=ID:group1,SM:P21213_121_S141,PL:illumina,LB:lib1,PU:unit1  --substitutionrate=0.08 -o $outdir'/P21213_121_S141_0.08.sam' -M reads.fastq &
stampy.py -g $ref -h $ref -t 5 --fast --overwrite --readgroup=ID:group1,SM:P21213_121_S141,PL:illumina,LB:lib1,PU:unit1  --substitutionrate=0.12 -o $outdir'/P21213_121_S141_0.12.sam' -M reads.fastq &

wait
