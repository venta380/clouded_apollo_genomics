#!/bin/sh
#PBS -S /bin/bash
#PBS -q batch
#PBS -l nodes=1:ppn=20
#PBS -N RepeatModeler
#PBS -l walltime=80:00:00
#PBS -o RepeatModeler.out
#PBS -e RepeatModeler.err
#PBS -l mem=100gb


module load RepeatModeler

fasta_CA='/scratch/vt20265/clouded_apollo/Genome/pt_114_001.Clouded_apollo.ipa1.5.0.purged.primary.fasta'
fasta_A='/scratch/vt20265/clouded_apollo/Genome/GCA_907164705.1_Parnassius_apollo_genomic.fasta'

#/apps/eb/RepeatModeler/2.0.1-foss-2019b/BuildDatabase -name clouded_apollo $fasta_CA
#/apps/eb/RepeatModeler/2.0.1-foss-2019b/BuildDatabase -name apollo $fasta_A

nohup /apps/eb/RepeatModeler/2.0.1-foss-2019b/RepeatModeler -database clouded_apollo -pa 10 -LTRStruct >& run.out &
nohup /apps/eb/RepeatModeler/2.0.1-foss-2019b/RepeatModeler -database apollo -pa 10 -LTRStruct >& run_A.out &
wait


