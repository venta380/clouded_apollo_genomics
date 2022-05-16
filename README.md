# clouded apollo genomics. 
This is a dedicated repository for the conservation genomics of Clouded apollo butterfly. 
![alt text](https://upload.wikimedia.org/wikipedia/commons/thumb/1/10/Parnassius_mnemosyne_MHNT_CUT_2013_3_5_Le_Mont_Dore_Male_dos.jpg/1920px-Parnassius_mnemosyne_MHNT_CUT_2013_3_5_Le_Mont_Dore_Male_dos.jpg "Logo Title Text 1")
## Maaping to the Reference genome. 
Apollo_stampy.py is a python script that writes sbatch scripts. These scripts were submitted to uppmax to get our final bam files. 

## SNP calling
SNP calling was performed on all bam files using BCFtools. 

bcftoolsCommand=mpileup --threads 10 --skip-indels -Ou -f /scratch/vt20265/clouded_apollo/Genome/GCA_907164705.1_Parnassius_apollo_genomic.fasta /scratch/vt20265/clouded_apollo/mapping/S121.bam /scratch/vt20265/clouded_apollo/mapping/S122.bam /scratch/vt20265/clouded_apollo/mapping/S123.bam /scratch/vt20265/clouded_apollo/mapping/S124.bam /scratch/vt20265/clouded_apollo/mapping/S125.bam /scratch/vt20265/clouded_apollo/mapping/S126.bam /scratch/vt20265/clouded_apollo/mapping/S127.bam /scratch/vt20265/clouded_apollo/mapping/S128.bam /scratch/vt20265/clouded_apollo/mapping/S129.bam /scratch/vt20265/clouded_apollo/mapping/S130.bam /scratch/vt20265/clouded_apollo/mapping/S131.bam /scratch/vt20265/clouded_apollo/mapping/S132.bam /scratch/vt20265/clouded_apollo/mapping/S133.bam /scratch/vt20265/clouded_apollo/mapping/S134.bam /scratch/vt20265/clouded_apollo/mapping/S135.bam /scratch/vt20265/clouded_apollo/mapping/S136.bam /scratch/vt20265/clouded_apollo/mapping/S137.bam /scratch/vt20265/clouded_apollo/mapping/S138.bam /scratch/vt20265/clouded_apollo/mapping/S139.bam /scratch/vt20265/clouded_apollo/mapping/S140.bam /scratch/vt20265/clouded_apollo/mapping/S141.bam /scratch/vt20265/clouded_apollo/mapping/S142.bam /scratch/vt20265/clouded_apollo/mapping/S143.bam /scratch/vt20265/clouded_apollo/mapping/S144.bam /scratch/vt20265/clouded_apollo/mapping/S145.bam /scratch/vt20265/clouded_apollo/mapping/S146.bam /scratch/vt20265/clouded_apollo/mapping/S147.bam /scratch/vt20265/clouded_apollo/mapping/S148.bam /scratch/vt20265/clouded_apollo/mapping/S149.bam /scratch/vt20265/clouded_apollo/mapping/S151.bam /scratch/vt20265/clouded_apollo/mapping/S152.bam /scratch/vt20265/clouded_apollo/mapping/S153.bam /scratch/vt20265/clouded_apollo/mapping/S154.bam /scratch/vt20265/clouded_apollo/mapping/S155.bam /scratch/vt20265/clouded_apollo/mapping/S156.bam /scratch/vt20265/clouded_apollo/mapping/S157.bam /scratch/vt20265/clouded_apollo/mapping/S158.bam /scratch/vt20265/clouded_apollo/mapping/S159.bam


## PCA


## LEA/Admixture

## IBD (Identity by descent) analysis and Plots

## FST


