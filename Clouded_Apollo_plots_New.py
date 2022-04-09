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
import re
import random 
import scipy
import threading 
import multiprocessing
from multiprocessing import  Pool



def join_raw_data_base(lists):
	for i in range(1,len(lists)):
		j=i-1
		if i == 1:
			new_2=pandas.merge(lists[j], lists[i], on=['CHROM','BIN_START'], how='outer')
		else:
			new_2=pandas.merge(nextone, lists[i], on=['CHROM','BIN_START'], how='outer')
		nextone=new_2
	return nextone


def join_raw_data_base_FST(lists):
	for i in range(1,len(lists)):
		j=i-1
		if i == 1:
			new_2=pandas.merge(lists[j], lists[i], on=['CHROM','BIN_START','BIN_END','N_sites'], how='outer')
		else:
			new_2=pandas.merge(nextone, lists[i], on=['CHROM','BIN_START','BIN_END','N_sites'], how='outer')
		nextone=new_2
	return nextone



output_4D_VA=pandas.read_csv('Vasternorrland_4fold_1_X__stats.csv', skiprows=1, names=['UI', 'CHROM','BIN_START', 'Pi_4D_Vasternorrland', 'Tw_4D_Vasternorrland', 'Td_4D_Vasternorrland', 'N_sites_4D_Vasternorrland'])
output_4D_BL=pandas.read_csv('Blekinge_4fold_1_X__stats.csv', skiprows=1, names=['UI', 'CHROM','BIN_START', 'Pi_4D_Blekinge', 'Tw_4D_Blekinge', 'Td_4D_Blekinge', 'N_sites_4D_Blekinge'])
output_4D_NA=pandas.read_csv('NordensArk_4fold_1_X__stats.csv', skiprows=1, names=['UI', 'CHROM','BIN_START', 'Pi_4D_NordensArk', 'Tw_4D_NordensArk', 'Td_4D_NordensArk', 'N_sites_4D_NordensArk'])
output_4D_UP=pandas.read_csv('Uppland_4fold_1_X__stats.csv', skiprows=1, names=['UI', 'CHROM','BIN_START', 'Pi_4D_Uppland', 'Tw_4D_Uppland', 'Td_4D_Uppland', 'N_sites_4D_Uppland'])

output_4D_Br=pandas.read_csv('Brudskär_4fold_1_X__stats.csv', skiprows=1, names=['UI', 'CHROM','BIN_START', 'Pi_4D_Brudskär', 'Tw_4D_Brudskär', 'Td_4D_Brudskär', 'N_sites_4D_Brudskär'])
output_4D_Lö=pandas.read_csv('Lötaholmen_4fold_1_X__stats.csv', skiprows=1, names=['UI', 'CHROM','BIN_START', 'Pi_4D_Lötaholmen', 'Tw_4D_Lötaholmen', 'Td_4D_Lötaholmen', 'N_sites_4D_Lötaholmen'])


output_3D_VA=pandas.read_csv('Vasternorrland_codon3_1_X__stats.csv', skiprows=1, names=['UI', 'CHROM','BIN_START', 'Pi_3D_Vasternorrland', 'Tw_3D_Vasternorrland', 'Td_3D_Vasternorrland', 'N_sites_3D_Vasternorrland'])
output_3D_BL=pandas.read_csv('Blekinge_codon3_1_X__stats.csv', skiprows=1, names=['UI', 'CHROM','BIN_START', 'Pi_3D_Blekinge', 'Tw_3D_Blekinge', 'Td_3D_Blekinge', 'N_sites_3D_Blekinge'])
output_3D_NA=pandas.read_csv('NordensArk_codon3_1_X__stats.csv', skiprows=1, names=['UI', 'CHROM','BIN_START', 'Pi_3D_NordensArk', 'Tw_3D_NordensArk', 'Td_3D_NordensArk', 'N_sites_3D_NordensArk'])
output_3D_UP=pandas.read_csv('Uppland_codon3_1_X__stats.csv', skiprows=1, names=['UI', 'CHROM','BIN_START', 'Pi_3D_Uppland', 'Tw_3D_Uppland', 'Td_3D_Uppland', 'N_sites_3D_Uppland'])

output_3D_Br=pandas.read_csv('Brudskär_codon3_1_X__stats.csv', skiprows=1, names=['UI', 'CHROM','BIN_START', 'Pi_3D_Brudskär', 'Tw_3D_Brudskär', 'Td_3D_Brudskär', 'N_sites_3D_Brudskär'])
output_3D_Lö=pandas.read_csv('Lötaholmen_codon3_1_X__stats.csv', skiprows=1, names=['UI', 'CHROM','BIN_START', 'Pi_3D_Lötaholmen', 'Tw_3D_Lötaholmen', 'Td_3D_Lötaholmen', 'N_sites_3D_Lötaholmen'])


output_2D_VA=pandas.read_csv('Vasternorrland_codon2_1_X__stats.csv', skiprows=1, names=['UI', 'CHROM','BIN_START', 'Pi_2D_Vasternorrland', 'Tw_2D_Vasternorrland', 'Td_2D_Vasternorrland', 'N_sites_2D_Vasternorrland'])
output_2D_BL=pandas.read_csv('Blekinge_codon2_1_X__stats.csv', skiprows=1, names=['UI', 'CHROM','BIN_START', 'Pi_2D_Blekinge', 'Tw_2D_Blekinge', 'Td_2D_Blekinge', 'N_sites_2D_Blekinge'])
output_2D_NA=pandas.read_csv('NordensArk_codon2_1_X__stats.csv', skiprows=1, names=['UI', 'CHROM','BIN_START', 'Pi_2D_NordensArk', 'Tw_2D_NordensArk', 'Td_2D_NordensArk', 'N_sites_2D_NordensArk'])
output_2D_UP=pandas.read_csv('Uppland_codon2_1_X__stats.csv', skiprows=1, names=['UI', 'CHROM','BIN_START', 'Pi_2D_Uppland', 'Tw_2D_Uppland', 'Td_2D_Uppland', 'N_sites_2D_Uppland'])


output_2D_Br=pandas.read_csv('Brudskär_codon2_1_X__stats.csv', skiprows=1, names=['UI', 'CHROM','BIN_START', 'Pi_2D_Brudskär', 'Tw_2D_Brudskär', 'Td_2D_Brudskär', 'N_sites_2D_Brudskär'])
output_2D_Lö=pandas.read_csv('Lötaholmen_codon2_1_X__stats.csv', skiprows=1, names=['UI', 'CHROM','BIN_START', 'Pi_2D_Lötaholmen', 'Tw_2D_Lötaholmen', 'Td_2D_Lötaholmen', 'N_sites_2D_Lötaholmen'])


output_1D_VA=pandas.read_csv('Vasternorrland_codon1_1_X__stats.csv', skiprows=1, names=['UI', 'CHROM','BIN_START', 'Pi_1D_Vasternorrland', 'Tw_1D_Vasternorrland', 'Td_1D_Vasternorrland', 'N_sites_1D_Vasternorrland'])
output_1D_BL=pandas.read_csv('Blekinge_codon1_1_X__stats.csv', skiprows=1, names=['UI', 'CHROM','BIN_START', 'Pi_1D_Blekinge', 'Tw_1D_Blekinge', 'Td_1D_Blekinge', 'N_sites_1D_Blekinge'])
output_1D_NA=pandas.read_csv('NordensArk_codon1_1_X__stats.csv', skiprows=1, names=['UI', 'CHROM','BIN_START', 'Pi_1D_NordensArk', 'Tw_1D_NordensArk', 'Td_1D_NordensArk', 'N_sites_1D_NordensArk'])
output_1D_UP=pandas.read_csv('Uppland_codon1_1_X__stats.csv', skiprows=1, names=['UI', 'CHROM','BIN_START', 'Pi_1D_Uppland', 'Tw_1D_Uppland', 'Td_1D_Uppland', 'N_sites_1D_Uppland'])


output_1D_Br=pandas.read_csv('Brudskär_codon1_1_X__stats.csv', skiprows=1, names=['UI', 'CHROM','BIN_START', 'Pi_1D_Brudskär', 'Tw_1D_Brudskär', 'Td_1D_Brudskär', 'N_sites_1D_Brudskär'])
output_1D_Lö=pandas.read_csv('Lötaholmen_codon1_1_X__stats.csv', skiprows=1, names=['UI', 'CHROM','BIN_START', 'Pi_1D_Lötaholmen', 'Tw_1D_Lötaholmen', 'Td_1D_Lötaholmen', 'N_sites_1D_Lötaholmen'])


Blekinge_NordensArk_Fst=pandas.read_csv('Blekinge_NordensArk_dxyFst_.csv')
Blekinge_Uppland_Fst=pandas.read_csv('Blekinge_Uppland_dxyFst_.csv')
Uppland_NordensArk_Fst=pandas.read_csv('Uppland_NordensArk_dxyFst_.csv')
Vasternorrland_Blekinge_Fst=pandas.read_csv('Vasternorrland_Blekinge_dxyFst_.csv')
Vasternorrland_NordensArk_Fst=pandas.read_csv('Vasternorrland_NordensArk_dxyFst_.csv')
Vasternorrland_Uppland_Fst=pandas.read_csv('Vasternorrland_Uppland_dxyFst_.csv')
Brudskär_4D_Lötaholmen_4D_Fst=pandas.read_csv('Brudskär_4D_Lötaholmen_4D_dxyFst_.csv')


Blekinge_NordensArk_Fst['Blekinge_NordensArk_Temp_Nsites']=Blekinge_NordensArk_Fst['Blekinge_NordensArk_fixed']+Blekinge_NordensArk_Fst['Blekinge_NordensArk_private_a']+Blekinge_NordensArk_Fst['Blekinge_NordensArk_private_b']+Blekinge_NordensArk_Fst['Blekinge_NordensArk_shared']
Blekinge_Uppland_Fst['Blekinge_Uppland_Temp_Nsites']=Blekinge_Uppland_Fst['Blekinge_Uppland_fixed']+Blekinge_Uppland_Fst['Blekinge_Uppland_private_a']+Blekinge_Uppland_Fst['Blekinge_Uppland_private_b']+Blekinge_Uppland_Fst['Blekinge_Uppland_shared']
Uppland_NordensArk_Fst['Uppland_NordensArk_Temp_Nsites']=Uppland_NordensArk_Fst['Uppland_NordensArk_fixed']+Uppland_NordensArk_Fst['Uppland_NordensArk_private_a']+Uppland_NordensArk_Fst['Uppland_NordensArk_private_b']+Uppland_NordensArk_Fst['Uppland_NordensArk_shared']
Vasternorrland_Blekinge_Fst['Vasternorrland_Blekinge_Temp_Nsites']=Vasternorrland_Blekinge_Fst['Vasternorrland_Blekinge_fixed']+Vasternorrland_Blekinge_Fst['Vasternorrland_Blekinge_private_a']+Vasternorrland_Blekinge_Fst['Vasternorrland_Blekinge_private_b']+Vasternorrland_Blekinge_Fst['Vasternorrland_Blekinge_shared']
Vasternorrland_NordensArk_Fst['Vasternorrland_NordensArk_Temp_Nsites']=Vasternorrland_NordensArk_Fst['Vasternorrland_NordensArk_fixed']+Vasternorrland_NordensArk_Fst['Vasternorrland_NordensArk_private_a']+Vasternorrland_NordensArk_Fst['Vasternorrland_NordensArk_private_b']+Vasternorrland_NordensArk_Fst['Vasternorrland_NordensArk_shared']
Vasternorrland_Uppland_Fst['Vasternorrland_Uppland_Temp_Nsites']=Vasternorrland_Uppland_Fst['Vasternorrland_Uppland_fixed']+Vasternorrland_Uppland_Fst['Vasternorrland_Uppland_private_a']+Vasternorrland_Uppland_Fst['Vasternorrland_Uppland_private_b']+Vasternorrland_Uppland_Fst['Vasternorrland_Uppland_shared']
Brudskär_4D_Lötaholmen_4D_Fst['Brudskär_4D_Lötaholmen_4D_Temp_Nsites']=Brudskär_4D_Lötaholmen_4D_Fst['Brudskär_4D_Lötaholmen_4D_fixed']+Brudskär_4D_Lötaholmen_4D_Fst['Brudskär_4D_Lötaholmen_4D_private_a']+Brudskär_4D_Lötaholmen_4D_Fst['Brudskär_4D_Lötaholmen_4D_private_b']+Brudskär_4D_Lötaholmen_4D_Fst['Brudskär_4D_Lötaholmen_4D_shared']



Blekinge_NordensArk_Fst=Blekinge_NordensArk_Fst[Blekinge_NordensArk_Fst['Blekinge_NordensArk_Temp_Nsites']>5]
Blekinge_Uppland_Fst=Blekinge_Uppland_Fst[Blekinge_Uppland_Fst['Blekinge_Uppland_Temp_Nsites']>5]
Uppland_NordensArk_Fst=Uppland_NordensArk_Fst[Uppland_NordensArk_Fst['Uppland_NordensArk_Temp_Nsites']>5]
Vasternorrland_Blekinge_Fst=Vasternorrland_Blekinge_Fst[Vasternorrland_Blekinge_Fst['Vasternorrland_Blekinge_Temp_Nsites']>5]
Vasternorrland_NordensArk_Fst=Vasternorrland_NordensArk_Fst[Vasternorrland_NordensArk_Fst['Vasternorrland_NordensArk_Temp_Nsites']>5]
Vasternorrland_Uppland_Fst=Vasternorrland_Uppland_Fst[Vasternorrland_Uppland_Fst['Vasternorrland_Uppland_Temp_Nsites']>5]
Brudskär_4D_Lötaholmen_4D_Fst=Brudskär_4D_Lötaholmen_4D_Fst[Brudskär_4D_Lötaholmen_4D_Fst['Brudskär_4D_Lötaholmen_4D_Temp_Nsites']>5]




output_1=[Blekinge_NordensArk_Fst,Blekinge_Uppland_Fst,Uppland_NordensArk_Fst,Vasternorrland_Blekinge_Fst,Vasternorrland_NordensArk_Fst,Vasternorrland_Uppland_Fst,Brudskär_4D_Lötaholmen_4D_Fst]
list_2_Fst=join_raw_data_base_FST(output_1)


output_1=[output_4D_VA,output_4D_BL,output_4D_NA,output_4D_UP,output_3D_VA,output_3D_BL,output_3D_NA,output_3D_UP,output_2D_VA,output_2D_BL,output_2D_NA,output_2D_UP,output_1D_VA,output_1D_BL,output_1D_NA,output_1D_UP,output_4D_Br, output_4D_Lö, output_3D_Br, output_3D_Lö, output_2D_Br, output_2D_Lö, output_1D_Br, output_1D_Lö]
list_2=join_raw_data_base(output_1)


##0 SNPS
#list_2[(list_2.Tw_2D_Vasternorrland==0.000000)&(list_2.Tw_1D_Vasternorrland==0.000000)&(list_2.Tw_3D_Vasternorrland==0.000000)&(list_2.Tw_4D_Vasternorrland==0.000000)]
#list_2[(list_2.Tw_2D_Vasternorrland==0.000000)&(list_2.Tw_1D_Vasternorrland==0.000000)&(list_2.Tw_3D_Vasternorrland==0.000000)&(list_2.Tw_4D_Vasternorrland==0.000000)&(list_2.Chr_number=='Chr_Z')]
##259 all 22 Z
#list_2[(list_2.Tw_2D_Blekinge==0.000000)&(list_2.Tw_1D_Blekinge==0.000000)&(list_2.Tw_3D_Blekinge==0.000000)&(list_2.Tw_4D_Blekinge==0.000000)]
#list_2[(list_2.Tw_2D_Blekinge==0.000000)&(list_2.Tw_1D_Blekinge==0.000000)&(list_2.Tw_3D_Blekinge==0.000000)&(list_2.Tw_4D_Blekinge==0.000000)&(list_2.Chr_number=='Chr_Z')]
##321 all 30 Z
#list_2[(list_2.Tw_2D_NordensArk==0.000000)&(list_2.Tw_1D_NordensArk==0.000000)&(list_2.Tw_3D_NordensArk==0.000000)&(list_2.Tw_4D_NordensArk==0.000000)]
#list_2[(list_2.Tw_2D_NordensArk==0.000000)&(list_2.Tw_1D_NordensArk==0.000000)&(list_2.Tw_3D_NordensArk==0.000000)&(list_2.Tw_4D_NordensArk==0.000000)&(list_2.Chr_number=='Chr_Z')]
##257 all 19 Z
#list_2[(list_2.Tw_2D_Uppland==0.000000)&(list_2.Tw_1D_Uppland==0.000000)&(list_2.Tw_3D_Uppland==0.000000)&(list_2.Tw_4D_Uppland==0.000000)]
#list_2[(list_2.Tw_2D_Uppland==0.000000)&(list_2.Tw_1D_Uppland==0.000000)&(list_2.Tw_3D_Uppland==0.000000)&(list_2.Tw_4D_Uppland==0.000000)&(list_2.Chr_number=='Chr_Z')]
##276 all 25 Z






Chrom=pandas.read_csv("/scratch/vt20265/clouded_apollo/Genome/Apollo_chr.txt", skiprows=1, names=['CHROM', 'Chr_number', 'Prop'])

list_2=list_2.merge(Chrom[['CHROM','Chr_number']], on='CHROM')
list_2=list_2[list_2.N_sites_4D_Blekinge>100]


list_2_Fst=list_2_Fst.merge(Chrom[['CHROM','Chr_number']], on='CHROM')
list_2_Fst=list_2_Fst[list_2_Fst.N_sites>100]


CHROM=["Chr_2","Chr_3","Chr_4","Chr_5","Chr_6","Chr_7","Chr_8","Chr_9","Chr_10","Chr_11","Chr_12","Chr_13","Chr_14","Chr_15","Chr_16","Chr_17","Chr_18","Chr_19","Chr_20","Chr_21","Chr_22","Chr_23","Chr_24","Chr_25","Chr_26","Chr_27","Chr_28","Chr_29","Chr_30","Chr_Z"]
CHROM_dict={}

for i in CHROM:
	CH=list(set(list_2[list_2['Chr_number']==i]['CHROM']))
	CH.sort()
	CHROM_dict[i]=CH




Z_chr=list_2[list_2.Chr_number=='Chr_Z']
Auto_chr=list_2[list_2.Chr_number!='Chr_Z']
Auto_chr=Auto_chr.dropna(subset=['Tw_4D_Vasternorrland','Tw_4D_Blekinge','Tw_4D_NordensArk','Tw_4D_Uppland','Pi_4D_Vasternorrland','Pi_4D_Blekinge','Pi_4D_NordensArk','Pi_4D_Uppland','Tw_3D_Vasternorrland','Tw_3D_Blekinge','Tw_3D_NordensArk','Tw_3D_Uppland','Pi_3D_Vasternorrland','Pi_3D_Blekinge','Pi_3D_NordensArk','Pi_3D_Uppland','Tw_2D_Vasternorrland','Tw_2D_Blekinge','Tw_2D_NordensArk','Tw_2D_Uppland','Pi_2D_Vasternorrland','Pi_2D_Blekinge','Pi_2D_NordensArk','Pi_2D_Uppland','Tw_1D_Vasternorrland','Tw_1D_Blekinge','Tw_1D_NordensArk','Tw_1D_Uppland','Pi_1D_Vasternorrland','Pi_1D_Blekinge','Pi_1D_NordensArk','Pi_1D_Uppland'])
Z_chr=Z_chr.dropna(subset=['Tw_4D_Vasternorrland','Tw_4D_Blekinge','Tw_4D_NordensArk','Tw_4D_Uppland','Pi_4D_Vasternorrland','Pi_4D_Blekinge','Pi_4D_NordensArk','Pi_4D_Uppland','Tw_3D_Vasternorrland','Tw_3D_Blekinge','Tw_3D_NordensArk','Tw_3D_Uppland','Pi_3D_Vasternorrland','Pi_3D_Blekinge','Pi_3D_NordensArk','Pi_3D_Uppland','Tw_2D_Vasternorrland','Tw_2D_Blekinge','Tw_2D_NordensArk','Tw_2D_Uppland','Pi_2D_Vasternorrland','Pi_2D_Blekinge','Pi_2D_NordensArk','Pi_2D_Uppland','Tw_1D_Vasternorrland','Tw_1D_Blekinge','Tw_1D_NordensArk','Tw_1D_Uppland','Pi_1D_Vasternorrland','Pi_1D_Blekinge','Pi_1D_NordensArk','Pi_1D_Uppland'])


list_2.loc[(list_2.Chr_number!='Chr_Z'), 'CHROM_Type']='Auto'
list_2.loc[(list_2.Chr_number=='Chr_Z'), 'CHROM_Type']='Chr_Z'


list_2_Fst.loc[(list_2_Fst.Chr_number!='Chr_Z'), 'CHROM_Type']='Auto'
list_2_Fst.loc[(list_2_Fst.Chr_number=='Chr_Z'), 'CHROM_Type']='Chr_Z'


temp={}
temp['Pi']=list(Auto_chr['Pi_1D_Vasternorrland'])+list(Auto_chr['Pi_1D_Blekinge'])+list(Auto_chr['Pi_1D_NordensArk'])+list(Auto_chr['Pi_1D_Uppland'])
temp['Tw']=list(Auto_chr['Tw_1D_Vasternorrland'])+list(Auto_chr['Tw_1D_Blekinge'])+list(Auto_chr['Tw_1D_NordensArk'])+list(Auto_chr['Tw_1D_Uppland'])
temp['Td']=list(Auto_chr['Td_1D_Vasternorrland'])+list(Auto_chr['Td_1D_Blekinge'])+list(Auto_chr['Td_1D_NordensArk'])+list(Auto_chr['Td_1D_Uppland'])
temp['population']=['Vasternorrland']*len(Auto_chr['Td_1D_Vasternorrland'])+['Blekinge']*len(Auto_chr['Td_1D_Blekinge'])+['NordensArk']*len(Auto_chr['Td_1D_NordensArk'])+['Uppland']*len(Auto_chr['Td_1D_Uppland'])
temp['CHR']=['Auto']*len(Auto_chr['Td_1D_Vasternorrland'])+['Auto']*len(Auto_chr['Td_1D_Blekinge'])+['Auto']*len(Auto_chr['Td_1D_NordensArk'])+['Auto']*len(Auto_chr['Td_1D_Uppland'])
temp['Sites']=['Codon1']*len(Auto_chr['Td_1D_Vasternorrland'])+['Codon1']*len(Auto_chr['Td_1D_Blekinge'])+['Codon1']*len(Auto_chr['Td_1D_NordensArk'])+['Codon1']*len(Auto_chr['Td_1D_Uppland'])
temp_df_1D=pandas.DataFrame(temp)

temp_Z={}
temp_Z['Pi']=list(Z_chr['Pi_1D_Vasternorrland'])+list(Z_chr['Pi_1D_Blekinge'])+list(Z_chr['Pi_1D_NordensArk'])+list(Z_chr['Pi_1D_Uppland'])
temp_Z['Tw']=list(Z_chr['Tw_1D_Vasternorrland'])+list(Z_chr['Tw_1D_Blekinge'])+list(Z_chr['Tw_1D_NordensArk'])+list(Z_chr['Tw_1D_Uppland'])
temp_Z['Td']=list(Z_chr['Td_1D_Vasternorrland'])+list(Z_chr['Td_1D_Blekinge'])+list(Z_chr['Td_1D_NordensArk'])+list(Z_chr['Td_1D_Uppland'])
temp_Z['population']=['Vasternorrland']*len(Z_chr['Td_1D_Vasternorrland'])+['Blekinge']*len(Z_chr['Td_1D_Blekinge'])+['NordensArk']*len(Z_chr['Td_1D_NordensArk'])+['Uppland']*len(Z_chr['Td_1D_Uppland'])
temp_Z['CHR']=['Z']*len(Z_chr['Td_1D_Vasternorrland'])+['Z']*len(Z_chr['Td_1D_Blekinge'])+['Z']*len(Z_chr['Td_1D_NordensArk'])+['Z']*len(Z_chr['Td_1D_Uppland'])
temp_Z['Sites']=['Codon1']*len(Z_chr['Td_1D_Vasternorrland'])+['Codon1']*len(Z_chr['Td_1D_Blekinge'])+['Codon1']*len(Z_chr['Td_1D_NordensArk'])+['Codon1']*len(Z_chr['Td_1D_Uppland'])
temp_df_Z_1D=pandas.DataFrame(temp_Z)

temp={}
temp['Pi']=list(Auto_chr['Pi_2D_Vasternorrland'])+list(Auto_chr['Pi_2D_Blekinge'])+list(Auto_chr['Pi_2D_NordensArk'])+list(Auto_chr['Pi_2D_Uppland'])
temp['Tw']=list(Auto_chr['Tw_2D_Vasternorrland'])+list(Auto_chr['Tw_2D_Blekinge'])+list(Auto_chr['Tw_2D_NordensArk'])+list(Auto_chr['Tw_2D_Uppland'])
temp['Td']=list(Auto_chr['Td_2D_Vasternorrland'])+list(Auto_chr['Td_2D_Blekinge'])+list(Auto_chr['Td_2D_NordensArk'])+list(Auto_chr['Td_2D_Uppland'])
temp['population']=['Vasternorrland']*len(Auto_chr['Td_2D_Vasternorrland'])+['Blekinge']*len(Auto_chr['Td_2D_Blekinge'])+['NordensArk']*len(Auto_chr['Td_2D_NordensArk'])+['Uppland']*len(Auto_chr['Td_2D_Uppland'])
temp['CHR']=['Auto']*len(Auto_chr['Td_2D_Vasternorrland'])+['Auto']*len(Auto_chr['Td_2D_Blekinge'])+['Auto']*len(Auto_chr['Td_2D_NordensArk'])+['Auto']*len(Auto_chr['Td_2D_Uppland'])
temp['Sites']=['Codon2']*len(Auto_chr['Td_2D_Vasternorrland'])+['Codon2']*len(Auto_chr['Td_2D_Blekinge'])+['Codon2']*len(Auto_chr['Td_2D_NordensArk'])+['Codon2']*len(Auto_chr['Td_2D_Uppland'])
temp_df_2D=pandas.DataFrame(temp)

temp_Z={}
temp_Z['Pi']=list(Z_chr['Pi_2D_Vasternorrland'])+list(Z_chr['Pi_2D_Blekinge'])+list(Z_chr['Pi_2D_NordensArk'])+list(Z_chr['Pi_2D_Uppland'])
temp_Z['Tw']=list(Z_chr['Tw_2D_Vasternorrland'])+list(Z_chr['Tw_2D_Blekinge'])+list(Z_chr['Tw_2D_NordensArk'])+list(Z_chr['Tw_2D_Uppland'])
temp_Z['Td']=list(Z_chr['Td_2D_Vasternorrland'])+list(Z_chr['Td_2D_Blekinge'])+list(Z_chr['Td_2D_NordensArk'])+list(Z_chr['Td_2D_Uppland'])
temp_Z['population']=['Vasternorrland']*len(Z_chr['Td_2D_Vasternorrland'])+['Blekinge']*len(Z_chr['Td_2D_Blekinge'])+['NordensArk']*len(Z_chr['Td_2D_NordensArk'])+['Uppland']*len(Z_chr['Td_2D_Uppland'])
temp_Z['CHR']=['Z']*len(Z_chr['Td_2D_Vasternorrland'])+['Z']*len(Z_chr['Td_2D_Blekinge'])+['Z']*len(Z_chr['Td_2D_NordensArk'])+['Z']*len(Z_chr['Td_2D_Uppland'])
temp_Z['Sites']=['Codon2']*len(Z_chr['Td_2D_Vasternorrland'])+['Codon2']*len(Z_chr['Td_2D_Blekinge'])+['Codon2']*len(Z_chr['Td_2D_NordensArk'])+['Codon2']*len(Z_chr['Td_2D_Uppland'])
temp_df_Z_2D=pandas.DataFrame(temp_Z)



temp={}
temp['Pi']=list(Auto_chr['Pi_3D_Vasternorrland'])+list(Auto_chr['Pi_3D_Blekinge'])+list(Auto_chr['Pi_3D_NordensArk'])+list(Auto_chr['Pi_3D_Uppland'])
temp['Tw']=list(Auto_chr['Tw_3D_Vasternorrland'])+list(Auto_chr['Tw_3D_Blekinge'])+list(Auto_chr['Tw_3D_NordensArk'])+list(Auto_chr['Tw_3D_Uppland'])
temp['Td']=list(Auto_chr['Td_3D_Vasternorrland'])+list(Auto_chr['Td_3D_Blekinge'])+list(Auto_chr['Td_3D_NordensArk'])+list(Auto_chr['Td_3D_Uppland'])
temp['population']=['Vasternorrland']*len(Auto_chr['Td_3D_Vasternorrland'])+['Blekinge']*len(Auto_chr['Td_3D_Blekinge'])+['NordensArk']*len(Auto_chr['Td_3D_NordensArk'])+['Uppland']*len(Auto_chr['Td_3D_Uppland'])
temp['CHR']=['Auto']*len(Auto_chr['Td_3D_Vasternorrland'])+['Auto']*len(Auto_chr['Td_3D_Blekinge'])+['Auto']*len(Auto_chr['Td_3D_NordensArk'])+['Auto']*len(Auto_chr['Td_3D_Uppland'])
temp['Sites']=['Codon3']*len(Auto_chr['Td_3D_Vasternorrland'])+['Codon3']*len(Auto_chr['Td_3D_Blekinge'])+['Codon3']*len(Auto_chr['Td_3D_NordensArk'])+['Codon3']*len(Auto_chr['Td_3D_Uppland'])
temp_df_3D=pandas.DataFrame(temp)



temp_Z={}
temp_Z['Pi']=list(Z_chr['Pi_3D_Vasternorrland'])+list(Z_chr['Pi_3D_Blekinge'])+list(Z_chr['Pi_3D_NordensArk'])+list(Z_chr['Pi_3D_Uppland'])
temp_Z['Tw']=list(Z_chr['Tw_3D_Vasternorrland'])+list(Z_chr['Tw_3D_Blekinge'])+list(Z_chr['Tw_3D_NordensArk'])+list(Z_chr['Tw_3D_Uppland'])
temp_Z['Td']=list(Z_chr['Td_3D_Vasternorrland'])+list(Z_chr['Td_3D_Blekinge'])+list(Z_chr['Td_3D_NordensArk'])+list(Z_chr['Td_3D_Uppland'])
temp_Z['population']=['Vasternorrland']*len(Z_chr['Td_3D_Vasternorrland'])+['Blekinge']*len(Z_chr['Td_3D_Blekinge'])+['NordensArk']*len(Z_chr['Td_3D_NordensArk'])+['Uppland']*len(Z_chr['Td_3D_Uppland'])
temp_Z['CHR']=['Z']*len(Z_chr['Td_3D_Vasternorrland'])+['Z']*len(Z_chr['Td_3D_Blekinge'])+['Z']*len(Z_chr['Td_3D_NordensArk'])+['Z']*len(Z_chr['Td_3D_Uppland'])
temp_Z['Sites']=['Codon3']*len(Z_chr['Td_3D_Vasternorrland'])+['Codon3']*len(Z_chr['Td_3D_Blekinge'])+['Codon3']*len(Z_chr['Td_3D_NordensArk'])+['Codon3']*len(Z_chr['Td_3D_Uppland'])
temp_df_Z_3D=pandas.DataFrame(temp_Z)



temp={}
temp['Pi']=list(Auto_chr['Pi_4D_Vasternorrland'])+list(Auto_chr['Pi_4D_Blekinge'])+list(Auto_chr['Pi_4D_NordensArk'])+list(Auto_chr['Pi_4D_Uppland'])
temp['Tw']=list(Auto_chr['Tw_4D_Vasternorrland'])+list(Auto_chr['Tw_4D_Blekinge'])+list(Auto_chr['Tw_4D_NordensArk'])+list(Auto_chr['Tw_4D_Uppland'])
temp['Td']=list(Auto_chr['Td_4D_Vasternorrland'])+list(Auto_chr['Td_4D_Blekinge'])+list(Auto_chr['Td_4D_NordensArk'])+list(Auto_chr['Td_4D_Uppland'])
temp['population']=['Vasternorrland']*len(Auto_chr['Td_4D_Vasternorrland'])+['Blekinge']*len(Auto_chr['Td_4D_Blekinge'])+['NordensArk']*len(Auto_chr['Td_4D_NordensArk'])+['Uppland']*len(Auto_chr['Td_4D_Uppland'])
temp['CHR']=['Auto']*len(Auto_chr['Td_4D_Vasternorrland'])+['Auto']*len(Auto_chr['Td_4D_Blekinge'])+['Auto']*len(Auto_chr['Td_4D_NordensArk'])+['Auto']*len(Auto_chr['Td_4D_Uppland'])
temp['Sites']=['4fold']*len(Auto_chr['Td_4D_Vasternorrland'])+['4fold']*len(Auto_chr['Td_4D_Blekinge'])+['4fold']*len(Auto_chr['Td_4D_NordensArk'])+['4fold']*len(Auto_chr['Td_4D_Uppland'])
temp_df_4D=pandas.DataFrame(temp)



temp_Z={}
temp_Z['Pi']=list(Z_chr['Pi_4D_Vasternorrland'])+list(Z_chr['Pi_4D_Blekinge'])+list(Z_chr['Pi_4D_NordensArk'])+list(Z_chr['Pi_4D_Uppland'])
temp_Z['Tw']=list(Z_chr['Tw_4D_Vasternorrland'])+list(Z_chr['Tw_4D_Blekinge'])+list(Z_chr['Tw_4D_NordensArk'])+list(Z_chr['Tw_4D_Uppland'])
temp_Z['Td']=list(Z_chr['Td_4D_Vasternorrland'])+list(Z_chr['Td_4D_Blekinge'])+list(Z_chr['Td_4D_NordensArk'])+list(Z_chr['Td_4D_Uppland'])
temp_Z['population']=['Vasternorrland']*len(Z_chr['Td_4D_Vasternorrland'])+['Blekinge']*len(Z_chr['Td_4D_Blekinge'])+['NordensArk']*len(Z_chr['Td_4D_NordensArk'])+['Uppland']*len(Z_chr['Td_4D_Uppland'])
temp_Z['CHR']=['Z']*len(Z_chr['Td_4D_Vasternorrland'])+['Z']*len(Z_chr['Td_4D_Blekinge'])+['Z']*len(Z_chr['Td_4D_NordensArk'])+['Z']*len(Z_chr['Td_4D_Uppland'])
temp_Z['Sites']=['4fold']*len(Z_chr['Td_4D_Vasternorrland'])+['4fold']*len(Z_chr['Td_4D_Blekinge'])+['4fold']*len(Z_chr['Td_4D_NordensArk'])+['4fold']*len(Z_chr['Td_4D_Uppland'])
temp_df_Z_4D=pandas.DataFrame(temp_Z)



temp_df_ALL=pandas.concat([temp_df_1D, temp_df_Z_1D, temp_df_2D, temp_df_Z_2D, temp_df_3D, temp_df_Z_3D, temp_df_4D, temp_df_Z_4D], ignore_index=True)


ax=sns.catplot(x="Sites", y="Pi",  hue="population", showmeans=True,meanprops={"marker": "+","markeredgecolor": "white","markersize": "10"}, data=temp_df_ALL[temp_df_ALL.CHR=='Auto'], showfliers = False,  kind="box", palette=sns.color_palette(['#D35400','#F4D03F','#2C3E50','#27AE60']), hue_order=['Blekinge','NordensArk','Uppland','Vasternorrland'])
#ax.fig.suptitle('')
ax.set(ylabel=r'$\pi$')
plt.savefig('pi_Auto.png', dpi=1200)
plt.close()

ax=sns.catplot(x="Sites", y="Tw",  hue="population", showmeans=True,meanprops={"marker": "+","markeredgecolor": "white","markersize": "10"}, data=temp_df_ALL[temp_df_ALL.CHR=='Auto'], showfliers = False,  kind="box", palette=sns.color_palette(['#D35400','#F4D03F','#2C3E50','#27AE60']), hue_order=['Blekinge','NordensArk','Uppland','Vasternorrland'])
#ax.fig.suptitle('')
ax.set(ylabel=r'$\Theta$')
plt.savefig('Tw_Auto.png', dpi=1200)
plt.close()

ax=sns.catplot(x="Sites", y="Td",  hue="population", showmeans=True,meanprops={"marker": "+","markeredgecolor": "white","markersize": "10"}, data=temp_df_ALL[temp_df_ALL.CHR=='Auto'], showfliers = False,  kind="box", palette=sns.color_palette(['#D35400','#F4D03F','#2C3E50','#27AE60']), hue_order=['Blekinge','NordensArk','Uppland','Vasternorrland'])
#ax.fig.suptitle('')
ax.set(ylabel='TajimaD')
plt.savefig('Td_Auto.png', dpi=1200)
plt.close()



from pandas.plotting import scatter_matrix
axes = scatter_matrix(Auto_chr[['Tw_4D_Vasternorrland','Tw_4D_Blekinge','Tw_4D_NordensArk','Tw_4D_Uppland']], alpha=0.5, diagonal='kde')

corr = list_2[['Tw_4D_Vasternorrland','Tw_4D_Blekinge','Tw_4D_NordensArk','Tw_4D_Uppland']].corr()
corr = np.asmatrix(corr)
for i, j in zip(*plt.np.triu_indices_from(axes, k=1)):
    axes[i, j].annotate("%.3f" %corr[i,j], (0.8, 0.8), xycoords='axes fraction', ha='center', va='center')


for ax in axes.ravel():
    ax.set_xlabel(ax.get_xlabel(), fontsize = 10, rotation = 90)
    ax.set_ylabel(ax.get_ylabel(), fontsize = 10, rotation = 0)

plt.savefig('scatter_matrix_Tw.png', dpi=1200)
plt.close()




#################


print (str(round(Auto_chr['Pi_4D_Vasternorrland'].median(), 4))+' ± '+str(round(Auto_chr['Pi_4D_Vasternorrland'].std(), 4))+' Vasternorrland Pi 4D Auto'+'  '+str(round(Z_chr['Pi_4D_Vasternorrland'].median(), 4))+' ± '+str(round(Z_chr['Pi_4D_Vasternorrland'].std(), 4))+' Vasternorrland Pi 4D Z-Chr') 
print (str(round(Auto_chr['Tw_4D_Vasternorrland'].median(), 4))+' ± '+str(round(Auto_chr['Tw_4D_Vasternorrland'].std(), 4))+' Vasternorrland Tw 4D Auto'+'  '+str(round(Z_chr['Tw_4D_Vasternorrland'].median(), 4))+' ± '+str(round(Z_chr['Tw_4D_Vasternorrland'].std(), 4))+' Vasternorrland Tw 4D Z-Chr') 
print (str(round(Auto_chr['Td_4D_Vasternorrland'].median(), 4))+' ± '+str(round(Auto_chr['Td_4D_Vasternorrland'].std(), 4))+' Vasternorrland Td 4D Auto'+'  '+str(round(Z_chr['Td_4D_Vasternorrland'].median(), 4))+' ± '+str(round(Z_chr['Td_4D_Vasternorrland'].std(), 4))+' Vasternorrland Td 4D Z-Chr') 
print (str(round(Auto_chr['Pi_4D_Blekinge'].median(), 4))+' ± '+str(round(Auto_chr['Pi_4D_Blekinge'].std(), 4))+' Blekinge Pi 4D Auto'+'  '+str(round(Z_chr['Pi_4D_Blekinge'].median(), 4))+' ± '+str(round(Z_chr['Pi_4D_Blekinge'].std(), 4))+' Blekinge Pi 4D Z-Chr') 
print (str(round(Auto_chr['Tw_4D_Blekinge'].median(), 4))+' ± '+str(round(Auto_chr['Tw_4D_Blekinge'].std(), 4))+' Blekinge Tw 4D Auto'+'  '+str(round(Z_chr['Tw_4D_Blekinge'].median(), 4))+' ± '+str(round(Z_chr['Tw_4D_Blekinge'].std(), 4))+' Blekinge Tw 4D Z-Chr') 
print (str(round(Auto_chr['Td_4D_Blekinge'].median(), 4))+' ± '+str(round(Auto_chr['Td_4D_Blekinge'].std(), 4))+' Blekinge Td 4D Auto'+'  '+str(round(Z_chr['Td_4D_Blekinge'].median(), 4))+' ± '+str(round(Z_chr['Td_4D_Blekinge'].std(), 4))+' Blekinge Td 4D Z-Chr') 
print (str(round(Auto_chr['Pi_4D_NordensArk'].median(), 4))+' ± '+str(round(Auto_chr['Pi_4D_NordensArk'].std(), 4))+' NordensArk Pi 4D Auto'+'  '+str(round(Z_chr['Pi_4D_NordensArk'].median(), 4))+' ± '+str(round(Z_chr['Pi_4D_NordensArk'].std(), 4))+' NordensArk Pi 4D Z-Chr') 
print (str(round(Auto_chr['Tw_4D_NordensArk'].median(), 4))+' ± '+str(round(Auto_chr['Tw_4D_NordensArk'].std(), 4))+' NordensArk Tw 4D Auto'+'  '+str(round(Z_chr['Tw_4D_NordensArk'].median(), 4))+' ± '+str(round(Z_chr['Tw_4D_NordensArk'].std(), 4))+' NordensArk Tw 4D Z-Chr') 
print (str(round(Auto_chr['Td_4D_NordensArk'].median(), 4))+' ± '+str(round(Auto_chr['Td_4D_NordensArk'].std(), 4))+' NordensArk Td 4D Auto'+'  '+str(round(Z_chr['Td_4D_NordensArk'].median(), 4))+' ± '+str(round(Z_chr['Td_4D_NordensArk'].std(), 4))+' NordensArk Td 4D Z-Chr') 
print (str(round(Auto_chr['Pi_4D_Uppland'].median(), 4))+' ± '+str(round(Auto_chr['Pi_4D_Uppland'].std(), 4))+' Uppland Pi 4D Auto'+'  '+str(round(Z_chr['Pi_4D_Uppland'].median(), 4))+' ± '+str(round(Z_chr['Pi_4D_Uppland'].std(), 4))+' Uppland Pi 4D Z-Chr') 
print (str(round(Auto_chr['Tw_4D_Uppland'].median(), 4))+' ± '+str(round(Auto_chr['Tw_4D_Uppland'].std(), 4))+' Uppland Tw 4D Auto'+'  '+str(round(Z_chr['Tw_4D_Uppland'].median(), 4))+' ± '+str(round(Z_chr['Tw_4D_Uppland'].std(), 4))+' Uppland Tw 4D Z-Chr') 
print (str(round(Auto_chr['Td_4D_Uppland'].median(), 4))+' ± '+str(round(Auto_chr['Td_4D_Uppland'].std(), 4))+' Uppland Td 4D Auto'+'  '+str(round(Z_chr['Td_4D_Uppland'].median(), 4))+' ± '+str(round(Z_chr['Td_4D_Uppland'].std(), 4))+' Uppland Td 4D Z-Chr') 
print (str(round(Auto_chr['Pi_3D_Vasternorrland'].median(), 4))+' ± '+str(round(Auto_chr['Pi_3D_Vasternorrland'].std(), 4))+' Vasternorrland Pi 3D Auto'+'  '+str(round(Z_chr['Pi_3D_Vasternorrland'].median(), 4))+' ± '+str(round(Z_chr['Pi_3D_Vasternorrland'].std(), 4))+' Vasternorrland Pi 3D Z-Chr') 
print (str(round(Auto_chr['Tw_3D_Vasternorrland'].median(), 4))+' ± '+str(round(Auto_chr['Tw_3D_Vasternorrland'].std(), 4))+' Vasternorrland Tw 3D Auto'+'  '+str(round(Z_chr['Tw_3D_Vasternorrland'].median(), 4))+' ± '+str(round(Z_chr['Tw_3D_Vasternorrland'].std(), 4))+' Vasternorrland Tw 3D Z-Chr') 
print (str(round(Auto_chr['Td_3D_Vasternorrland'].median(), 4))+' ± '+str(round(Auto_chr['Td_3D_Vasternorrland'].std(), 4))+' Vasternorrland Td 3D Auto'+'  '+str(round(Z_chr['Td_3D_Vasternorrland'].median(), 4))+' ± '+str(round(Z_chr['Td_3D_Vasternorrland'].std(), 4))+' Vasternorrland Td 3D Z-Chr') 
print (str(round(Auto_chr['Pi_3D_Blekinge'].median(), 4))+' ± '+str(round(Auto_chr['Pi_3D_Blekinge'].std(), 4))+' Blekinge Pi 3D Auto'+'  '+str(round(Z_chr['Pi_3D_Blekinge'].median(), 4))+' ± '+str(round(Z_chr['Pi_3D_Blekinge'].std(), 4))+' Blekinge Pi 3D Z-Chr') 
print (str(round(Auto_chr['Tw_3D_Blekinge'].median(), 4))+' ± '+str(round(Auto_chr['Tw_3D_Blekinge'].std(), 4))+' Blekinge Tw 3D Auto'+'  '+str(round(Z_chr['Tw_3D_Blekinge'].median(), 4))+' ± '+str(round(Z_chr['Tw_3D_Blekinge'].std(), 4))+' Blekinge Tw 3D Z-Chr') 
print (str(round(Auto_chr['Td_3D_Blekinge'].median(), 4))+' ± '+str(round(Auto_chr['Td_3D_Blekinge'].std(), 4))+' Blekinge Td 3D Auto'+'  '+str(round(Z_chr['Td_3D_Blekinge'].median(), 4))+' ± '+str(round(Z_chr['Td_3D_Blekinge'].std(), 4))+' Blekinge Td 3D Z-Chr') 
print (str(round(Auto_chr['Pi_3D_NordensArk'].median(), 4))+' ± '+str(round(Auto_chr['Pi_3D_NordensArk'].std(), 4))+' NordensArk Pi 3D Auto'+'  '+str(round(Z_chr['Pi_3D_NordensArk'].median(), 4))+' ± '+str(round(Z_chr['Pi_3D_NordensArk'].std(), 4))+' NordensArk Pi 3D Z-Chr') 
print (str(round(Auto_chr['Tw_3D_NordensArk'].median(), 4))+' ± '+str(round(Auto_chr['Tw_3D_NordensArk'].std(), 4))+' NordensArk Tw 3D Auto'+'  '+str(round(Z_chr['Tw_3D_NordensArk'].median(), 4))+' ± '+str(round(Z_chr['Tw_3D_NordensArk'].std(), 4))+' NordensArk Tw 3D Z-Chr') 
print (str(round(Auto_chr['Td_3D_NordensArk'].median(), 4))+' ± '+str(round(Auto_chr['Td_3D_NordensArk'].std(), 4))+' NordensArk Td 3D Auto'+'  '+str(round(Z_chr['Td_3D_NordensArk'].median(), 4))+' ± '+str(round(Z_chr['Td_3D_NordensArk'].std(), 4))+' NordensArk Td 3D Z-Chr') 
print (str(round(Auto_chr['Pi_3D_Uppland'].median(), 4))+' ± '+str(round(Auto_chr['Pi_3D_Uppland'].std(), 4))+' Uppland Pi 3D Auto'+'  '+str(round(Z_chr['Pi_3D_Uppland'].median(), 4))+' ± '+str(round(Z_chr['Pi_3D_Uppland'].std(), 4))+' Uppland Pi 3D Z-Chr') 
print (str(round(Auto_chr['Tw_3D_Uppland'].median(), 4))+' ± '+str(round(Auto_chr['Tw_3D_Uppland'].std(), 4))+' Uppland Tw 3D Auto'+'  '+str(round(Z_chr['Tw_3D_Uppland'].median(), 4))+' ± '+str(round(Z_chr['Tw_3D_Uppland'].std(), 4))+' Uppland Tw 3D Z-Chr') 
print (str(round(Auto_chr['Td_3D_Uppland'].median(), 4))+' ± '+str(round(Auto_chr['Td_3D_Uppland'].std(), 4))+' Uppland Td 3D Auto'+'  '+str(round(Z_chr['Td_3D_Uppland'].median(), 4))+' ± '+str(round(Z_chr['Td_3D_Uppland'].std(), 4))+' Uppland Td 3D Z-Chr') 
print (str(round(Auto_chr['Pi_2D_Vasternorrland'].median(), 4))+' ± '+str(round(Auto_chr['Pi_2D_Vasternorrland'].std(), 4))+' Vasternorrland Pi 2D Auto'+'  '+str(round(Z_chr['Pi_2D_Vasternorrland'].median(), 4))+' ± '+str(round(Z_chr['Pi_2D_Vasternorrland'].std(), 4))+' Vasternorrland Pi 2D Z-Chr') 
print (str(round(Auto_chr['Tw_2D_Vasternorrland'].median(), 4))+' ± '+str(round(Auto_chr['Tw_2D_Vasternorrland'].std(), 4))+' Vasternorrland Tw 2D Auto'+'  '+str(round(Z_chr['Tw_2D_Vasternorrland'].median(), 4))+' ± '+str(round(Z_chr['Tw_2D_Vasternorrland'].std(), 4))+' Vasternorrland Tw 2D Z-Chr') 
print (str(round(Auto_chr['Td_2D_Vasternorrland'].median(), 4))+' ± '+str(round(Auto_chr['Td_2D_Vasternorrland'].std(), 4))+' Vasternorrland Td 2D Auto'+'  '+str(round(Z_chr['Td_2D_Vasternorrland'].median(), 4))+' ± '+str(round(Z_chr['Td_2D_Vasternorrland'].std(), 4))+' Vasternorrland Td 2D Z-Chr') 
print (str(round(Auto_chr['Pi_2D_Blekinge'].median(), 4))+' ± '+str(round(Auto_chr['Pi_2D_Blekinge'].std(), 4))+' Blekinge Pi 2D Auto'+'  '+str(round(Z_chr['Pi_2D_Blekinge'].median(), 4))+' ± '+str(round(Z_chr['Pi_2D_Blekinge'].std(), 4))+' Blekinge Pi 2D Z-Chr') 
print (str(round(Auto_chr['Tw_2D_Blekinge'].median(), 4))+' ± '+str(round(Auto_chr['Tw_2D_Blekinge'].std(), 4))+' Blekinge Tw 2D Auto'+'  '+str(round(Z_chr['Tw_2D_Blekinge'].median(), 4))+' ± '+str(round(Z_chr['Tw_2D_Blekinge'].std(), 4))+' Blekinge Tw 2D Z-Chr') 
print (str(round(Auto_chr['Td_2D_Blekinge'].median(), 4))+' ± '+str(round(Auto_chr['Td_2D_Blekinge'].std(), 4))+' Blekinge Td 2D Auto'+'  '+str(round(Z_chr['Td_2D_Blekinge'].median(), 4))+' ± '+str(round(Z_chr['Td_2D_Blekinge'].std(), 4))+' Blekinge Td 2D Z-Chr') 
print (str(round(Auto_chr['Pi_2D_NordensArk'].median(), 4))+' ± '+str(round(Auto_chr['Pi_2D_NordensArk'].std(), 4))+' NordensArk Pi 2D Auto'+'  '+str(round(Z_chr['Pi_2D_NordensArk'].median(), 4))+' ± '+str(round(Z_chr['Pi_2D_NordensArk'].std(), 4))+' NordensArk Pi 2D Z-Chr') 
print (str(round(Auto_chr['Tw_2D_NordensArk'].median(), 4))+' ± '+str(round(Auto_chr['Tw_2D_NordensArk'].std(), 4))+' NordensArk Tw 2D Auto'+'  '+str(round(Z_chr['Tw_2D_NordensArk'].median(), 4))+' ± '+str(round(Z_chr['Tw_2D_NordensArk'].std(), 4))+' NordensArk Tw 2D Z-Chr') 
print (str(round(Auto_chr['Td_2D_NordensArk'].median(), 4))+' ± '+str(round(Auto_chr['Td_2D_NordensArk'].std(), 4))+' NordensArk Td 2D Auto'+'  '+str(round(Z_chr['Td_2D_NordensArk'].median(), 4))+' ± '+str(round(Z_chr['Td_2D_NordensArk'].std(), 4))+' NordensArk Td 2D Z-Chr') 
print (str(round(Auto_chr['Pi_2D_Uppland'].median(), 4))+' ± '+str(round(Auto_chr['Pi_2D_Uppland'].std(), 4))+' Uppland Pi 2D Auto'+'  '+str(round(Z_chr['Pi_2D_Uppland'].median(), 4))+' ± '+str(round(Z_chr['Pi_2D_Uppland'].std(), 4))+' Uppland Pi 2D Z-Chr') 
print (str(round(Auto_chr['Tw_2D_Uppland'].median(), 4))+' ± '+str(round(Auto_chr['Tw_2D_Uppland'].std(), 4))+' Uppland Tw 2D Auto'+'  '+str(round(Z_chr['Tw_2D_Uppland'].median(), 4))+' ± '+str(round(Z_chr['Tw_2D_Uppland'].std(), 4))+' Uppland Tw 2D Z-Chr') 
print (str(round(Auto_chr['Td_2D_Uppland'].median(), 4))+' ± '+str(round(Auto_chr['Td_2D_Uppland'].std(), 4))+' Uppland Td 2D Auto'+'  '+str(round(Z_chr['Td_2D_Uppland'].median(), 4))+' ± '+str(round(Z_chr['Td_2D_Uppland'].std(), 4))+' Uppland Td 2D Z-Chr') 
print (str(round(Auto_chr['Pi_1D_Vasternorrland'].median(), 4))+' ± '+str(round(Auto_chr['Pi_1D_Vasternorrland'].std(), 4))+' Vasternorrland Pi 1D Auto'+'  '+str(round(Z_chr['Pi_1D_Vasternorrland'].median(), 4))+' ± '+str(round(Z_chr['Pi_1D_Vasternorrland'].std(), 4))+' Vasternorrland Pi 1D Z-Chr') 
print (str(round(Auto_chr['Tw_1D_Vasternorrland'].median(), 4))+' ± '+str(round(Auto_chr['Tw_1D_Vasternorrland'].std(), 4))+' Vasternorrland Tw 1D Auto'+'  '+str(round(Z_chr['Tw_1D_Vasternorrland'].median(), 4))+' ± '+str(round(Z_chr['Tw_1D_Vasternorrland'].std(), 4))+' Vasternorrland Tw 1D Z-Chr') 
print (str(round(Auto_chr['Td_1D_Vasternorrland'].median(), 4))+' ± '+str(round(Auto_chr['Td_1D_Vasternorrland'].std(), 4))+' Vasternorrland Td 1D Auto'+'  '+str(round(Z_chr['Td_1D_Vasternorrland'].median(), 4))+' ± '+str(round(Z_chr['Td_1D_Vasternorrland'].std(), 4))+' Vasternorrland Td 1D Z-Chr') 
print (str(round(Auto_chr['Pi_1D_Blekinge'].median(), 4))+' ± '+str(round(Auto_chr['Pi_1D_Blekinge'].std(), 4))+' Blekinge Pi 1D Auto'+'  '+str(round(Z_chr['Pi_1D_Blekinge'].median(), 4))+' ± '+str(round(Z_chr['Pi_1D_Blekinge'].std(), 4))+' Blekinge Pi 1D Z-Chr') 
print (str(round(Auto_chr['Tw_1D_Blekinge'].median(), 4))+' ± '+str(round(Auto_chr['Tw_1D_Blekinge'].std(), 4))+' Blekinge Tw 1D Auto'+'  '+str(round(Z_chr['Tw_1D_Blekinge'].median(), 4))+' ± '+str(round(Z_chr['Tw_1D_Blekinge'].std(), 4))+' Blekinge Tw 1D Z-Chr') 
print (str(round(Auto_chr['Td_1D_Blekinge'].median(), 4))+' ± '+str(round(Auto_chr['Td_1D_Blekinge'].std(), 4))+' Blekinge Td 1D Auto'+'  '+str(round(Z_chr['Td_1D_Blekinge'].median(), 4))+' ± '+str(round(Z_chr['Td_1D_Blekinge'].std(), 4))+' Blekinge Td 1D Z-Chr') 
print (str(round(Auto_chr['Pi_1D_NordensArk'].median(), 4))+' ± '+str(round(Auto_chr['Pi_1D_NordensArk'].std(), 4))+' NordensArk Pi 1D Auto'+'  '+str(round(Z_chr['Pi_1D_NordensArk'].median(), 4))+' ± '+str(round(Z_chr['Pi_1D_NordensArk'].std(), 4))+' NordensArk Pi 1D Z-Chr') 
print (str(round(Auto_chr['Tw_1D_NordensArk'].median(), 4))+' ± '+str(round(Auto_chr['Tw_1D_NordensArk'].std(), 4))+' NordensArk Tw 1D Auto'+'  '+str(round(Z_chr['Tw_1D_NordensArk'].median(), 4))+' ± '+str(round(Z_chr['Tw_1D_NordensArk'].std(), 4))+' NordensArk Tw 1D Z-Chr') 
print (str(round(Auto_chr['Td_1D_NordensArk'].median(), 4))+' ± '+str(round(Auto_chr['Td_1D_NordensArk'].std(), 4))+' NordensArk Td 1D Auto'+'  '+str(round(Z_chr['Td_1D_NordensArk'].median(), 4))+' ± '+str(round(Z_chr['Td_1D_NordensArk'].std(), 4))+' NordensArk Td 1D Z-Chr') 
print (str(round(Auto_chr['Pi_1D_Uppland'].median(), 4))+' ± '+str(round(Auto_chr['Pi_1D_Uppland'].std(), 4))+' Uppland Pi 1D Auto'+'  '+str(round(Z_chr['Pi_1D_Uppland'].median(), 4))+' ± '+str(round(Z_chr['Pi_1D_Uppland'].std(), 4))+' Uppland Pi 1D Z-Chr') 
print (str(round(Auto_chr['Tw_1D_Uppland'].median(), 4))+' ± '+str(round(Auto_chr['Tw_1D_Uppland'].std(), 4))+' Uppland Tw 1D Auto'+'  '+str(round(Z_chr['Tw_1D_Uppland'].median(), 4))+' ± '+str(round(Z_chr['Tw_1D_Uppland'].std(), 4))+' Uppland Tw 1D Z-Chr') 
print (str(round(Auto_chr['Td_1D_Uppland'].median(), 4))+' ± '+str(round(Auto_chr['Td_1D_Uppland'].std(), 4))+' Uppland Td 1D Auto'+'  '+str(round(Z_chr['Td_1D_Uppland'].median(), 4))+' ± '+str(round(Z_chr['Td_1D_Uppland'].std(), 4))+' Uppland Td 1D Z-Chr') 












































list_2_Fst=list_2_Fst.merge(Chrom[['CHROM','Chr_number']], on='CHROM')


Z_chr_fst=list_2_Fst[list_2_Fst.Chr_number=='Chr_Z']
Auto_chr_fst=list_2_Fst[list_2_Fst.Chr_number!='Chr_Z']


list_2_Fst.loc[(list_2_Fst.Chr_number!='Chr_Z'), 'CHROM_Type']='Auto'
list_2_Fst.loc[(list_2_Fst.Chr_number=='Chr_Z'), 'CHROM_Type']='Chr_Z'





print (str(round(Auto_chr_fst['Blekinge_Uppland_fst'].median(), 4))+' ± '+str(round(Auto_chr_fst['Blekinge_Uppland_fst'].std(), 4))+' Blekinge_Uppland_fst 4D Auto' + '  '+str(round(Z_chr_fst['Blekinge_Uppland_fst'].median(), 4))+' ± '+str(round(Z_chr_fst['Blekinge_Uppland_fst'].std(), 4))+' Blekinge_Uppland_fst 4D Z-Chr')
print (str(round(Auto_chr_fst['Uppland_NordensArk_fst'].median(), 4))+' ± '+str(round(Auto_chr_fst['Uppland_NordensArk_fst'].std(), 4))+' Uppland_NordensArk_fst 4D Auto' + '  '+str(round(Z_chr_fst['Uppland_NordensArk_fst'].median(), 4))+' ± '+str(round(Z_chr_fst['Uppland_NordensArk_fst'].std(), 4))+' Uppland_NordensArk_fst 4D Z-Chr')
print (str(round(Auto_chr_fst['Vasternorrland_Uppland_fst'].median(), 4))+' ± '+str(round(Auto_chr_fst['Vasternorrland_Uppland_fst'].std(), 4))+' Vasternorrland_Uppland_fst 4D Auto' + '  '+str(round(Z_chr_fst['Vasternorrland_Uppland_fst'].median(), 4))+' ± '+str(round(Z_chr_fst['Vasternorrland_Uppland_fst'].std(), 4))+' Vasternorrland_Uppland_fst 4D Z-Chr')
print (str(round(Auto_chr_fst['Blekinge_NordensArk_fst'].median(), 4))+' ± '+str(round(Auto_chr_fst['Blekinge_NordensArk_fst'].std(), 4))+' Blekinge_NordensArk_fst 4D Auto' + '  '+str(round(Z_chr_fst['Blekinge_NordensArk_fst'].median(), 4))+' ± '+str(round(Z_chr_fst['Blekinge_NordensArk_fst'].std(), 4))+' Blekinge_NordensArk_fst 4D Z-Chr')
print (str(round(Auto_chr_fst['Vasternorrland_Blekinge_fst'].median(), 4))+' ± '+str(round(Auto_chr_fst['Vasternorrland_Blekinge_fst'].std(), 4))+' Vasternorrland_Blekinge_fst 4D Auto' + '  '+str(round(Z_chr_fst['Vasternorrland_Blekinge_fst'].median(), 4))+' ± '+str(round(Z_chr_fst['Vasternorrland_Blekinge_fst'].std(), 4))+' Vasternorrland_Blekinge_fst 4D Z-Chr')
print (str(round(Auto_chr_fst['Vasternorrland_NordensArk_fst'].median(), 4))+' ± '+str(round(Auto_chr_fst['Vasternorrland_NordensArk_fst'].std(), 4))+' Vasternorrland_NordensArk_fst 4D Auto' + '  '+str(round(Z_chr_fst['Vasternorrland_NordensArk_fst'].median(), 4))+' ± '+str(round(Z_chr_fst['Vasternorrland_NordensArk_fst'].std(), 4))+' Vasternorrland_NordensArk_fst 4D Z-Chr')





print (str(round(Auto_chr_fst['Blekinge_Uppland_dxy'].median(), 4))+' ± '+str(round(Auto_chr_fst['Blekinge_Uppland_dxy'].std(), 4))+' Blekinge_Uppland_dxy 4D Auto' + '  '+str(round(Z_chr_fst['Blekinge_Uppland_dxy'].median(), 4))+' ± '+str(round(Z_chr_fst['Blekinge_Uppland_dxy'].std(), 4))+' Blekinge_Uppland_dxy 4D Z-Chr')
print (str(round(Auto_chr_fst['Uppland_NordensArk_dxy'].median(), 4))+' ± '+str(round(Auto_chr_fst['Uppland_NordensArk_dxy'].std(), 4))+' Uppland_NordensArk_dxy 4D Auto' + '  '+str(round(Z_chr_fst['Uppland_NordensArk_dxy'].median(), 4))+' ± '+str(round(Z_chr_fst['Uppland_NordensArk_dxy'].std(), 4))+' Uppland_NordensArk_dxy 4D Z-Chr')
print (str(round(Auto_chr_fst['Vasternorrland_Uppland_dxy'].median(), 4))+' ± '+str(round(Auto_chr_fst['Vasternorrland_Uppland_dxy'].std(), 4))+' Vasternorrland_Uppland_dxy 4D Auto' + '  '+str(round(Z_chr_fst['Vasternorrland_Uppland_dxy'].median(), 4))+' ± '+str(round(Z_chr_fst['Vasternorrland_Uppland_dxy'].std(), 4))+' Vasternorrland_Uppland_dxy 4D Z-Chr')
print (str(round(Auto_chr_fst['Blekinge_NordensArk_dxy'].median(), 4))+' ± '+str(round(Auto_chr_fst['Blekinge_NordensArk_dxy'].std(), 4))+' Blekinge_NordensArk_dxy 4D Auto' + '  '+str(round(Z_chr_fst['Blekinge_NordensArk_dxy'].median(), 4))+' ± '+str(round(Z_chr_fst['Blekinge_NordensArk_dxy'].std(), 4))+' Blekinge_NordensArk_dxy 4D Z-Chr')
print (str(round(Auto_chr_fst['Vasternorrland_Blekinge_dxy'].median(), 4))+' ± '+str(round(Auto_chr_fst['Vasternorrland_Blekinge_dxy'].std(), 4))+' Vasternorrland_Blekinge_dxy 4D Auto' + '  '+str(round(Z_chr_fst['Vasternorrland_Blekinge_dxy'].median(), 4))+' ± '+str(round(Z_chr_fst['Vasternorrland_Blekinge_dxy'].std(), 4))+' Vasternorrland_Blekinge_dxy 4D Z-Chr')
print (str(round(Auto_chr_fst['Vasternorrland_NordensArk_dxy'].median(), 4))+' ± '+str(round(Auto_chr_fst['Vasternorrland_NordensArk_dxy'].std(), 4))+' Vasternorrland_NordensArk_dxy 4D Auto' + '  '+str(round(Z_chr_fst['Vasternorrland_NordensArk_dxy'].median(), 4))+' ± '+str(round(Z_chr_fst['Vasternorrland_NordensArk_dxy'].std(), 4))+' Vasternorrland_NordensArk_dxy 4D Z-Chr')








#################

def manhattan_plot_fst_dxy_pi(sites_filter_0, fst, dxy,pi,taj,fixed,chromosomes_dict,colour):
	species_colour=sns.color_palette("Set2", 8)[:3]
	chromosome=["Chr_2","Chr_3","Chr_4","Chr_5","Chr_6","Chr_7","Chr_8","Chr_9","Chr_10","Chr_11","Chr_12","Chr_13","Chr_14","Chr_15","Chr_16","Chr_17","Chr_18","Chr_19","Chr_20","Chr_21","Chr_22","Chr_23","Chr_24","Chr_25","Chr_26","Chr_27","Chr_28","Chr_29","Chr_30","Chr_Z"]
	pi_selection=pi
	fst_selection=fst
	dxy_selection=dxy
	taj_selection=taj
	fst_set=[]
	dxy_set=[]
	fixed_set=[]
	pi_set=[]
	taj_set=[]
	chrom={}
	for j in chromosome:
		chrom[j]=0
		for i in chromosomes_dict[str(j)]:
			pi_list=sites_filter_0[(sites_filter_0['CHROM']==i)][pi_selection]
			fixed_1=sites_filter_0[(sites_filter_0['CHROM']==i)][fixed]
			chrom[j]+=len(sites_filter_0[(sites_filter_0['CHROM']==i)])
			fst_list=sites_filter_0[(sites_filter_0['CHROM']==i)][fst_selection]
			dxy_list=sites_filter_0[(sites_filter_0['CHROM']==i)][dxy_selection]
			taj_list=sites_filter_0[(sites_filter_0['CHROM']==i)][taj_selection]
			if len(fst_list)>=5:
				fst_=fst_list[fst].rolling(window=5).mean()
				fst_set.append(fst_)
			if len(dxy_list)>=5:
				dxy_=dxy_list[dxy].rolling(window=5).mean()
				dxy_set.append(dxy_)
			if len(pi_list)>=5:
				pi_=pi_list.rolling(window=5).mean()
				pi_set.append(pi_)
				#pi_.mean().plot(fontsize=10,legend=True)
			if len(fixed_1)>=1:
				fixed_=fixed_1.rolling(window=5).mean()
				fixed_set.append(fixed_)
			if len(taj_list)>=5:
				taj_=taj_list.rolling(window=5).mean()
				taj_set.append(taj_)
				#pi_.mean().plot(fontsize=10,legend=True)
	fst_set2=pandas.DataFrame(columns=fst_selection)
	for j in fst_set:
		list_j=[fst_set2,j]
		fst_set2=pandas.concat(list_j,ignore_index=True)
	
		
	dxy_set2=pandas.DataFrame(columns=dxy_selection)
	for n in dxy_set:
		list_n=[dxy_set2,n]
		dxy_set2=pandas.concat(list_n,ignore_index=True)

	fixed_set2=pandas.DataFrame(columns=fixed)
	for f in fixed_set:
		list_f=[fixed_set2,f]
		fixed_set2=pandas.concat(list_f,ignore_index=True)

	pi_set2=pandas.DataFrame(columns=pi_selection)
	for p in pi_set:
		list_p=[pi_set2,p]
		pi_set2=pandas.concat(list_p,ignore_index=True)

	taj_set2=pandas.DataFrame(columns=taj_selection)
	for p in taj_set:
		list_p=[taj_set2,p]
		taj_set2=pandas.concat(list_p,ignore_index=True)

	fig, axes = plt.subplots(nrows=5, sharex=True)
	fig.subplots_adjust(hspace=0.1)
	fixed_set2.plot.area(ax=axes[0],legend=False,color=sns.color_palette(['#27AE60','#D35400','#F4D03F','#2C3E50']),stacked=False, alpha=0.5).set_ylabel('Taj 4D', fontsize=15)
	fst_set2.plot(ax=axes[1],legend=False, alpha=0.6,color=sns.color_palette(['#27AE60','#D35400','#F4D03F','#2C3E50'])).set_ylabel(r'$\pi$ 4D', fontsize=15)
	dxy_set2.plot(ax=axes[2], legend=False, alpha=0.6,color=sns.color_palette(['#27AE60','#D35400','#F4D03F','#2C3E50'])).set_ylabel(r'$\Theta$ 4D', fontsize=15)
	pi_set2.plot(ax=axes[3], legend=False, alpha=0.6 ,color=sns.color_palette()).set_ylabel(r'Fst', fontsize=15)
	taj_set2.plot(ax=axes[4], legend=False, alpha=0.6 ,color=sns.color_palette()).set_ylabel(r'Dxy', fontsize=15)
	li=[]
	for c in chromosome:
		if c == 'Chr_2':
			start=0
			end=int(chrom[str(c)])
		else:
			start=end
			end=start+int(chrom[str(c)])
		li.append([start, end])
	for l in li[0::2]:
		axes[0].axvspan(l[0], l[1], facecolor='#dbd9d9', alpha=0.5)
		axes[1].axvspan(l[0], l[1], facecolor='#dbd9d9', alpha=0.5)
		axes[2].axvspan(l[0], l[1], facecolor='#dbd9d9', alpha=0.5)
		axes[3].axvspan(l[0], l[1], facecolor='#dbd9d9', alpha=0.5)
		axes[4].axvspan(l[0], l[1], facecolor='#dbd9d9', alpha=0.5)
	ticks = axes[4].get_xticks()/100
	axes[4].set_xticklabels(ticks)
	axes[4].set_xlabel("Mega bases")
	return fig





list_2[['Pi_4D_Vasternorrland','Tw_4D_Vasternorrland','Td_4D_Vasternorrland','N_sites_4D_Vasternorrland','Pi_4D_Blekinge','Tw_4D_Blekinge','Td_4D_Blekinge','N_sites_4D_Blekinge','Pi_4D_NordensArk','Tw_4D_NordensArk','Td_4D_NordensArk','N_sites_4D_NordensArk','Pi_4D_Uppland','Tw_4D_Uppland','Td_4D_Uppland','N_sites_4D_Uppland','Pi_3D_Vasternorrland','Tw_3D_Vasternorrland','Td_3D_Vasternorrland','N_sites_3D_Vasternorrland','Pi_3D_Blekinge','Tw_3D_Blekinge','Td_3D_Blekinge','N_sites_3D_Blekinge','Pi_3D_NordensArk','Tw_3D_NordensArk','Td_3D_NordensArk','N_sites_3D_NordensArk','Pi_3D_Uppland','Tw_3D_Uppland','Td_3D_Uppland','N_sites_3D_Uppland','Pi_2D_Vasternorrland','Tw_2D_Vasternorrland','Td_2D_Vasternorrland','N_sites_2D_Vasternorrland','Pi_2D_Blekinge','Tw_2D_Blekinge','Td_2D_Blekinge','N_sites_2D_Blekinge','Pi_2D_NordensArk','Tw_2D_NordensArk','Td_2D_NordensArk','N_sites_2D_NordensArk','Pi_2D_Uppland','Tw_2D_Uppland','Td_2D_Uppland','N_sites_2D_Uppland','Pi_1D_Vasternorrland','Tw_1D_Vasternorrland','Td_1D_Vasternorrland','N_sites_1D_Vasternorrland','Pi_1D_Blekinge','Tw_1D_Blekinge','Td_1D_Blekinge','N_sites_1D_Blekinge','Pi_1D_NordensArk','Tw_1D_NordensArk','Td_1D_NordensArk','N_sites_1D_NordensArk','Pi_1D_Uppland','Tw_1D_Uppland','Td_1D_Uppland','N_sites_1D_Uppland']]=list_2[['Pi_4D_Vasternorrland','Tw_4D_Vasternorrland','Td_4D_Vasternorrland','N_sites_4D_Vasternorrland','Pi_4D_Blekinge','Tw_4D_Blekinge','Td_4D_Blekinge','N_sites_4D_Blekinge','Pi_4D_NordensArk','Tw_4D_NordensArk','Td_4D_NordensArk','N_sites_4D_NordensArk','Pi_4D_Uppland','Tw_4D_Uppland','Td_4D_Uppland','N_sites_4D_Uppland','Pi_3D_Vasternorrland','Tw_3D_Vasternorrland','Td_3D_Vasternorrland','N_sites_3D_Vasternorrland','Pi_3D_Blekinge','Tw_3D_Blekinge','Td_3D_Blekinge','N_sites_3D_Blekinge','Pi_3D_NordensArk','Tw_3D_NordensArk','Td_3D_NordensArk','N_sites_3D_NordensArk','Pi_3D_Uppland','Tw_3D_Uppland','Td_3D_Uppland','N_sites_3D_Uppland','Pi_2D_Vasternorrland','Tw_2D_Vasternorrland','Td_2D_Vasternorrland','N_sites_2D_Vasternorrland','Pi_2D_Blekinge','Tw_2D_Blekinge','Td_2D_Blekinge','N_sites_2D_Blekinge','Pi_2D_NordensArk','Tw_2D_NordensArk','Td_2D_NordensArk','N_sites_2D_NordensArk','Pi_2D_Uppland','Tw_2D_Uppland','Td_2D_Uppland','N_sites_2D_Uppland','Pi_1D_Vasternorrland','Tw_1D_Vasternorrland','Td_1D_Vasternorrland','N_sites_1D_Vasternorrland','Pi_1D_Blekinge','Tw_1D_Blekinge','Td_1D_Blekinge','N_sites_1D_Blekinge','Pi_1D_NordensArk','Tw_1D_NordensArk','Td_1D_NordensArk','N_sites_1D_NordensArk','Pi_1D_Uppland','Tw_1D_Uppland','Td_1D_Uppland','N_sites_1D_Uppland']].astype(float)


manhattan_plot_fst_dxy_pi(list_2, ['Pi_4D_Vasternorrland','Pi_4D_Blekinge','Pi_4D_NordensArk','Pi_4D_Uppland'], ['Tw_4D_Vasternorrland','Tw_4D_Blekinge','Tw_4D_NordensArk','Tw_4D_Uppland'], ['Blekinge_Uppland_fst','Uppland_NordensArk_fst','Vasternorrland_Uppland_fst','Blekinge_NordensArk_fst','Vasternorrland_Blekinge_fst','Vasternorrland_NordensArk_fst'], ['Blekinge_Uppland_dxy','Uppland_NordensArk_dxy','Vasternorrland_Uppland_dxy','Blekinge_NordensArk_dxy','Vasternorrland_Blekinge_dxy','Vasternorrland_NordensArk_dxy'], ['Td_4D_Vasternorrland','Td_4D_Blekinge','Td_4D_NordensArk','Td_4D_Uppland'],CHROM_dict,'#050505')
plt.savefig('scan_FST_dxy.png', dpi=1200)
plt.close()


sns.boxplot(y="variable", x="value", showmeans=True,meanprops={"marker": "+","markeredgecolor": "white","markersize": "10"},showfliers = False, data=pandas.melt(list_2_Fst[list_2_Fst.CHROM_Type=='Auto'][['Blekinge_NordensArk_fst','Blekinge_Uppland_fst', 'Vasternorrland_Blekinge_fst','Uppland_NordensArk_fst','Vasternorrland_NordensArk_fst','Vasternorrland_Uppland_fst','CHROM_Type']], id_vars='CHROM_Type'))
plt.savefig('FST_box.png', dpi=1200)
plt.close()





sns.boxplot(y="variable", x="value", showmeans=True,meanprops={"marker": "+","markeredgecolor": "white","markersize": "10"},showfliers = False, data=pandas.melt(list_2_Fst[list_2_Fst.CHROM_Type=='Auto'][['Blekinge_NordensArk_dxy','Blekinge_Uppland_dxy', 'Vasternorrland_Blekinge_dxy','Uppland_NordensArk_dxy','Vasternorrland_NordensArk_dxy','Vasternorrland_Uppland_dxy','CHROM_Type']], id_vars='CHROM_Type'))
plt.savefig('Dxy_box.png', dpi=1200)
plt.close()





axes = sns.pairplot(list_2_Fst[['Vasternorrland_Blekinge_fst','Uppland_NordensArk_fst','Vasternorrland_Uppland_fst','Vasternorrland_NordensArk_fst','Blekinge_Uppland_fst','Blekinge_NordensArk_fst','CHROM_Type']], hue="CHROM_Type", plot_kws = {'alpha': 0.3, 'edgecolor': 'k'}, size = 4)

grid = sns.PairGrid(data= list_2_Fst[['Vasternorrland_Blekinge_fst','Uppland_NordensArk_fst','Vasternorrland_Uppland_fst','Vasternorrland_NordensArk_fst','Blekinge_Uppland_fst','Blekinge_NordensArk_fst','CHROM_Type']], hue="CHROM_Type", size = 4)

grid = grid.map_upper(plt.scatter)

grid = grid.map_lower(sns.kdeplot)

grid = grid.map_diag(sns.kdeplot)





axes = sns.pairplot(list_2[['Tw_4D_Vasternorrland','Tw_4D_Blekinge','Tw_4D_NordensArk','Tw_4D_Uppland']], plot_kws = {'alpha': 0.6, 'edgecolor': 'k'}, size = 4)

grid = sns.PairGrid(data= list_2[['Tw_4D_Vasternorrland','Tw_4D_Blekinge','Tw_4D_NordensArk','Tw_4D_Uppland']], size = 4)

grid = grid.map_upper(plt.scatter)

grid = grid.map_lower(sns.kdeplot)

grid = grid.map_diag(plt.hist, bins = 10, edgecolor = 'k')





corr = list_2[['Vasternorrland_Blekinge_fst','Uppland_NordensArk_fst','Vasternorrland_Uppland_fst','Vasternorrland_NordensArk_fst','Blekinge_Uppland_fst','Blekinge_NordensArk_fst']].corr()

sns.jointplot(x='Vasternorrland_Blekinge_fst', y='N_sites', data=list_2)


list_2[['Vasternorrland_Blekinge_fst','Uppland_NordensArk_fst','Vasternorrland_Uppland_fst','Vasternorrland_NordensArk_fst','Blekinge_Uppland_fst','Blekinge_NordensArk_fst', 'N_sites']].corr()





print (Blekinge_NordensArk_Fst[Blekinge_NordensArk_Fst['CHROM'].isin(set(Auto_chr.CHROM))][['Blekinge_NordensArk_fixed','Blekinge_NordensArk_private_a','Blekinge_NordensArk_private_b','Blekinge_NordensArk_shared']].sum())
print (Blekinge_Uppland_Fst[Blekinge_Uppland_Fst['CHROM'].isin(set(Auto_chr.CHROM))][['Blekinge_Uppland_fixed','Blekinge_Uppland_private_a','Blekinge_Uppland_private_b','Blekinge_Uppland_shared']].sum())
print (Uppland_NordensArk_Fst[Uppland_NordensArk_Fst['CHROM'].isin(set(Auto_chr.CHROM))][['Uppland_NordensArk_fixed','Uppland_NordensArk_private_a','Uppland_NordensArk_private_b','Uppland_NordensArk_shared']].sum())
print (Vasternorrland_Blekinge_Fst[Vasternorrland_Blekinge_Fst['CHROM'].isin(set(Auto_chr.CHROM))][['Vasternorrland_Blekinge_fixed','Vasternorrland_Blekinge_private_a','Vasternorrland_Blekinge_private_b','Vasternorrland_Blekinge_shared']].sum())
print (Vasternorrland_NordensArk_Fst[Vasternorrland_NordensArk_Fst['CHROM'].isin(set(Auto_chr.CHROM))][['Vasternorrland_NordensArk_fixed','Vasternorrland_NordensArk_private_a','Vasternorrland_NordensArk_private_b','Vasternorrland_NordensArk_shared']].sum())
print (Vasternorrland_Uppland_Fst[Vasternorrland_Uppland_Fst['CHROM'].isin(set(Auto_chr.CHROM))][['Vasternorrland_Uppland_fixed','Vasternorrland_Uppland_private_a','Vasternorrland_Uppland_private_b','Vasternorrland_Uppland_shared']].sum())





print (Blekinge_NordensArk_Fst[['Blekinge_NordensArk_fixed','Blekinge_NordensArk_private_a','Blekinge_NordensArk_private_b','Blekinge_NordensArk_shared']].sum())
print (Blekinge_Uppland_Fst[['Blekinge_Uppland_fixed','Blekinge_Uppland_private_a','Blekinge_Uppland_private_b','Blekinge_Uppland_shared']].sum())
print (Vasternorrland_Blekinge_Fst[['Vasternorrland_Blekinge_fixed','Vasternorrland_Blekinge_private_a','Vasternorrland_Blekinge_private_b','Vasternorrland_Blekinge_shared']].sum())
print (Uppland_NordensArk_Fst[['Uppland_NordensArk_fixed','Uppland_NordensArk_private_a','Uppland_NordensArk_private_b','Uppland_NordensArk_shared']].sum())
print (Vasternorrland_NordensArk_Fst[['Vasternorrland_NordensArk_fixed','Vasternorrland_NordensArk_private_a','Vasternorrland_NordensArk_private_b','Vasternorrland_NordensArk_shared']].sum())
print (Vasternorrland_Uppland_Fst[['Vasternorrland_Uppland_fixed','Vasternorrland_Uppland_private_a','Vasternorrland_Uppland_private_b','Vasternorrland_Uppland_shared']].sum())

