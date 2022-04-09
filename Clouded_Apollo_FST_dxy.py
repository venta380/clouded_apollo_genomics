import sys
import os
import string
import pandas
import personal_popgen
import itertools
import numpy as np
import argparse
import threading
import multiprocessing
from multiprocessing import  Pool



#python Monarch_east_west_pi_DXY.py -P1 East -P2 1977 -o East_West_Dxy_0 -N east_1977_N_sites -F east_1977_all_0_N_sites

pwd='/scratch/vt20265/clouded_apollo/mapping/'
os.chdir(pwd)

freq_list=["Brudskär_4D","Lötaholmen_4D"]
sites_col='N_sites'



lists_1=pandas.read_csv('windows.csv')
CHROMs=sorted(set(lists_1.CHROM))

primary_keys=lists_1

all_sites=pandas.read_csv('codon_df_4d.csv',names=['CHROM','POS'], sep='\t', header=None)
all_sites['POS']=all_sites['POS'].astype('int')
all_sites['BIN_START']=0
all_sites['BIN_START']=(np.floor(all_sites['POS']/100000)*100000)+1
all_sites['BIN_END']=all_sites['BIN_START']+(100000-1)
all_sites_2=all_sites.groupby(['CHROM','BIN_START','BIN_END'])['POS'].count()
all_sites_2=all_sites_2.reset_index()
all_sites_2=all_sites_2.rename(index=str, columns={'POS': 'N_sites' })




###########in process###########

#sites_filter=nextone[(nextone['irish_juvernica_N_sites'] >= 3000) & (nextone['kazak_juvernica_N_sites'] >= 3000) & (nextone['kaz_sin_N_sites'] >= 3000) & (nextone['spanish_reali_N_sites'] >= 3000) & (nextone['spanish_sinapis_N_sites'] >= 3000)& (nextone['swe_sin_allele_N_sites'] >= 3000)& (nextone['sinapis_N_sites'] >= 3000)& (nextone['reali_N_sites'] >= 3000)& (nextone['juvernica_N_sites'] >= 3000)]

##population script

#sites_col=args.sites
#freq_list=[args.POP1, args.POP2]
dxy_comb={str(i[0]+'_'+i[1]+'_dxy'): [pwd+(i[0]+'.frq'), (pwd+i[1]+'.frq')] for i in list(itertools.combinations(freq_list, 2))}
columns={str(i[0]+'_'+i[1]+'_dxy'): [(i[0]), (i[1])] for i in list(itertools.combinations(freq_list, 2))}
for i in dxy_comb.keys():
        filea=dxy_comb[i][0]
        fileb=dxy_comb[i][1]
        output_col=i
        fstcol=i[:-3]+'fst'
        fixedcol=i[:-3]+'fixed'
        privatea=i[:-3]+'private_a'
        privateb=i[:-3]+'private_b'
        sharedcol=i[:-3]+'shared'
        samecol=i[:-3]+'same_sites'
        Pi_a=i[:-3]+'Pi_a'
        Pi_b=i[:-3]+'Pi_b'
        primary_keys[str(output_col)]=0.0
        primary_keys[str(fstcol)]=0.0
        primary_keys[str(fixedcol)]=0.0
        primary_keys[str(privatea)]=0.0
        primary_keys[str(privateb)]=0.0
        primary_keys[str(sharedcol)]=0.0
        primary_keys[str(samecol)]=0.0
        primary_keys[str(Pi_a)]=0.0
        primary_keys[str(Pi_b)]=0.0


def main(filea,fileb,primary_keys, sites_col, part):
        primary_keys_1=np.array_split(primary_keys, 4)[part]
        for j in personal_popgen.dxy_window_function(filea,fileb,primary_keys_1, sites_col):
            #output_list [window_index, np.nanmean(avg_dxy), np.nanmean(avg_fst),fixed_diff, private_pop1, private_pop2, shared, fixed_same]
            primary_keys_1.at[j[0], str(output_col)] = j[1]
            primary_keys_1.at[j[0], str(fstcol)] = j[2]
            primary_keys_1.at[j[0], str(fixedcol)] = j[3]
            primary_keys_1.at[j[0], str(privatea)] = j[4]
            primary_keys_1.at[j[0], str(privateb)] = j[5]
            primary_keys_1.at[j[0], str(sharedcol)] = j[6]
            primary_keys_1.at[j[0], str(samecol)] = j[7]
            primary_keys_1.at[j[0], str(Pi_a)] = j[8]
            primary_keys_1.at[j[0], str(Pi_b)] = j[9]
            #print (j)
            #sys.stdout.flush()
        return primary_keys_1

#0      window_index
#1      np.nansum(avg_dxy)/(fixed_diff+private_pop1+private_pop2+shared)
#2      np.nanmean(avg_fst)
#3      fixed_diff
#4      private_pop1
#5      private_pop2
#6      shared
#7      fixed_same
#8      np.nansum(avg_pi_a)/(fixed_diff+private_pop1+private_pop2+shared)
#9      np.nansum(avg_pi_b)/(fixed_diff+private_pop1+private_pop2+shared)

#multiprocessing using 20 cores

primary_keys_a=all_sites_2

for i in dxy_comb:
    p = Pool(processes=4)
    filea=dxy_comb[i][0]
    fileb=dxy_comb[i][1]
    output_col=i
    fstcol=i[:-3]+'fst'
    fixedcol=i[:-3]+'fixed'
    privatea=i[:-3]+'private_a'
    privateb=i[:-3]+'private_b'
    sharedcol=i[:-3]+'shared'
    samecol=i[:-3]+'same_sites'
    Pi_a=i[:-3]+'Pi_a'
    Pi_b=i[:-3]+'Pi_b'
    primary_keys_a[str(output_col)]=0.0
    primary_keys_a[str(fstcol)]=0.0
    primary_keys_a[str(fixedcol)]=0.0
    primary_keys_a[str(privatea)]=0.0
    primary_keys_a[str(privateb)]=0.0
    primary_keys_a[str(sharedcol)]=0.0
    primary_keys_a[str(samecol)]=0.0
    primary_keys_a[str(Pi_a)]=0.0
    primary_keys_a[str(Pi_b)]=0.0
    data = p.starmap(main, zip([filea]*4,[fileb]*4,[primary_keys_a]*4, [sites_col]*4, range(4)))
    FINAL=pandas.concat(data)
    primary_keys_a=FINAL
    p.close()


FINAL.to_csv(i+'Fst_.csv')
#merging subprocess



