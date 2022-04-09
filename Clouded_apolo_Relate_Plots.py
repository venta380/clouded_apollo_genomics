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


Relate=pandas.read_csv('Relate.csv')
POP1=pandas.read_csv('pop.group', sep=' ',names=['ID1','POP1','Col1'])
Relate=Relate.merge(POP1, on='ID1')
POP2=pandas.read_csv('pop.group', sep=' ',names=['ID2','POP2','Col2'])
Relate=Relate.merge(POP2, on='ID2')

Different=Relate[Relate.POP1 != Relate.POP2]
Same=Relate[Relate.POP1 == Relate.POP2]



Uppland=["S121","S122","S123","S124","S125","S126","S127","S128","S129","S130"]
Blekinge=["S131","S132","S133","S134","S135","S136","S137","S138","S139","S140"]
NordensArk=["P21213_121_S141","S142","S143","S144","S145","S146","S147","S148","S149"]
Västernorrland=["S151","S152","S153","S154","S155","S156","S157","S158","S159"]

Uppland_DF=Different[(Different.ID1.isin(Uppland)) | (Different.ID2.isin(Uppland))]
Blekinge_DF=Different[(Different.ID1.isin(Blekinge)) | (Different.ID2.isin(Blekinge))]
NordensArk_DF=Different[(Different.ID1.isin(NordensArk)) | (Different.ID2.isin(NordensArk))]
Västernorrland_DF=Different[(Different.ID1.isin(Västernorrland)) | (Different.ID2.isin(Västernorrland))]

Uppland_DF['Main_Pop']='Uppland'
Blekinge_DF['Main_Pop']='Blekinge'
NordensArk_DF['Main_Pop']='NordensArk'
Västernorrland_DF['Main_Pop']='Västernorrland'


Uppland_DF.loc[Uppland_DF.POP1 != 'Uppland', 'Other_POP']=Uppland_DF[(Uppland_DF.POP1 != 'Uppland')]['POP1']
Blekinge_DF.loc[Blekinge_DF.POP1 != 'Blekinge', 'Other_POP']=Blekinge_DF[(Blekinge_DF.POP1 != 'Blekinge')]['POP1']
NordensArk_DF.loc[NordensArk_DF.POP1 != 'NordensArk', 'Other_POP']=NordensArk_DF[(NordensArk_DF.POP1 != 'NordensArk')]['POP1']
Västernorrland_DF.loc[Västernorrland_DF.POP1 != 'Västernorrland', 'Other_POP']=Västernorrland_DF[(Västernorrland_DF.POP1 != 'Västernorrland')]['POP1']

Uppland_DF.loc[Uppland_DF.POP1 == 'Uppland', 'Other_POP']=Uppland_DF[(Uppland_DF.POP1 == 'Uppland')]['POP2']
Blekinge_DF.loc[Blekinge_DF.POP1 == 'Blekinge', 'Other_POP']=Blekinge_DF[(Blekinge_DF.POP1 == 'Blekinge')]['POP2']
NordensArk_DF.loc[NordensArk_DF.POP1 == 'NordensArk', 'Other_POP']=NordensArk_DF[(NordensArk_DF.POP1 == 'NordensArk')]['POP2']
Västernorrland_DF.loc[Västernorrland_DF.POP1 == 'Västernorrland', 'Other_POP']=Västernorrland_DF[(Västernorrland_DF.POP1 == 'Västernorrland')]['POP2']


Uppland_DF['Main_Pop']='Uppland'
Blekinge_DF['Main_Pop']='Blekinge'
NordensArk_DF['Main_Pop']='NordensArk'
Västernorrland_DF['Main_Pop']='Västernorrland'


all_DF=pandas.concat([Uppland_DF,Blekinge_DF,NordensArk_DF,Västernorrland_DF],ignore_index=True)
all_DF['Main_Pop']=all_DF['Main_Pop']+'_'+all_DF['Other_POP']


sns.scatterplot(data=Uppland_DF, x="k0", y="k1", hue="Other_POP", palette=sns.color_palette(['#D35400','#F4D03F', '#27AE60']), hue_order=['Blekinge','NordensArk','Västernorrland'])
plt.ylim(0.0, 1.0)
plt.xlim(0.0, 1.0)
plt.title('Uppland')
plt.savefig('Uppland_Relate.png', dpi=1200)
plt.close()

sns.scatterplot(data=Blekinge_DF, x="k0", y="k1", hue="Other_POP", palette=sns.color_palette(['#2C3E50','#F4D03F', '#27AE60']), hue_order=['Uppland','NordensArk','Västernorrland'])
plt.ylim(0.0, 1.0)
plt.xlim(0.0, 1.0)
plt.title('Blekinge')
plt.savefig('Blekinge_Relate.png', dpi=1200)
plt.close()

sns.scatterplot(data=NordensArk_DF, x="k0", y="k1", hue="Other_POP", palette=sns.color_palette(['#2C3E50','#D35400','#27AE60']), hue_order=['Uppland','Blekinge','Västernorrland'])
plt.ylim(0.0, 1.0)
plt.xlim(0.0, 1.0)
plt.title('NordensArk')
plt.savefig('NordensArk_Relate.png', dpi=1200)
plt.close()

sns.scatterplot(data=Västernorrland_DF, x="k0", y="k1", hue="Other_POP", palette=sns.color_palette(['#2C3E50','#F4D03F','#D35400']), hue_order=['Uppland','NordensArk','Blekinge'])
plt.ylim(0.0, 1.0)
plt.xlim(0.0, 1.0)
plt.title('Västernorrland')
plt.savefig('Västernorrland_Relate.png', dpi=1200)
plt.close()



sns.boxplot(y="Main_Pop", x="kinship", data=all_DF[all_DF.Main_Pop.isin(['Blekinge_NordensArk','Uppland_Blekinge','Blekinge_Västernorrland','Uppland_NordensArk','NordensArk_Västernorrland','Uppland_Västernorrland'])], order=['Blekinge_NordensArk','Uppland_Blekinge','Blekinge_Västernorrland','Uppland_NordensArk','NordensArk_Västernorrland','Uppland_Västernorrland'])



sns.boxplot(y="POP2", x="kinship", data=Same, palette=sns.color_palette(['#D35400','#F4D03F','#2C3E50','#27AE60']), order=['Blekinge','NordensArk','Uppland','Västernorrland'])



 


