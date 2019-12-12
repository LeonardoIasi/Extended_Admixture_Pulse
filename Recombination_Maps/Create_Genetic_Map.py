import msprime
import pandas as pd
import numpy as np
from random import randint
from scipy import stats
import os
import os.path
from pathlib import Path
import matplotlib as mpl
### Not used
def Genetic_Map(Genetic_Map,length,start,Folder_name):
	Map_names=['position', 'COMBINED_rateGenetic_Map','Genetic_Map(cM)']
	Map=pd.read_table('{0}'.format(Genetic_Map),header=None,delim_whitespace=True,names=Map_names)
	Cleaved_Map=Map.loc[Map['position']>= start]
	Starting_Position=Cleaved_Map.iloc[0][0]
	Starting_distance=Cleaved_Map.iloc[0][2]
	Cleaved_Map.iloc[:, 0] = Cleaved_Map.iloc[:,0]- Starting_Position
	Cleaved_Map.iloc[:, 2] = Cleaved_Map.iloc[:,2]- Starting_distance
	Scaled_Map=Cleaved_Map.loc[Cleaved_Map['position'] <= float(length)]
	Scaled_Map.iloc[0][0]=1
	Last_entrie_Genetic_Map=(Scaled_Map['COMBINED_rate(cM/Mb)'].iloc[-1]*((float(length)-Scaled_Map['position'].iloc[-1])/1e06))+Scaled_Map['Genetic_Map(cM)'].iloc[-1]
	Last_entrie=pd.DataFrame([[float(length),0,Last_entrie_Genetic_Map]],columns=Map_names)
	Scaled_Map_2=Scaled_Map.append(Last_entrie,ignore_index=True)
	First_Col = pd.DataFrame(['chr1']*len(Scaled_Map_2))
	Scaled_Map_2['CHR']=First_Col
	cols=['CHR','position', 'COMBINED_rate(cM/Mb)','Genetic_Map(cM)']
	Scaled_Map_2 = Scaled_Map_2[cols]
	Scaled_Map_2.to_csv('{0}/Scaled_Genetic_AAMap.txt'.format(Folder_name),sep='\t',index=False)
### This is used to create Scaled_Genetic_AAMap.txt
def Genetic_Map(Genetic_Map_Name,start,length,Folder_name):
	Map=pd.read_csv('{0}'.format(Genetic_Map_Name),header=0,delim_whitespace=True)
	Rate = []
	for i in range(0,len(Map),1):
		if i == 0:
			Rate.append( (Map.iloc[i,1]-Map.iloc[i,1])/((Map.iloc[i,0])/1e6))
		else:
			Rate.append( (Map.iloc[i,1]-Map.iloc[(i-1),1])/((Map.iloc[i,0]-Map.iloc[(i-1),0])/1e6))
	Map_new=pd.DataFrame(Map[Map.columns[0]])
	Map_new['Rate(cM/Mb)']=Rate
	Map_new['Genetic_Map(cM)']=Map[Map.columns[1]]
	Map_names=['Physical_Position_Build36(hg18)', 'Rate(cM/Mb)','Genetic_Map(cM)']
	Cleaved_Map=Map_new.loc[Map_new['Physical_Position_Build36(hg18)']>= start]
	Starting_Position=Cleaved_Map.iloc[0][0]
	Starting_distance=Cleaved_Map.iloc[0][2]
	Cleaved_Map.iloc[:, 0] = Cleaved_Map.iloc[:,0]- Starting_Position
	Cleaved_Map.iloc[:, 2] = Cleaved_Map.iloc[:,2]- Starting_distance
	Scaled_Map=Cleaved_Map.loc[Cleaved_Map['Physical_Position_Build36(hg18)'] <= float(length)]
	Scaled_Map.iloc[0][0]=1
	Last_entrie_Genetic_Map=(Scaled_Map['Rate(cM/Mb)'].iloc[-1]*((float(length)-Scaled_Map['Physical_Position_Build36(hg18)'].iloc[-1])/1e06))+Scaled_Map['Genetic_Map(cM)'].iloc[-1]
	Last_entrie=pd.DataFrame([[float(length),0,Last_entrie_Genetic_Map]],columns=Map_names)
	Scaled_Map_2=Scaled_Map.append(Last_entrie,ignore_index=True)
	First_Col = pd.DataFrame(['chr1']*len(Scaled_Map_2))
	Scaled_Map_2['CHR']=First_Col
	cols=['CHR','Physical_Position_Build36(hg18)', 'Rate(cM/Mb)','Genetic_Map(cM)']
	Scaled_Map_2 = Scaled_Map_2[cols]
	Scaled_Map_2=Scaled_Map_2.rename(index=str, columns={"CHR": "Chromosome", "Physical_Position_Build36(hg18)": "Position(bp)", "Rate(cM/Mb)": "Rate(cM/Mb)", "Genetic_Map(cM)": "Map(cM)"})
	Scaled_Map_2.to_csv('{0}/Scaled_Genetic_AAMap.txt'.format(Folder_name),sep='\t',index=False)
### This is used to create a scaled map for each chromosome
def Genetic_Map(Genetic_Map_Name,start,stop,Folder_name,Chr):
	Map=pd.read_csv('{0}'.format(Genetic_Map_Name),header=0,delim_whitespace=True)
	Rate = []
	for i in range(0,len(Map),1):
		if i == 0:
			Rate.append( (Map.iloc[i,1]-Map.iloc[i,1])/((Map.iloc[i,0])/1e6))
		else:
			Rate.append( (Map.iloc[i,1]-Map.iloc[(i-1),1])/((Map.iloc[i,0]-Map.iloc[(i-1),0])/1e6))
	Map_new=pd.DataFrame(Map[Map.columns[0]])
	Map_new['Rate(cM/Mb)']=Rate
	Map_new['Genetic_Map(cM)']=Map[Map.columns[1]]
	Map_names=['Physical_Position_Build36(hg18)', 'Rate(cM/Mb)','Genetic_Map(cM)']
	Cleaved_Map=Map_new.loc[Map_new['Physical_Position_Build36(hg18)']>= start]
	Starting_Position=Cleaved_Map.iloc[0][0]
	Starting_distance=Cleaved_Map.iloc[0][2]
	Cleaved_Map.iloc[:, 0] = Cleaved_Map.iloc[:,0]- Starting_Position
	Cleaved_Map.iloc[:, 2] = Cleaved_Map.iloc[:,2]- Starting_distance
	Scaled_Map=Cleaved_Map.loc[Cleaved_Map['Physical_Position_Build36(hg18)'] <= (Cleaved_Map.iloc[-1,0]-stop)]
	Scaled_Map.iloc[0][0]=1
	Scaled_Map.iloc[-1,1]=0
	Scaled_Map_2=Scaled_Map.reset_index()  
	First_Col = pd.DataFrame(['chr{0}'.format(Chr)]*len(Scaled_Map_2))
	Scaled_Map_2['CHR']=First_Col
	cols=['CHR','Physical_Position_Build36(hg18)', 'Rate(cM/Mb)','Genetic_Map(cM)']
	Scaled_Map_2 = Scaled_Map_2[cols]
	Scaled_Map_2=Scaled_Map_2.rename(index=str, columns={"CHR": "Chromosome", "Physical_Position_Build36(hg18)": "Position(bp)", "Rate(cM/Mb)": "Rate(cM/Mb)", "Genetic_Map(cM)": "Map(cM)"})
	Scaled_Map_2.to_csv('{0}/Scaled_Genetic_AAMap_Chr{1}.txt'.format(Folder_name,Chr),sep='\t',index=False)

for i in range(1,21,1):
	Genetic_Map('Recombination_Maps/AAmap/AAmap.chr{0}.txt'.format(i),5e6,5e6,'Recombination_Maps',i)

