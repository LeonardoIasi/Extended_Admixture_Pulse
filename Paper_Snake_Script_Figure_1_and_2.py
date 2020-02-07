configfile:	"config/control.yaml"
configfile:	"config/master.yaml"
configfile:	"config/sim_parameters_new.yaml"
configfile:	"config/ALDER.yaml"
import msprime
import pandas as pd
import numpy as np
from random import randint
from scipy import stats
import os
import os.path
from pathlib import Path
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
import matplotlib.pyplot as plt

def load_config_master(config,wildcard):
	master = config['_Default_Scenario_'].copy()
	master.update(config[wildcard])
	return master

def load_config_sim_parameters(config,master):
	params=config['_MS_default_'].copy()
	params.update(config[master])
	return params

def load_ALDER_parameters(config,master):
	alder =config['_ALDER_default_'].copy()
	alder.update(config[master])
	return alder



def create_parfile_for_ALDER(alder,master_name,number,Folder_name,GF_Model,Ascertainment):
	file=open("{0}/parfile_alder-{1}-run{2}-{3}-ascertainment-{4}.txt".format(Folder_name,master_name,number,GF_Model,Ascertainment),"w")
	file.write("""genotypename:\t{0}/Ascertained_Simulation_merged-{1}-run{2}-{3}-ascertainment-{4}.eigenstratgeno
snpname:\t{0}/Simulation_Genetic_Distance_set-{1}-run{2}-{3}-ascertainment-{4}.snp
indivname:\t{0}/Simulation-{1}-Chr0-run{2}-{3}-ascertainment-{4}.ind
admixpop:\t{5}
refpops:\t{6}
checkmap:\t{7}
binsize:\t{8}
mindis:\t{9}
maxdis:\t{10}
fast_snp_read:\t{11}
num_threads:\t{12}
raw_outname:\t{0}/Raw_ALDER_output-{1}-run{2}-{3}-ascertainment-{4}.txt""".format(Folder_name,master_name,number,GF_Model,Ascertainment,str(alder['admixpop']),str(alder['refpops']),str(alder['checkmap']),float(alder['binsize']),float(alder['mindis']),float(alder['maxdis']),str(alder['fast_snp_read']),int(alder['num_threads'])))
	file.close()

def GF_Model_III(a,d,GF_start,demographic_events):
	# 1 minus is for changing the direction of the decline
	m=[1-x/d for x in range(d+1)]
	#m=[x/d for x in range(d+1)]
	m_sum=sum(m)
	m= [i*a/m_sum for i in m]
	for x in range(d+1):
		Gene_Flow_Rate_Change=msprime.MigrationRateChange(time=(GF_start+x),rate= m[x],matrix_index=(1,2) )
		demographic_events.append(Gene_Flow_Rate_Change)
	Gene_Flow_end=msprime.MigrationRateChange(time=(GF_start+(d+1)),rate= 0,matrix_index=(1,2) )
	demographic_events.append(Gene_Flow_end)
	return demographic_events

def GF_Model_IV(TimeSpan,GF_start,GF_stop,demographic_events,GF_rate,loc,Folder_name,master_name,number,GF_Model):
	check_m2_1=list()
	check_m2_2=list()
	event_length=list()
	event_length_2=list([0])
	EX= ((GF_stop - GF_start)/2 +GF_start)
	VarX= ((GF_stop - GF_start)/4)**2
	b= EX/VarX
	a=EX*b
	m=[stats.gamma.pdf(x=range(TimeSpan+1),a=a+1,loc=loc,scale=1/b)]
	m=m[0]
	m[abs(m) < 1e-6] = 0
	m2 = [i * GF_rate/ sum(m) for i in m]
	#fig = plt.figure(figsize=(6, 5))
	#plt.plot(m2)
	#plt.xlabel('Generations')
	#plt.title("shape={0}, 1/scale={1}".format(a,b))
	#plt.savefig('{0}/Gene_Flow_Model-{1}-{2}.png'.format(Folder_name,master_name,GF_Model))
	#plt.close(fig)

	for x in range(TimeSpan+1):
		if x >= GF_start and x <= GF_stop:
			Gene_Flow_Rate_Change=msprime.MigrationRateChange(time=(x),rate= m2[x],matrix_index=(1,2) )
			demographic_events.append(Gene_Flow_Rate_Change)
			event_length.append(x)
			check_m2_1.append(m2[x])
	for x in range(TimeSpan+1):
		if x < GF_start or x > GF_stop:
			if m2[x] > 1e-5:
				Gene_Flow_Rate_Change_2=msprime.MigrationRateChange(time=(x),rate= m2[x],matrix_index=(1,2) )
				demographic_events.append(Gene_Flow_Rate_Change_2)
				event_length_2.append(x)
				check_m2_2.append(m2[x])
	if event_length[-1] > event_length_2[-1]:
		Neandertal_Gene_Flow_end=msprime.MigrationRateChange(time=(event_length[-1]+1),rate=0)
		demographic_events.append(Neandertal_Gene_Flow_end)
	else:
		Neandertal_Gene_Flow_end=msprime.MigrationRateChange(time=(event_length_2[-1]+1),rate=0)
		demographic_events.append(Neandertal_Gene_Flow_end)

	return demographic_events

def simulation(params,Folder_name,master_name,number, master,GF_Model):

	samples = [msprime.Sample(population=0, time=int(params['sample_Africans_at']))] * int(params['sample_Africans']) + [msprime.Sample(population=1, time=int(params['sample_non_Africans_at']))]*int(params['sample_non_Africans']) + [msprime.Sample(population=2, time=int(params['sample_Neandertals_at']))]*int(params['sample_Neandertals'])


	population_configuration=[
		msprime.PopulationConfiguration(
			initial_size=int(params['N_Africans'])),
		msprime.PopulationConfiguration(
			initial_size=int(params['N_non_Africans'])),
		msprime.PopulationConfiguration(
			initial_size=int(params['N_Neandertals']))
		]


	if str(master['GF_Model'][0]) == 'GF_Model_IV':
		Neandertal_Gene_Flow_absolute_start=msprime.MigrationRateChange(time=int(params['sample_Neandertals']),rate=0)
		Neandertal_Gene_Flow_absolute_end=msprime.MigrationRateChange(time=int(params['split_time_non_Africans']),rate=0)
		EUR_AFR_Gene_Flow_start=msprime.MigrationRateChange(time=int(params['EUR_AFR_migration_start']),rate=float(params['EUR_AFR_migration_rate_GF'])/(int(params['EUR_AFR_migration_end'])-int(params['EUR_AFR_migration_start'])),matrix_index=(0,1) )
		AFR_EUR_Gene_Flow_start=msprime.MigrationRateChange(time=int(params['EUR_AFR_migration_start']),rate=float(params['EUR_AFR_migration_rate_GF'])/(int(params['EUR_AFR_migration_end'])-int(params['EUR_AFR_migration_start'])),matrix_index=(1,0) )
		EUR_AFR_Gene_Flow_end=msprime.MigrationRateChange(time=int(params['EUR_AFR_migration_end']),rate=0,matrix_index=(0,1) )
		AFR_EUR_Gene_Flow_end=msprime.MigrationRateChange(time=int(params['EUR_AFR_migration_end']),rate=0,matrix_index=(1,0) )
		Split_Time_Neandertals=msprime.MassMigration(time=int(params['split_time_Neandertals']),source=0,destination=2,proportion=1.0)
		Split_Time_non_Africans=msprime.MassMigration(time=int(params['split_time_non_Africans']),source=1,destination=0,proportion=1.0)

		demographic_events=[EUR_AFR_Gene_Flow_start,AFR_EUR_Gene_Flow_start,EUR_AFR_Gene_Flow_end,AFR_EUR_Gene_Flow_end,Neandertal_Gene_Flow_absolute_end,Neandertal_Gene_Flow_absolute_start,Split_Time_Neandertals,Split_Time_non_Africans]

		demographic_events=GF_Model_IV(TimeSpan=int(params['split_time_non_Africans']),GF_start=int(params['GF_start']),GF_stop=int(params['GF_stop']),demographic_events=demographic_events,GF_rate=float(params['migration_rate_GF']),loc=int(master['GF_Model'][1]),Folder_name=Folder_name,master_name=master_name,number=number,GF_Model=GF_Model)


	if str(master['Inferred_Demographies'][0]) == 'True':
		Initial_Neandertal_Size=msprime.PopulationParametersChange(time=int(params['split_time_Neandertals']),initial_size=int(params['N_Africans']), growth_rate=0,population_id=2)
		Initial_YRI_Size=msprime.PopulationParametersChange(time=int(params['split_time_Neandertals']),initial_size=int(params['N_Africans']), growth_rate=0,population_id=0)
		Initial_CEU_Size=msprime.PopulationParametersChange(time=int(params['split_time_Neandertals']),initial_size=int(params['N_Africans']), growth_rate=0,population_id=1)
		demographic_events.append(Initial_Neandertal_Size)
		demographic_events.append(Initial_YRI_Size)
		demographic_events.append(Initial_CEU_Size)

		Vindjia=pd.DataFrame(data=pd.read_csv("{0}".format(str(master['Inferred_Demographies'][1])), sep= " ",index_col=0))
		for index,row in Vindjia.iterrows():
			if float(row[0]) <= float(int(params['split_time_Neandertals'])*float(params['mutation_rate'])):
				Neandertal_Size_Change=msprime.PopulationParametersChange(time=(int(row[0]/float(params['mutation_rate']))),initial_size=(int(row[1]/float(params['mutation_rate']))), growth_rate=0,population_id=2)
				demographic_events.append(Neandertal_Size_Change)
		YRI=pd.DataFrame(data=pd.read_csv("{0}".format(str(master['Inferred_Demographies'][2])), sep= " ",index_col=0))
		for index,row in YRI.iterrows():
			if float(row[0]) <= float(int(params['split_time_Neandertals'])*float(params['mutation_rate'])):
				YRI_Size_Change=msprime.PopulationParametersChange(time=(int(row[0]/float(params['mutation_rate']))),initial_size=(int(row[1]/float(params['mutation_rate']))), growth_rate=0,population_id=0)
				demographic_events.append(YRI_Size_Change)
		CEU=pd.DataFrame(data=pd.read_csv("{0}".format(str(master['Inferred_Demographies'][3])), sep= " ",index_col=0))
		for index,row in CEU.iterrows():
			if float(row[0]) <= float(int(params['split_time_non_Africans'])*float(params['mutation_rate'])):
				CEU_Size_Change=msprime.PopulationParametersChange(time=(int(row[0]/float(params['mutation_rate']))),initial_size=(int(row[1]/float(params['mutation_rate']))), growth_rate=0,population_id=1)
				demographic_events.append(CEU_Size_Change)

	if str(master['Complex_Demography'][0]) == 'True':
		Ancient_Substructure_end=msprime.MigrationRateChange(time=int(params['Substructure_end']),rate=float(params['Substructure_mig_rate']),matrix_index=(0,1) )
		Ancient_Substructure_start=msprime.MigrationRateChange(time=int(params['Substructure_mig_start']),rate=0,matrix_index=(0,1) )
		Ancient_Substructure_end2=msprime.MigrationRateChange(time=int(params['Substructure_end']),rate=float(params['Substructure_mig_rate']),matrix_index=(1,0) )
		Ancient_Substructure_start2=msprime.MigrationRateChange(time=int(params['Substructure_mig_start']),rate=0,matrix_index=(1,0) )
		Initial_Neandertal_Size=msprime.PopulationParametersChange(time=int(params['split_time_Neandertals']),initial_size=int(params['N_Africans']), growth_rate=0,population_id=2)
		Initial_YRI_Size=msprime.PopulationParametersChange(time=int(params['split_time_Neandertals']),initial_size=int(params['N_Africans']), growth_rate=0,population_id=0)
		Initial_CEU_Size=msprime.PopulationParametersChange(time=int(params['split_time_Neandertals']),initial_size=int(params['N_Africans']), growth_rate=0,population_id=1)
		demographic_events.append(Ancient_Substructure_end)
		demographic_events.append(Ancient_Substructure_start)
		demographic_events.append(Ancient_Substructure_end2)
		demographic_events.append(Ancient_Substructure_start2)
		demographic_events.append(Initial_Neandertal_Size)
		demographic_events.append(Initial_YRI_Size)
		demographic_events.append(Initial_CEU_Size)

		Vindjia=pd.DataFrame(data=pd.read_csv("{0}".format(str(master['Complex_Demography'][1])), sep= " ",index_col=0))
		for index,row in Vindjia.iterrows():
			if float(row[0]) <= float(int(params['split_time_Neandertals'])*float(params['mutation_rate'])):
				Neandertal_Size_Change=msprime.PopulationParametersChange(time=(int(row[0]/float(params['mutation_rate']))),initial_size=(int(row[1]/float(params['mutation_rate']))), growth_rate=0,population_id=2)
				demographic_events.append(Neandertal_Size_Change)
		YRI=pd.DataFrame(data=pd.read_csv("{0}".format(str(master['Complex_Demography'][2])), sep= " ",index_col=0))
		for index,row in YRI.iterrows():
			if float(row[0]) <= float(int(params['split_time_Neandertals'])*float(params['mutation_rate'])):
				YRI_Size_Change=msprime.PopulationParametersChange(time=(int(row[0]/float(params['mutation_rate']))),initial_size=(int(row[1]/float(params['mutation_rate']))), growth_rate=0,population_id=0)
				demographic_events.append(YRI_Size_Change)
		CEU=pd.DataFrame(data=pd.read_csv("{0}".format(str(master['Complex_Demography'][3])), sep= " ",index_col=0))
		for index,row in CEU.iterrows():
			if float(row[0]) <= float(int(params['split_time_non_Africans'])*float(params['mutation_rate'])):
				CEU_Size_Change=msprime.PopulationParametersChange(time=(int(row[0]/float(params['mutation_rate']))),initial_size=(int(row[1]/float(params['mutation_rate']))), growth_rate=0,population_id=1)
				demographic_events.append(CEU_Size_Change)

	demographic_events.sort(key = lambda x: x.time)
	migration_matrix=[[0,0,0],
			[0,0,0],
			[0,0,0]]

	check=msprime.DemographyDebugger(
	population_configurations=population_configuration,
	migration_matrix=migration_matrix,
	demographic_events=demographic_events)
	check.print_history(output=open('{0}/Demography_history-{1}-{2}.txt'.format(Folder_name,master_name,GF_Model), 'w'))

	if str(master['Recombination_Map'][0]) == 'False':
		simulation=msprime.simulate(
			samples=samples,
			mutation_rate=float(params['mutation_rate']),
			length=float(params['chr_length']),
			recombination_rate=float(params['recombination_rate']),
			population_configurations=population_configuration,
			migration_matrix=migration_matrix,
			demographic_events=demographic_events
			)
	elif str(master['Recombination_Map'][0]) == 'True':
		simulation=msprime.simulate(
			samples=samples,
			mutation_rate=float(params['mutation_rate']),
			recombination_map=msprime.RecombinationMap.read_hapmap('{0}'.format(str(master['Recombination_Map'][1]))),
			population_configurations=population_configuration,
			migration_matrix=migration_matrix,
			demographic_events=demographic_events)



	return simulation




############################################################################################
# Choose which Set to process by giving Set_name the name of the Set to be processed #
Set_name='Fig_2_D_Complex'

# Choose number of replicates #
replicates=100

# Choose Folder #
Folder_name='../Fig_2_D_Complex'

# Choose Result Folder Name #
Result_Folder='../Fig_2_D_Complex/Result'
############################################################################################


master_name=config['CONTROL']['%s' % Set_name]
master_Scenario= list(master_name.keys())
Sim_name=config['MASTER']['%s' % master_Scenario[0]]['simulation']


rule alle:
	input:
		#expand("{Folder_name}/Ascertained_Simulation_merged-{master_name}-run{number}-{GF_Model}-ascertainment-{Ascertainment}.eigenstratgeno",master_name=config['CONTROL']['%s' % Set_name],number=range(replicates),Set_name=Set_name,GF_Model=config['MASTER']['%s' % master_Scenario[0]]['GF_Model'][0],Result_Folder=Result_Folder,Folder_name=Folder_name,Ascertainment=config['MASTER']['%s' % master_Scenario[0]]['ascertainment']),
		#expand("{Folder_name}/Simulation_Genetic_Distance_set-{master_name}-run{number}-{GF_Model}-ascertainment-{Ascertainment}.snp",master_name=config['CONTROL']['%s' % Set_name],number=range(replicates),Set_name=Set_name,GF_Model=config['MASTER']['%s' % master_Scenario[0]]['GF_Model'][0],Result_Folder=Result_Folder,Folder_name=Folder_name,Ascertainment=config['MASTER']['%s' % master_Scenario[0]]['ascertainment']),
		#expand("{Folder_name}/Simulation-{master_name}-Chr0-run{number}-{GF_Model}-ascertainment-{Ascertainment}.ind",master_name=config['CONTROL']['%s' % Set_name],number=range(replicates),Set_name=Set_name,GF_Model=config['MASTER']['%s' % master_Scenario[0]]['GF_Model'][0],Result_Folder=Result_Folder,Folder_name=Folder_name,Ascertainment=config['MASTER']['%s' % master_Scenario[0]]['ascertainment']),
		#expand("{Folder_name}/Simulation-{master_name}-Chr{nChr}-run{number}-{GF_Model}-ascertainment-{Ascertainment}.vcf",master_name=config['CONTROL']['%s' % Set_name],number=range(replicates),Set_name=Set_name,GF_Model=config['MASTER']['%s' % master_Scenario[0]]['GF_Model'][0],Folder_name=Folder_name,nChr=range(config['sim_parameters']['_MS_default_']['num_chr']),Ascertainment=config['MASTER']['%s' % master_Scenario[0]]['ascertainment'])
		expand("{Result_Folder}/Result_file_SIM_Raw_ALDER-Fit-{Set_name}-{GF_Model}-min_dist_Fit-{min_dist_Fit}-ascertainment-{Ascertainment}.txt",master_name=config['CONTROL']['%s' % Set_name],number=range(replicates),Set_name=Set_name,GF_Model=config['MASTER']['%s' % master_Scenario[0]]['GF_Model'][0],Result_Folder=Result_Folder,Folder_name=Folder_name,min_dist_Fit=config['MASTER']['%s' % master_Scenario[0]]['min_dist_Fit'],Ascertainment=config['MASTER']['%s' % master_Scenario[0]]['ascertainment']),

rule basic_simulation:
	input:
	output:
		temp("{Folder_name}/Simulation-{master_name}-Chr{nChr}-run{number}-{GF_Model}-ascertainment-{Ascertainment}.eigenstratgeno"),
		temp("{Folder_name}/Simulation-{master_name}-Chr{nChr}-run{number}-{GF_Model}-ascertainment-{Ascertainment}.snp"),
		temp("{Folder_name}/Simulation-{master_name}-Chr{nChr}-run{number}-{GF_Model}-ascertainment-{Ascertainment}.ind"),
		#"{Folder_name}/Simulation-{master_name}-Chr{nChr}-run{number}-{GF_Model}-ascertainment-{Ascertainment}.vcf"
	run:
		master=load_config_master(config['MASTER'],wildcards.master_name)
		params=load_config_sim_parameters(config['sim_parameters'],master['simulation'])
		simulation_results=simulation(params=params,Folder_name=wildcards.Folder_name,master_name=wildcards.master_name,number=wildcards.number,master=master,GF_Model=wildcards.GF_Model);
		#with open('{}'.format(output[3]),'w') as vcf_file:
		#	simulation_results.write_vcf(vcf_file,2,str(int(wildcards.nChr)+1))
		AFR = simulation_results.get_samples(0)
		EUR = simulation_results.get_samples(1)
		NEA = simulation_results.get_samples(2)
		# Write geno and snp files
		with open('{}'.format(output[0]),'w') as outgeno , open('{}'.format(output[1]),'w') as snpfile:
			for variant in simulation_results.variants():
				AFR_geno = ''
				EUR_geno = ''
				NEA_geno = ''
				# This line ascertaines the snps !!!
				if master['ascertainment'] == int(0):
					if np.sum(variant.genotypes[AFR]) == 0 and np.sum(variant.genotypes[NEA]) != 0:
						pos = int(variant.site.position)
						if variant.site.id == 0:
							pos_genetic = 0
						else:
							pos_genetic = pos*float(params['recombination_rate'])
						snpfile.write('rs{0}\t{1}\t{2}\t{3}\n'.format(pos,str(int(wildcards.nChr)+1),pos_genetic,pos))
						for i in range(0,len(AFR),2):
							AFR_geno += str(np.sum(variant.genotypes[AFR][i:i+2]))
						for i in range(0,len(EUR),2):
							EUR_geno += str(np.sum(variant.genotypes[EUR][i:i+2]))
						for i in range(0,len(NEA),2):
							NEA_geno += str(np.sum(variant.genotypes[NEA][i:i+2]))
						outgeno.write('{0}{1}{2}\n'.format(AFR_geno,EUR_geno,NEA_geno))
				if master['ascertainment'] == int(1):
					if np.sum(variant.genotypes[AFR]) == 0 and np.sum(variant.genotypes[NEA]) != 0 and np.sum(variant.genotypes[EUR]) != 0:
						pos = int(variant.site.position)
						if variant.site.id == 0:
							pos_genetic = 0
						else:
							pos_genetic = pos*float(params['recombination_rate'])
						snpfile.write('rs{0}\t{1}\t{2}\t{3}\n'.format(pos,str(int(wildcards.nChr)+1),pos_genetic,pos))
						for i in range(0,len(AFR),2):
							AFR_geno += str(np.sum(variant.genotypes[AFR][i:i+2]))
						for i in range(0,len(EUR),2):
							EUR_geno += str(np.sum(variant.genotypes[EUR][i:i+2]))
						for i in range(0,len(NEA),2):
							NEA_geno += str(np.sum(variant.genotypes[NEA][i:i+2]))
						outgeno.write('{0}{1}{2}\n'.format(AFR_geno,EUR_geno,NEA_geno))
		# Write the ind file
		with open('{}'.format(output[2]),'w') as outind:
			for i in range(0,int(len(AFR)/2)):
				outind.write('Ind_{0}\tU\tAfricans\n'.format(i))
			for i in range(0,int(len(EUR)/2)):
				outind.write('Ind_{0}\tU\tNon_Africans\n'.format(i))
			for i in range(0,int(len(NEA)/2)):
				outind.write('Ind_{0}\tU\tNeandertals\n'.format(i))
			outind.close()


def merge_geno_files(wildcards):
	files = expand("{Folder_name}/Simulation-{master_name}-Chr{nChr}-run{number}-{GF_Model}-ascertainment-{Ascertainment}.eigenstratgeno",master_name=wildcards.master_name,number=wildcards.number,Set_name=Set_name,GF_Model=config['MASTER']['%s' % master_Scenario[0]]['GF_Model'][0],Folder_name=Folder_name,nChr=range(config['sim_parameters']['_MS_default_']['num_chr']),Ascertainment=config['MASTER']['%s' % master_Scenario[0]]['ascertainment'])
	return files

rule merge_geno_files:
	input:
		merge_geno_files
	output:
		temp("{Folder_name}/Ascertained_Simulation_merged-{master_name}-run{number}-{GF_Model}-ascertainment-{Ascertainment}.eigenstratgeno")
	run:
		with open('{}'.format(output[0]),'w') as out_merged_geno:
			for i, file in enumerate(input):
				with open(file) as data:
					for line in data:
						out_merged_geno.write(line)

def merge_snp_files(wildcards):
	files = expand("{Folder_name}/Simulation-{master_name}-Chr{nChr}-run{number}-{GF_Model}-ascertainment-{Ascertainment}.snp",master_name=wildcards.master_name,number=wildcards.number,Set_name=Set_name,GF_Model=config['MASTER']['%s' % master_Scenario[0]]['GF_Model'][0],Folder_name=Folder_name,nChr=range(config['sim_parameters']['_MS_default_']['num_chr']),Ascertainment=config['MASTER']['%s' % master_Scenario[0]]['ascertainment'])
	return files

rule merge_snp_files:
	input:
		merge_snp_files
	output:
		temp("{Folder_name}/Ascertained_Simulation_merged-{master_name}-run{number}-{GF_Model}-ascertainment-{Ascertainment}.snp")
	run:
		with open('{}'.format(output[0]),'w') as out_merged_geno:
			for i, file in enumerate(input):
				with open(file) as data:
					for line in data:
						out_merged_geno.write(line)

rule Set_Genetic_Distance:
	input:
		"{Folder_name}/Ascertained_Simulation_merged-{master_name}-run{number}-{GF_Model}-ascertainment-{Ascertainment}.snp"
	output:
		temp("{Folder_name}/Simulation_Genetic_Distance_set-{master_name}-run{number}-{GF_Model}-ascertainment-{Ascertainment}.snp")
	run:
		master=load_config_master(config['MASTER'],wildcards.master_name)
		params=load_config_sim_parameters(config['sim_parameters'],master['simulation'])
		if str(master['Recombination_Map'][0]) == 'False':
			x=float(params['recombination_rate'])
			shell("""awk '{{print $1,$2, $3=$4*{x}, $4}}' {input} | sed '{{s/ /\t/g}}' > {output} """)
		elif str(master['Recombination_Map'][0]) == 'True':
			x=float(1.843e-8)
			shell("""awk '{{print $1,$2, $3=$4*{x}, $4}}' {input} | sed '{{s/ /\t/g}}' > {output} """)

rule run_ALDER:
	input:
		"{Folder_name}/Ascertained_Simulation_merged-{master_name}-run{number}-{GF_Model}-ascertainment-{Ascertainment}.eigenstratgeno",
		"{Folder_name}/Simulation_Genetic_Distance_set-{master_name}-run{number}-{GF_Model}-ascertainment-{Ascertainment}.snp",
		"{Folder_name}/Simulation-{master_name}-Chr0-run{number}-{GF_Model}-ascertainment-{Ascertainment}.ind"


	output:
		"{Folder_name}/Raw_ALDER_output-{master_name}-run{number}-{GF_Model}-ascertainment-{Ascertainment}.txt"

	threads: 8
	run:
		master=load_config_master(config['MASTER'],wildcards.master_name)
		alder=load_ALDER_parameters(config['ALDER'],master['alder'])
		create_parfile_for_ALDER(alder=alder, master_name=wildcards.master_name,number=wildcards.number,Folder_name=wildcards.Folder_name,GF_Model=wildcards.GF_Model,Ascertainment=wildcards.Ascertainment)
		shell("""set +e
~leonardo_iasi/alder/alder -p {wildcards.Folder_name}/parfile_alder-{wildcards.master_name}-run{wildcards.number}-{wildcards.GF_Model}-ascertainment-{wildcards.Ascertainment}.txt
exitcode=$?
if [ $exitcode -eq 1 ]
then
		exit 0
else
	exit 0
fi""")

rule fit_Curve:
	input:
		"{Folder_name}/Raw_ALDER_output-{master_name}-run{number}-{GF_Model}-ascertainment-{Ascertainment}.txt"
	output:
		"{Folder_name}/Model_Fit/Summary-Fit-{master_name}-run{number}-{GF_Model}-min_dist_Fit-{min_dist_Fit}-ascertainment-{Ascertainment}.log",
		#"{Folder_name}/Model_Fit/Plot-Fit-{master_name}-run{number}-{GF_Model}-min_dist_Fit-{min_dist_Fit}-ascertainment-{Ascertainment}.pdf",


	run:
		master=load_config_master(config['MASTER'],wildcards.master_name)
		params=load_config_sim_parameters(config['sim_parameters'],master['simulation'])
		#shell("""Rscript /r1/people/leonardo_iasi/Desktop/Neandertal_Human_Introgression_Project/Paper/Paper_Scripts/Fit_Exponential.R {input[0]} {output[0]} {output[1]} {wildcards.min_dist_Fit}""")
		shell("""Rscript /r1/people/leonardo_iasi/Desktop/Neandertal_Human_Introgression_Project/Paper/Paper_Scripts/Fit_Exponential.R {input[0]} {output[0]} {output[0]} {wildcards.min_dist_Fit}""")

rule grep_results_and_merge_with_Scenarios:
	input:
		"{Folder_name}/Model_Fit/Summary-Fit-{master_name}-run{number}-{GF_Model}-min_dist_Fit-{min_dist_Fit}-ascertainment-{Ascertainment}.log"
	output:
		temp("{Folder_name}/Result-Fit_Output-{master_name}-run{number}-{GF_Model}-min_dist_Fit-{min_dist_Fit}-ascertainment-{Ascertainment}.txt")
	run:
		master=load_config_master(config['MASTER'],wildcards.master_name)
		params=load_config_sim_parameters(config['sim_parameters'],master['simulation'])
		shell("""grep "A, m, c, RSS_Expo, nls status:" {input} |  sed 's/ /\t/g' |
awk '{{print $7,$8,$9,$10,$11,$12="{wildcards.master_name}",$13="{params[GF_start]}",$14="{params[GF_stop]}",$15="{master[ascertainment]}",$16="{master[min_dist_Fit]}",$17="{master[GF_Model][0]}"}}' > {output} """)




rule merge_all_files:
	input:
		expand("{Folder_name}/Result-Fit_Output-{master_name}-run{number}-{GF_Model}-min_dist_Fit-{min_dist_Fit}-ascertainment-{Ascertainment}.txt",master_name=config['CONTROL']['%s' % Set_name],
number=range(replicates),Set_name=Set_name,Folder_name=Folder_name,GF_Model=config['MASTER']['%s' % master_Scenario[0]]['GF_Model'][0],min_dist_Fit=config['MASTER']['%s' % master_Scenario[0]]['min_dist_Fit'],
Ascertainment=config['MASTER']['%s' % master_Scenario[0]]['ascertainment'])
	output:
		"{Result_Folder}/Result_file_SIM_Raw_ALDER-Fit-{Set_name}-{GF_Model}-min_dist_Fit-{min_dist_Fit}-ascertainment-{Ascertainment}.txt"
	run:
		results=list()
		for i in input:
			with open(i) as line:
				lines= line.readlines()
			content= lines[0]
			results.append(content)
			file=open(output[0],'w')
			for n in results:
				file.write("%s" % n)
