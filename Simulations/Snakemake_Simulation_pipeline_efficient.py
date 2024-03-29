configfile:    "config/control_efficient.yaml"
configfile:    "config/sim_parameters_efficient.yaml"
configfile:    "config/master_efficient.yaml"
configfile:    "config/ALDER.yaml"
import sys
import pandas as pd
import msprime
import numpy as np
import random
import tskit
import gzip
from scipy import stats
import os
import os.path

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



def create_parfile_for_ALDER(alder,master_name,number,Folder_name,GF_Model,Ascertainment,downsample,Recomb_correction):
    file=open("{0}/parfile_alder-{1}-run{2}-{3}-ascertainment{4}-downsampled{5}-Recc{6}.txt".format(Folder_name,master_name,number,GF_Model,Ascertainment,downsample,Recomb_correction),"w")
    file.write("""genotypename:\t{0}/Ascertained_Simulation_merged-{1}-run{2}-{3}-ascertainment{4}-downsampled{5}.eigenstratgeno
snpname:\t{0}/Simulation_Genetic_Distance_set-{1}-run{2}-{3}-ascertainment{4}-downsampled{5}-Recc{6}.snp
indivname:\t{0}/Simulation-{1}-Chr0-run{2}-{3}-ascertainment{4}-downsampled{5}.ind
admixpop:\t{7}
refpops:\t{8}
checkmap:\t{9}
binsize:\t{10}
mindis:\t{11}
maxdis:\t{12}
fast_snp_read:\t{13}
num_threads:\t{14}
raw_outname:\t{0}/Raw_ALDER_output-{1}-run{2}-{3}-ascertainment{4}-downsampled{5}-Recc{6}.txt""".format(Folder_name,master_name,number,GF_Model,Ascertainment,downsample,Recomb_correction,str(alder['admixpop']),str(alder['refpops']),str(alder['checkmap']),float(alder['binsize']),float(alder['mindis']),float(alder['maxdis']),str(alder['fast_snp_read']),int(alder['num_threads'])))
    file.close()

def get_all_snps(ts,pop_names=['AFR','EUR','NEA']):
    sample_names = []
    for p in range(len(pop_names)):
        s_name=pop_names[p]
        n_pop = len(ts.get_samples(p))
        sample_names.extend([f"{s_name}{i}" for i in range(n_pop)])
    return pd.concat(
        pd.DataFrame(
            t.genotype_matrix(),
            columns=sample_names,
            index=(v.position for v in t.variants())
        ) for t in [ts]
    )

def write_all_snps(file_name, ts, pop_names=['AFR','EUR','NEA'], chrom='1'):
    """Write genotype matrix to file.

    More memory efficient than get_all_snps as it iterates over the
    tree sequence directly.

    """
    sample_names = []
    for p in range(len(pop_names)):
        s_name=pop_names[p]
        n_pop = len(ts.get_samples(p))
        sample_names.extend([f"{s_name}{i}" for i in range(n_pop)])
    with gzip.open(file_name, 'wt') as f:
        f.write("chrom\tpos\t")
        f.write("\t".join(sample_names))
        f.write("\n")
        for v in ts.variants():
            gt = '\t'.join(str(gt) for gt in v.genotypes)
            f.write(f"{chrom}\t{int(v.site.position)}\t{gt}\n")

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
    m=[stats.gamma.pdf(x=range(TimeSpan+1),a=a,loc=loc,scale=1/b)]
    m=m[0]
    m[abs(m) < 1e-6] = 0
    m2 = [i * GF_rate/ sum(m) for i in m]

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

def simulation(params,Folder_name,output_prefix,master_name,number, master,GF_Model):
    os.makedirs(os.path.dirname(output_prefix), exist_ok=True)
    samples = [msprime.Sample(population=0, time=int(params['sample_Africans_at']))] * int(params['sample_Africans']) + [msprime.Sample(population=1, time=int(params['sample_non_Africans_at']))]*int(params['sample_non_Africans']) + [msprime.Sample(population=2, time=int(params['sample_Neandertals_at']))]*int(params['sample_Neandertals'])


    population_configuration=[
        msprime.PopulationConfiguration(
            initial_size=int(params['N_Africans'])),
        msprime.PopulationConfiguration(
            initial_size=int(params['N_non_Africans'])),
        msprime.PopulationConfiguration(
            initial_size=int(params['N_Neandertals']))
        ]


    #if str(master['GF_Model'][0]) == 'GF_Model_IV':
    #change_sim_model=msprime.SimulationModelChange(time=int(params['sim_time_change']), model="hudson")
    Neandertal_Gene_Flow_absolute_start=msprime.MigrationRateChange(time=int(params['sample_Neandertals_at']),rate=0)
    Neandertal_Gene_Flow_absolute_end=msprime.MigrationRateChange(time=int(params['split_time_non_Africans']),rate=0)
    EUR_AFR_Gene_Flow_start=msprime.MigrationRateChange(time=params['EUR_AFR_migration_start'],rate=float(params['EUR_AFR_migration_rate_GF']/(params['EUR_AFR_migration_end']-params['EUR_AFR_migration_start'])),matrix_index=(0,1) )
    AFR_EUR_Gene_Flow_start=msprime.MigrationRateChange(time=params['EUR_AFR_migration_start'],rate=float(params['EUR_AFR_migration_rate_GF']/(params['EUR_AFR_migration_end']-params['EUR_AFR_migration_start'])),matrix_index=(1,0) )
    EUR_AFR_Gene_Flow_end=msprime.MigrationRateChange(time=int(params['EUR_AFR_migration_end']),rate=0,matrix_index=(0,1) )
    AFR_EUR_Gene_Flow_end=msprime.MigrationRateChange(time=int(params['EUR_AFR_migration_end']),rate=0,matrix_index=(1,0) )
    Substructure_start1=msprime.MigrationRateChange(time=params['Substructure_mig_start'],rate=float(params['Substructure_mig_rate']),matrix_index=(0,1) )
    Substructure_start2=msprime.MigrationRateChange(time=params['Substructure_mig_start'],rate=float(params['Substructure_mig_rate']),matrix_index=(1,0) )
    Substructure_end1=msprime.MigrationRateChange(time=params['Substructure_end'],rate=0,matrix_index=(0,1) )
    Substructure_end2=msprime.MigrationRateChange(time=params['Substructure_end'],rate=0,matrix_index=(1,0) )
    Split_Time_Neandertals=msprime.MassMigration(time=int(params['split_time_Neandertals']),source=0,destination=2,proportion=1.0)
    Split_Time_non_Africans=msprime.MassMigration(time=int(params['split_time_non_Africans']),source=1,destination=0,proportion=1.0)

    #demographic_events=[EUR_AFR_Gene_Flow_start,AFR_EUR_Gene_Flow_start,EUR_AFR_Gene_Flow_end,AFR_EUR_Gene_Flow_end,
    #Substructure_start1,Substructure_start2,Substructure_end1,Substructure_end2,Neandertal_Gene_Flow_absolute_end,Neandertal_Gene_Flow_absolute_start,Split_Time_Neandertals,Split_Time_non_Africans]
    demographic_events=[EUR_AFR_Gene_Flow_start,AFR_EUR_Gene_Flow_start,EUR_AFR_Gene_Flow_end,AFR_EUR_Gene_Flow_end,
    Substructure_start1,Substructure_start2,Substructure_end1,Substructure_end2,Neandertal_Gene_Flow_absolute_end,Neandertal_Gene_Flow_absolute_start,Split_Time_Neandertals,Split_Time_non_Africans]

    demographic_events=GF_Model_IV(TimeSpan=int(params['split_time_non_Africans']),GF_start=int(params['GF_start']),GF_stop=int(params['GF_stop']),demographic_events=demographic_events,GF_rate=float(params['migration_rate_GF']),loc=int(master['GF_Model'][1]),Folder_name=Folder_name,master_name=master_name,number=number,GF_Model=GF_Model)


    if str(master['Complex_Demography'][0]) == 'True':
        print("Simulating using inferred demographies")
        Initial_Neandertal_Size=msprime.PopulationParametersChange(time=int(params['split_time_Neandertals']),initial_size=int(params['N_Africans']), growth_rate=0,population_id=2)
        Initial_YRI_Size=msprime.PopulationParametersChange(time=int(params['split_time_Neandertals']),initial_size=int(params['N_Africans']), growth_rate=0,population_id=0)
        Initial_CEU_Size=msprime.PopulationParametersChange(time=int(params['split_time_Neandertals']),initial_size=int(params['N_Africans']), growth_rate=0,population_id=1)
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
    check.print_history(output=open(f"{Folder_name}_{master_name}_Demography_history.txt", 'w'))

    if str(master['Recombination_Map'][0]) == 'False':
        simulation=msprime.simulate(
            samples=samples,
            #model="dtwf",
            mutation_rate=float(params['mutation_rate']),
            length=float(params['chr_length']),
            recombination_rate=float(params['recombination_rate']),
            population_configurations=population_configuration,
            migration_matrix=migration_matrix,
            demographic_events=demographic_events
            )
    elif str(master['Recombination_Map'][0]) == 'True':
        print("Simulating using empirical recombination maps")
        simulation=msprime.simulate(
            samples=samples,
            #model="dtwf",
            mutation_rate=float(params['mutation_rate']),
            recombination_map=msprime.RecombinationMap.read_hapmap('{0}'.format(str(master['Recombination_Map'][1]))),
            population_configurations=population_configuration,
            migration_matrix=migration_matrix,
            demographic_events=demographic_events)



    return simulation




############################################################################################
# Choose which Set to process by giving Set_name the name of the Set to be processed #
Set_name='Fig_2_B_complex_revision'

# Choose number of replicates #
replicates=25

# Choose Folder #
Folder_name='../../Fig_2_B_complex_revision'

# Choose Result Folder Name #
Result_Folder='../../Fig_2_B_complex_revision/Result_classic_Exp'

# Choose if lambda chould be fixed or variable while fitting the AIC_Lomax
Fix_lambda=False
Simple_pulse_only=True
############################################################################################


master_name=config['CONTROL']['%s' % Set_name]
master_Scenario= list(master_name.keys())
Sim_name=config['MASTER']['%s' % master_Scenario[0]]['simulation']


rule alle:
    input:
        expand("{Result_Folder}/Result_file_SIM_Raw_ALDER-Fit-{Set_name}-{GF_Model}-min_dist_Fit-{min_dist_Fit}-ascertainment{Ascertainment}-downsampled{downsample}-Recc{Recomb_correction}.txt",master_name=config['CONTROL']['%s' % Set_name],number=range(replicates),Set_name=Set_name,GF_Model=config['MASTER']['%s' % master_Scenario[0]]['GF_Model'][0],Result_Folder=Result_Folder,Folder_name=Folder_name,min_dist_Fit=config['MASTER']['%s' % master_Scenario[0]]['min_dist_Fit']
,Ascertainment=config['MASTER']['%s' % master_Scenario[0]]['ascertainment'],downsample=config['MASTER']['%s' % master_Scenario[0]]['downsample'],Recomb_correction=config['MASTER']['%s' % master_Scenario[0]]['Recombination_Correction']),



rule basic_simulation:
    input:
    output:
        sim_out="{Folder_name}/Simulation-{master_name}-Chr{nChr}-run{number}-{GF_Model}.trees"
        #sim_out="{Folder_name}/sim/Simulation-{master_name}-Chr{nChr}-run{number}-{GF_Model}_snps.tsv.gz"
    run:
        master=load_config_master(config['MASTER'],wildcards.master_name)
        params=load_config_sim_parameters(config['sim_parameters'],master['simulation'])
        output_prefix=f"{wildcards.Folder_name}/sim/Simulation-{wildcards.master_name}-Chr{wildcards.nChr}-run{wildcards.number}-{wildcards.GF_Model}"
        simulation_results=simulation(params=params,Folder_name=wildcards.Folder_name,output_prefix=output_prefix,master_name=wildcards.master_name,number=wildcards.number,master=master,GF_Model=wildcards.GF_Model);
        # write whole tree sequence using tskit
        simulation_results.dump('{}'.format(output[0]))
        #write_all_snps(f"{output_prefix}_snps.tsv.gz",
        #                   simulation_results,
        #                  chrom=wildcards.nChr)


rule ascertainment:
    input:
        #sim_in="{Folder_name}/sim/Simulation-{master_name}-Chr{nChr}-run{number}-{GF_Model}_snps.tsv.gz"
        sim_in="{Folder_name}/Simulation-{master_name}-Chr{nChr}-run{number}-{GF_Model}.trees"
    output:
        temp("{Folder_name}/Simulation-{master_name}-Chr{nChr}-run{number}-{GF_Model}-ascertainment{Ascertainment}-f.eigenstratgeno"),
        temp("{Folder_name}/Simulation-{master_name}-Chr{nChr}-run{number}-{GF_Model}-ascertainment{Ascertainment}-f.snp"),
        temp("{Folder_name}/Simulation-{master_name}-Chr{nChr}-run{number}-{GF_Model}-ascertainment{Ascertainment}-f.ind")
        #snp_data="{Folder_name}/sim/Simulation-{master_name}-Chr{nChr}-run{number}-{GF_Model}_snps.tsv.gz"
        #"{Folder_name}/Simulation-{master_name}-Chr{nChr}-run{number}-{GF_Model}-ascertainment{Ascertainment}-f.eigenstratgeno",
        #"{Folder_name}/Simulation-{master_name}-Chr{nChr}-run{number}-{GF_Model}-ascertainment{Ascertainment}-f.snp",
        #"{Folder_name}/Simulation-{master_name}-Chr{nChr}-run{number}-{GF_Model}-ascertainment{Ascertainment}-f.ind"
    run:
        master=load_config_master(config['MASTER'],wildcards.master_name)
        params=load_config_sim_parameters(config['sim_parameters'],master['simulation'])
        print(wildcards.Ascertainment)
        simulation_results=tskit.load(input[0])
        snps2=get_all_snps(simulation_results)
        #snps = pd.read_csv(input.sim_in, sep="\t")
        #snps.pos = snps.pos.astype(int)
        #snps2 = snps.drop_duplicates(subset=['pos'])

        # Write geno and snp files


        afr_cols = [col for col in snps2.columns if col.startswith("AFR")]
        eur_cols = [col for col in snps2.columns if col.startswith("EUR")]
        nea_cols = [col for col in snps2.columns if col.startswith("NEA")]
        n_afr, n_eur, n_nea = len(afr_cols), len(eur_cols), len(nea_cols)

        D = dict()
        #D['rs'] = snps2.inde.astype(int)
        D['rs'] = snps2.index.astype(int)
        D['chrom'] = str(int(wildcards.nChr)+1)
        #D['pos'] = snps2.pos.astype(int)
        #D['map'] = snps2.pos.astype(int) * float(params['recombination_rate'])
        D['pos'] = snps2.index.astype(int)
        D['map'] = snps2.index.astype(int) * float(params['recombination_rate'])

        D["AFR_alt"] = np.sum(snps2[afr_cols], 1)
        D["EUR_alt"] = np.sum(snps2[eur_cols], 1)
        D["NEA_alt"] = np.sum(snps2[nea_cols], 1)

        D["AFR_ref"] = n_afr - D['AFR_alt']
        D["EUR_ref"] = n_eur - D['EUR_alt']
        D["NEA_ref"] = n_nea - D['NEA_alt']

        geno_snp = pd.DataFrame.from_dict(D)

        pd.options.mode.chained_assignment = None
        if wildcards.Ascertainment == "I" :
            filter_ = (geno_snp.AFR_alt == 0) & (geno_snp.NEA_alt != 0)
            snp_file = geno_snp[filter_]
            snp_file = snp_file[['rs', 'chrom','map','pos']]
            snp_file.iat[0,2] = 0
            snp_file['rs'] = 'rs' + snp_file['rs'].astype(str)
            geno_file = snps2[filter_]
            geno_file = pd.DataFrame(np.add.reduceat(geno_file.values, np.arange(len(geno_file.columns))[::2], axis=1))
        elif wildcards.Ascertainment == "II":
            filter_ = (geno_snp.AFR_alt == 0) & (geno_snp.NEA_alt != 0) & (geno_snp.EUR_alt != 0)
            snp_file = geno_snp[filter_]
            snp_file = snp_file[['rs', 'chrom','map','pos']]
            snp_file.iat[0,2] = 0
            snp_file['rs'] = 'rs' + snp_file['rs'].astype(str)
            geno_file = snps2[filter_]
            geno_file = pd.DataFrame(np.add.reduceat(geno_file.values, np.arange(len(geno_file.columns))[::2], axis=1))
        else:
            raise ValueError("ascertainment not known")

        np.savetxt('{}'.format(output[0]), geno_file.astype(int), fmt='%i', delimiter="")
        snp_file.to_csv('{}'.format(output[1]), header=None, index=None, sep='\t', mode='a')
        with open('{}'.format(output[2]),'w') as outind:
            for i in range(0,int(n_afr/2)):
                outind.write('Ind_{0}\tU\tAfricans\n'.format(i))
            for i in range(0,int(n_eur/2)):
                outind.write('Ind_{0}\tU\tNon_Africans\n'.format(i))
            for i in range(0,int(n_nea/2)):
                outind.write('Ind_{0}\tU\tNeandertals\n'.format(i))
            outind.close()

rule downsample:
    input:
        ingeno="{Folder_name}/Simulation-{master_name}-Chr{nChr}-run{number}-{GF_Model}-ascertainment{Ascertainment}-f.eigenstratgeno",
        insnp="{Folder_name}/Simulation-{master_name}-Chr{nChr}-run{number}-{GF_Model}-ascertainment{Ascertainment}-f.snp",
        inind="{Folder_name}/Simulation-{master_name}-Chr{nChr}-run{number}-{GF_Model}-ascertainment{Ascertainment}-f.ind"
    output:
        genoout=temp("{Folder_name}/Simulation-{master_name}-Chr{nChr}-run{number}-{GF_Model}-ascertainment{Ascertainment}-downsampled{downsample}.eigenstratgeno"),
        snpout=temp("{Folder_name}/Simulation-{master_name}-Chr{nChr}-run{number}-{GF_Model}-ascertainment{Ascertainment}-downsampled{downsample}.snp"),
        indout=temp("{Folder_name}/Simulation-{master_name}-Chr{nChr}-run{number}-{GF_Model}-ascertainment{Ascertainment}-downsampled{downsample}.ind"),
    run:
        downsample=wildcards.downsample
        downsample=int(downsample)/100
        genofile=pd.read_csv('{}'.format(input.ingeno),header=None)
        snpfile=pd.read_csv('{}'.format(input.insnp),header=None, delimiter="\t")
        with open('{}'.format(input.inind),'r')  as inind:
            indfile = inind.readlines()
            inind.close()
        if len(snpfile) != len(genofile):
            raise ValueError("snp and genofile must have same number of rows")
        with open('{}'.format(output.indout),'w')  as outind_d:
            for row in indfile:
              outind_d.write(str(row))
        if(downsample==1):
            snpfile.to_csv('{}'.format(output.snpout), header=None, index=None, sep='\t')
            with open('{}'.format(output.genoout),'w')  as outgeno_d:
                for row in genofile[0]:
                  outgeno_d.write('{}\n'.format(str(row)))

        else:
            k=len(snpfile)
            drop_indices = np.random.choice(snpfile.index[1:], int(k*(1-downsample)), replace=False)
            snpfile_d = snpfile.drop(drop_indices)
            genofile_d = genofile.drop(drop_indices)
            snpfile_d.to_csv('{}'.format(output.snpout), header=None, index=None, sep='\t')
            with open('{}'.format(output.genoout),'w')  as outgeno_d:
                for row in genofile_d[0]:
                  outgeno_d.write('{}\n'.format(str(row)))

def merge_geno_files(wildcards):
    files = expand("{Folder_name}/Simulation-{master_name}-Chr{nChr}-run{number}-{GF_Model}-ascertainment{Ascertainment}-downsampled{downsample}.eigenstratgeno",
master_name=wildcards.master_name,number=wildcards.number,Set_name=Set_name,GF_Model=wildcards.GF_Model,Folder_name=Folder_name,
nChr=range(config['sim_parameters']['_MS_default_']['num_chr']),Ascertainment=wildcards.Ascertainment,downsample=wildcards.downsample)
    return files

rule merge_geno_files:
    input:
        merge_geno_files
    output:
        temp("{Folder_name}/Ascertained_Simulation_merged-{master_name}-run{number}-{GF_Model}-ascertainment{Ascertainment}-downsampled{downsample}.eigenstratgeno")
        #"{Folder_name}/Ascertained_Simulation_merged-{master_name}-run{number}-{GF_Model}-ascertainment{Ascertainment}-downsampled{downsample}.eigenstratgeno"
    run:
        with open('{}'.format(output[0]),'w') as out_merged_geno:
            for i, file in enumerate(input):
                with open(file) as data:
                    for line in data:
                        out_merged_geno.write(line)

def merge_snp_files(wildcards):
    files = expand("{Folder_name}/Simulation-{master_name}-Chr{nChr}-run{number}-{GF_Model}-ascertainment{Ascertainment}-downsampled{downsample}.snp",
master_name=wildcards.master_name,number=wildcards.number,Set_name=Set_name,GF_Model=wildcards.GF_Model,
Folder_name=Folder_name,nChr=range(config['sim_parameters']['_MS_default_']['num_chr']),Ascertainment=wildcards.Ascertainment,downsample=wildcards.downsample)
    return files

rule merge_snp_files:
    input:
        merge_snp_files
    output:
        temp("{Folder_name}/Ascertained_Simulation_merged-{master_name}-run{number}-{GF_Model}-ascertainment{Ascertainment}-downsampled{downsample}.snp")
        #"{Folder_name}/Ascertained_Simulation_merged-{master_name}-run{number}-{GF_Model}-ascertainment{Ascertainment}-downsampled{downsample}.snp"
    run:
        with open('{}'.format(output[0]),'w') as out_merged_geno:
            for i, file in enumerate(input):
                with open(file) as data:
                    for line in data:
                        out_merged_geno.write(line)

def get_recomb_rate_for_correctio(master_name):
    master_name=master_name
    cor_rr=config['MASTER']['%s' % master_name]['correction_recomb_rate']

    return cor_rr

rule Set_Genetic_Distance:
    input:
        "{Folder_name}/Ascertained_Simulation_merged-{master_name}-run{number}-{GF_Model}-ascertainment{Ascertainment}-downsampled{downsample}.snp"
    output:
        temp("{Folder_name}/Simulation_Genetic_Distance_set-{master_name}-run{number}-{GF_Model}-ascertainment{Ascertainment}-downsampled{downsample}-Recc{Recomb_correction}.snp")
    run:
        master=load_config_master(config['MASTER'],wildcards.master_name)
        params=load_config_sim_parameters(config['sim_parameters'],master['simulation'])
        if str(master['Recombination_Map'][0]) == 'False':
            x=float(params['recombination_rate'])
            print(x)
            shell("""awk '{{print $1,$2, $3=$4*{x}, $4}}' {input} | sed '{{s/ /\t/g}}' > {output} """)
        elif str(master['Recombination_Map'][0]) == 'True' and wildcards.Recomb_correction == 'No_correction':
            cor_rr=get_recomb_rate_for_correctio(wildcards.master_name)
            print('Correcting genetic distance using: ',cor_rr, ' M/bp')
            x=cor_rr
            shell("""awk '{{print $1,$2, $3=$4*{x}, $4}}' {input} | sed '{{s/ /\t/g}}' > {output} """)
        elif str(master['Recombination_Map'][0]) == 'True' and wildcards.Recomb_correction == 'AAMap_correction':
            print('Correcting genetic distance using: the African-American-Map')
            shell("""Rscript /r1/people/leonardo_iasi/Desktop/Neandertal_Human_Introgression_Project/Paper/Paper_Scripts/AAMap_interpolation.R {input[0]} "/home/leonardo_iasi/Desktop/Neandertal_Human_Introgression_Project/Check_Cox_2019_Paper/Recombination_Maps/Scaled_Genetic_AAMap.txt" {output[0]}""")
        elif str(master['Recombination_Map'][0]) == 'True' and wildcards.Recomb_correction == 'HapMap_correction':
            print('Correcting genetic distance using: the HapMap-Map')
            shell("""Rscript /r1/people/leonardo_iasi/Desktop/Neandertal_Human_Introgression_Project/Paper/Paper_Scripts/HapMap_interpolation.R {input[0]} "/home/leonardo_iasi/Desktop/Neandertal_Human_Introgression_Project/Paper/Recombination_Maps/Scaled_HapMap.txt" {output[0]}""")


rule run_ALDER:
    input:
        "{Folder_name}/Ascertained_Simulation_merged-{master_name}-run{number}-{GF_Model}-ascertainment{Ascertainment}-downsampled{downsample}.eigenstratgeno",
        "{Folder_name}/Simulation_Genetic_Distance_set-{master_name}-run{number}-{GF_Model}-ascertainment{Ascertainment}-downsampled{downsample}-Recc{Recomb_correction}.snp",
        "{Folder_name}/Simulation-{master_name}-Chr0-run{number}-{GF_Model}-ascertainment{Ascertainment}-downsampled{downsample}.ind"

    output:
        "{Folder_name}/Raw_ALDER_output-{master_name}-run{number}-{GF_Model}-ascertainment{Ascertainment}-downsampled{downsample}-Recc{Recomb_correction}.txt"

    threads: 8
    run:
        master=load_config_master(config['MASTER'],wildcards.master_name)
        alder=load_ALDER_parameters(config['ALDER'],master['alder'])
        create_parfile_for_ALDER(alder=alder, master_name=wildcards.master_name,number=wildcards.number,Folder_name=wildcards.Folder_name,GF_Model=wildcards.GF_Model,Ascertainment=wildcards.Ascertainment,downsample=wildcards.downsample,Recomb_correction=wildcards.Recomb_correction)
        shell("""set +e
~leonardo_iasi/alder/alder -p {wildcards.Folder_name}/parfile_alder-{wildcards.master_name}-run{wildcards.number}-{wildcards.GF_Model}-ascertainment{wildcards.Ascertainment}-downsampled{wildcards.downsample}-Recc{wildcards.Recomb_correction}.txt
exitcode=$?
if [ $exitcode -eq 1 ]
then
        exit 0
else
    exit 0
fi""")

def get_max_length(wildcards):
    master_name=wildcards.master_name
    max_dist=config['MASTER']['%s' % master_name]['max_dist_Fit']

    return max_dist


rule fit_Curve:
    input:
        "{Folder_name}/Raw_ALDER_output-{master_name}-run{number}-{GF_Model}-ascertainment{Ascertainment}-downsampled{downsample}-Recc{Recomb_correction}.txt"
    params:
        max_dist=get_max_length,
        Fix_lambda=Fix_lambda,
        Only_Simple_Pulse_fit=Simple_pulse_only
    output:
        #"{Folder_name}/Model_Fit/Summary-Fit-{master_name}-run{number}-{GF_Model}-min_dist_Fit-{min_dist_Fit}-ascertainment{Ascertainment}-Recc{Recomb_correction}.log"
        "{Folder_name}/Model_Fit_classic_Exp/Summary-Fit-{master_name}-run{number}-{GF_Model}-min_dist_Fit-{min_dist_Fit}-ascertainment{Ascertainment}-downsampled{downsample}-Recc{Recomb_correction}.log"
        #"{Folder_name}/Model_Fit/Plot-Fit-{master_name}-run{number}-{GF_Model}-min_dist_Fit-{min_dist_Fit}-ascertainment{Ascertainment}.pdf",
    script:
        #"/r1/people/leonardo_iasi/Desktop/Neandertal_Human_Introgression_Project/Paper/Paper_Scripts/Fit_Exponential_and_Lomax_Snake.R"
        "/r1/people/leonardo_iasi/Desktop/Neandertal_Human_Introgression_Project/Paper/Paper_Scripts/Fit_Exponential_new.R"

rule grep_results_and_merge_with_Scenarios:
    input:
        #"{Folder_name}/Model_Fit/Summary-Fit-{master_name}-run{number}-{GF_Model}-min_dist_Fit-{min_dist_Fit}-ascertainment{Ascertainment}-Recc{Recomb_correction}.log"
        "{Folder_name}/Model_Fit_classic_Exp/Summary-Fit-{master_name}-run{number}-{GF_Model}-min_dist_Fit-{min_dist_Fit}-ascertainment{Ascertainment}-downsampled{downsample}-Recc{Recomb_correction}.log"
    output:
        temp("{Folder_name}/Result-Fit_Output-{master_name}-run{number}-{GF_Model}-min_dist_Fit-{min_dist_Fit}-ascertainment{Ascertainment}-downsampled{downsample}-Recc{Recomb_correction}.txt")

    run:
        master=load_config_master(config['MASTER'],wildcards.master_name)
        params=load_config_sim_parameters(config['sim_parameters'],master['simulation'])
        #shell("""grep "A, s, c, RSS_Expo, AIC_Expo, A, s, w,c, RSS_Lomax, AIC_Lomax, F_Test:" {input} |  sed 's/ /\t/g' |
#awk '{{print $12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24="{wildcards.master_name}",$25="{params[GF_start]}",$26="{params[GF_stop]}",$27="{wildcards.Ascertainment}",$28="{wildcards.min_dist_Fit}",$29="{wildcards.downsample}"}}' > {output} """)
        shell("""grep "A, m, c, RSS_Expo, nls status:" {input} |  sed 's/ /\t/g' |
awk '{{print $7,$8,$9,$10,$11,$12="{wildcards.master_name}",$13="{params[GF_start]}",$14="{params[GF_stop]}",$15="{wildcards.Ascertainment}",$16="{wildcards.min_dist_Fit}",$17="{wildcards.downsample}"}}' > {output} """)


def merge_all_files(wildcards):
    files = expand("{Folder_name}/Result-Fit_Output-{master_name}-run{number}-{GF_Model}-min_dist_Fit-{min_dist_Fit}-ascertainment{Ascertainment}-downsampled{downsample}-Recc{Recomb_correction}.txt",master_name=config['CONTROL']['%s' % Set_name],number=range(replicates),Set_name=Set_name,Folder_name=Folder_name,GF_Model=wildcards.GF_Model,min_dist_Fit= wildcards.min_dist_Fit,Ascertainment=wildcards.Ascertainment,downsample=wildcards.downsample,Recomb_correction=wildcards.Recomb_correction)
    return files

rule merge_all_files:
    input:
        merge_all_files

    output:
        "{Result_Folder}/Result_file_SIM_Raw_ALDER-Fit-{Set_name}-{GF_Model}-min_dist_Fit-{min_dist_Fit}-ascertainment{Ascertainment}-downsampled{downsample}-Recc{Recomb_correction}.txt"
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
