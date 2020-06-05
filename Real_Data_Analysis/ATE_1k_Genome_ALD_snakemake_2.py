#######################################
CHR=range(1,23)
Recomb_correction='CEU_specific_map'
#Recomb_correction='constant'
#######################################
def create_parfile_for_ALDER(output_file):
	file=open("parfile_alder-ALL_chr_hcNea_YRI_CEU_1k_G_LES_ascertained_corrected_rr.txt","w")
	file.write("""genotypename:\tALL_chr_hcNea_YRI_CEU_1k_G_LES_ascertained.geno
snpname:\tALL_chr_hcNea_YRI_CEU_1k_G_LES_ascertained_corrected_constant_rr.snp
indivname:\tATE_samples.ind
admixpop:\tEUR
refpops:\tAFR;Neandertal
checkmap:\tFalse
binsize:\t1e-05
mindis:\t1e-5
maxdis:\t0.1
fast_snp_read:\tTrue
num_threads:\t4
raw_outname:\t{0}""".format(output_file))
	file.close()

rule alle:
	input:
		expand("Raw_ALDER_output-ALL_chr_hcNea_YRI_CEU_1k_G_LES_ascertained_corrected_rr_{rrMap}.txt",rrMap=Recomb_correction)

rule Subsample_VCF:
	input:
		vcf="/mnt/sequencedb/gendivdata/2_genotypes/giantVcfs/merged_var_nosing_sites_arch_apes_sgdp1_g1000_chr{chr}.vcf.gz",
		samples="ATE_samples_filter.txt"
	output:
		"1kG_AFR_EUR_NEA_vcf_chr_{chr}.vcf.gz"
	run:
		shell("""bcftools view  -Oz -S {input.samples} {input.vcf} -o {output} """)

rule Ascertain_VCF:
	input:
		"1kG_AFR_EUR_NEA_vcf_chr_{chr}.vcf.gz"
	output:
		#temp("1kG_AFR_EUR_NEA_vcf_chr_{chr}_LES_ascertained.vcf")
		"1kG_AFR_EUR_NEA_vcf_chr_{chr}_LES_ascertained.vcf"
	run:
		shell("""zcat {input}  | awk 'NR>10{{for (i=10;i<=117;i++) {{if ($i!="0|0") next}}}}{{print}}'  | awk 'NR>10{{for (i=217;i<=220;i++) {{if ($i=="0/0") next}}}}{{print}}' > {output}""")

rule merge_VCF_1:
	input:
		"1kG_AFR_EUR_NEA_vcf_chr_1_LES_ascertained.vcf"
	output:
		temp("VCF_Header.vcf")
	run:
		shell("""grep '^#' {input} > {output}""")

def merge_vcf_files(wildcards):
	files = expand("1kG_AFR_EUR_NEA_vcf_chr_{chr}_LES_ascertained.vcf",chr=CHR)
	return files

rule merge_VCF_2:
	input:
		merge_vcf_files,
	output:
		temp("ALL_chr_hcNea_YRI_CEU_1k_G_LES_ascertained_0.vcf")
	run:
		shell("""for i in {input};
		do
		  grep -v '^#' $i >> {output}
		done""")


rule merge_VCF_3:
	input:
		vcf="ALL_chr_hcNea_YRI_CEU_1k_G_LES_ascertained_0.vcf",
		header="VCF_Header.vcf"
	output:
		"ALL_chr_hcNea_YRI_CEU_1k_G_LES_ascertained.vcf"
	run:
		shell("""cat {input.header} {input.vcf} >> {output}""")

rule vcf2eigen:
	input:
		"ALL_chr_hcNea_YRI_CEU_1k_G_LES_ascertained.vcf",
	output:
		"ALL_chr_hcNea_YRI_CEU_1k_G_LES_ascertained.geno",
		temp("ALL_chr_hcNea_YRI_CEU_1k_G_LES_ascertained.snp"),
		temp("ALL_chr_hcNea_YRI_CEU_1k_G_LES_ascertained.ind")
	run:
		shell("""cd ~/
		./vcf2eigenstrat/bin/vcf2eigenstrat /mnt/diversity/leonardo_iasi/Adm_Time_Dating_Sim/Real_Data_Analysis/1k_Genomes_ALD/ALL_chr_hcNea_YRI_CEU_1k_G_LES_ascertained.vcf /mnt/diversity/leonardo_iasi/Adm_Time_Dating_Sim/Real_Data_Analysis/1k_Genomes_ALD/ALL_chr_hcNea_YRI_CEU_1k_G_LES_ascertained
		cd /mnt/diversity/leonardo_iasi/Adm_Time_Dating_Sim/Real_Data_Analysis/1k_Genomes_ALD/""")

rule Set_Genetic_Distance:
    input:
        "ALL_chr_hcNea_YRI_CEU_1k_G_LES_ascertained.snp"
    output:
        "ALL_chr_hcNea_YRI_CEU_1k_G_LES_ascertained_corrected_rr_{rrMap}.snp"
    run:
		if wildcards.rrMap == 'constant':
			x=float(1.8e-8)
			shell("""awk '{{print "rs"$4,$2, $3=$4*{x}, $4, $5, $6}}' {input} | sed '{{s/ /\t/g}}' > {output}""")
		elif wildcards.rrMap == 'CEU_specific_map':
            shell("""Rscript /mnt/diversity/leonardo_iasi/Adm_Time_Dating_Sim/Real_Data_Analysis/1k_Genomes_ALD/CEU_RecomMap_Interpolation.R {input[0]} "/mnt/diversity/leonardo_iasi/Adm_Time_Dating_Sim/Real_Data_Analysis/1k_Genomes_ALD/All_chr_CEU_recombination_map_hapmap_format_hg19.txt" {output[0]}""")


rule run_ALDER:
	input:
		"ALL_chr_hcNea_YRI_CEU_1k_G_LES_ascertained.geno",
        "ALL_chr_hcNea_YRI_CEU_1k_G_LES_ascertained_corrected_rr_{rrMap}.snp",
		"ATE_samples.ind"

	output:
		"Raw_ALDER_output-ALL_chr_hcNea_YRI_CEU_1k_G_LES_ascertained_corrected_rr_{rrMap}.txt"

	threads: 8
	run:
		create_parfile_for_ALDER(output_file='{}'.format(output))
		shell("""set +e
~leonardo_iasi/alder/alder -p parfile_alder-ALL_chr_hcNea_YRI_CEU_1k_G_LES_ascertained_corrected_rr.txt
exitcode=$?
if [ $exitcode -eq 1 ]
then
		exit 0
else
	exit 0
fi""")
