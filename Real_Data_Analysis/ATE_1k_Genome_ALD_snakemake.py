import numpy as np

CHR=range(22)

rule alle:
	input:
		"merged_var_nosing_sites_arch_apes_sgdp1_g1000_chr{chr}_red_ascertained.vcf.gz",chr=CHR

rule Subsample_VCF:
	input:
		vcf="/mnt/sequencedb/gendivdata/2_genotypes/giantVcfs/merged_var_nosing_sites_arch_apes_sgdp1_g1000_chr{chr}.vcf.gz",
		samples="ATE_samples_red.txt"
	output:
		temp("merged_var_nosing_sites_arch_apes_sgdp1_g1000_chr{chr}_red.vcf.gz")
	run:
		shell("""bcftools view  -S {input.samples} {input.vcf} -o {output} """)

rule Ascertaine_VCF:
	input:
		temp("merged_var_nosing_sites_arch_apes_sgdp1_g1000_chr{chr}_red.vcf.gz")
	output:
		"merged_var_nosing_sites_arch_apes_sgdp1_g1000_chr{chr}_red_ascertained.vcf.gz"
	run:
		shell("""cat {input} | awk 'NR>10{{for (i=10;i<=117;i++) {{if ($i!="0|0") next}}}}{{print}}' | awk 'NR>10{{for (i=217;i<=220;i++) {{if ($i=="0/0") next}}}}{{print}}' > {output}}
 """)

 rule merge_VCF_1:
 	input:
		"merged_var_nosing_sites_arch_apes_sgdp1_g1000_chr1_red_ascertained.vcf.gz"
 	output:
 		"merge.vcf"
 	run:
 		shell("""grep '^#' {input} > merge.vcf""")

def merge_vcf_files(wildcards):
	files = expand(merged_var_nosing_sites_arch_apes_sgdp1_g1000_chr{chr}_red_ascertained.vcf.gz",chr=chr)
	return files

 rule merge_VCF_2:
 	input:
		merge_vcf_files,
		"merge.vcf"
 	output:
 		"merged_var_nosing_sites_arch_apes_sgdp1_g1000_ALL_chr_red_ascertained.vcf.gz"
 	run:
 		shell("""grep -v '^#' {input[0]}   >> {input[1]} | mv {input[1]} {output}""")



	rule reshape:
		input:
			"{Folder_name}/Simulation-{master_name}-run{number}-{GF_Model}.vcf"
		output:
			temp("{Folder_name}/Reshaped_Simulation-{master_name}-run{number}-{GF_Model}.vcf")

		run:
			shell("""sed '7,${{/^#/d;}}' {input} > {output} """)

	rule ascertainment_4:
		input:
			temp("{Folder_name}/Reshaped_and_ascerted_2_Simulation-{master_name}-run{number}-{GF_Model}.vcf"),
			temp("{Folder_name}/Vcf_Header-{master_name}-run{number}-{GF_Model}.vcf")
		output:
			"{Folder_name}/Reshaped_and_ascerted_Simulation-{master_name}-run{number}-{GF_Model}.vcf"
		shell:
			"""cat {input[1]} {input[0]} > {output}"""

	rule vcftools_conversion:
		input:
			"{Folder_name}/Reshaped_and_ascerted_Simulation-{master_name}-run{number}-{GF_Model}.vcf"

		output:
			temp("{Folder_name}/Simulation-{master_name}-run{number}-{GF_Model}.ped"),
			temp("{Folder_name}/Simulation-{master_name}-run{number}-{GF_Model}.map")

		shell:
			"vcftools --vcf {input} --plink --out {wildcards.Folder_name}/Simulation-{wildcards.master_name}-run{wildcards.number}-{wildcards.GF_Model}"

	rule remodel_pad_file:
		input:
			"{Folder_name}/Simulation-{master_name}-run{number}-{GF_Model}.ped"

		output:
			temp("{Folder_name}/Remodeled_Simulation-{master_name}-run{number}-{GF_Model}.ped")

		shell:
			"""awk '$6="1"' {input} | sed 's/ /\t/g' > {output}"""

	rule remodel_map_file:
		input:
			"{Folder_name}/Simulation-{master_name}-run{number}-{GF_Model}.map"
		output:
			#temp("{Folder_name}/Remodeled_Simulation-{master_name}-run{number}-{GF_Model}.map")
			"{Folder_name}/Remodeled_Simulation-{master_name}-run{number}-{GF_Model}.map"
		shell:
			"""sed 's/://g' {input} | awk '{{print $1,"rs"$2,$3,$4}}' | sed 's/ /\t/g'  > {output}"""

	rule EIGEN_conversion:
		input:
			temp("{Folder_name}/Remodeled_Simulation-{master_name}-run{number}-{GF_Model}.ped"),
			temp("{Folder_name}/Remodeled_Simulation-{master_name}-run{number}-{GF_Model}.map")

		output:
			temp("{Folder_name}/Simulation-{master_name}-run{number}-{GF_Model}.ind"),
			temp("{Folder_name}/Simulation-{master_name}-run{number}-{GF_Model}.snp"),
			"{Folder_name}/Simulation-{master_name}-run{number}-{GF_Model}.eigenstratgeno"
		run:
			create_parfile_for_convertf(wildcards.Folder_name,wildcards.master_name,wildcards.number,wildcards.GF_Model)
			shell("~leonardo_iasi/EIG-master/bin/convertf -p {wildcards.Folder_name}/parfile_convertf-{wildcards.master_name}-run{wildcards.number}-{wildcards.GF_Model}.txt ")
