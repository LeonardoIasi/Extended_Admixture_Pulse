mkdir ../autosnake_Variyng_Time_of_Recent_sampling_GF_l_fixed_Recomn_Map_Hap_Map
snakemake -w 120 --jobname "_{rule}_{jobid}" -j500 \
    --cluster-config config/cluster.yaml \
    --cluster "qsub -q all.q -l h_vmem={cluster.mem} -o ../autosnake_Variyng_Time_of_Recent_sampling_GF_l_fixed_Recomn_Map_Hap_Map  -e ../autosnake_Variyng_Time_of_Recent_sampling_GF_l_fixed_Recomn_Map_Hap_Map -cwd -V -S /bin/bash" $*  

