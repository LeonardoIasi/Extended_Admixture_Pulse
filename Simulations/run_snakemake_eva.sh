<<<<<<< HEAD
mkdir ../../autosnake_Fig_2_C_corrected_revision
snakemake -s Snakemake_Simulation_pipeline_efficient.py --keep-going -w 120 --jobname "_{rule}_{jobid}" -j800 \
    --cluster-config config/cluster.yaml \
    --cluster "qsub -q all.q -l h_vmem={cluster.mem} -o ../../autosnake_Fig_2_C_corrected_revision  -e ../../autosnake_Fig_2_C_corrected_revision -cwd -V -S /bin/bash" $*  
=======
mkdir ../autosnake_Fig_2_B_complex_revision
snakemake -w 120 --jobname "_{rule}_{jobid}" -j800 \
    --cluster-config config/cluster.yaml \
    --cluster "qsub -q all.q -l h_vmem={cluster.mem} -o ../../autosnake_Fig_2_B_complex_revision  -e ../../autosnake_Fig_2_B_complex_revision -cwd -V -S /bin/bash" $*  
>>>>>>> 2375b3eed395c8debef4b08da0f411d92039204f

