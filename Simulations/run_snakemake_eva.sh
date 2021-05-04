mkdir ../../autosnake_Fig_2_D_Complex_corrected_revision
snakemake -s Snakemake_Simulation_pipeline_efficient_trees_in_memory.py --keep-going -w 120 --jobname "_{rule}_{jobid}" -j800 \
    --cluster-config config/cluster.yaml \
    --cluster "qsub -q all.q -l h_vmem={cluster.mem} -o ../../autosnake_Fig_2_D_Complex_corrected_revision  -e ../../autosnake_Fig_2_D_Complex_corrected_revision -cwd -V -S /bin/bash" $*
