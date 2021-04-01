mkdir ../autosnake_Fig_2_B_complex_revision
snakemake -w 120 --jobname "_{rule}_{jobid}" -j800 \
    --cluster-config config/cluster.yaml \
    --cluster "qsub -q all.q -l h_vmem={cluster.mem} -o ../../autosnake_Fig_2_B_complex_revision  -e ../../autosnake_Fig_2_B_complex_revision -cwd -V -S /bin/bash" $*  

