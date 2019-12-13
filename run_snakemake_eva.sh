mkdir ../autosnake_Recent_Sampling
snakemake -w 120 --jobname "_{rule}_{jobid}" -j500 \
    --cluster-config config/cluster.yaml \
    --cluster "qsub -q all.q -l h_vmem={cluster.mem} -o ../autosnake_Recent_Sampling/ -e ../autosnake_Recent_Sampling/ -cwd -V -S /bin/bash" $*  

