mkdir ../autosnake_Fig_2_B_complex
snakemake -w 120 --jobname "_{rule}_{jobid}" -j500 \
    --cluster-config config/cluster.yaml \
    --cluster "qsub -q all.q -l h_vmem={cluster.mem} -o ../autosnake_Fig_2_B_complex  -e ../autosnake_Fig_2_B_complex -cwd -V -S /bin/bash" $*  

