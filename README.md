# Extended Admixture Pulse

This repository contains the scripts, pipelines and analysis steps used in this project. 

## Extended Admixture Pulse Inferrence
This folder contains the a short description of the simple and extended pulse model and the R scripts to fit the model to data.
A full descriptions of the models can be found here (https://www.biorxiv.org/content/10.1101/2021.04.04.438357v1)
The script [Fitting_simple_and_extended_pulse_to_ALD.R](Extended_Admixture_Pulse_inferrence/Fitting_simple_and_extended_pulse_to_ALD.R) contains the functions to fit the two models to ALD data.
The script [Fitting_simple_and_extended_pulse_to_segment_length.R](Extended_Admixture_Pulse_inferrence/Fitting_simple_and_extended_pulse_to_segment_length.R) contains the functions to fit the two models to segment data.


## Simulation pipeline
The folder [Simulations](Simulations) contains 
a Snakemake script that executes msprime (Kelleher et al. 2016) simulations for Neandertal admixture into non-Africans. The parameters for the simulations are given
in the [Simulations/config](Simulations/config) folder. The Recombination maps and effective population size changes are read in from [Simulations/Recombination_Maps](Simulations/Recombination_Maps) and 
[Simulations/Inferred_Pop_History](Simulations/Inferred_Pop_History). The ALDER program (Loh et al. 2013) is used to compute the weighted LD for Neandertal introgressed variants (using the raw_outname option to get the raw weighted LD per genetic distance file). The LD is then
fitted to a simple and/or extended pulse as described before.

## Paper

This folder contains the Rmd script calculating all statistics and producing all the figures in the paper.
