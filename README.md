# Extended Admixture Pulse

This repository contains the scripts, pipelines and analysis steps used in this project. 

## Fitting the simple and extended pulse model
The descriptions of the models can be found here (..put biorivx link ...). In brief, one approach to learn about admixture dates from genetic data uses a recombination clock model: Conceptually, admixture segments are the result of the introduced chromosomes being broken down by recombination; the offspring of an archaic and a modern human parent will have one whole chromosome each of either ancestry. Thus, the offsprings' markers are in full ancestry linkage disequilibrium (ALD); all archaic variants are present on one DNA molecule, and all modern human one on the other one. In each generation meiotic recombination will reshuffle the chromosomes, progressively breaking down the ancestral chromosome down into shorter segments of archaic ancestry and the Admixture induced Linkage Disequilibrium (ALD) similarly decreases with each generation after gene flow. One can fit the decline of segment length or ALD to resolve for the time since the admixture event. The simple pulse model assumes that all admixture happens in one generation and uses an exponential distribution to fit the segment length (using the exponential PDF) or ALD (using the tail function of the exponential) to resolve for the mean time of gene flow (tm). 

The extended pulse relaxes the one generation assumption by modeling the admixture process over a certain duration (td) as a gamma distribution with shape = k and scale = tm/k.
Where, td= 4 tm k^(-1/2). 

This results in a heavier tailed length distribution of segments/ALD (Lomax pdf/Lomax tail function). 

The admixture estimates both for ALD or directly inferred segments rely on precise genetic distances. We do not recommend to use a constant recombination rate to convert physical distances into genetic distances, rather population specific recombination maps.

### Fitting using ALD

To fit the simple and extended pulse using ALD, the raw output from the ALDER program (Loh et al. 2013) can be used. The data has 3 columns: col 1 = distance between bins of SNPs the LD is computed for in centiMorgan, col 2 = the LD per bin, col 3 (optional) = bin count.
The script [Fitting_simple_and_extended_pulse_to_ALD.R](Extended_Admixture_Pulse_inferrence/Fitting_simple_and_extended_pulse_to_ALD.R) is based on Moorjani et al. 2016. It uses an residual sum of squares optimization function (DEoptim) to get good starting parameters for the nls function. Three parameters must be provided lval = lower value of distance between SNPs, hval = maximum value of distance and 
affine = (logical) if a parameter modeling background LD should be used. The lower and upper boundaries for the optimization are currently hard coded and are reasonable boundaries for Neandertal introgression, but one might want to change them for other scenarios.

### Fitting using segment length

The script [Fitting_simple_and_extended_pulse_to_segment_length.R](Extended_Admixture_Pulse_inferrence/Fitting_simple_and_extended_pulse_to_segment_length.R) to fit the simple/extended pulses to directly inferred segments (e.g. using Skov et al. 2018) uses the R optim function (method="L-BFGS-B") to fit the data to an exponential or lomax pdf (optionally a truncated pdf which is recommended). The input data must at least contain one column (length_cM) giving the length of unique fragments in centiMorgan. There are two additional parameters lower and upper truncation which indicates the bounds in which fragments can be relaiably be called (very small fragments will be missing, long segments might be called as multiple seperated segments).

## Simulation pipeline
The folder [Simulations](Simulations) contains 
a Snakemake script that executes msprime (Kelleher et al. 2016) simulations for Neandertal admixture into non-Africans. The parameters for the simulations are given
in the [Simulations/config](Simulations/config) folder. The Recombination maps and effective population size changes are read in from [Simulations/Recombination_Maps](Simulations/Recombination_Maps) and 
[Simulations/Inferred_Pop_History](Simulations/Inferred_Pop_History). The ALDER program (Loh et al. 2013) is used to compute the weighted LD for Neandertal introgressed variants (using the raw_outname option to get the raw weighted LD per genetic distance file). The LD is then
fitted to a simple and/or extended pulse as described before.

## Paper

This folder contains the Rmd script calculating all statistics and producing all the figures in the paper.
