# Extended Admixture Pulse

This repository contains the scripts, pipelines and analysis steps used in this project. 

## Fitting the simple and extended pulse model
The descriptions of the models can be found here (..put biorivx link ...). In brief, one approach to learn about admixture dates from genetic data uses a recombination clock model: Conceptually, admixture segments are the result of the introduced chromosomes being broken down by recombination; the offspring of an archaic and a modern human parent will have one whole chromosome each of either ancestry. Thus, the offsprings' markers are in full ancestry linkage disequilibrium (ALD); all archaic variants are present on one DNA molecule, and all modern human one on the other one. In each generation meiotic recombination will reshuffle the chromosomes, progressively breaking down the ancestral chromosome down into shorter segments of archaic ancestry and the Admixture induced Linkage Disequilibrium (ALD) similarly decreases with each generation after gene flow. One can fit the decline of segment length or ALD to resolve for the time since the admixture event. The simple pulse model assumes that all admixture happens in one generation and uses an exponential distribution to fit the segment length (using the exponential PDF) or ALD (using the tail function of the exponential) to resolve for the mean time of gene flow (tm). 

The extended pulse relaxes the one generation assumption by modeling the admixture process over a certain duration (td) as a gamma distribution with shape = k and scale = tm/k.
Where, td= 4 tm k^(-1/2). 

This results in a havier tailed length distribution of segments/ALD (Lomax pdf/Lomax tail function). 

### Fitting using ALD

To fit the simple and extended pulse using ALD, the raw output from the ALDER program (Loh et al. 2013) can be used. The data has 3 columns: col 1 = distance between bins of SNPs the LD is computed for in centiMorgan, col 2 = the LD per bin, col 3 (optional) = bin count.
The script is based on Moorjani et al. 2016 and can be found in [Extended_Admixture_Pulse](Extended_Admixture_Pulse_inferrence/). Three parameters must be provided lval = lower value of distance between SNPs, hval = maximum value of distance and 
affine = (logical) if a parameyter modelin background LD should be used.

### Fitting using segment length

Comming soon

## Simulation pipeline
The folder [Extended_Admixture_Pulse](Simulations) contains 
a Snakemake script that executes msprime (Kelleher et al. 2016) simulations for Neandertal admixture into non-Africans. The parameters for the simulations are given
in the [Extended_Admixture_Pulse](Simulations\config) folder. The Recombination maps and effective population size changes are read in from [Extended_Admixture_Pulse](Simulations\Recombination_Maps) and 
[Extended_Admixture_Pulse](Simulations\Inferred_Pop_History). The ALDER program (Loh et al. 2013) is used to compute the weighted LD for Neandertal introgressed variants (using the raw_outname option to get the raw weighted LD per genetic distance file). The LD is then
fitted to a simple and/or extended pulse as described before.
