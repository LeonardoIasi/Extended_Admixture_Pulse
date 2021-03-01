# Extended Admixture Pulse

This repository contains the scripts, pipelines and analysis steps used in this project. 

## Fitting the simple and extended pulse model
The descriptions of the models can be found here (..put biorivx link ...). In brife, one approach to learn about admixture dates from genetic data uses a \emph{recombination clock} model: Conceptually, admixture segments are the result of the introduced chromosomes being broken down by recombination; the offspring of an archaic and a modern human parent will have one whole chromosome each of either ancestry. Thus, the offsprings' markers are in full ancestry linkage disequilibrium (ALD); all archaic variants are present on one DNA molecule, and all modern human one on the other one. In each generation meiotic recombination will reshuffle the chromosomes, progressively breaking down the ancestral chromosome down into shorter segments of archaic ancestry and the Admixture induced Linkage Disequilibrium (ALD) similarly decreases with each generation after gene flow. One can fit the decline of segment length or ALD to resolve for the time since the admixture event. The simple pulse model assumes that all admixture happens in one generation and uses an exponential distribution to fit the segment length (using the exponential PDF) or ALD (using the tail function of the exponential) to resolve for the mean time of gene flow ($t_m$). 

The extended pulse relaxes the one generation assumption by modeling the admixture process over a ceartain duration ($t_d$) as a gamma distribution with shape = k and scale = tm/k.
Where, $t_d=4t_m k^{-\frac{1}{2}} $

$$  P(T_i=t)=\frac{1}{\Gamma(k)(\frac{t_m}{k})^k}t^{k-1}e^{-t\frac{k}{t_m}} $$

This results in a havier tailed length distribution of segments/ALD (Lomax pdf/Lomax tail function). 

$$  P(L=l | k, t_m) = t_{m}^{-k} \ \Bigg( \frac{k}{l+\frac{k}{t_{m}}}\Bigg)^{k+1} $$

### Fitting using ALD

To fit the simple and extended pulse using ALD, the raw output from the ALDER program (Loh et al. 2013) can be used. 



## Simulation pipeline
The folder [Extended_Admixture_Pulse](Simulations) contains 
a Snakemake script that executes msprime (Kelleher et al. 2016) simlations for Neandertal admixture into non-Africans. The parameters for the simulations are given
in the [Extended_Admixture_Pulse](Simulations\config) folder. The Recombination maps and effective population size changes are read in from [Extended_Admixture_Pulse](Simulations\Recombination_Maps) and 
[Extended_Admixture_Pulse](Simulations\Inferred_Pop_History). The ALDER program (Loh et al. 2013) is used to compute the weighted LD for Neadertal introgressed variants (using the raw_outname option to get the raw weighted LD per genetic distance file). The LD is then
fitted to a simple and/or extended pulse as described before.
