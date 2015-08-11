# LDhat
LDhat is a package written in the C and C++ languages for the analysis of recombination rates from population genetic data. 
The package is available either as C/C++ source code. <br><br>
The following programs are included:

- <b>convert</b>: Simple manipulation of DNA sequence data
- <b>pairwise</b>: Parametric and nonparametric analyses of recombination
- <b>interval</b>: Estimation of variable recombination rates
- <b>rhomap</b>: Estimation of variable recombination rates in the presence of hotspots
- <b>fin</b>: Simulation of genetic variation data in LDhat format

Other programs in the package can be used to generate lookup tables and summarise the output of analyses. 
Programs are run from the command line. 
A [PDF manual](manual.pdf) containing information about installing, running and interpreting results from the package should be read before proceeding. 
An [example data set](Example) - LPL sequences from a Finnish population ([Nickerson <i>et al</i>., 1998](http://www.ncbi.nlm.nih.gov/pubmed/9662394)) - is included in the package.

Installation
------------

A makefile is included, so it should be possible to compile LDhat by just typing `make`. The resulting executables are run from the command line.

Likelihood Lookup Tables
------------------------

A number of precomputed likelihood lookup tables are available for [download](lk_files). 
Using these tables in conjunction with <b>lkgen</b> will speed up analyses considerably, as the calculation of 2-locus coalescent likelihoods is the most computationally intesive aspect of the core LDhat algorithms. 
These tables all use the same grid for 4N<sub>e</sub>r (rho), 101 points evenly spaced between 0 and 100. 
They differ only in the number of sequences (n) and the value of population mutation rate parameter (theta). 
Note that <b>lkgen</b> can be used to generated lookup tables from these for smaller numbers of chromosomes, and that minor differences in theta do not appear to strongly influence the results. 
<br><br>
All of the likelihood lookup files are gzip compressed, and should be decompressed before use.
<br>
- [n=192, theta=0.001 per site](lk_files/lk_n192_t0.001.gz)
- [n=120, theta=0.001 per site](lk_files/lk_n120_t0.001.gz)
- [n=100, theta=0.001 per site](lk_files/lk_n100_t0.001.gz)
- [n=100, theta=0.01 per site](lk_files/lk_n100_t0.01.gz)
- [n=50, theta=0.001 per site](lk_files/lk_n50_t0.001.gz)
- [n=50, theta=0.1 per site](lk_files/lk_n50_t0.1.gz)
- [n=50, theta=0.5 per site](lk_files/lk_n50_t0.5.gz)

Getting help
------------

Send questions, comments, or suggestions to ldhat-help@lists.sourceforge.net

