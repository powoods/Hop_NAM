### README File For Hop NAM Mapping Project ###
## Author: Patrick Woods ##

The code contained in this repo are those used to process and analyze all data generated generated for the Hop NAM experiment conducted by the USDA - ARS. 

Two script files are present, one written in R (version 4.2.0) and the other written in Python (version 3.11.3).

Both scripts achieve the same results by estimating broad sense heritability, calculating estimated marginal means and converting these means to binary (0, 1 format).

The R script focuses more on a line by line approach while the python script uses a single function to perform all required tasks.

The python script contains two versions of the main function: 1) the initial version (est_nam_H2()) which only estimates heritability and marginal means, and 2) the most recent version 
(est_nam_H2_means()) which estimates 
heritability, marginal means, and can convert marginal means to binary (0, 1) format.

Currently only the R script possesses the code for quantile normalization of susceptibility phenotypes. This functionality will be available in the python script in the next releases.
