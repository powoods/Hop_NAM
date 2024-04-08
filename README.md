### README File For Hop NAM Mapping Project ###
## Author: Patrick Woods ##

The code contained in this repo are those used to process and analyze all data generated generated for the Hop NAM experiment conducted by the USDA - ARS. 

Four script files are present, one written in R (version 4.2.0) another written in Python (version 3.11.8)and the two more written in shell.

The R and Python scripts achieve the same results by estimating broad sense heritability, calculating estimated marginal means and converting these means to binary (0, 1 format).

The R and Python scripts also contain the code used for quantilnormalization of the phenotype data.

The R script uses a line by line approach while the python script uses a single function to perform all required tasks.

The R script also has additional code to calculate statistics across all NAM families and within family phenotype segregation tests.

The shell script titled "NAM_parent_fastq_processing.sh" contains all code used to process raw paired end fastq.gz files for NAM parents from start to finish (read trimmed to 
variant calling).
