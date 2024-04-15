### README File For Hop NAM Mapping Project ###
## Author: Patrick Woods ##

The code contained in this repo are those used to process and analyze all data generated generated for the Hop NAM experiment conducted by the USDA - ARS. 

Four script files are present, one written in R (version 4.2.0) another written in Python (version 3.11.8)and the two more written in shell.

The R and Python scripts called "NAM_phenotype_data_processing.R/.py" include all scripts used to process the raw NAM phenotype data.

The shell script titled "NAM_parent_fastq_processing.sh" contains all code used to process raw paired end fastq.gz files for NAM parents from start to finish (read trimmed to 
variant calling).

The shell scripts titled "NAM_parent_VCF_filtering_1/2/3.sh" contain all code (in order) used to filter the NAM parent VCF from start to finish.