### First script used to filter the NAM parent raw Variant Call File (VCF): filtering for SNP quality, bi-allelic SNPs, and missing data ###

# Keep only loci with a quality score greater than or equal to 30.0, no missing data, no INDELs, and sites that are bi-allelic
vcftools --gzvcf vcftools --gzvcf NAM_parents_all.vcf.gz --minQ 30.0 --max-missing 1.0 --max-alleles 2 --min-alleles 2 --remove-indels --recode --recode-INFO-all --out NAM_parents_all_q30_nomiss_biSNPs

# Compress the lightly filtered NAM parent VCF file for input into vcfR for the second filtering script #
bgzip NAM_parents_all_q30_nomiss_biSNPs.recode.vcf