### Third and final script for filtering the NAM parent VCF file: Filtering for the 10 assembled chromosomes max read depth, minor allele count and missing data ###

# Extract all loci on chromosomes 1 to 10
vcftools --gzvcf NAM_parents_q30_biSNPs_mndp8_90CIdp_33_66hetdp.vcf.gz --chr Scaffold_1531 --chr Scaffold_19 --chr Scaffold_76 --chr Scaffold_24 --chr Scaffold_172 --chr Scaffold_77 --chr Scaffold_73 --chr Scaffold_49 --chr Scaffold_191 --chr Scaffold_1533 --recode --recode-INFO-all --out NAM_parents_q30_biSNPs_mndp8_90CIdp_33_66hetdp_nomiss_chr1_10

# Impose final filters of max read depth of 40, minor allele count of 1, and no missing data
vcftools --gzvcf NAM_parents_q30_biSNPs_mndp8_90CIdp_33_66hetdp_nomiss_chr1_10.recode.vcf.gz --maxDP 40.0 --mac 1 --max-missing 1.0 --recode --recode-INFO-all --out NAM_parents_q30_biSNPs_mndp8_90CIdp_33_66hetdp_mxdp40_mac1_nomiss_chr1_10

# BGzip and index the final VCF using Tabix
bgzip NAM_parents_q30_biSNPs_mndp8_90CIdp_33_66hetdp_mxdp40_mac1_nomiss_chr1_10.recode.vcf

tabix -p vcf NAM_parents_q30_biSNPs_mndp8_90CIdp_33_66hetdp_mxdp40_mac1_nomiss_chr1_10.recode.vcf.gz
