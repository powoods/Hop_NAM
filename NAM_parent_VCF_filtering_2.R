### Second script for filtering the NAM parent VCF: filtering for coverage ###

library(vcfR) #R package for basic coverage filtering
library(SNPfiltR) #R package for heterozygote coverage filtering

#Read VCF file from filtering script #1, assess the pre-filtering coverage and apply a 90% confidence interval + minimum depth filter of 8

vcf <- read.vcfR('NAM_parents_all_q30_nomiss_biSNPs.recode.vcf.gz', verbose = F) #read in the lightly filtered VCF
dp <- extract.gt(vcf, element = "DP", as.numeric = T) #calculate the raw read depth of every variant per sample
boxplot(dp, ylab = "Depth") #make a boxplot to visualize the raw read depth per locus per sample
sums <- apply(dp, MARGIN=2, quantile, probs=c(0.05,0.95), na.rm=TRUE) #calculate a 90% confidence interval for sample read depth
dp2 <- sweep(dp, MARGIN=2, FUN = "-", sums[1,])
dp[dp2 < 0] <- NA
dp2 <- sweep(dp, MARGIN =2, FUN = "-", sums[2,])
dp[dp2 > 0] <- NA
dp[dp < 8] <- NA #set a minimum read depth of 8 for each genotype call
boxplot(dp, ylab = "Depth") #visualize the filtered read depth of loci within the 90% confidence interval
is.na( vcf@gt[,-1][is.na(dp) ] ) <- TRUE #update the original vcfR object so that all genotype calls not within the 90% confidence interval are changed to "no call"
vcf #observe the increased percent of missing data
write.vcf(vcf, file = "NAM_parents_q30_biSNPs_mndp8_90CIdp.vcf.gz")


# Read in the VCF file we just produced for additional coverage filtering for heterozygotes
vcf2 <- read.vcfR('NAM_parents_q30_biSNPs_mndp8_90CIdp.vcf.gz') # read in the coverage filtered VCF file we just produced
vcf2 # make sure the VCF qualities are the same as we left off
vcf2<-filter_allele_balance(vcfR, min.ratio = 0.33, max.ratio = 0.66) #Impose a 33% minor allele and 66% major allele count read depth filter to heterozygous genotype calls. Set genotypes with allele balances outside this ratio to missing.
vcf2
write.vcf(vcf2, file = 'NAM_parents_q30_biSNPs_mndp8_90CIdp_33_66hetdp.vcf.gz')






