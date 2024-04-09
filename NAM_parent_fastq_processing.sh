# Script used for processing raw fastq.gz files containing the whole 
genome sequencing data for each hop NAM parent. #

# Rename files by shortening their length #
for f in *.fastq.gz
do
echo mv -v "$f" "${f#*P001_}"
done

# Inspect raw reads for quality and adapter contamination using 'fastqc' #
for f in *.fastqc.gz
do
fastqc -t7 $f
done

# Trim sequences with poor quality and/or adapter contamination using 
'fastp' #
for f in *_L003-L004_R1_001.fastq.gz
do
r1=$f
r2=${f/_L003-L004_R1_001.fastq.gz}_L003-L004_R2_001.fastq.gz
r3=${f/_L003-L004_R1_001.fastq.gz}_R1_trim.fastq.gz
r4=${f/_L003-L004_R1_001.fastq.gz}_r2_trim.fastq.gz
r5=${f/_L003-L004_R1_001.fastq.gz}_report.html
echo  "fastp --thread 8 -p -l 36 -g --cut_front --cut_tail -y --detect_adapter_for_pe --adapter_fasta adapters.fa -h $r5 -i $f -I $r2 -o $r3 -O $r4"
done

# Use fastqc as a secondary check to make sure parent fastq.gz files look good for alignment #
for f in *.fastqc.gz
do
fastqc -t7 $f
done

# Align quality + adapter trimmed paired end reads to the indexed and masked reference Cascade genome using BWA MEM and pipe output to samtools view to convert to a .bam file #
for f in *_R1_trim.fastq.gz
do
r1=$f
r2=4{f/_R1_trim.fastq.gz}_R2_trim.fastq.gz
r3=${f/_R1_trim.fastq.gz}.bam
bwa mem -t 10 dovetailCascadeFullAssemblyMasked.fasta $r1 $r2 | samtools view -@10  -bS - > $r3
done

# Sort each read within the .bam files by name (required for marking PCR duplicates) using samtools sort #
for f in *.bam
do
r1=$f
r2=${f/.bam}_nsort.bam
samtools sort -@ 11 -n $r1 -o $r2
done

# Fill in various coordinates and flags from name-sorted .bam files with samtools fixmate #
for f in *_nsort.bam
do
r1=$f
r2=${f/_nsort.bam}_fx_nsort.bam
samtools fixmate -@ 11 -m $r1 $r2
done

# Sort alignments in .bam files by physical position using samtools sort #
for f in *_fx_nsort.bam
do
r1=$f
r2=${f/_fx_nsort.bam}_pos_fx_nsort.bam
samtools sort -@ 11 $r1 -o $r2
done

# Mark and remove PCR duplicates using samtools markdup #
for f in *_pos_fx_nsort.bam
do
r1=$f
r2=${f/_pos_fx_nsort.bam}_md_pos_fx_nsort.bam
samtools markdup -@ 11 $r1 $r2
done

# Apply a final filter to .bam alignment to remove alignements with a mapping quality less than 10 and keep only properly paired alignments using samtools view #
for f in *_md_pos_fx_nsort.bam
do
r1=$f
r2=${f/_md_pos_fx_nsort.bam}prop.bam
samtools view -q 10 -f 0x002 -b $r1 > $r2
done

# Index all final filtered alignment .bam files using samtools index #
for f in *prop.bam
do
samtools index $f
done

# Now call variants separately for each Hop chromosome (1-10) across all filtered .bam files using bcftools mpileup and pipe the output to bcftools call #
# Start a separate screen session for each Hop chromosome to call variants on them individually. #
# Each Cascade_chrN.txt file contains 3 columns that specify the Chromosome ID, Chromosome Start Position, and Chromosome End Position #

ls *prop.bam > NAM_parents.txt #create a simple text file with a single column containing the name of all final filtered .bam files for each parent.


screen

bcftools mpileup -R Cascade_chr1.txt -a FORMAT/DP -a FORMAT/SP -a FORMAT/AD -a INFO/AD -Ou -A -f dovetailCascadeFullAssemblyMasked.fasta -b NAM_parents.txt | bcftools call -R 
Cascade_chr1.txt  -vmO z -f GQ,GP -o 
NAM_parents_chr1.vcf.gz

screen

bcftools mpileup -R Cascade_chr2.txt -a FORMAT/DP -a FORMAT/SP -a FORMAT/AD -a INFO/AD -Ou -A -f dovetailCascadeFullAssemblyMasked.fasta -b NAM_parents.txt | bcftools call -R Cascade_chr2.txt -vmO z -f GQ,GP -o 
NAM_parents_chr2.vcf.gz

screen

bcftools mpileup -R Cascade_chr3.txt -a FORMAT/DP -a FORMAT/SP -a FORMAT/AD -a INFO/AD -Ou -A -f dovetailCascadeFullAssemblyMasked.fasta -b NAM_parents.txt | bcftools call -R Cascade_chr3.txt -vmO z -f GQ,GP -o 
NAM_parents_chr3.vcf.gz

screen

bcftools mpileup -R Cascade_chr4.txt -a FORMAT/DP -a FORMAT/SP -a FORMAT/AD -a INFO/AD -Ou -A -f dovetailCascadeFullAssemblyMasked.fasta -b NAM_parents.txt | bcftools call -R Cascade_chr4.txt -vmO z -f GQ,GP -o 
NAM_parents_chr4.vcf.gz

screen

bcftools mpileup -R Cascade_chr5.txt -a FORMAT/DP -a FORMAT/SP -a FORMAT/AD -a INFO/AD -Ou -A -f dovetailCascadeFullAssemblyMasked.fasta -b NAM_parents.txt | bcftools call -R Cascade_chr5.txt -vmO z -f GQ,GP -o 
NAM_parents_chr5.vcf.gz

screen

bcftools mpileup -R Cascade_chr6.txt -a FORMAT/DP -a FORMAT/SP -a FORMAT/AD -a INFO/AD -Ou -A -f dovetailCascadeFullAssemblyMasked.fasta -b NAM_parents.txt | bcftools call -R Cascade_chr6.txt -vmO z -f GQ,GP -o 
NAM_parents_chr6.vcf.gz

screen

bcftools mpileup -R Cascade_chr7.txt -a FORMAT/DP -a FORMAT/SP -a FORMAT/AD -a INFO/AD -Ou -A -f dovetailCascadeFullAssemblyMasked.fasta -b NAM_parents.txt | bcftools call -R Cascade_chr7.txt -vmO z -f GQ,GP -o 
NAM_parents_chr7.vcf.gz

screen

bcftools mpileup -R Cascade_chr8.txt -a FORMAT/DP -a FORMAT/SP -a FORMAT/AD -a INFO/AD -Ou -A -f dovetailCascadeFullAssemblyMasked.fasta -b NAM_parents.txt | bcftools call -R Cascade_chr8.txt -vmO z -f GQ,GP -o 
NAM_parents_chr8.vcf.gz

screen

bcftools mpileup -R Cascade_chr9.txt -a FORMAT/DP -a FORMAT/SP -a FORMAT/AD -a INFO/AD -Ou -A -f dovetailCascadeFullAssemblyMasked.fasta -b NAM_parents.txt | bcftools call -R Cascade_chr9.txt -vmO z -f GQ,GP -o 
NAM_parents_chr9.vcf.gz

screen

bcftools mpileup -R Cascade_chrX.txt -a FORMAT/DP -a FORMAT/SP -a FORMAT/AD -a INFO/AD -Ou -A -f dovetailCascadeFullAssemblyMasked.fasta -b NAM_parents.txt | bcftools call -R Cascade_chrX.txt -vmO z -f GQ,GP -o 
NAM_parents_chrX.vcf.gz

# Concatenate all raw VCF files to create one final raw VCF for Chromosomes 1-10 using bcftools concat #
bcftools concat -Ovz NAM_parents_chr1.vcf.gz NAM_parents_chr2.vcf.gz NAM_parents_chr3.vcf.gz NAM_parents_chr4.vcf.gz NAM_parents_chr5.vcf.gz NAM_parents_chr6.vcf.gz NAM_parents_chr7.vcf.gz NAM_parents_chr8.vcf.gz 
NAM_parents_chr9.vcf.gz NAM_parents_chrX.vcf.gz -o NAM_parents_chr1_X.vcf.gz






