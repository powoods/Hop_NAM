# Script for preparing the Cascade reference genome for alignment #

# Index the reference genome using BWA #
bwa index dovetailCascadeFullAssemblyMasked.fasta

# Align quality trimmed reads to the index reference genome #
for f in *_R1_trim.fastqc.gz
do
r1=$f
r2=${f/_R1_trim.fastq.gz}_R2_trim.fastq.gz
r3=${f/_R1_trim.fastq.gz}.sam
bwa mem -t 10 -M $r1 $r2 > $r3
done
