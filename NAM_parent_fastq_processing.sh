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


# Align quality + adapter trimmed paired end reads to the indexed and masked reference Cascade genome using BWA MEM and pipe output to samtools view to convert to a .bam file#
for f in *_R1_trim.fastq.gz
r1=$f
r2=4{f/_R1_trim.fastq.gz}_R2_trim.fastq.gz
r3=${f/_R1_trim.fastq.gz}.bam
bwa mem -t 10 dovetailCascadeFullAssemblyMasked.fasta $r1 $r2 | samtools view -@10  -bS - > $r3
done

# Sort each read within the .bam files by name (required for marking PCR duplicates) using samtools #
for f in *.bam
r1=$f
r2=${f/.bam}_nsort.bam
do
samtools sort -@ 11 -n $r1 -o $r2
done
