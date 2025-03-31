#!/bin/bash

set -e # Exit immediately if a command exits with a non-zero status

#------------------CONFIGURATION-----------------
THREADS=8
RAW_READS="/home/Desktop/test/raw_reads"
TRIMMED="/home/Desktop/test/trimmed_reads"
BAM_DIR="/home/Desktop/test/bam/"
DEDUP_BAM="/home/Desktop/test/picard/"
REFERENCE="/home/Desktop/test/reference"
ANNOTATION="/home/Desktop/test/snpEff"
LOGFILE="analysis.log"
exec > "$LOGFILE" 2>&1


#-------------Check if tools are installed-------
for tool in cutadapt bowtie2 samtools gatk java bcftools ; do
  if ! command -v $tool &> /dev/null; then
    echo "$tool could not be found. Exiting."
    exit 1
  fi
done

#-------create the necessary directories--------- 
mkdir -p "$TRIMMED" "$BAM_DIR" "$DEDUP_BAM"

#-----------------FILES--------------------------
Tumor_Raw_R1="$RAW_READS/Tumor_R1.fastq.gz"
Tumor_Raw_R2="$RAW_READS/Tumor_R2.fastq.gz"
Normal_Raw_R1="$RAW_READS/Norm_R1.fastq.gz"
Normal_Raw_R2="$RAW_READS/Norm_R2.fastq.gz"

ref="$REFERENCE/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

#---------Pre-processing of raw reads------------
## Please change according to the adapters present in you data and set the read length filter parameter accordingly based on the quality check
## The raw reads passed the read length filter and trimmed data sets have a distribution of read lengths. The adapters present in the dataset were Trueseq adapters.

echo "Trimming tumor sample " | tee -a "$LOGFILE"
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -m 50 -j "$THREADS" \
-o "$TRIMMED/Tumor_R1.fastq" -p "$TRIMMED/Tumor_R2.fastq" \
"$Tumor_Raw_R1" "$Tumor_Raw_R2"
echo "completed preprocessing of tumor sample"

echo "Trimming normal sample " | tee -a "$LOGFILE"
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -m 50 -j "$THREADS" \
-o "$TRIMMED/Normal_R1.fastq" -p "$TRIMMED/Normal_R2.fastq" \
"$Normal_Raw_R1" "$Normal_Raw_R2"
echo "completed preprocessing of normal sample"

#-------------reference genome indexing----------
bowtie2-build --threads $THREADS "$ref" hg38

#--------------------alignment-------------------
# replase with correct lane_id, sample_name, library, platform and instrument ID and platform unit 

bowtie2 --threads $THREADS --rg-id tumor_MHG12345-1789_lane1 --rg SM:tumor --rg LB:lib2 --rg PL:ILLUMINA --rg PU:DS10000000_11_MHG12345-1789_1_tumor \
-x /home/Desktop/test/hg38 -1 "$TRIMMED/Tumor_R1.fastq" -2 "$TRIMMED/Tumor_R2.fastq" \
-S "$BAM_DIR/Tum_var.sam"
echo "completed assembly of tumor sample"

bowtie2 --threads $THREADS --rg-id normal_MHG12345-1789_lane1 --rg SM:normal --rg LB:lib2 --rg PL:ILLUMINA --rg PU:DS10000000_11_MHG12345-1789_1_normal \
-x /home/Desktop/test/hg38 -1 "$TRIMMED/Normal_R1.fastq" -2 "$TRIMMED/Normal_R2.fastq" \
-S "$BAM_DIR/Norm_var.sam"
echo "completed assembly of normal sample"


#------Bam creation, sorting and indexing--------
for sample in Tum Norm; do
  samtools view -bS "$BAM_DIR/${sample}_var.sam" -o "$BAM_DIR/${sample}_var.bam"
  samtools sort "$BAM_DIR/${sample}_var.bam" -o "$BAM_DIR/${sample}_var.sorted.bam"
  samtools index "$BAM_DIR/${sample}_var.sorted.bam"
done

#---------Duplicate removal using picard---------
for sample in Tum Norm; do
 java -jar picard.jar MarkDuplicates \
  I="$BAM_DIR/${sample}_var.sorted.bam" \
  O="$DEDUP_BAM/${sample}_var_markdup.bam" \
  M="$DEDUP_BAM/marked_dup_metrics_${sample}_var.txt" \
  REMOVE_DUPLICATES=true \
  REMOVE_SEQUENCING_DUPLICATES=true
 samtools index "$DEDUP_BAM/${sample}_var_markdup.bam"
done 

#----------reference_Dictionary------------------
if [ ! -f $ref.fai ]; then
 samtools fqidx "$ref" 
fi

gatk CreateSequenceDictionary -R "$ref"

#----------------variant calling-----------------
gatk HaplotypeCaller  -R "$ref" -I "$DEDUP_BAM/Norm_var_markdup.bam" -O Normal.vcf.gz
gatk VariantFiltration -R "$ref" -V Normal.vcf.gz -O Germline.vcf.gz -L fil.bed --variant-output-filtering CONTAINED --filter-name "LowQual" --filter-expression "QUAL < 30.0 || DP < 10"
bcftools view -f PASS Germline.vcf.gz -o Germline_PASS.vcf.gz -O z
gatk IndexFeatureFile -I Germline_PASS.vcf.gz
gatk VariantAnnotator -R "$ref" -V Germline_PASS.vcf.gz -O germline_resource.vcf.gz --annotation AlleleFraction
gatk IndexFeatureFile -I germline_resource.vcf.gz
 
if [ ! -f fil.bed ]; then
 echo "bed file not found. Exiting.." 
exit 1
fi

gatk Mutect2 -R "$ref" -I "$DEDUP_BAM/Norm_var_markdup.bam" -tumor normal -L fil.bed --germline-resource germline_resource.vcf.gz -O normal_pon.vcf.gz --max-mnp-distance 0
gatk CreateSomaticPanelOfNormals -V normal_pon.vcf.gz -O pon.vcf.gz
gatk Mutect2 -R "$REFERENCE/Homo_sapiens.GRCh38.dna.primary_assembly.fa" -I "$DEDUP_BAM/Tum_var_markdup.bam" -I "$DEDUP_BAM/Norm_var_markdup.bam" -tumor tumor -normal normal -L fil.bed --germline-resource germline_resource.vcf.gz --panel-of-normals pon.vcf.gz -O somatic.vcf.gz
gatk FilterMutectCalls -R "$ref" -V somatic.vcf.gz --max-events-in-region 3 --max-alt-allele-count 2 -O somatic_filtered.vcf.gz
bcftools view -f PASS somatic_filtered.vcf.gz -o somatic_only.vcf.gz -O z

#-------------variant annotation-----------------
# update java to latest version

java -jar $ANNOTATION/snpEff.jar download GRCh38.99
gunzip somatic_only.vcf.gz
java -Xmx8g -jar $ANNOTATION/snpEff.jar "$ANNOTATION/GRCh38.99" somatic_only.vcf > somatic_ann.vcf

#-------------mutation rate calc-----------------

# Tumor mutation burden was calculated across total bases of targeted regions.
total_bases=$(cat fil.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2-1 }END{print SUM}')
total_kb=$(echo "$total_bases / 1000" | bc -l)
total_mutations=$(bcftools view somatic_only.vcf| grep -v "^#" | wc -l)
somatic_mutation_burden=$(echo "scale=6; $total_mutations / $total_kb" | bc)

echo "Total Somatic Mutations : $total_mutations"
echo "Tumor Mutation Burden (TMB): $somatic_mutation_burden mutations/kb
