#!/bin/bash -l
#SBATCH --output=/path/to/output/%u/%j.out
#SBATCH --job-name=RNAseq_pipeline_BOH
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=72:00:00

##########################################################################################################
# RNA-seq Data Analysis Workflow: 50-1208954829 GeneWiz (Cardiomyocytes + SARS-COV2 alpha variant)
# Author: Béibhinn O'Hora (2025)
# Illumina PE 2x150bp
#
# Steps:
#   1. QC 
#   2. Trimming 
#   3. Reference genome setup (GENCODE / GRCh38)
#   4. STAR genome index building
#   5. STAR alignment
#   6. Read quantification 
#   7. Merge all count tables
#   8. Add gene names to final matrix
#
# Modules required:
#   graalvm, fastqc, MultiQC, Trimmomatic, STAR, Subread
##########################################################################################################

# Note: replace full pathnames based on your preferred area of data storage

set -e  
set -x  

#--------------------------------------------#
# 0. Load modules
#--------------------------------------------#

# Modules installed on Maestro cluster hosted by Institut Pasteur
module load graalvm
module load fastqc/0.12.1
module load MultiQC/1.12
module load Trimmomatic
module load STAR/2.7.11b
module load Subread/2.0.6

#--------------------------------------------#
# 1. Quality control (FastQC + MultiQC)
#--------------------------------------------#

RAW_FASTQ_DIR="/path/to/raw_fastq"
FASTQC_DIR="/path/to/output/00_fastqc"
MULTIQC_DIR="/path/to/output/01_multiqc"

mkdir -p "$FASTQC_DIR" "$MULTIQC_DIR"
cd "$RAW_FASTQ_DIR"

fastqc *.fastq.gz --outdir "$FASTQC_DIR"
multiqc "$FASTQC_DIR" -o "$MULTIQC_DIR"

#--------------------------------------------#
# 2. Trimming (Trimmomatic)
#--------------------------------------------#

TRIMMED_DIR="/path/to/output/02_trimmed"
mkdir -p "$TRIMMED_DIR"
cd "$RAW_FASTQ_DIR"

for R1 in *_R1_001.fastq.gz; do
    SAMPLE="${R1%_R1_001.fastq.gz}"
    R2="${SAMPLE}_R2_001.fastq.gz"

    if [[ ! -f "$R2" ]]; then
        echo "Skipping $SAMPLE: R2 file not found."
        continue
    fi

    R1_PAIRED="${TRIMMED_DIR}/${SAMPLE}_R1_paired.fastq.gz"
    R1_UNPAIRED="${TRIMMED_DIR}/${SAMPLE}_R1_unpaired.fastq.gz"
    R2_PAIRED="${TRIMMED_DIR}/${SAMPLE}_R2_paired.fastq.gz"
    R2_UNPAIRED="${TRIMMED_DIR}/${SAMPLE}_R2_unpaired.fastq.gz"

    Trimmomatic PE -threads 8 \
        "$R1" "$R2" \
        "$R1_PAIRED" "$R1_UNPAIRED" \
        "$R2_PAIRED" "$R2_UNPAIRED" \
        ILLUMINACLIP:"$TRIMMOMATIC_ADAPTERS/TruSeq3-PE.fa":2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
done

#--------------------------------------------#
# 3. Reference genome setup (GENCODE)
#--------------------------------------------#

REF_DIR="/path/to/reference/GENCODE_GRCh38"
mkdir -p "$REF_DIR"
cd "$REF_DIR"

if [[ ! -f "GRCh38.primary_assembly.genome.fa" ]]; then
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
    gunzip GRCh38.primary_assembly.genome.fa.gz
fi
if [[ ! -f "gencode.v44.annotation.gtf" ]]; then
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
    gunzip gencode.v44.annotation.gtf.gz
fi

#--------------------------------------------#
# 4. Build STAR index
#--------------------------------------------#

STAR_INDEX_DIR="/path/to/reference/STAR_index_GRCh38"
mkdir -p "$STAR_INDEX_DIR"

if [[ ! -f "$STAR_INDEX_DIR/SA" ]]; then
    STAR --runThreadN 16 \
         --runMode genomeGenerate \
         --genomeDir "$STAR_INDEX_DIR" \
         --genomeFastaFiles "$REF_DIR/GRCh38.primary_assembly.genome.fa" \
         --sjdbGTFfile "$REF_DIR/gencode.v44.annotation.gtf" \
         --sjdbOverhang 149
fi

#--------------------------------------------#
# 5. STAR alignment
#--------------------------------------------#

ALIGN_DIR="/path/to/output/03_star_align"
mkdir -p "$ALIGN_DIR"
cd "$TRIMMED_DIR"

for R1 in *_R1_paired.fastq.gz; do
    SAMPLE=$(basename "$R1" _R1_paired.fastq.gz)
    R2="${SAMPLE}_R2_paired.fastq.gz"

    STAR --runThreadN 8 \
         --genomeDir "$STAR_INDEX_DIR" \
         --readFilesIn "$R1" "$R2" \
         --readFilesCommand zcat \
         --outFileNamePrefix "${ALIGN_DIR}/${SAMPLE}_" \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts
done

#--------------------------------------------#
# 6. Quantification (featureCounts)
#--------------------------------------------#

COUNT_DIR="/path/to/output/04_featureCounts"
mkdir -p "$COUNT_DIR"
cd "$ALIGN_DIR"

featureCounts -T 8 -p -s 0 \
  -a "$REF_DIR/gencode.v44.annotation.gtf" \
  -o "${COUNT_DIR}/gene_counts_raw.txt" *_Aligned.sortedByCoord.out.bam

#--------------------------------------------#
# 7. Merge per-sample count files
#--------------------------------------------#

MERGE_DIR="/path/to/output/05_merged_counts"
mkdir -p "$MERGE_DIR"
cd "$COUNT_DIR"

files=$(ls *_counts.txt | sort)
first_file=$(echo "$files" | head -n1)
cut -f1 "$first_file" > "$MERGE_DIR/merged_counts_matrix.tsv"

for f in $files; do
    cut -f2 "$f" > "$MERGE_DIR/tmp_counts.txt"
    paste "$MERGE_DIR/merged_counts_matrix.tsv" "$MERGE_DIR/tmp_counts.txt" > "$MERGE_DIR/tmp_merged.tsv"
    mv "$MERGE_DIR/tmp_merged.tsv" "$MERGE_DIR/merged_counts_matrix.tsv"
done

# Add header
header="gene_id"
for f in $files; do
    sample=$(basename "$f" _counts.txt)
    header="${header}\t${sample}"
done
echo -e "$header" | cat - "$MERGE_DIR/merged_counts_matrix.tsv" > "$MERGE_DIR/merged_counts_matrix_with_header.tsv"

rm -f "$MERGE_DIR/tmp_counts.txt" "$MERGE_DIR/merged_counts_matrix.tsv"

#--------------------------------------------#
# 8. Add gene names (from GTF)
#--------------------------------------------#

awk -F'\t' '$3=="gene" {
    match($9, /gene_id "([^"]+)";/, id)
    match($9, /gene_name "([^"]+)";/, name)
    if (id[1] != "" && name[1] != "")
        print id[1] "\t" name[1]
}' "$REF_DIR/gencode.v44.annotation.gtf" > "$MERGE_DIR/gene_id_to_name.tsv"

sort -k1,1 "$MERGE_DIR/merged_counts_matrix_with_header.tsv" > "$MERGE_DIR/merged_counts_sorted.tsv"
sort -k1,1 "$MERGE_DIR/gene_id_to_name.tsv" > "$MERGE_DIR/gene_id_to_name_sorted.tsv"

tail -n +2 "$MERGE_DIR/merged_counts_sorted.tsv" > "$MERGE_DIR/tmp_body.tsv"
head -1 "$MERGE_DIR/merged_counts_sorted.tsv" > "$MERGE_DIR/header.tmp"

join -t $'\t' -1 1 -2 1 "$MERGE_DIR/tmp_body.tsv" "$MERGE_DIR/gene_id_to_name_sorted.tsv" > "$MERGE_DIR/counts_with_symbols_tmp.tsv"

echo -e "$(cat "$MERGE_DIR/header.tmp")\tgene_name" > "$MERGE_DIR/header_final.tmp"
cat "$MERGE_DIR/header_final.tmp" "$MERGE_DIR/counts_with_symbols_tmp.tsv" > "$MERGE_DIR/final_counts_matrix_with_symbols.tsv"

rm -f "$MERGE_DIR/tmp_body.tsv" "$MERGE_DIR/header.tmp" "$MERGE_DIR/header_final.tmp" \
      "$MERGE_DIR/counts_with_symbols_tmp.tsv" "$MERGE_DIR/gene_id_to_name.tsv" \
      "$MERGE_DIR/gene_id_to_name_sorted.tsv" "$MERGE_DIR/merged_counts_sorted.tsv"

#--------------------------------------------#
# 9. Final QC summary
#--------------------------------------------#

FINAL_QC_DIR="/path/to/output/06_final_multiqc"
mkdir -p "$FINAL_QC_DIR"
multiqc "$ALIGN_DIR" "$COUNT_DIR" -o "$FINAL_QC_DIR"

