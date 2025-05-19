#!/bin/bash
# methylation_analysis.sh

# Comprehensive methylation analysis pipeline for epibac-heritome
# Here are the steps for analyzing the methylated CpGs from the aligned BAM files obtained from Nanopore.

set -euo pipefail

# Configuration
REF_ASSEMBLY="acropora_ref_assembly.fa"
REF_GFF="acropora_ref_assembly.gff"

# Generate traditional pileup BED files using a standradized threshold for each sample
# This will generate the bedMethyl files tabulate the counts of base modifications from every sequencing
# read over each aligned reference genomic position.
# modkit v0.4.4 (https://nanoporetech.github.io/modkit)
for BAM_FILE in *.aligned.bam; do
    SAMPLE_NAME=$(basename "${BAM_FILE}" .aligned.bam)
    PILEUP_BED="${SAMPLE_NAME}.traditional.pileup.bed"
    LOG_FILE="${SAMPLE_NAME}_traditional_pileup.log"
    modkit pileup "${BAM_FILE}" "${PILEUP_BED}" --ref "${REF_ASSEMBLY}" --log-filepath "${LOG_FILE}" --preset traditional --filter-threshold 0.95 --only-tabs
done

# Filter pileup files for coverage >= 10
# Minimum read coverage parameters were s based on Dimond et al., 2021 (https://doi.org/10.1093/g3journal/jkab148)
for PILEUP in *.traditional.pileup.bed; do
    SAMPLE_NAME=$(basename "${PILEUP}" .traditional.pileup.bed)
    awk '$5 >= 10' "${PILEUP}" > "${SAMPLE_NAME}.10x.pileup.bed"
done

# Sort filtered pileup files
# bedtools v2.30.0 (https://bedtools.readthedocs.io/en/latest/index.html)
for BED in *.10x.pileup.bed; do
    sort -k1,1 -k2,2n "${BED}" > "${BED}.sorted"
done

# Filter so just the CpG positions with 10x coverage in across all samples remain
bedtools intersect -a "${BED_FILES[0]}" "${BED_FILES[@]:1/#/-b }" -c > common_positions_w_counts.bed
awk -v num_files=11 '$NF == num_files' common_positions_w_counts.bed > common_positions.bed

# ---- Use this data to generate the files for NanoMethViz ----

# Extract methylation data at common positions
# The output from this step calculates the methylation probabilities for each position in each read
# and will be used for analyzing methylation using NanoMethViz
for BAM_FILE in *.bam; do
    SAMPLE_NAME=$(basename "${BAM_FILE}" .bam)
    modkit extract full --include-bed common_positions.bed --ref "${REF_ASSEMBLY}" "${BAM_FILE}" "${SAMPLE_NAME}.modkit_extract.tsv"
done
gzip *.modkit_extract.tsv

# Add sample name column to each file (because this will be lost when combining in the next step)
for FILE in *.modkit_extract.tsv.gz; do
    SAMPLE=$(basename $FILE .modkit_extract.tsv.gz)
    zcat $FILE | awk -v sample=$SAMPLE 'BEGIN {OFS="\t"}
        NR == 1 {print "sample", $0}
        NR > 1 {print sample, $0}' > ${FILE%.tsv.gz}_samp.tsv
done

# Combine all sample files (without headers)
for FILE in *_samp.tsv; do
    tail -n +2 "$FILE" >> combined_data_tmp.tsv
done

# Sort and create tabix file for NanoMethViz
# samtools v1.16.1 (https://github.com/samtools/samtools)
awk -F'\t' -v OFS='\t' '{ print $1, $5, $4, $6, $14, $2 }' combined_data_tmp.tsv > combined_data.tsv
sort -k2,2 -k3,3n combined_data.tsv > combined_data.sorted.tsv
bgzip combined_data.sorted.tsv
tabix -s 2 -b 3 -e 3 combined_data.sorted.tsv.gz

# ---- Return to methylation pipeline ----

# Filter individual pileup files by common positions
for BED in *.10x.pileup.bed.sorted; do
    SAMPLE_NAME=$(basename "${BED}" .10x.pileup.bed.sorted)
    bedtools intersect -a ${BED} -b common_positions.bed > ${SAMPLE_NAME}.filtered.pileup.bed
done

# ---- The following steps use python scripts found here: https://github.com/arbarno/epibac_heritome ----

# Merge filtered pileup files to create a single file with compiled methylation data
# python v3.11.0
ALL_FILTERED=$(echo *.filtered.pileup.bed)
merge_bed_files.py ${ALL_FILTERED} -v > merged_bed_output.bed

# Annotate merged BED file
# python v3.11.0
annotate_merged_bed.py ${REF_ASSEMBLY} ${REF_GFF} merged_bed_output.bed > merged_bed_annotated.bed
cut -f7 merged_bed_annotated.bed | sort -u > gene_universe.txt

# Add annotation to filtered pileup files
cut -f 7- merged_bed_annotated.bed > tmp
for BED in *filtered.pileup.bed; do
    SAMPLE_NAME=$(basename "${BED}" .filtered.pileup.bed)
    cut -f1,2,3,11,12,10 ${BED} | paste - tmp > ${SAMPLE_NAME}_filtered_annotated.bed
done

# Calculate CpG bias
# python v3.11.0
calc_cpg_bias.py ${REF_ASSEMBLY} ${REF_GFF} > /ibex/project/c2208/nanopore/output/acropora_ref_cpg_bias.tsv

# Calculate mean methylation per sample
# python v3.11.0
for BED in *_filtered_annotated.bed; do
    SAMPLE_NAME=$(basename "${BED}" _filtered_annotated.bed)
    calc_mean_methylation.py ${REF_ASSEMBLY} ${REF_GFF} ${BED} > ${SAMPLE_NAME}_mean_meths.tsv
done

# Tabulate all mean methylation files
# python v3.11.0
ALL_FILTERED=$(echo *filtered_annotated.bed)
tabulate_tsvs.py ${ALL_FILTERED} -k 0 1 -c 3 -v > all_mean_meths.tsv

# Tabulate filtered context percentages
# python v3.11.0
ALL_FILTERED=$(echo *filtered_annotated.bed)
tabulate_tsvs.py ${ALL_FILTERED} -k 0 1 6 10 -c 5 -v > all_filt_pct_context.tsv
sed -i '' 's/^\t\t/scaffold\tpos\tgene\tcontext/' all_filt_pct_context.tsv
sed -i '' 's/\.bed//g' all_filt_pct_context.tsv
gzip all_filt_pct_context.tsv
