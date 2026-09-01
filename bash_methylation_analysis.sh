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

# Change the second position of the bed file to look for conserved dinucleotides (CG) in the next step
for PILEUP in *.10x.pileup.bed; do
    SAMPLE_NAME=$(basename "${PILEUP}" .10x.pileup.bed)
    awk -F'\t' 'BEGIN{OFS="\t"} {$3=$2+2; print}' "${PILEUP}" > "${SAMPLE_NAME}.10x.cpgs.bed"
done

# Filter out positions that are not conserved with clair3
# clair3 v2.0.2 (https://github.com/HKU-BAL/clair3)
conda activate clair3
for BED in *.10x.cpgs.bed; do
    SAMPLE_NAME=$(basename "${BED}" .10x.cpgs.bed)
    run_clair3.sh \
        --bam_fn="${SAMPLE_NAME}.aligned.bam" \
        --ref_fn="${REF_ASSEMBLY}" \
        --threads=16 \
        --platform=ont \
        --model_path="${CONDA_PREFIX}/bin/models/r1041_e82_400bps_sup_v420" \
        --include_all_ctgs \
        --bed_fn="${BED}" \
        --output="${SAMPLE_NAME}_clair3"
done
conda deactivate

# Filter out CpG variants -> keep the CpG positions that don't PASS as variants
# bedtools/2.30.0 bcftools/1.16
for BED in *.10x.cpgs.bed; do
    SAMPLE_NAME=$(basename "${BED}" .10x.cpgs.bed)
    bcftools view -f PASS "${SAMPLE_NAME}_clair3"/merge_output.vcf.gz > clair3_pass/"${SAMPLE_NAME}_pass.vcf"
    bedtools intersect -a "${BED}" -b clair3_pass/"${SAMPLE_NAME}_pass.vcf" -v > 10x_clair3_filt/"${SAMPLE_NAME}.10x.clair3.filt.bed"
done

# Go back to just the C positions (from CpGs)
cd 10x_clair3_filt
for BED in *.10x.clair3.filt.bed; do
    SAMPLE_NAME=$(basename "${BED}" .10x.clair3.filt.bed)
    awk -F'\t' 'BEGIN{OFS="\t"} {$3=$2+1; print}' "${BED}" > ${SAMPLE_NAME}.10x.clair3.filt.Cs.bed
done

# Sort filtered clair3 passing files
# bedtools v2.30.0
cd 10x_clair3_filt
for BED in *.10x.clair3.filt.Cs.bed; do
    sort -k1,1 -k2,2n "${BED}" > "${BED}.sorted"
done

# Filter each sample to have conserved CpG positions across all samples
SORTED_BEDS=(*.10x.clair3.filt.Cs.bed.sorted)
NUM_SAMPLES=${#SORTED_BEDS[@]}
bedtools multiinter -i "${SORTED_BEDS[@]}" > multiinter_output.bed
awk -v n="${NUM_SAMPLES}" '$4 == n {print $1"\t"$2"\t"$3}' multiinter_output.bed | sort -k1,1 -k2,2n -u > common_positions_clair3.bed

# ---- Use this data to generate the files for NanoMethViz ----

# Extract methylation data at common positions
# The output from this step calculates the methylation probabilities for each position in each read
# and will be used for analyzing methylation using NanoMethViz
for BAM_FILE in *.bam; do
    SAMPLE_NAME=$(basename "${BAM_FILE}" .bam)
    modkit extract full --include-bed common_positions_clair3.bed --ref "${REF_ASSEMBLY}" "${BAM_FILE}" "${SAMPLE_NAME}.modkit_extract.tsv"
done
gzip *.modkit_extract.tsv

# Add sample name column to each file (because this will be lost when combining in the next step)
for FILE in *.modkit_extract.tsv.gz; do
    SAMPLE=$(basename $FILE .modkit_extract.tsv.gz)
    zcat $FILE | awk -v sample=$SAMPLE 'BEGIN {OFS="\t"}
        NR == 1 {print "sample", $0}
        NR > 1 {print sample, $0}' > ${FILE%.tsv.gz}_samp.tsv
done

# Extract the header from the first file and save it to the combined file
head -n 1 $(ls *_samp.tsv | head -n 1) > combined_data_tmp.tsv

# Combine all sample files (without headers)
for FILE in *_samp.tsv; do
    tail -n +2 "$FILE" >> combined_data_tmp.tsv
done

# Sort and create tabix file for NanoMethViz
# samtools v1.16.1 (https://github.com/samtools/samtools)
awk -F'\t' -v OFS='\t' '{ print $1, $5, $4, $6, $14, $2 }' combined_data_tmp.tsv > combined_data_clair3.tsv
sort -k2,2 -k3,3n combined_data_clair3.tsv > combined_data_clair3.sorted.tsv
bgzip combined_data_clair3.sorted.tsv
tabix -s 2 -b 3 -e 3 combined_data_clair3.sorted.tsv.gz

# ---- Return to methylation pipeline ----

# Filter individual pileup files by common positions
cd 10x_clair3_filt
for BED in *.clair3.filt.Cs.bed.sorted; do
    SAMPLE_NAME=$(basename "${BED}" .clair3.filt.Cs.bed.sorted)
    bedtools intersect -a ${BED} -b common_positions_clair3.bed > ${SAMPLE_NAME}.clair3.filt.conservedCpGs.bed
done

# ---- The following steps use python scripts found here: https://github.com/arbarno/epibac_heritome/tree/main/python_scripts ----

# Merge filtered pileup files to create a single file with compiled methylation data
# python v3.11.0
cd 10x_clair3_filt
ALL_FILTERED=$(echo *.clair3.filt.conservedCpGs.bed)
merge_bed_files.py ${ALL_FILTERED} -v > merged_clair3_bed_output.bed

# Annotate merged BED file
# python v3.11.0
cd 10x_clair3_filt
annotate_merged_bed.py ${REF_ASSEMBLY} ${REF_GFF} merged_clair3_bed_output.bed > merged_clair3_bed_annotated.bed
cut -f7 merged_clair3_bed_annotated.bed | sort -u > clair3_gene_universe.txt

# Add annotation to filtered pileup files
cd 10x_clair3_filt
cut -f 7- merged_clair3_bed_annotated.bed > tmp
for BED in *.clair3.filt.conservedCpGs.bed; do
    SAMPLE_NAME=$(basename "${BED}" .clair3.filt.conservedCpGs.bed)
    cut -f1,2,3,11,12,10 ${BED} | paste - tmp > ${SAMPLE_NAME}.clair3.filt.conservedCpGs.annotated.bed
done

# Calculate CpG bias
# python v3.11.0
calc_cpg_bias.py ${REF_ASSEMBLY} ${REF_GFF} > /ibex/project/c2208/nanopore/output/acropora_ref_cpg_bias.tsv

# Aggregate methylation data
# python v3.11.0
cd 10x_clair3_filt
ALL_FILTERED=$(echo *.10x.clair3.filt.conservedCpGs.annotated.bed)
python3 aggregate_gene_context_counts.py ${ALL_FILTERED} --out-gene-tsv gene_level_counts.tsv --out-context-tsv context_level_counts.tsv

# Tabulate filtered context percentages
# python v3.11.0
ALL_FILTERED=$(echo *.clair3.filt.conservedCpGs.annotated.bed)
tabulate_tsvs.py ${ALL_FILTERED} -k 0 1 6 10 -c 4 -v > all_clair3_pct_context.tsv
sed -i 's/^\t\t\t/scaffold\tpos\tgene\tcontext/' all_clair3_pct_context.tsv && sed -i 's/\.bed//g' all_clair3_pct_context.tsv
gzip all_clair3_pct_context.tsv
