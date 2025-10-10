#!/usr/bin/env bash
#SBATCH --partition=batch
#SBATCH --time=96:00:00
#SBATCH --job-name=htseq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --output=htseq.out
#SBATCH --error=htseq.err

ANNOT=/work/nclab/lizzy/hawcofrag/genome/Fracy1_GeneModels_FilteredModels2.gtf 
INPUT=/work/nclab/lizzy/hawcofrag/genome/star_align/

# Loop through all sample folders
for dir in "${INPUT}"/25099FL-06-01-*_S*_L001; do
    SAMPLE=$(basename "$dir")
    BAM="${dir}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
    OUT="${SAMPLE}.htseq.model2.txt"

    # Run htseq-count
    htseq-count \
      -f bam \
      -r pos \
      -s reverse \
      -t exon \
      -i gene_id \
      --minaqual 10 \
      --secondary-alignments ignore \
      --supplementary-alignments ignore \
      "$BAM" "$ANNOT" > "$OUT"
done
