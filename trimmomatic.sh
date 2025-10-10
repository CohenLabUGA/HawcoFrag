#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=10G
#SBATCH --job-name=trimmomatic
#SBATCH --output=%x.out
#SBATCH --error=%x.err

#load module trimmomatic
module load Trimmomatic

#name indir and outdir
indir="/work/nclab/data/Hawco/25099FL_06_01_01"
outdir="/work/nclab/lizzy/hawcofrag/trimmed-reads"

# Define adapter file path
adapter_file="/work/nclab/lizzy/hawcofrag/scripts/adapters.fa"

for r1 in "$indir"/*_R1_001.fastq.gz; do
    # Get the prefix
    base=$(basename "$r1" _R1_001.fastq.gz)
    r2="$indir/${base}_R2_001.fastq.gz"

 if [[ -f "$outdir/${base}_R1_trimmed_paired.fastq.gz" && -f "$outdir/${base}_R2_trimmed_paired.fastq.gz" ]]; then
        echo "Skipping $base — already processed"
        continue
    fi

    echo "Processing $base at $(date)"

    trimmomatic PE \
        "$r1" "$r2" \
        "$outdir/${base}_R1_trimmed_paired.fastq.gz" "$outdir/${base}_R1_trimmed_unpaired.fastq.gz" \
        "$outdir/${base}_R2_trimmed_paired.fastq.gz" "$outdir/${base}_R2_trimmed_unpaired.fastq.gz" \
        ILLUMINACLIP:"$adapter_file":2:30:7 \
        LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2 MINLEN:50
done
