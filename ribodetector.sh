#!/bin/bash
#SBATCH --partition=batch
#SBATCH --time=96:00:00
#SBATCH --job-name=ribodetector
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G		
#SBATCH --threads-per-core=1
#SBATCH --output=ribodetector.out
#SBATCH --error=ribodetector.err
#SBATCH --verbose
export OMP_NUM_THREADS=20
indir="/work/nclab/lizzy/hawcofrag/trimmed-reads"
outdir="/work/nclab/lizzy/hawcofrag/ribodetector"

for r1_file in "$indir"/*_R1_trimmed_paired.fastq.gz; do
    # Extract the sample prefix by removing _R1_trimmed_paired.fastq.gz
    sample_prefix=$(basename "$r1_file" | sed 's/_R1_trimmed_paired.fastq.gz//')
    # Construct corresponding R2 file name
    r2_file="$indir/${sample_prefix}_R2_trimmed_paired.fastq.gz"
    # Define expected output files
    out_r1="$outdir/${sample_prefix}.R1.fq"
    out_r2="$outdir/${sample_prefix}.R2.fq"
    # Check if output files already exist
    if [[ -f "$out_r1" && -f "$out_r2" ]]; then
        echo "Sample $sample_prefix already processed, skipping."
        continue
    fi
    # Check if the R2 file exists
    if [[ -f "$r2_file" ]]; then
        echo "Processing sample: $sample_prefix"
        echo "R1: $(basename "$r1_file")"
        echo "R2: $(basename "$r2_file")"
        ribodetector_cpu -t 20 -l 151 \
            -i "$r1_file" "$r2_file" \
            -e rrna \
            --chunk_size 128 \
            -o "$outdir/${sample_prefix}.R1.fq" "$outdir/${sample_prefix}.R2.fq"
    else
        echo "Warning: R2 file not found for $sample_prefix, skipping."
    fi
done

