#!/usr/bin/env bash
#SBATCH --partition=batch
#SBATCH --time=96:00:00
#SBATCH --job-name=STAR
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=180G
#SBATCH --output=STAR.out
#SBATCH --error=STAR.err

GENOME_DIR="/work/nclab/lizzy/hawcofrag/genome/index"
READDIR="/work/nclab/lizzy/hawcofrag/ribodetector"
OUTBASE="/work/nclab/lizzy/hawcofrag/genome/star_align"
THREADS="${SLURM_CPUS_PER_TASK:-12}"


# Map S1..S12 (filenames like 25099FL-06-01-01_S1_L001.R1.fq / .R2.fq)
for i in 4; do
#for i in $(seq 1 12); do
  I2=$(printf "%02d" "$i")

  R1="${READDIR}/25099FL-06-01-${I2}_S${i}_L001.R1.fq"
  R2="${READDIR}/25099FL-06-01-${I2}_S${i}_L001.R2.fq"

  SAMPLE=$(basename "$R1" | sed -E 's/\.R1\.fq$//')   # e.g., 25099FL-06-01-01_S1_L001
  OUTDIR="${OUTBASE}/${SAMPLE}"
  mkdir -p "$OUTDIR"
  PREFIX="${OUTDIR}/${SAMPLE}_"

  ./bin/Linux_x86_64_static/STAR --runThreadN "$THREADS" \
       --genomeDir "$GENOME_DIR" \
       --readFilesIn "$R1" "$R2" \
       --outFileNamePrefix "$PREFIX" \
       --outSAMtype BAM SortedByCoordinate \
       --limitBAMsortRAM 180000000000 \
       --outBAMsortingBinsN 200 \
       --quantMode GeneCounts \
       --outSAMattrRGline ID:${SAMPLE} SM:${SAMPLE} PL:ILLUMINA
done
