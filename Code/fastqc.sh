#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=10G
#SBATCH --job-name=fastqc
#SBATCH --output=%x.out
#SBATCH --error=%x.err

indir=/work/nclab/lizzy/hawcofrag/ribodetector/
outdir=/work/nclab/lizzy/hawcofrag/ribodetector/fastqc/


#on the command line, execute:
#module load FastQC

fastqc $indir/*_L001.R1.fq -o $outdir
