#!/bin/bash
#SBATCH --job-name=samtools_SNP
#SBATCH --partition=highmem_30d_p
#SBATCH --nodes=2 # use just 1
#SBATCH --ntasks=3
#SBATCH --mem=170gb # Per processor memory request
#SBATCH --export=NONE
#SBATCH --time=10-00:00:00
#SBATCH --mail-user=jcf87188@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=4_samtools_SNP_matchGATK_formatQ.out
#SBATCH --error=4_samtools_SNP_matchGATK_formatQ.err

cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#loading the software in the cluster
ml SAMtools/1.14-GCC-8.3.0-fixed
ml BCFtools/1.13-GCC-8.3.0

genome=/home/jcf87188/202107Pal_IscapGenome/GCF_016920785.2_ASM1692078v2_genomic.fna
dir=/scratch/jcf87188/reseq_iscap/v2_clean_all_outputs

# From https://www.biostars.org/p/377842/

bcftools mpileup \
    --redo-BAQ \
    --per-sample-mF \
    --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
    -f $genome \
    -b $dir/1_bamfiles_md_fx.txt | \
  bcftools call \
    --multiallelic-caller \
    --ploidy 2 \
    --variants-only \
    --annotate FORMAT/GQ \
    -Ov \
  > $dir/samtools_snpcalling/reseq_all_calls_formatQ.vcf