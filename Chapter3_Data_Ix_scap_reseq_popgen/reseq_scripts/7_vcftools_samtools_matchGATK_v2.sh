#!/bin/bash
#SBATCH --job-name=samtools_VCF
#SBATCH --partition=batch
#SBATCH --nodes=2 # use just 1
#SBATCH --ntasks=3
#SBATCH --mem=50gb # Per processor memory request
#SBATCH --export=NONE
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=jcf87188@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=7_samtools_VCF_matchsamtools_formatQ_v2.out
#SBATCH --error=7_samtools_VCF_matchsamtools_formatQ_v2.err

cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

dir=/scratch/jcf87188/reseq_iscap/v2_clean_all_outputs/samtools_snpcalling

cd $dir/reseq_all_samtools_formatQ

ml VCFtools/0.1.16-GCC-8.3.0-Perl-5.30.0
VCF=$dir/reseq_all_samtools_formatQ/reseq_all_calls_formatQ.vcf.gz
OUT=$dir/reseq_all_samtools_formatQ/reseq_all_calls_formatQ


# Filtering for samtools
vcftools --gzvcf $VCF --max-missing 0.95 --minDP 4 --maxDP 50 --minGQ 25 --minQ 30 --maf 0.05 --max-alleles 2 --recode --recode-INFO-all --out $dir/reseq_all_samtools_formatQ/samtools_formatQ_m95_4x50x_GQ25_Q30_maf05_ma2
vcftools --gzvcf $VCF --remove-indels --max-missing 0.95 --minDP 4 --maxDP 50 --minGQ 25 --minQ 30 --maf 0.05 --max-alleles 2 --recode --recode-INFO-all --out $dir/reseq_all_samtools_formatQ/samtools_formatQ_rmIndels_m95_4x50x_GQ25_Q30_maf05_ma2