#!/bin/bash
#SBATCH --job-name=3RAD_VCF
#SBATCH --partition=batch
#SBATCH --nodes=2 # use just 1
#SBATCH --ntasks=3
#SBATCH --mem=50gb # Per processor memory request
#SBATCH --export=NONE
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=jcf87188@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=2_3RAD_vcftools.out
#SBATCH --error=2_3RAD_vcftools.err

cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

dir=/scratch/jcf87188/reseq_iscap/v2_clean_all_outputs/vcf_compare

mkdir $dir/vcftools/3RAD

ml VCFtools/0.1.16-GCC-8.3.0-Perl-5.30.0
VCF_all=/scratch/jcf87188/reseq_iscap/3RAD_analysis/stacks_reseq_mapped_bam/vcftools/populations.snps_sorted.vcf.gz
VCF=$dir/filtered_vcfs/3RAD_m90_6x150x_GQ30_maf05_ma2.recode.vcf.gz
OUT=$dir/vcftools/compare_filters

# Filtering for 3RAD
vcftools --gzvcf $VCF_all --max-missing 0.90 --recode --recode-INFO-all --out $OUT/3RAD_m90
vcftools --gzvcf $VCF_all --max-missing 0.90 --minDP 6 --recode --recode-INFO-all --out $OUT/3RAD_m90_6x
vcftools --gzvcf $VCF_all --max-missing 0.90 --minDP 6 --maxDP 150 --recode --recode-INFO-all --out $OUT/3RAD_m90_6x150x
vcftools --gzvcf $VCF_all --max-missing 0.90 --minDP 6 --maxDP 150 --minGQ 30 --recode --recode-INFO-all --out $OUT/3RAD_m90_6x150x_GQ30
vcftools --gzvcf $VCF_all --max-missing 0.90 --minDP 6 --maxDP 150 --minGQ 30 --maf 0.05 --recode --recode-INFO-all --out $OUT/3RAD_m90_6x150x_GQ30_maf05
vcftools --gzvcf $VCF_all --max-missing 0.90 --minDP 6 --maxDP 150 --minGQ 30 --maf 0.05 --max-alleles 2 --recode --recode-INFO-all $OUT/3RAD_m90_6x150x_GQ30_maf05_ma2
