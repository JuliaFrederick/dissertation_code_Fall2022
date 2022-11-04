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
#SBATCH --output=2_samtools_vcftools_v2.out
#SBATCH --error=2_samtools_vcftools_v2.err

cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

dir=/scratch/jcf87188/reseq_iscap/v2_clean_all_outputs/vcf_compare

#mkdir $dir/vcftools/samtools

ml VCFtools/0.1.16-GCC-8.3.0-Perl-5.30.0
VCF_all=/scratch/jcf87188/reseq_iscap/v2_clean_all_outputs/samtools_snpcalling/reseq_all_samtools_formatQ/reseq_all_calls_formatQ.vcf.gz
VCF=$dir/filtered_vcfs/samtools_formatQ_m90_6x50x_GQ30_Q30_maf05_ma2.recode.vcf.gz
OUT=$dir/vcftools/compare_filters

# Filtering for Samtools 
vcftools --gzvcf $VCF_all --max-missing 0.90 --recode --recode-INFO-all --out $OUT/samtools_m90
vcftools --gzvcf $VCF_all --max-missing 0.90 --minDP 6 --recode --recode-INFO-all --out $OUT/samtools_m90_6x
vcftools --gzvcf $VCF_all --max-missing 0.90 --minDP 6 --maxDP 50 --recode --recode-INFO-all --out $OUT/samtools_m90_6x50x
vcftools --gzvcf $VCF_all --max-missing 0.90 --minDP 6 --maxDP 50 --minGQ 30 --recode --recode-INFO-all --out $OUT/samtools_m90_6x50x_GQ30
vcftools --gzvcf $VCF_all --max-missing 0.90 --minDP 6 --maxDP 50 --minGQ 30 --minQ 30 --recode --recode-INFO-all --out $OUT/samtools_m90_6x50x_GQ30_Q30
vcftools --gzvcf $VCF_all --max-missing 0.90 --minDP 6 --maxDP 50 --minGQ 30 --minQ 30 --maf 0.05 --recode --recode-INFO-all --out $OUT/samtools_m90_6x50x_GQ30_Q30_maf05
vcftools --gzvcf $VCF_all --max-missing 0.90 --minDP 6 --maxDP 50 --minGQ 30 --minQ 30 --maf 0.05 --max-alleles 2 --recode --recode-INFO-all --out $OUT/samtools_m90_6x50x_GQ30_Q30_maf05_ma2
vcftools --gzvcf $VCF_all --remove-indels --max-missing 0.90 --minDP 6 --maxDP 50 --minGQ 30 --minQ 30 --maf 0.05 --max-alleles 2 --recode --recode-INFO-all --out $OUT/samtools_m90_6x50x_GQ30_Q30_maf05_ma2_rmInd
vcftools --gzvcf $VCF_all --keep-only-indels --max-missing 0.90 --minDP 6 --maxDP 50 --minGQ 30 --minQ 30 --maf 0.05 --max-alleles 2 --recode --recode-INFO-all --out $OUT/samtools_m90_6x50x_GQ30_Q30_maf05_ma2_onlyInd
