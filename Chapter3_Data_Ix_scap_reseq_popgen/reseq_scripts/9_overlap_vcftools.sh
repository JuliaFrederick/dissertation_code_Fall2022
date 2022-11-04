#!/bin/bash
#SBATCH --job-name=overlap_VCF
#SBATCH --partition=batch
#SBATCH --nodes=2 # use just 1
#SBATCH --ntasks=3
#SBATCH --mem=50gb # Per processor memory request
#SBATCH --export=NONE
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=jcf87188@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=2_overlap_vcftools.out
#SBATCH --error=2_overlap_vcftools.err

cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

dir=/scratch/jcf87188/reseq_iscap/v2_clean_all_outputs/vcf_compare

mkdir $dir/vcftools/overlap

ml VCFtools/0.1.16-GCC-8.3.0-Perl-5.30.0
VCF=$dir/filtered_vcfs/dir/0000.vcf.gz
OUT=$dir/vcftools/overlap/overlap0

vcftools --gzvcf $VCF --freq2 --out $OUT --max-alleles 2
vcftools --gzvcf $VCF --depth --out $OUT
vcftools --gzvcf $VCF --site-mean-depth --out $OUT
vcftools --gzvcf $VCF --site-quality --out $OUT
vcftools --gzvcf $VCF --missing-site --out $OUT
vcftools --gzvcf $VCF --missing-indv --out $OUT
vcftools --gzvcf $VCF --het --out $OUT

VCF=$dir/filtered_vcfs/dir/0001.vcf.gz
OUT=$dir/vcftools/overlap/overlap1

vcftools --gzvcf $VCF --freq2 --out $OUT --max-alleles 2
vcftools --gzvcf $VCF --depth --out $OUT
vcftools --gzvcf $VCF --site-mean-depth --out $OUT
vcftools --gzvcf $VCF --site-quality --out $OUT
vcftools --gzvcf $VCF --missing-site --out $OUT
vcftools --gzvcf $VCF --missing-indv --out $OUT
vcftools --gzvcf $VCF --het --out $OUT