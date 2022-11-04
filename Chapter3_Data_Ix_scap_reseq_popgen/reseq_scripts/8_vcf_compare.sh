#!/bin/bash
#SBATCH --job-name=vcf_compare
#SBATCH --partition=highmem_p
#SBATCH --nodes=2 # use just 1
#SBATCH --ntasks=3
#SBATCH --mem=80gb # Per processor memory request
#SBATCH --export=NONE
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=jcf87188@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=1_vcf_compare.out
#SBATCH --error=1_vcf_compare.err

cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#loading the software in the cluster
ml VCFtools/0.1.16-GCC-8.3.0-Perl-5.30.0
ml SAMtools/1.14-GCC-8.3.0-fixed

vcf_3RAD=/scratch/jcf87188/reseq_iscap/3RAD_analysis/stacks_reseq_mapped_bam/vcftools/populations.snps_sorted.vcf.gz
vcf_samtools=/scratch/jcf87188/reseq_iscap/v2_clean_all_outputs/samtools_snpcalling/reseq_all_samtools_formatQ/reseq_all_calls_formatQ.vcf.gz
vcf_GATK=/scratch/jcf87188/reseq_iscap/v2_clean_all_outputs/gatk_snpcalling/6_gatk_output.vcf.gz

genome=/home/jcf87188/202107Pal_IscapGenome/GCF_016920785.2_ASM1692078v2_genomic.fna

vcf-compare $vcf_GATK $vcf_samtools $vcf_3RAD 

#-R $genome -w 50

