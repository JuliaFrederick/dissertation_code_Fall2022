#!/bin/bash
#SBATCH --job-name=vcf_compare_GS
#SBATCH --partition=highmem_p
#SBATCH --nodes=2 # use just 1
#SBATCH --ntasks=3
#SBATCH --mem=20gb # Per processor memory request
#SBATCH --export=NONE
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=jcf87188@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=1_vcf_compare_GS.out
#SBATCH --error=1_vcf_compare_GS.err

cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#loading the software in the cluster
ml VCFtools/0.1.16-GCC-8.3.0-Perl-5.30.0
ml SAMtools/1.14-GCC-8.3.0-fixed

vcf_3RAD=3RAD_m90_6x150x_GQ30_maf05_ma2.recode.vcf.gz
vcf_samtools=samtools_formatQ_m90_6x50x_GQ30_Q30_maf05_ma2.recode.vcf.gz
vcf_GATK=gatk_m90_6x50x_GQ30_Q30_maf05_ma2.recode.vcf.gz

genome=/home/jcf87188/202107Pal_IscapGenome/GCF_016920785.2_ASM1692078v2_genomic.fna

vcf-compare $vcf_GATK $vcf_samtools

#-R $genome -w 50

