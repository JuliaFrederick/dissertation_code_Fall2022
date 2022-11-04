#!/bin/bash
#SBATCH --job-name=2_vcftools
#SBATCH --partition=batch
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --ntasks=4                    # Number of MPI ranks
#SBATCH --ntasks-per-node=4           # Number of MPI ranks per node
#SBATCH --cpus-per-task=4             # Number of OpenMP threads for each MPI process/rank
#SBATCH --mem-per-cpu=2000mb          # Per processor memory request
#SBATCH --export=NONE
#SBATCH --time=1-00:00:00
#SBATCH --output=2_vcftools.out
#SBATCH --error=2_vcftools.err
#SBATCH --mail-user=jcf87188@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL

cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#loading the software in the cluster
ml VCFtools/0.1.16-GCC-8.3.0-Perl-5.30.0

hmdir=/lustre2/scratch/jcf87188/3RAD_popgen_allPOPS/stacks_ams_bam_pops_bestSamples

## Getting Data 
vcftools --vcf $hmdir/populations.snps.vcf --out $hmdir/vcftools/rmInd_p50_5x200x --remove-indels --max-missing 0.50 --minDP 5 --maxDP 200 --recode --recode-INFO-all
vcftools --vcf $hmdir/populations.snps.vcf --out $hmdir/vcftools/rmInd_p10_5x200x --remove-indels --max-missing 0.90 --minDP 5 --maxDP 200 --recode --recode-INFO-all

vcftools --vcf $hmdir/populations.snps.vcf --out $hmdir/vcftools/rmInd_p50_5x200x_maf05 --remove-indels --max-missing 0.50 --minDP 5 --maxDP 200 --maf 0.05 --recode --recode-INFO-all
vcftools --vcf $hmdir/populations.snps.vcf --out $hmdir/vcftools/rmInd_p10_5x200x_maf05 --remove-indels --max-missing 0.90 --minDP 5 --maxDP 200 --maf 0.05 --recode --recode-INFO-all


## All
vcftools --vcf $hmdir/populations.snps.vcf --freq2 --out $hmdir/vcftools/all --max-alleles 2
vcftools --vcf $hmdir/populations.snps.vcf --depth --out $hmdir/vcftools/all
vcftools --vcf $hmdir/populations.snps.vcf --site-mean-depth --out $hmdir/vcftools/all
vcftools --vcf $hmdir/populations.snps.vcf --site-quality --out $hmdir/vcftools/all
vcftools --vcf $hmdir/populations.snps.vcf --missing-indv --out $hmdir/vcftools/all
vcftools --vcf $hmdir/populations.snps.vcf --missing-site --out $hmdir/vcftools/all
vcftools --vcf $hmdir/populations.snps.vcf --het --out $hmdir/vcftools/all

## 50% missing
vcftools --vcf $hmdir/vcftools/rmInd_p50_5x200x.recode.vcf --freq2 --out $hmdir/vcftools/rmInd_p50_5x200x --max-alleles 2
vcftools --vcf $hmdir/vcftools/rmInd_p50_5x200x.recode.vcf --depth --out $hmdir/vcftools/rmInd_p50_5x200x
vcftools --vcf $hmdir/vcftools/rmInd_p50_5x200x.recode.vcf --site-mean-depth --out $hmdir/vcftools/rmInd_p50_5x200x
vcftools --vcf $hmdir/vcftools/rmInd_p50_5x200x.recode.vcf --site-quality --out $hmdir/vcftools/rmInd_p50_5x200x
vcftools --vcf $hmdir/vcftools/rmInd_p50_5x200x.recode.vcf --missing-indv --out $hmdir/vcftools/rmInd_p50_5x200x
vcftools --vcf $hmdir/vcftools/rmInd_p50_5x200x.recode.vcf --missing-site --out $hmdir/vcftools/rmInd_p50_5x200x
vcftools --vcf $hmdir/vcftools/rmInd_p50_5x200x.recode.vcf --het --out $hmdir/vcftools/rmInd_p50_5x200x

## 80% missing
vcftools --vcf $hmdir/vcftools/rmInd_p20_5x200x.recode.vcf --freq2 --out $hmdir/vcftools/rmInd_p20_5x200x --max-alleles 2
vcftools --vcf $hmdir/vcftools/rmInd_p20_5x200x.recode.vcf --depth --out $hmdir/vcftools/rmInd_p20_5x200x
vcftools --vcf $hmdir/vcftools/rmInd_p20_5x200x.recode.vcf --site-mean-depth --out $hmdir/vcftools/rmInd_p20_5x200x
vcftools --vcf $hmdir/vcftools/rmInd_p20_5x200x.recode.vcf --site-quality --out $hmdir/vcftools/rmInd_p20_5x200x
vcftools --vcf $hmdir/vcftools/rmInd_p20_5x200x.recode.vcf --missing-indv --out $hmdir/vcftools/rmInd_p20_5x200x
vcftools --vcf $hmdir/vcftools/rmInd_p20_5x200x.recode.vcf --missing-site --out $hmdir/vcftools/rmInd_p20_5x200x
vcftools --vcf $hmdir/vcftools/rmInd_p20_5x200x.recode.vcf --het --out $hmdir/vcftools/rmInd_p20_5x200x

## 90% missing
vcftools --vcf $hmdir/vcftools/rmInd_p10_5x200x.recode.vcf --freq2 --out $hmdir/vcftools/rmInd_p10_5x200x --max-alleles 2
vcftools --vcf $hmdir/vcftools/rmInd_p10_5x200x.recode.vcf --depth --out $hmdir/vcftools/rmInd_p10_5x200x
vcftools --vcf $hmdir/vcftools/rmInd_p10_5x200x.recode.vcf --site-mean-depth --out $hmdir/vcftools/rmInd_p10_5x200x
vcftools --vcf $hmdir/vcftools/rmInd_p10_5x200x.recode.vcf --site-quality --out $hmdir/vcftools/rmInd_p10_5x200x
vcftools --vcf $hmdir/vcftools/rmInd_p10_5x200x.recode.vcf --missing-indv --out $hmdir/vcftools/rmInd_p10_5x200x
vcftools --vcf $hmdir/vcftools/rmInd_p10_5x200x.recode.vcf --missing-site --out $hmdir/vcftools/rmInd_p10_5x200x
vcftools --vcf $hmdir/vcftools/rmInd_p10_5x200x.recode.vcf --het --out $hmdir/vcftools/rmInd_p10_5x200x


## 50% missing with maf 5%
vcftools --vcf $hmdir/vcftools/rmInd_p50_5x200x_maf05.recode.vcf --freq2 --out $hmdir/vcftools/rmInd_p50_5x200x_maf05 --max-alleles 2
vcftools --vcf $hmdir/vcftools/rmInd_p50_5x200x_maf05.recode.vcf --depth --out $hmdir/vcftools/rmInd_p50_5x200x_maf05
vcftools --vcf $hmdir/vcftools/rmInd_p50_5x200x_maf05.recode.vcf --site-mean-depth --out $hmdir/vcftools/rmInd_p50_5x200x_maf05
vcftools --vcf $hmdir/vcftools/rmInd_p50_5x200x_maf05.recode.vcf --site-quality --out $hmdir/vcftools/rmInd_p50_5x200x_maf05
vcftools --vcf $hmdir/vcftools/rmInd_p50_5x200x_maf05.recode.vcf --missing-indv --out $hmdir/vcftools/rmInd_p50_5x200x_maf05
vcftools --vcf $hmdir/vcftools/rmInd_p50_5x200x_maf05.recode.vcf --missing-site --out $hmdir/vcftools/rmInd_p50_5x200x_maf05
vcftools --vcf $hmdir/vcftools/rmInd_p50_5x200x_maf05.recode.vcf --het --out $hmdir/vcftools/rmInd_p50_5x200x_maf05

## 80% missing with maf 5%
vcftools --vcf $hmdir/vcftools/rmInd_p20_5x200x_maf05.recode.vcf --freq2 --out $hmdir/vcftools/rmInd_p20_5x200x_maf05 --max-alleles 2
vcftools --vcf $hmdir/vcftools/rmInd_p20_5x200x_maf05.recode.vcf --depth --out $hmdir/vcftools/rmInd_p20_5x200x_maf05
vcftools --vcf $hmdir/vcftools/rmInd_p20_5x200x_maf05.recode.vcf --site-mean-depth --out $hmdir/vcftools/rmInd_p20_5x200x_maf05
vcftools --vcf $hmdir/vcftools/rmInd_p20_5x200x_maf05.recode.vcf --site-quality --out $hmdir/vcftools/rmInd_p20_5x200x_maf05
vcftools --vcf $hmdir/vcftools/rmInd_p20_5x200x_maf05.recode.vcf --missing-indv --out $hmdir/vcftools/rmInd_p20_5x200x_maf05
vcftools --vcf $hmdir/vcftools/rmInd_p20_5x200x_maf05.recode.vcf --missing-site --out $hmdir/vcftools/rmInd_p20_5x200x_maf05
vcftools --vcf $hmdir/vcftools/rmInd_p20_5x200x_maf05.recode.vcf --het --out $hmdir/vcftools/rmInd_p20_5x200x_maf05

## 90% missing with maf 5%
vcftools --vcf $hmdir/vcftools/rmInd_p10_5x200x_maf05.recode.vcf --freq2 --out $hmdir/vcftools/rmInd_p10_5x200x_maf05 --max-alleles 2
vcftools --vcf $hmdir/vcftools/rmInd_p10_5x200x_maf05.recode.vcf --depth --out $hmdir/vcftools/rmInd_p10_5x200x_maf05
vcftools --vcf $hmdir/vcftools/rmInd_p10_5x200x_maf05.recode.vcf --site-mean-depth --out $hmdir/vcftools/rmInd_p10_5x200x_maf05
vcftools --vcf $hmdir/vcftools/rmInd_p10_5x200x_maf05.recode.vcf --site-quality --out $hmdir/vcftools/rmInd_p10_5x200x_maf05
vcftools --vcf $hmdir/vcftools/rmInd_p10_5x200x_maf05.recode.vcf --missing-indv --out $hmdir/vcftools/rmInd_p10_5x200x_maf05
vcftools --vcf $hmdir/vcftools/rmInd_p10_5x200x_maf05.recode.vcf --missing-site --out $hmdir/vcftools/rmInd_p10_5x200x_maf05
vcftools --vcf $hmdir/vcftools/rmInd_p10_5x200x_maf05.recode.vcf --het --out $hmdir/vcftools/rmInd_p10_5x200x_maf05

