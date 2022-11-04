#!/bin/bash
#SBATCH --job-name=gatk_gVCF
#SBATCH --partition=highmem_30d_p
#SBATCH --nodes=1 # use just 1
#SBATCH --ntasks=4
#SBATCH --mem=450gb # Per processor memory request
#SBATCH --export=NONE
#SBATCH --time=30-00:00:00
#SBATCH --mail-user=jcf87188@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=6_gatk_gVCF.out
#SBATCH --error=6_gatk_gVCF.err

#loading the software in the cluster
ml SAMtools/1.10-iccifort-2019.5.281
ml GATK/4.1.8.1-GCCcore-8.3.0-Java-1.8

cd $SLURM_SUBMIT_DIR

# set reference
REFERENCE=/home/jcf87188/202107Pal_IscapGenome/GCF_016920785.2_ASM1692078v2_genomic.fna

# now go ahead and run and set batch size equal to cores on node
# also set max RAM a little low, because there is additional overhead
# involved
gatk --java-options "-Xmx250g" \
    GenotypeGVCFs \
    -R $REFERENCE \
    -V gendb://my_database \
    -O 6_gatk_output.vcf.gz
