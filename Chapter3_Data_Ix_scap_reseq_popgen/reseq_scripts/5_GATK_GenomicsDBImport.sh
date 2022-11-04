#!/bin/bash
#SBATCH --job-name=gatk_DBI
#SBATCH --partition=highmem_30d_p
#SBATCH --nodes=1 # use just 1
#SBATCH --ntasks=4
#SBATCH --mem=450gb # Per processor memory request
#SBATCH --export=NONE
#SBATCH --time=30-00:00:00
#SBATCH --mail-user=jcf87188@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=5_gatk_DBI.out
#SBATCH --error=5_gatk_DBI.err

cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#loading the software in the cluster
ml SAMtools/1.10-iccifort-2019.5.281
ml GATK/4.1.8.1-GCCcore-8.3.0-Java-1.8

cd $SLURM_SUBMIT_DIR
# we need a tmp dir that is large for the program to use

mkdir -p $SLURM_SUBMIT_DIR/tmp

# set batch size equal to cores on node also set max RAM a
# little low, because there is additional overhead
# involved per GATK website
gatk --java-options "-Xmx250g" GenomicsDBImport \
   --genomicsdb-workspace-path my_database \
   --sample-name-map 5.1_gvcfs-for-db-import.sample_map_v6 \
   --tmp-dir $SLURM_SUBMIT_DIR/tmp \
   --batch-size 32 \
   -L /home/jcf87188/202107Pal_IscapGenome/GCF_016920785.2_ASM1692078v2_genomic.fna.bed