#!/bin/bash
#SBATCH --job-name=indexgenome
#SBATCH --partition=highmem_p
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=1gb
#SBATCH --export=NONE
#SBATCH --time=24:00:00
#SBATCH --output=indexgenome.out
#SBATCH --error=indexgenome.err
#SBATCH --mail-user=jcf87188@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL

cd $SLURM_SUBMIT_DIR

#loading the software in the cluster
#ml Stacks/2.5-iccifort-2019.5.281
#ml Python/3.8.2-GCCcore-8.3.0
ml BWA/0.7.17-GCC-8.3.0
#ml SAMtools/1.10-iccifort-2019.5.281

gunzip /lustre2/scratch/jcf87188/ricinus/genomes/cmidichloria/GCF_000219355.1_ASM21935v1_genomic.fna.gz
bwa index /lustre2/scratch/jcf87188/ricinus/genomes/cmidichloria/GCF_000219355.1_ASM21935v1_genomic.fna
