#!/bin/bash
#SBATCH --job-name=Full_Enzyme_Stacks_1
#SBATCH --partition=highmem_p
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=12
#SBATCH --mem=400gb
#SBATCH --export=NONE
#SBATCH --time=96:00:00
#SBATCH --output=Full_Enzyme_Stacks_1.out
#SBATCH --error=Full_Enzyme_Stacks_1.err
#SBATCH --mail-user=jcf87188@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL

cd $SLURM_SUBMIT_DIR

#loading the software in the cluster
#ml Stacks/2.5-iccifort-2019.5.281


datafolder_new=/lustre2/scratch/jcf87188/oct2021_plate4_PG/data_new

mkdir $datafolder_new/process_output_final_v2
cp $datafolder_new/process_output_second_*/*.fq.gz $datafolder_new/process_output_final_v2/
cd $datafolder_new/process_output_final_v2/
rename .1.1. .1. *
rename .2.2. .2. *
gunzip *.fq.gz
sed -Ei 's/\/1\/1/\/1/' *.1.fq
sed -Ei 's/\/2\/2/\/2/' *.2.fq
gzip *.fq
