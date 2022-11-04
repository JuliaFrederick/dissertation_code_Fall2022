#!/bin/bash
#SBATCH --job-name=cleanP4new
#SBATCH --partition=batch
#SBATCH --nodes=1 # use just 1
#SBATCH --ntasks=1 # use just 1
#SBATCH --cpus-per-task=4 # Number of OpenMP threads per job element
#SBATCH --mem-per-cpu=200mb # Per processor memory request
#SBATCH --export=NONE
#SBATCH --time=96:00:00
#SBATCH --mail-user=jcf87188@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=1-133 # Here is where you set the number of job elements

cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#loading the software in the cluster
ml Stacks/2.5-iccifort-2019.5.281
#ml Python/3.8.2-GCCcore-8.3.0
ml BWA/0.7.17-GCC-8.3.0
ml SAMtools/1.10-iccifort-2019.5.281

mkdir /lustre2/scratch/jcf87188/oct2021_plate4_PG/data_new/

datafolder=/lustre2/scratch/jcf87188/oct2021_plate4_PG/data_new
barcodes=/lustre2/scratch/jcf87188/oct2021_plate4_PG/barcodes_newsamples
reads=/lustre2/scratch/jcf87188/oct2021_plate4_PG/rawdata_newsamples

i=$(cat input_${SLURM_ARRAY_TASK_ID})

mkdir $datafolder/process_output_first_$i
process_radtags -P -1 $reads/$i/$i\_1.fq.gz -2 $reads/$i/$i\_2.fq.gz -b $barcodes/$i\.txt --inline_inline -c -q -r -t 140 --disable_rad_check -o $datafolder/process_output_first_$i/ -i gzfastq

mkdir $datafolder/process_output_second_$i
process_radtags -P -1 $datafolder/process_output_first_$i/$i\.1.fq.gz -2 $datafolder/process_output_first_$i/$i\.2.fq.gz -c -q -t 140 --renz_1 xbaI --renz_2 ecoRI -o $datafolder/process_output_second_$i/ -i gzfastq

## Needs to be run separately now because of the array
#mkdir $datafolder/process_output_final
#cp $datafolder/process_output_second_*/*.fq.gz $datafolder/process_output_final/
#cd $datafolder/process_output_final/
#rename .1.1. .1. *
#rename .2.2. .2. *
#gunzip *.fq.gz
#sed -Ei 's/\/1\/1/\/1/' *.1.fq
#sed -Ei 's/\/2\/2/\/2/' *.2.fq
#gzip *.fq
