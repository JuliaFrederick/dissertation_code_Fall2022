#!/bin/bash
#SBATCH --job-name=mito_mapped
#SBATCH --partition=batch
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --ntasks=1                    # Number of MPI ranks
#SBATCH --ntasks-per-node=2           # Number of MPI ranks per node
#SBATCH --cpus-per-task=2             # Number of OpenMP threads for each MPI process/rank
#SBATCH --mem-per-cpu=45gb          # Per processor memory request
#SBATCH --export=NONE
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=jcf87188@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=2,21 # Here is where you set the number of job elements

cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

i=$(cat input_${SLURM_ARRAY_TASK_ID})

fastq=/scratch/jcf87188/reseq_iscap/v2_clean_reads_2022/paired_Illumina_trim_reads
mtDNA=/scratch/jcf87188/reseq_iscap/v2_clean_all_outputs/mt_DNA/mitofinder_mapped/Ixscap_mtDNA_MZ645749.1.gb
out=/scratch/jcf87188/reseq_iscap/v2_clean_all_outputs/mt_DNA/mitofinder_mapped

cd $out
singularity run /apps/singularity-images/mitofinder_v1.4.1.sif -j $i -1 $fastq/$i\_R1_pair_trim.fastq.gz -2 $fastq/$i\_R2_pair_trim.fastq.gz -r $mtDNA -o 5
