#!/bin/bash
#SBATCH --job-name=denovo_amer
#SBATCH --partition=batch
#SBATCH --nodes=2                     # Number of nodes
#SBATCH --ntasks=8                    # Number of MPI ranks
#SBATCH --ntasks-per-node=4           # Number of MPI ranks per node
#SBATCH --cpus-per-task=4             # Number of OpenMP threads for each MPI process/rank
#SBATCH --mem-per-cpu=2000mb          # Per processor memory request
#SBATCH --export=NONE
#SBATCH --time=3-00:00:00
#SBATCH --output=denovo_amer.out
#SBATCH --error=denovo_amer.err
#SBATCH --mail-user=jcf87188@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL

cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#loading the software in the cluster
#ml Stacks/2.55-foss-2019b
ml Stacks/2.61-foss-2019b
#ml Python/3.8.2-GCCcore-8.3.0
#ml BWA/0.7.17-GCC-8.3.0
#ml SAMtools/1.10-iccifort-2019.5.281


populations -P /lustre2/scratch/jcf87188/amblyomma_3rad/denovo_output/denovo_n_8 -M /lustre2/scratch/jcf87188/amblyomma_3rad/AamerJuly2019_populationmap.txt -O /lustre2/scratch/jcf87188/amblyomma_3rad/populations_M7n8/ -r 0.20 -p 1 --fasta-loci --fasta-samples --fasta-samples-raw --vcf --fstats