#!/bin/bash
#SBATCH --job-name=3_denovo_n
#SBATCH --partition=batch
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --ntasks=1                    # Number of MPI ranks
#SBATCH --ntasks-per-node=2           # Number of MPI ranks per node
#SBATCH --cpus-per-task=2             # Number of OpenMP threads for each MPI process/rank
#SBATCH --mem-per-cpu=30gb          # Per processor memory request
#SBATCH --export=NONE
#SBATCH --time=4-10:00:00
#SBATCH --mail-user=jcf87188@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=1-10 # Here is where you set the number of job elements

cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#loading the software in the cluster
ml Stacks/2.61-foss-2019b
ml Perl/5.30.0-GCCcore-8.3.0

#mkdir  /lustre2/scratch/jcf87188/3RAD_popgen_allPOPS/NE_3RAD_malefemale/denovo_output_eqMF
dir=/lustre2/scratch/jcf87188/3RAD_popgen_allPOPS/NE_3RAD_malefemale

i=$(cat input_${SLURM_ARRAY_TASK_ID})

mkdir $dir/denovo_output_eqMF_v3/denovo_n_$i
denovo_map.pl --samples $dir/clean_combine/ --popmap $dir/popmap_male_female_equal.txt -o $dir/denovo_output_eqMF_v3/denovo_n_$i --paired -m 3 -n $i -M 4 -T 12 -X "populations: -r 0.7 -p 1 --fstats --fasta-loci --fasta-samples --vcf --genepop"