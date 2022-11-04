#!/bin/bash
#SBATCH --job-name=refmapstacks_v4
#SBATCH --partition=batch
#SBATCH --nodes=2                     # Number of nodes
#SBATCH --ntasks=8                    # Number of MPI ranks
#SBATCH --ntasks-per-node=4           # Number of MPI ranks per node
#SBATCH --cpus-per-task=4             # Number of OpenMP threads for each MPI process/rank
#SBATCH --mem-per-cpu=2000mb          # Per processor memory request
#SBATCH --export=NONE
#SBATCH --time=7-00:00:00
#SBATCH --output=refmapstacks_v4.out
#SBATCH --error=refmapstacks_v4.err
#SBATCH --mail-user=jcf87188@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL

cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#loading the software in the cluster
#ml Stacks/2.55-foss-2019b
ml Stacks/2.5-iccifort-2019.5.281
#ml Python/3.8.2-GCCcore-8.3.0
#ml BWA/0.7.17-GCC-8.3.0
#ml SAMtools/1.10-iccifort-2019.5.281


ref_map.pl --samples /lustre2/scratch/jcf87188/3RAD_popgen_allPOPS/ams_all_bam_files --popmap /lustre2/scratch/jcf87188/3RAD_popgen_allPOPS/populationmaps_best_noOKC.txt --out-path /lustre2/scratch/jcf87188/3RAD_popgen_allPOPS/stacks_ams_bam_pops_best_noOKC -X "populations: -r 0.6 -p 2 --hwe --fasta-loci --fasta-samples --fasta-samples-raw --vcf --fstats --genepop --vcf --phylip --phylip-var --structure --treemix"