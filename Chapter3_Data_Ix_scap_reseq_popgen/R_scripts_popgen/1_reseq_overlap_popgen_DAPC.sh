#!/bin/bash
#SBATCH --job-name=1_R_overlap_popgen_DAPC
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=75gb
#SBATCH --export=NONE
#SBATCH --time=7-00:00:00
#SBATCH --output=1_R_overlap_popgen_DAPC.out
#SBATCH --error=1_R_overlap_popgen_DAPC.err
#SBATCH --mail-user=jcf87188@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL

cd $SLURM_SUBMIT_DIR

#loading the software in the cluster
ml R/4.2.1-foss-2020b
ml Rstudio/2022.07.1-554 
ml rgdal/1.5-27-foss-2019b-R-4.0.0 
#ml rgdal/1.4-8-foss-2019b-R-4.0.0

R CMD BATCH 1_reseq_overlap_popgen_DAPC.R 