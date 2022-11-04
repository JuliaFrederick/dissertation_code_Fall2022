#!/bin/bash
#SBATCH --job-name=gatk_
#SBATCH --partition=batch
#SBATCH --nodes=1 # use just 1
#SBATCH --ntasks=1 # use just 1
#SBATCH --cpus-per-task=4 # Number of OpenMP threads per job element
#SBATCH --mem-per-cpu=40gb # Per processor memory request
#SBATCH --export=NONE
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=jcf87188@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=1-30 # Here is where you set the number of job elements

cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#loading the software in the cluster
ml SAMtools/1.10-iccifort-2019.5.281
ml GATK/4.1.8.1-GCCcore-8.3.0-Java-1.8

genome=/home/jcf87188/202107Pal_IscapGenome/GCF_016920785.2_ASM1692078v2_genomic.fna
bam_file=/lustre2/scratch/jcf87188/reseq_iscap/v2_clean_all_outputs/bwa_mapped/properly_mapped_reads
output=/lustre2/scratch/jcf87188/reseq_iscap/v2_clean_all_outputs/gatk_duplicates

i=$(cat input_${SLURM_ARRAY_TASK_ID})
# run duplicate marking using Spark (in local mode) and setting the cores appropriately
# we are assuming the number of threads here will be 5
gatk --java-options "-Xmx16G" MarkDuplicatesSpark --spark-runner LOCAL --input $bam_file/$i\.sorted.bam --output $output/$i\.md.bam --conf 'spark.executor.cores=4' &&
gatk --java-options "-Xmx16G" SetNmMdAndUqTags --INPUT $output/$i\.md.bam --OUTPUT $output/$i\.md.fx.bam --REFERENCE_SEQUENCE $genome &&
gatk --java-options "-Xmx16G" BuildBamIndex --INPUT $output/$i\.md.fx.bam