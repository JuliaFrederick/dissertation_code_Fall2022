#!/bin/bash
#SBATCH --job-name=map2ref_new
#SBATCH --partition=batch
#SBATCH --nodes=1 # use just 1
#SBATCH --ntasks=1 # use just 1
#SBATCH --cpus-per-task=4 # Number of OpenMP threads per job element
#SBATCH --mem-per-cpu=2000mb # Per processor memory request
#SBATCH --export=NONE
#SBATCH --time=96:00:00
#SBATCH --mail-user=jcf87188@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=1-133 # Here is where you set the number of job elements

cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#loading the software in the cluster
#ml Stacks/2.5-iccifort-2019.5.281
#ml Python/3.8.2-GCCcore-8.3.0
ml BWA/0.7.17-GCC-8.3.0
ml SAMtools/1.10-iccifort-2019.5.281

mkdir /lustre2/scratch/jcf87188/oct2021_plate4_PG/aligned_new_v2
mkdir /lustre2/scratch/jcf87188/oct2021_plate4_PG/aligned_new_v2/sorted

cleanreads=/lustre2/scratch/jcf87188/oct2021_plate4_PG/data_new/process_output_final_v2
bamoutput=/lustre2/scratch/jcf87188/oct2021_plate4_PG/aligned_new_v2

i=$(cat input_${SLURM_ARRAY_TASK_ID})

# run bwa and output BAM
bwa mem -t 4 /home/jcf87188/202107Pal_IscapGenome/GCF_016920785.2_ASM1692078v2_genomic.fna $cleanreads/$i\.1.fq.gz $cleanreads/$i\.2.fq.gz | samtools view -bS - > $bamoutput/$i\.bam

# filter BAM for primary mapping reads, imperfect matches, and >5 SNPs per read (NM:i:[0-5], below)
samtools view -h -q 25 -F 4 -F 256 $bamoutput/$i\.bam | grep -v XA:Z | grep -v SA:Z | awk '{if($0 ~ /^@/ || $6 ~ /140M/) {print $0}}' | grep -E '^@|NM:i:[0-5]\s' | samtools view -bS - > $bamoutput/$i\.q25.unique.perfect.bam

samtools sort $bamoutput/$i\.q25.unique.perfect.bam -o $bamoutput/sorted/$i\.bam