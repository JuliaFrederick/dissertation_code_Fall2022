#!/bin/bash
#SBATCH --job-name=clean_map
#SBATCH --partition=batch
#SBATCH --nodes=1 # use just 1
#SBATCH --ntasks=1 # use just 1
#SBATCH --cpus-per-task=4 # Number of OpenMP threads per job element
#SBATCH --mem-per-cpu=20gb # Per processor memory request
#SBATCH --export=NONE
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=jcf87188@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=1-30 # Here is where you set the number of job elements

cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#loading the software in the cluster
ml Trimmomatic/0.39-Java-1.8.0_144

mkdir /scratch/jcf87188/reseq_iscap/clean_reads_2022

inputpath_r1=/scratch/jcf87188/reseq_iscap/raw_novaseq_data_aug2021/FastqFiles_NVS108B_Glenn_R1
inputpath_r2=/scratch/jcf87188/reseq_iscap/raw_novaseq_data_aug2021/FastqFiles_NVS108B_Glenn_R2
cleanpath=/scratch/jcf87188/reseq_iscap/v2_clean_reads_2022

mkdir /scratch/jcf87188/reseq_iscap/v2_clean_reads_2022/paired_Illumina_trim_reads/
mkdir /scratch/jcf87188/reseq_iscap/v2_clean_reads_2022/unpaired_Illumina_trim_reads/

i=$(cat input_${SLURM_ARRAY_TASK_ID})
#clean and trim
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 12 $inputpath_r1/$i\_R1.fastq.gz $inputpath_r2/$i\_R2.fastq.gz \
$cleanpath/paired_Illumina_trim_reads/$i\_R1_pair_trim.fastq.gz $cleanpath/unpaired_Illumina_trim_reads/$i\_R1_unpair_trim.fastq.gz  \
$cleanpath/paired_Illumina_trim_reads/$i\_R2_pair_trim.fastq.gz $cleanpath/unpaired_Illumina_trim_reads/$i\_R2_unpair_trim.fastq.gz  \
ILLUMINACLIP:/apps/eb/Trimmomatic/0.39-Java-1.8.0_144/adapters/TruSeq3-PE-2.fa:2:30:10:2:TRUE SLIDINGWINDOW:5:20 MINLEN:50
