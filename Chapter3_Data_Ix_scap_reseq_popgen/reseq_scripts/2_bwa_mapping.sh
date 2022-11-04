#!/bin/bash
#SBATCH --job-name=bwa_mapping
#SBATCH --partition=batch
#SBATCH --nodes=1 # use just 1
#SBATCH --ntasks=1 # use just 1
#SBATCH --cpus-per-task=4 # Number of OpenMP threads per job element
#SBATCH --mem-per-cpu=50gb # Per processor memory request
#SBATCH --export=NONE
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=jcf87188@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=1-30 # Here is where you set the number of job elements

cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#loading the software in the cluster
ml BWA/0.7.17-GCC-8.3.0
ml SAMtools/1.10-iccifort-2019.5.281

clean_reads=/lustre2/scratch/jcf87188/reseq_iscap/v2_clean_reads_2022/paired_Illumina_trim_reads
output=/lustre2/scratch/jcf87188/reseq_iscap/v2_clean_all_outputs/bwa_mapped

genome=/home/jcf87188/202107Pal_IscapGenome/GCF_016920785.2_ASM1692078v2_genomic.fna

i=$(cat input_${SLURM_ARRAY_TASK_ID})
# create the RG header for each sample to it gets in there while mapping
# this assume all data are from the same read-group (library)
HEADER=`printf @RG%sID:%s%sSM:%s%sPL:ILLUMINA '\\t' $i '\\t' $i '\\t'`

# run bwa and output BAM, sort that BAM and index it
bwa mem -t 4 -R "${HEADER}" $genome $clean_reads/$i\_R1_pair_trim.fastq.gz $clean_reads/$i\_R2_pair_trim.fastq.gz | samtools view -bS - > $output/$i\.bam
samtools sort $output/$i\.bam -o $output/properly_mapped_reads/$i\.sorted.bam
samtools index $output/properly_mapped_reads/$i\.sorted.bam

# filter for reads that did not map
samtools view -h -f 4 $output/$i\.bam | samtools view -bS - > $output/unmapped_reads/$i\.unmapped.bam