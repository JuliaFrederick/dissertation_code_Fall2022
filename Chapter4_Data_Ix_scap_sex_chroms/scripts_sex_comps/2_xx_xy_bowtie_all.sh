#!/bin/bash
#SBATCH --job-name=bt_xx_xy
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
ml SAMtools/1.10-iccifort-2019.5.281
ml Bowtie2/2.4.1-GCC-8.3.0 
#ml seqtk/1.3-GCC-8.3.0
ml BBMap/38.93-GCC-8.3.0 
ml deepTools/3.5.1-intel-2020b-Python-3.8.6

fqdir=/lustre2/scratch/jcf87188/reseq_iscap/v2_clean_reads_2022/paired_Illumina_trim_reads
subselect=/lustre2/scratch/jcf87188/reseq_iscap/xx_xy_cov/fq_downselected
outdir=/lustre2/scratch/jcf87188/reseq_iscap/xx_xy_cov/v2_clean_outputs/bt_mapped
piledir=/lustre2/scratch/jcf87188/reseq_iscap/xx_xy_cov/v2_clean_outputs/bbmap_pileup
bwdir=/lustre2/scratch/jcf87188/reseq_iscap/xx_xy_cov/v2_clean_outputs/bigwig

genome=/home/jcf87188/202107Pal_IscapGenome/GCF_016920785.2_ASM1692078v2_genomic

#bowtie2-build $genome /home/jcf87188/202107Pal_IscapGenome/GCF_016920785.2_ASM1692078v2_genomic

i=$(cat input_${SLURM_ARRAY_TASK_ID})
# subsample to 61 million reads
#reformat.sh in=$fqdir/$i\_R1_pair_trim.fastq.gz out=$subselect/$i\_R1_sub.fq samplereadstarget=14500000 sampleseed=444
#reformat.sh in=$fqdir/$i\_R2_pair_trim.fastq.gz out=$subselect/$i\_R2_sub.fq samplereadstarget=14500000 sampleseed=444
#seqtk sample -s100 $fqdir/$i\_R1_pair_trim.fastq.gz 14500000 > $subselect/$i\_R1_sub.fq.gz
#seqtk sample -s100 $fqdir/$i\_R2_pair_trim.fastq.gz 14500000 > $subselect/$i\_R2_sub.fq.gz

# map to genome using bowtie
bowtie2 -x $genome -1 $fqdir/$i\_R1_pair_trim.fastq.gz -2 $fqdir/$i\_R2_pair_trim.fastq.gz --very-sensitive | samtools view -bS - > $outdir/$i\.bam
samtools view -h -F 4 -F 256 $outdir/$i\.bam | samtools view -bS - > $outdir/$i\.propmapped.bam
samtools sort $outdir/$i\.propmapped.bam -o $outdir/$i\.sorted.bam
samtools index $outdir/$i\.sorted.bam

# calculated stats
pileup.sh in=$outdir/$i\.sorted.bam out=$piledir/$i\_all.pileup.txt

# bam coverage
bamCoverage -b $outdir/$i\.sorted.bam -o $bwdir/$i\_all.bw -of bigwig --binSize 1000