#!/bin/bash
#SBATCH --job-name=bt_5x_RI
#SBATCH --partition=highmem_p
#SBATCH --nodes=2 # use just 1
#SBATCH --ntasks=3
#SBATCH --mem=100gb # Per processor memory request
#SBATCH --export=NONE
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=jcf87188@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=3_bt_5x_RI.out
#SBATCH --error=3_bt_5x_RI.err

cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#loading the software in the cluster
ml SAMtools/1.10-iccifort-2019.5.281
#ml Bowtie2/2.4.1-GCC-8.3.0 
#ml seqtk/1.3-GCC-8.3.0
ml BBMap/38.93-GCC-8.3.0 
ml deepTools/3.5.1-intel-2020b-Python-3.8.6

btdir=/lustre2/scratch/jcf87188/reseq_iscap/xx_xy_cov/v2_clean_outputs/bt_mapped_new
outdir=/lustre2/scratch/jcf87188/reseq_iscap/xx_xy_cov/v2_clean_outputs/bt_mapped_5x
mergedir=/lustre2/scratch/jcf87188/reseq_iscap/xx_xy_cov/v2_clean_outputs/bt_mapped_5x/merged

##Louisiana Males
samtools view -bh -s 444.760225027 $btdir/RI00906.sorted.bam > $outdir/RI00906.5x.bam
samtools view -bh -s 444.823045267 $btdir/RI00907.sorted.bam > $outdir/RI00907.5x.bam
samtools view -bh -s 444.833611204 $btdir/RI00909.sorted.bam > $outdir/RI00909.5x.bam

samtools merge $mergedir/RI009_male.bam $outdir/RI00906.5x.bam $outdir/RI00907.5x.bam $outdir/RI00909.5x.bam
samtools sort $mergedir/RI009_male.bam -o $mergedir/RI009_male_sorted.bam
samtools index $mergedir/RI009_male_sorted.bam

##Louisiana Females
samtools view -bh -s 444.572016932 $btdir/RI00911.sorted.bam > $outdir/RI00911.5x.bam
samtools view -bh -s 444.551389502 $btdir/RI00916.sorted.bam > $outdir/RI00916.5x.bam
samtools view -bh -s 444.558098002 $btdir/RI00926.sorted.bam > $outdir/RI00926.5x.bam

samtools merge $mergedir/RI009_female.bam $outdir/RI00911.5x.bam $outdir/RI00916.5x.bam $outdir/RI00926.5x.bam
samtools sort $mergedir/RI009_female.bam -o $mergedir/RI009_female_sorted.bam
samtools index $mergedir/RI009_female_sorted.bam

## calculated stats
pileup.sh in=$mergedir/RI009_male_sorted.bam out=$mergedir/pile/RI009_male_15x.pileup.txt
pileup.sh in=$mergedir/RI009_female_sorted.bam out=$mergedir/pile/RI009_female_15x.pileup.txt

## bam coverage
bamCoverage -b $mergedir/RI009_male_sorted.bam -o $mergedir/bigwig/RI009_male_15x.bw -of bigwig --binSize 1000
bamCoverage -b $mergedir/RI009_female_sorted.bam -o $mergedir/bigwig/RI009_female_15x.bw -of bigwig --binSize 1000

## Compare Females to Male
bigwigCompare -b1 $mergedir/bigwig/RI009_female_15x.bw -b2 $mergedir/bigwig/RI009_male_15x.bw --binSize 1000 -o $mergedir/bigwig/RI009_FM_15x.bw -of bigwig