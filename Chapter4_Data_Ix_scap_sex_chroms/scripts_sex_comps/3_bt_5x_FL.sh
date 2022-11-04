#!/bin/bash
#SBATCH --job-name=bt_5x_FL
#SBATCH --partition=highmem_p
#SBATCH --nodes=2 # use just 1
#SBATCH --ntasks=3
#SBATCH --mem=100gb # Per processor memory request
#SBATCH --export=NONE
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=jcf87188@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=3_bt_5x_FL.out
#SBATCH --error=3_bt_5x_FL.err

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

##Florida Males
cp $btdir/FL09702.sorted.bam $outdir/FL09702.5x.bam
cp $btdir/FL09708.sorted.bam $outdir/FL09708.5x.bam
samtools view -bh --subsample-seed 444 --subsample 0.957297865 $btdir/FL09710.sorted.bam > $outdir/FL09710.5x.bam

samtools merge -o $mergedir/FL097_male.bam $outdir/FL09702.5x.bam $outdir/FL09708.5x.bam $outdir/FL09710.5x.bam
samtools sort $mergedir/FL097_male.bam -o $mergedir/FL097_male_sorted.bam
samtools index $mergedir/FL097_male_sorted.bam

##Florida Females
samtools view -bh --subsample-seed 444 --subsample 0.729927007 $btdir/FL09711.sorted.bam > $outdir/FL09711.5x.bam
samtools view -bh --subsample-seed 444 --subsample 0.364219114 $btdir/FL09712.sorted.bam > $outdir/FL09712.5x.bam
samtools view -bh --subsample-seed 444 --subsample 0.676132522 $btdir/FL09715.sorted.bam > $outdir/FL09715.5x.bam

samtools merge -o $mergedir/FL097_female.bam $outdir/FL09711.5x.bam $outdir/FL09712.5x.bam $outdir/FL09715.5x.bam
samtools sort $mergedir/FL097_female.bam -o $mergedir/FL097_female_sorted.bam
samtools index $mergedir/FL097_female_sorted.bam

## calculated stats
pileup.sh in=$mergedir/FL097_male_sorted.bam out=$mergedir/pile/FL097_male_15x.pileup.txt
pileup.sh in=$mergedir/FL097_female_sorted.bam out=$mergedir/pile/FL097_female_15x.pileup.txt

## bam coverage
bamCoverage -b $mergedir/FL097_male_sorted.bam -o $mergedir/bigwig/FL097_male_15x.bw -of bigwig --binSize 1000
bamCoverage -b $mergedir/FL097_female_sorted.bam -o $mergedir/bigwig/FL097_female_15x.bw -of bigwig --binSize 1000

## Compare Females to Male
bigwigCompare -b1 $mergedir/bigwig/FL097_female_15x.bw -b2 $mergedir/bigwig/FL097_male_15x.bw --binSize 1000 -o $mergedir/bigwig/FL097_FM_15x.bw --of bigwig