#!/bin/bash
#SBATCH --job-name=bt_5x_LA
#SBATCH --partition=highmem_p
#SBATCH --nodes=2 # use just 1
#SBATCH --ntasks=3
#SBATCH --mem=100gb # Per processor memory request
#SBATCH --export=NONE
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=jcf87188@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=3_bt_5x_LA.out
#SBATCH --error=3_bt_5x_LA.err

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
samtools view -bh -s 444.792267469 $btdir/LA07902.sorted.bam > $outdir/LA07902.5x.bam
samtools view -bh -s 444.936154278 $btdir/LA07904.sorted.bam > $outdir/LA07904.5x.bam
samtools view -bh -s 444.67558438 $btdir/LA07907.sorted.bam > $outdir/LA07907.5x.bam

samtools merge $mergedir/LA079_male.bam $outdir/LA07902.5x.bam $outdir/LA07904.5x.bam $outdir/LA07907.5x.bam
samtools sort $mergedir/LA079_male.bam -o $mergedir/LA079_male_sorted.bam
samtools index $mergedir/LA079_male_sorted.bam

##Louisiana Females
samtools view -bh -s 444.707313623 $btdir/LA07910.sorted.bam > $outdir/LA07910.5x.bam
samtools view -bh -s 444.790888959 $btdir/LA07911.sorted.bam > $outdir/LA07911.5x.bam
samtools view -bh -s 444.647668394 $btdir/LA07914.sorted.bam > $outdir/LA07914.5x.bam

samtools merge $mergedir/LA079_female.bam $outdir/LA07910.5x.bam $outdir/LA07911.5x.bam $outdir/LA07914.5x.bam
samtools sort $mergedir/LA079_female.bam -o $mergedir/LA079_female_sorted.bam
samtools index $mergedir/LA079_female_sorted.bam

## calculated stats
pileup.sh in=$mergedir/LA079_male_sorted.bam out=$mergedir/pile/LA079_male_15x.pileup.txt
pileup.sh in=$mergedir/LA079_female_sorted.bam out=$mergedir/pile/LA079_female_15x.pileup.txt

## bam coverage
bamCoverage -b $mergedir/LA079_male_sorted.bam -o $mergedir/bigwig/LA079_male_15x.bw -of bigwig --binSize 1000
bamCoverage -b $mergedir/LA079_female_sorted.bam -o $mergedir/bigwig/LA079_female_15x.bw -of bigwig --binSize 1000

## Compare Females to Male
bigwigCompare -b1 $mergedir/bigwig/LA079_female_15x.bw -b2 $mergedir/bigwig/LA079_male_15x.bw --binSize 1000 -o $mergedir/bigwig/LA079_FM_15x.bw -of bigwig