#!/bin/bash
#SBATCH --job-name=bt_5x_SC
#SBATCH --partition=highmem_p
#SBATCH --nodes=2 # use just 1
#SBATCH --ntasks=3
#SBATCH --mem=100gb # Per processor memory request
#SBATCH --export=NONE
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=jcf87188@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=3_bt_5x_SC.out
#SBATCH --error=3_bt_5x_SC.err

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
samtools view -bh -s 444.736160188 $btdir/SC00303.sorted.bam > $outdir/SC00303.5x.bam
samtools view -bh -s 444.876731545 $btdir/SC00304.sorted.bam > $outdir/SC00304.5x.bam
samtools view -bh -s 444.559722378 $btdir/SC00308.sorted.bam > $outdir/SC00308.5x.bam

samtools merge $mergedir/SC003_male.bam $outdir/SC00303.5x.bam $outdir/SC00304.5x.bam $outdir/SC00308.5x.bam
samtools sort $mergedir/SC003_male.bam -o $mergedir/SC003_male_sorted.bam
samtools index $mergedir/SC003_male_sorted.bam

##Louisiana Females
samtools view -bh -s 444.800897005 $btdir/SC00306.sorted.bam > $outdir/SC00306.5x.bam
samtools view -bh -s 444.579508577 $btdir/SC00314.sorted.bam > $outdir/SC00314.5x.bam
samtools view -bh -s 444.517116558 $btdir/SC00315.sorted.bam > $outdir/SC00315.5x.bam

samtools merge $mergedir/SC003_female.bam $outdir/SC00306.5x.bam $outdir/SC00314.5x.bam $outdir/SC00315.5x.bam
samtools sort $mergedir/SC003_female.bam -o $mergedir/SC003_female_sorted.bam
samtools index $mergedir/SC003_female_sorted.bam

## calculated stats
pileup.sh in=$mergedir/SC003_male_sorted.bam out=$mergedir/pile/SC003_male_15x.pileup.txt
pileup.sh in=$mergedir/SC003_female_sorted.bam out=$mergedir/pile/SC003_female_15x.pileup.txt

## bam coverage
bamCoverage -b $mergedir/SC003_male_sorted.bam -o $mergedir/bigwig/SC003_male_15x.bw -of bigwig --binSize 1000
bamCoverage -b $mergedir/SC003_female_sorted.bam -o $mergedir/bigwig/SC003_female_15x.bw -of bigwig --binSize 1000

## Compare Females to Male
bigwigCompare -b1 $mergedir/bigwig/SC003_female_15x.bw -b2 $mergedir/bigwig/SC003_male_15x.bw --binSize 1000 -o $mergedir/bigwig/SC003_FM_15x.bw -of bigwig