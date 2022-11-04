#!/bin/bash
#SBATCH --job-name=bt_5x_WI
#SBATCH --partition=highmem_p
#SBATCH --nodes=2 # use just 1
#SBATCH --ntasks=3
#SBATCH --mem=100gb # Per processor memory request
#SBATCH --export=NONE
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=jcf87188@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=3_bt_5x_WI.out
#SBATCH --error=3_bt_5x_WI.err

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
samtools view -bh -s 444.881057269 $btdir/WI08103.sorted.bam > $outdir/WI08103.5x.bam
samtools view -bh -s 444.777000777 $btdir/WI08110.sorted.bam > $outdir/WI08110.5x.bam
samtools view -bh -s 444.773156023 $btdir/WI08111.sorted.bam > $outdir/WI08111.5x.bam

samtools merge $mergedir/WI081_male.bam $outdir/WI08103.5x.bam $outdir/WI08110.5x.bam $outdir/WI08111.5x.bam
samtools sort $mergedir/WI081_male.bam -o $mergedir/WI081_male_sorted.bam
samtools index $mergedir/WI081_male_sorted.bam

##Louisiana Females
samtools view -bh -s 444.989297529 $btdir/WI08104.sorted.bam > $outdir/WI08104.5x.bam
cp $btdir/WI08105.sorted.bam $outdir/WI08105.5x.bam
samtools view -bh -s 444.871592662 $btdir/WI08106.sorted.bam > $outdir/WI08106.5x.bam

samtools merge $mergedir/WI081_female.bam $outdir/WI08104.5x.bam $outdir/WI08105.5x.bam $outdir/WI08106.5x.bam
samtools sort $mergedir/WI081_female.bam -o $mergedir/WI081_female_sorted.bam
samtools index $mergedir/WI081_female_sorted.bam

## calculated stats
pileup.sh in=$mergedir/WI081_male_sorted.bam out=$mergedir/pile/WI081_male_15x.pileup.txt
pileup.sh in=$mergedir/WI081_female_sorted.bam out=$mergedir/pile/WI081_female_15x.pileup.txt

## bam coverage
bamCoverage -b $mergedir/WI081_male_sorted.bam -o $mergedir/bigwig/WI081_male_15x.bw -of bigwig --binSize 1000
bamCoverage -b $mergedir/WI081_female_sorted.bam -o $mergedir/bigwig/WI081_female_15x.bw -of bigwig --binSize 1000

## Compare Females to Male
bigwigCompare -b1 $mergedir/bigwig/WI081_female_15x.bw -b2 $mergedir/bigwig/WI081_male_15x.bw --binSize 1000 -o $mergedir/bigwig/WI081_FM_15x.bw -of bigwig