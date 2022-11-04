#!/bin/bash
#SBATCH --job-name=1_info_more
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2gb
#SBATCH --export=NONE
#SBATCH --time=96:00:00
#SBATCH --output=1_info_more_ams.out
#SBATCH --error=1_info_more_ams.err
#SBATCH --mail-user=jcf87188@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL

cd $SLURM_SUBMIT_DIR

#loading the software in the cluster
#ml Stacks/2.5-iccifort-2019.5.281
#ml Python/3.8.2-GCCcore-8.3.0
#ml BWA/0.7.17-GCC-8.3.0
ml SAMtools/1.10-iccifort-2019.5.281


## Here are the specific commands we are running

set -ueo pipefail
SAMPLES="FL06501
FL06502
FL06506
FL09701
FL09702
FL09704
FL09705
FL09707
FL09708
FL09710
FL09711
FL09712
FL09713
FL09715
LA01503
LA01504
LA01509
LA01510
LA01513
LA04701
LA04702
LA07905
LA07906
LA07907
LA07908
LA07910
MI06501
MI06502
MI06503
MI06505
MI06506
MI06507
MI06508
MI06509
MI06510
MI06511
MI06512
MO15913
NC16901
NC16902
NC16903
NC16905
NC16906
NC16908
NC16911
NC18906
NC18909
NH00307
NH00308
NH00309
OH08301
OH08302
OKC11901
OKC11902
OKC11903
OKC11904
OKC11905
OKC11906
OKC11907
OKC11908
OKC11910
OKC11911
OKC11912
OKW11901
OKW11902
OKW11903
OKW11904
OKW11905
OKW11906
OKW11907
OKW11908
OKW11909
OKW11910
OKW11911
OKW11912
RI00902
RI00904
RI00906
RI00907
RI00909
WI07501
WI07502
WI07503
WI07504
WI07505
WI07506
WI07507
WI07508
WI07509
WI07510
WI07511
WI07512
WI08102
WI08103
WI08104
WI08105
WI08106
WI08107
WI08108
WI08109
WI08110
WI08112
"
i=1

for i in $SAMPLES
do
printf "$i ams q25 \n"
samtools flagstat /lustre2/scratch/jcf87188/3RAD_popgen_allPOPS/ams_all_bam_files/$i\.ms.bam
done