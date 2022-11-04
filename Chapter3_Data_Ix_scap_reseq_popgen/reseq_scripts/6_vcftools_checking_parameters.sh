#!/bin/bash

#mkdir ./reseq_all_calls_newPipe
cp reseq_all_calls_newPipe.vcf ./reseq_all_calls_newPipe
cd ./reseq_all_calls_newPipe

ml BCFtools/1.13-GCC-8.3.0
bgzip reseq_all_calls_newPipe.vcf
bcftools index reseq_all_calls_newPipe.vcf.gz
module unload BCFtools/1.13-GCC-8.3.0

ml VCFtools/0.1.16-GCC-8.3.0-Perl-5.30.0
VCF=./reseq_all_calls_newPipe.vcf.gz
OUT=./reseq_all_calls_newPipe

vcftools --gzvcf $VCF --freq2 --out $OUT --max-alleles 2
vcftools --gzvcf $VCF --depth --out $OUT
vcftools --gzvcf $VCF --site-mean-depth --out $OUT
vcftools --gzvcf $VCF --site-quality --out $OUT
vcftools --gzvcf $VCF --missing-site --out $OUT
vcftools --gzvcf $VCF --missing-indv --out $OUT
vcftools --gzvcf $VCF --het --out $OUT