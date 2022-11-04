library(tidyverse)
library(rstatix)

fst <- read_csv("/Users/juliafrederick/UGA_EHS Dropbox/Julia Frederick/Julia_Frederick_Projects/Ixscap_reseq_population_genetics/popgen/FST_comps.csv")
# fst_3RAD <- subset(fst, fst$Analysis == "3RAD")
# fst_gatk <- subset(fst, fst$Analysis == "GATK")
# fst_sam <- subset(fst, fst$Analysis == "Samtools")

cor(fst[,c('3RAD','GATK','Samtools')])
cor_test(fst[,c('3RAD','GATK','Samtools')])