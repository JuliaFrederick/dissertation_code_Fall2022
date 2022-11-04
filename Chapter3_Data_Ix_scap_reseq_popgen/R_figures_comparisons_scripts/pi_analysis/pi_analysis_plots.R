## vcftools --gzvcf ./dir/0000.vcf.gz --window-pi 1000 --out overlap_pi_1kb
library(ggplot2)

# setwd("~/Dropbox (UGA_EHS)/Julia_Frederick_Projects/Ixscap_reseq_population_genetics/popgen/R_scripts/pi_analysis")

top14 <- c("NW_024609835.1","NW_024609836.1","NW_024609837.1","NW_024609838.1","NW_024609839.1","NW_024609846.1",
           "NW_024609857.1","NW_024609868.1","NW_024609879.1","NW_024609880.1","NW_024609881.1","NW_024609882.1",
           "NW_024609883.1","NW_024609884.1")
top4 <- c("NW_024609857.1","NW_024609868.1","NW_024609879.1","NW_024609883.1")

pi.3RAD <- read.table("./3RAD_pi_1kb.windowed.pi",header=T)
pi.3RAD.chr14 <- subset(pi.3RAD, CHROM %in%  top14)

ggplot(pi.3RAD.chr14) +
  geom_point(aes(x=BIN_START, y=PI))+
  facet_wrap(.~CHROM)

pi.gatk <- read.table("./gatk_pi_1kb.windowed.pi",header=T)
pi.gatk.chr14 <- subset(pi.gatk, CHROM %in%  top14)

ggplot(pi.gatk.chr14) +
  geom_point(aes(x=BIN_START, y=PI))+
  facet_wrap(.~CHROM)

pi.sam <- read.table("./samtools_pi_1kb.windowed.pi",header=T)
pi.sam.chr14 <- subset(pi.sam, CHROM %in%  top14)

ggplot(pi.sam.chr14) +
  geom_point(aes(x=BIN_START, y=PI))+
  facet_wrap(.~CHROM)

pi.over <- read.table("./overlap_pi_1kb.windowed.pi",header=T)
pi.over.chr14 <- subset(pi.over, CHROM %in%  top14)

ggplot(pi.over.chr14) +
  geom_point(aes(x=BIN_START, y=PI))+
  facet_wrap(.~CHROM)

pi.3RAD.chr14$dataset <- "3RAD"
pi.gatk.chr14$dataset <- "GATK"
pi.sam.chr14$dataset <- "Samtools"
pi.over.chr14$dataset <- "Overlap"
pi.all.14 <- rbind(pi.3RAD.chr14, pi.gatk.chr14, pi.sam.chr14, pi.over.chr14)

ggplot(pi.all.14) +
  geom_point(aes(x=BIN_START, y=PI, col=dataset), alpha=0.7)+
  facet_wrap(.~CHROM, scales = "free_x") +
  scale_color_manual(breaks=c("3RAD","GATK","Samtools","Overlap"),
                     values=c("#7FB3EE","#332851","#CA3074","#F6C667")) +
  theme_bw() + 
  labs(x="Position on Scaffold using 1000bp Bins", y="Nucleotide Diversity (\u03C0)") +
  guides(color = guide_legend(override.aes = list(size=3))) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        strip.text.x = element_text(size=12),
        legend.position = "bottom", legend.box = "horizontal", 
        legend.text = element_text(size=10),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-10,-10,0,-10),
        legend.title = element_blank())
