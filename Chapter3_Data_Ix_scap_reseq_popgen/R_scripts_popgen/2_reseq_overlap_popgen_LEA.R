### R scripts of Ixodes scapularis resequencing popgen - 3RAD project
### Created by Julia Frederick, reference websites listed per section
### Last updated August 19 2022

# Load Libraries
#library(here)
#setwd("~/Dropbox (UGA_EHS)/Julia_Frederick_Projects/Ixscap_reseq_population_genetics/popgen/samtools/popgen_samtools_R")

##################################################################################
################################ Basic Statistics ################################
##################################################################################
# VCF file from samtools SNP calling and post-VCFtools filtering
# max-missing 95%
# 6x - 50x
# maf 05
# Q 30
# GQ 25
# remove indels
# maximum allele count 2
## using these to just test pipeline.
## has 775191 sites

# Load Libraries
library(vcfR)
library(hierfstat)
library(pegas)
library(adegenet)
library(poppr)

# From LEA 
library(LEA)
library(miscTools)
library(stringr)
library(reshape2)
library(dplyr)
library(grid)
library(ggpubr)


# Load Data
vcf_file <- ("./0000.vcf.gz")
vcf <- read.vcfR(vcf_file)
pop.data <- read.csv("../reseq_strata.csv")
head(pop.data)
names(pop.data)[1] <- 'sample' # rename the first column because why should it import normally
pop.data <- pop.data[match(colnames(vcf@gt)[-1], pop.data$sample),] # put pop.data file in same order as the VCF file

####################################################################################
############################## STRUCTURE-like Analysis ##############################
####################################################################################
# https://github.com/Tom-Jenkins/admixture_pie_chart_map_tutorial
# Load Libraries
# Load packages
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("LEA")

#library(LEA)
#library(miscTools)
#library(stringr)
#library(reshape2)
#library(dplyr)
#library(grid)
#library(ggpubr)

# Input Data
gi.lea <- vcfR2genind(vcf) # So there is no associated pop information
ploidy(gi.lea) <- 2
gi.lea
nLoc(gi.lea) # number of loci
nPop(gi.lea) # number of sites
nInd(gi.lea) # number of individuals
summary(gi.lea$pop) # sample size

# ----------------- #
# Admixture analysis using SNMF from LEA
# ----------------- #

# Export genotypes in STRUCTURE format
source("./TJ_genind2structure_function.R")
data=gi.lea
locN <- 1:nLoc(gi.lea) ## renaming loci, the plus in some of the names causes an error
locNames(data) <- locN
genind2structure(data, file="genotypes", pops = FALSE, markers=FALSE)

# Convert STRUCTURE file to .geno format
struct2geno("genotypes", ploidy = 2, FORMAT = 2, extra.column = 1)

# Run snmf algorithm
set.seed(123)
snmf1 = snmf("genotypes.geno",
             K = 1:10, # number of K ancestral populations to run
             repetitions = 10, # ten repetitions for each K
             entropy = TRUE, # calculate cross-entropy
             project = "new")

# Load snmf project
snmf1 = load.snmfProject("genotypes.snmfProject")

# Plot cross-entropy results to assess optimal number of K
# Smaller values of cross-entropy usually mean better runs
# A plateau usually represents the K that best fits the data
options(bitmapType='cairo')
png(file="./reseq_overlap_structure_K.png",width=600, height=500)
plot(snmf1, col = "blue", cex = 1.5, pch = 19)
dev.off()

# Extract the cross-entropy of all runs where K = 2
ce = cross.entropy(snmf1, K = 2)
ce

# Find the run with the lowest cross-entropy
lowest.ce = which.min(ce)
lowest.ce

# Extract Q-matrix for the best run
qmatrix = as.data.frame(Q(snmf1, K = 2, run = lowest.ce))
head(qmatrix)

# Label column names of qmatrix
ncol(qmatrix)
cluster_names = c()
for (i in 1:ncol(qmatrix)){
  cluster_names[i] = paste("Cluster", i)
}
cluster_names
colnames(qmatrix) = cluster_names
head(qmatrix)

# Add individual IDs
qmatrix$Ind = indNames(gi.lea)

# Add site IDs
pop(gi.lea) <- pop.data$location
qmatrix$Site = gi.lea$pop
head(qmatrix)

# Convert dataframe to long format
qlong = melt(qmatrix, id.vars=c("Ind","Site"))
head(qlong)

# Change order of sites by using the factor function
# site.order = c("Vig","Ios","Cor","Mul","She","Cro","Hel","Flo","Lys","Ber")
# qlong$Site_ord = factor(qlong$Site, levels = site.order)

# Adjust facet labels
levels(qlong$Site)
facet.labs = as.character(unique(pop.data$location))
levels(qlong$Site) = facet.labs
levels(qlong$Site)

# Define colour palette
#pal = colorRampPalette(c("blue","yellow","red"))
#cols = pal(length(unique(qlong$variable))) # blue, grey, yellow, orange, red
cols = c("#7F7F7F","#FF7F00","#FFFF00","#FF0000","#0000FF")

# Plot admixture barplot 
admix.bar = ggplot(data=qlong, aes(x=Ind, y=value, fill=variable))+
  geom_bar(stat = "identity")+
  scale_y_continuous(expand = c(0,0))+
  facet_wrap(~Site, scales = "free", ncol = 3)+
  scale_fill_manual(values = cols)+
  ylab("Admixture proportion")+
  # xlab("Individual")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(colour="black", size=10),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
admix.bar
ggsave(admix.bar,filename="reseq_overlapdata_LEA_admixture_barplot", device="png", width=10, height=6, dpi=300)

save.image(file="./LEA_reseq_overlap_m90_6x50x_GQ30_Q30_maf05_ma2.RData")
