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
vcf_file <- ("./samtools_formatQ_m90_6x50x_GQ30_Q30_maf05_ma2.recode.vcf.gz")
vcf <- read.vcfR(vcf_file)
pop.data <- read.csv("../reseq_strata.csv")
head(pop.data)
names(pop.data)[1] <- 'sample' # rename the first column because why should it import normally
pop.data <- pop.data[match(colnames(vcf@gt)[-1], pop.data$sample),] # put pop.data file in same order as the VCF file

###################################################################################
################################# Population DAPC #################################
###################################################################################
# Creating Input Variables
gl <- vcfR2genlight(vcf)
ploidy(gl) <- 2
pop(gl) <- pop.data$location # Assign Populations

pop(gl)
gl

# Running DAPC
dapc.pop <- adegenet::dapc(gl, n.da=100, n.pca=75, var.contrib=TRUE)
temp <- adegenet::optim.a.score(dapc.pop)
best<-as.integer(temp$best)
best

dapc.pop <- adegenet::dapc(gl, n.da=100, n.pca=best)
dapc.pop

# Plotting DAPC
options(bitmapType='cairo')
png(file="./reseq_samtoolsdata_DAPC_populations.png", width = 700, 600)
dapc.plot <- scatter(dapc.pop, cex = 2, clab= 0,
                     posi.pca = "topleft", scree.pca = TRUE,
                     posi.da = "bottomleft",
                     legend = TRUE,
                     #xlim=c(15,-10), ylim=c(25,-15),
                     inset.solid = 0.7)
dev.off()

png(file="./varcontrib_DAPC_samtools.png", width=1000,900)
contrib.s <- loadingplot(dapc.pop$var.contr, axis =2, lab.jitter=1)
dev.off()

contrib.s$var.names

####################################################################################
############################## Population DAPC w/o FL ##############################
####################################################################################
# Creating Input Variables
FL097 <- subset(pop.data, pop.data$state == "FL")
FL097 <- as.character(FL097$sample)
gl.rm <- gl[!(indNames(gl) %in% FL097)]

# Running DAPC
dapc.rm <- adegenet::dapc(gl.rm, n.da=100, n.pca=75, var.contrib=TRUE)
temp <- adegenet::optim.a.score(dapc.rm)
best<-as.integer(temp$best)
best

dapc.rm <- adegenet::dapc(gl.rm, n.da=100, n.pca=best)
dapc.rm

# Plot DAPC Clusters
png(file="./freseq_samtoolsdata_DAPC_rmFL097.png", width = 700,600)
dapc.plot.rm <- scatter(dapc.rm, cex = 2, clab= 0,
                        posi.pca = "topleft", scree.pca = TRUE,
                        posi.da = "bottomleft",
                        legend = TRUE, cleg = 1,
                        #xlim=c(10,-10), ylim=c(10,-12),
                        inset.solid = 0.7)
dev.off()

png(file="./varcontrib_DAPC_samtools_woFL097.png", width=1000,900)
contrib.rm <- loadingplot(dapc.rm$var.contr, axis =2, lab.jitter=1, thres=0.001)
dev.off()

contrib.rm$var.names
# 43 sites

####################################################################################
################################# Basic Statistics #################################
####################################################################################
# Setting up the genind 
gi <- vcfR2genind(vcf)
ploidy(gi) <- 2
pop(gi) <- pop.data$location

# Basic Stats from Hierfstat
bs.gi<-basic.stats(gi,diploid = TRUE)

## Heterozyogosity Obs per pop
Ho <- bs.gi[["Ho"]]
Ho_avg <- apply(Ho,2, mean, na.rm=TRUE)  # applies function 'mean' to 2nd dimension (columns)

## Heterozygosity Expected (mean gene diversity) per pop
Hs <- bs.gi[["Hs"]]
Hs_avg <- apply(Hs,2, mean, na.rm=TRUE)

## Inbreeding Coefficent per pop
Fis <- bs.gi[["Fis"]]
Fis_avg <- apply(Fis,2, mean, na.rm=TRUE)  

## Allelic richness per pop (hierfstat)
ar <- allelic.richness(gi, diploid=TRUE)
ar.2 <- ar[["Ar"]]
ar_avg <- apply(ar.2,2, mean, na.rm=TRUE)  

## Number of Alleles (hierfstat)
na<- nb.alleles(gi, diploid = TRUE)
na_sum <- apply(na,2, sum, na.rm=TRUE)

## Getting Total number of samples
totaln <- as.numeric(table(gi@pop))

## Getting Pop Names
pops <- as.character( unique(pop(gi)) )

sum_stats_pops <- data.frame("Populations"=pops, "N"=totaln, "Nalleles"=na_sum, 
                             "Ar" = as.numeric(ar_avg), "Ho"=as.numeric(Ho_avg),
                             "Hs"=as.numeric(Hs_avg), "Fis"=as.numeric(Fis_avg))
write.csv(sum_stats_pops,"./reseq_samtoolsdata_sum_stats_pops.csv")

# Pairwise Fst per population
hgi<- genind2hierfstat(gi)

pairnei <- pairwise.neifst(hgi) # Use this data
write.csv(pairnei,file="./reseq_samtoolsdata_pairnei.csv")
pairwc <- pairwise.WCfst(hgi) # This data is for comparison purposes
write.csv(pairwc,file="./reseq_samtoolsdata_pairwc.csv")

save.image(file="./DAPC_BStats_reseq_samtools_formatQ_m90_6x50x_GQ30_Q30_maf05_ma2.RData")
