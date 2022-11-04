### R scripts of Ixodes scapularis 3RAD project
### Created by Julia Frederick, reference websites listed per section
### Last updated February 7 2022

# Load Libraries
library(here)

##################################################################################
################################### Sample Map ###################################
##################################################################################
# This data set has 353 I. scapularis adult ticks

# Load libraries
library(scatterpie)
library(rgdal)
library(ggplot2)
library(tidyverse)

# Map of Samples
load(here("./data/DataForJulia_v02.RData")) # from Guha
locations <- read.csv(here("./data/353samples.csv")) # from my files

# Shape File of Endemic Levels
my_spdf <- readOGR(dsn=here("./data/SHAPES"), layer="LymeIncidence_Kernel95pc") # from Guha

# Set map variables
usa <- US.Map.Extent.LonLat
states <- US.State.Cropped
counties <- US.County.Cropped

# Plot map of I. scapularis sample locations with pie charts for sex breakdown
Map_3RAD.color <- ggplot(data = states, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color = "black", fill = "white",size=.1) + 
  geom_polygon(data = my_spdf, aes(x = long, y = lat, group = group), 
               colour = "#bdbdbd", fill = NA) +
  geom_polygon(color="black", fill=NA,size=.1) + #puts state borders back on top
  geom_scatterpie(data=locations, aes(x=longitude, y=latitude, group = location, r=.65), 
                  cols = colnames(locations[,c(7:8)]), alpha=.9) +
  scale_fill_manual(name="Sex of Samples",
                    values=c('#67a9cf','#ef8a62'),
                    labels=c('male','female')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.85,0.25),
        legend.title = element_text(size=11),
        legend.text = element_text(size=9),
        legend.key.size = unit(4.5,'mm')) +
  geom_hline(yintercept = PopGroup.hLine.Lat, size =.65, linetype = "dashed") +
  geom_vline(xintercept = PopGroup.vLine.Lon, size =.65, linetype = "dashed")
Map_3RAD.color
ggsave(here("./figures_formatted/Fig1a_3RADmap_pieCharts_color.pdf"), 
       width = 112, height = 112, units = "mm", dpi=600)

Map_3RAD.bw <- ggplot(data = states, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color = "black", fill = "white") + 
  geom_polygon(data = my_spdf, aes(x = long, y = lat, group = group), 
               colour = "#bdbdbd", fill = NA) +
  geom_polygon(color="black", fill=NA) + #puts state borders back on top
  geom_scatterpie(data=locations, aes(x=longitude, y=latitude, group = location, r=.65), 
                  cols = colnames(locations[,c(7:8)]), alpha=.9) +
  scale_fill_manual(name="Sex of Samples",
                    values=c('#252525','#f0f0f0'),
                    labels=c('male','female')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.85,0.25),
        legend.title = element_text(size=11),
        legend.text = element_text(size=9),
        legend.key.size = unit(6,'mm')) +
  geom_hline(yintercept = PopGroup.hLine.Lat, size =.65, linetype = "dashed") +
  geom_vline(xintercept = PopGroup.vLine.Lon, size =.65, linetype = "dashed")
Map_3RAD.bw
ggsave(here("./figures_formatted/Fig1a_3RADmap_pieCharts_color.pdf"), 
       width = 112, height = 112, units = "mm", dpi=600)

Map_3RAD.c <- ggplot(data = states, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color = "black", fill = "white") + 
  geom_polygon(data = my_spdf, aes(x = long, y = lat, group = group), 
               colour = "papayawhip", fill = NA) +
  geom_polygon(color="black", fill=NA) + #puts state borders back on top
  geom_scatterpie(data=locations, aes(x=longitude, y=latitude, group = location, r=.6), 
                  cols = colnames(locations[,c(7:8)]), alpha=.8) +
  scale_fill_manual(name="Sex of Samples",
                    values=c('blue','red'),
                    labels=c('male','female')) +
  xlab("Longitude")+
  ylab("Latitude")+
  theme(
    axis.text = element_text(colour="black", size=12),
    axis.title = element_text(colour="black", size=14),
    panel.background = element_rect(fill="lightsteelblue2"),
    panel.border = element_rect(fill = NA, colour = "black", size = 0.5),
    panel.grid.minor = element_line(colour="grey90", size=0.5),
    panel.grid.major = element_line(colour="grey90", size=0.5),
    legend.text = element_text(size=12),
    legend.title = element_blank(),
    legend.key.size = unit(0.7, "cm"),
    legend.position = "top") +
  geom_hline(yintercept = PopGroup.hLine.Lat, size =.75, linetype = "dashed") +
  geom_vline(xintercept = PopGroup.vLine.Lon, size =.75, linetype = "dashed")
Map_3RAD.c
ggsave(here("./figures/Fig1a_3RADmap_pieCharts_details.png"), width = 8, height = 8, dpi=400)

# Plot map of I. scapularis sample locations, no sex breakdown
Map_3RAD_sp <- ggplot(data = states, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color = "black", fill = "white") + 
  geom_polygon(data = my_spdf, aes(x = long, y = lat, group = group), 
               colour = "papayawhip", fill = NA) +
  geom_polygon(color="black", fill=NA) + #puts state borders back on top
  geom_point(data=locations,aes(x=longitude, y=latitude, group=location),
             size =3, color="red")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank()) +
  geom_hline(yintercept = PopGroup.hLine.Lat, size =.75, linetype = "dashed") +
  geom_vline(xintercept = PopGroup.vLine.Lon, size =.75, linetype = "dashed")
Map_3RAD_sp
ggsave(here("./figures/Fig1b_3RADmap_scatter.png"), width = 8, height = 8, dpi=400)

Map_3RAD_sp.c <- ggplot(data = states, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color = "black", fill = "white") + 
  geom_polygon(data = my_spdf, aes(x = long, y = lat, group = group), 
               colour = "papayawhip", fill = NA) +
  geom_polygon(color="black", fill=NA) + #puts state borders back on top
  geom_point(data=locations,aes(x=longitude, y=latitude, group=location),
             size =3, color="red")+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(
    axis.text = element_text(colour="black", size=12),
    axis.title = element_text(colour="black", size=14),
    panel.background = element_rect(fill="lightsteelblue2"),
    panel.border = element_rect(fill = NA, colour = "black", size = 0.5),
    panel.grid.minor = element_line(colour="grey90", size=0.5),
    panel.grid.major = element_line(colour="grey90", size=0.5),
    legend.text = element_text(size=12),
    legend.title = element_blank(),
    legend.key.size = unit(0.7, "cm"),
    legend.position = "top") +
  geom_hline(yintercept = PopGroup.hLine.Lat, size =.75, linetype = "dashed") +
  geom_vline(xintercept = PopGroup.vLine.Lon, size =.75, linetype = "dashed")
Map_3RAD_sp.c
ggsave(here("./figures/Fig1b_3RADmap_scatter_details.png"), width = 8, height = 8, dpi=400)

##################################################################################
################################ Basic Statistics ################################
##################################################################################
# VCF file from Stacks output and post-VCFtools filterins
  # 6x-200x
  # Minimum Allele Frequency = 5%
  # Max Missing data = 80% (80% of samples have to have that locus)
  # Indels Removed
  # Total sites = 7274

# Load Libraries
library(vcfR)
library(hierfstat)
library(pegas)

# Load Data
vcf_file <- here("./data/rmInd_p20_6x200x_maf05.recode.vcf")
vcf <- read.vcfR(vcf_file)
pop.data <- read.csv(here("./data/strat_info_353ams.csv"))
names(pop.data)[1] <- 'sample' # rename the first column because why should it import normally
pop.data <- pop.data[match(colnames(vcf@gt)[-1], pop.data$sample),] # put pop.data file in same order as the VCF file

# Setting up the genind 
gi <- vcfR2genind(vcf)
ploidy(gi) <- 2
pop(gi) <- pop.data$county

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
  pops <- colnames( bs.gi[["pop.freq"]][[1]])

sum_stats_pops <- data.frame("Populations"=pops, "N"=totaln, "Nalleles"=na_sum, 
                             "Ar" = as.numeric(ar_avg), "Ho"=as.numeric(Ho_avg),
                             "Hs"=as.numeric(Hs_avg), "Fis"=as.numeric(Fis_avg))
write.csv(sum_stats_pops,here("./tables/sum_stats_pops.csv"))

# Testing significance of stat differences between North Vs South
labelled_ssp <- sum_stats_pops %>% dplyr::mutate(NvS = dplyr::if_else(
  Populations %in% c("Aiken County, SC","Bossier County, LA","Claiborne County, TN","Iberville County, LA",
                     "Jefferson County, FL","Lauderdale County, AL","Osceola County, FL","Payne County (wild), OK",
                     "Rapides County, LA","Stokes County, NC","Walker County, GA","Watauga Count, NC"),
  "South", "North"))

temp <- subset(pop.data, endemic95=="In")
labelled_ssp <- labelled_ssp %>% dplyr::mutate(in95 = dplyr::if_else(
  Populations %in% levels(factor(temp$county)),
  "In", "Out"))

library(ggpubr)
Ar_bp <- ggboxplot(
  labelled_ssp, x = "NvS", y = "Ar",
  ylab = "Allelic Richness", xlab = "North vs South", add = "jitter"
)
png(file=here("./figures/Ar_bp.png"), width = 500, height = 500)
Ar_bp
dev.off()
Ar_tt <- t.test(Ar ~ in95, data=labelled_ssp) %>% broom::tidy()

Hs_bp <- ggboxplot(
  labelled_ssp, x = "NvS", y = "Hs",
  ylab = "Heterozygosity", xlab = "North vs South", add = "jitter"
)
png(file=here("./figures/Hs_bp.png"), width = 500, height = 500)
Hs_bp
dev.off()

# Testing Ho vs Hs - Bartlett's Test (hierfstat)
div <- summary(gi)
bartlett.test(list(div$Hexp, div$Hobs)) # a test : H0: Hexp = Hobs
    # Output - Bartlett's K-squared = 887.14, df = 1, p-value < 2.2e-16

# Per locus HW Test (pegas)
HWlocus <- hw.test(gi, B=1000)
HWdf <- as.data.frame(HWlocus)
nrow(HWdf)
  # 7274
HWgood <- subset(HWdf, Pr.exact > 0.005)
nrow(HWgood) 
  # 1272 in HW eq

###################################################################################
################################### AMOVA & FST ###################################
###################################################################################
# Load Libraries
library(poppr)

# AMOVA State/County
gis <- vcfR2genind(vcf)
ploidy(gis) <- 2
strata(gis)<- pop.data
nameStrata(gis) <- colnames(pop.data)
amova <- poppr.amova(gis, ~state/county)
amova
amova.test <- randtest(amova,nrepet=999)
amova.test

amova2 <- poppr.amova(gis, ~quadrant/county)
amova2

png(file=here("./figures/Supp_AMOVA.test_state_county.png"), width =600, height =700)
plot(amova.test)
dev.off()

# Pairwise Fst per population
pop(gi) <- pop.data$county # ensure pop is county
hgi<- genind2hierfstat(gi)

pairnei <- pairwise.neifst(hgi) # Use this data
write.csv(pairnei,file=here("./tables/pairnei.csv"))
pairwc <- pairwise.WCfst(hgi) # This data is for comparison purposes
write.csv(pairwc,file=here("./tables/pairwc.csv"))

# Plotting Fst Heatmap
compfst <- read.csv(here("./tables/pairnei.csv"))
rona <- compfst[,1]
row.names(compfst) <- rona
compfst$X <- NULL
dfpfst <- as.data.frame(compfst)
fst.hm.test <- dfpfst %>% rownames_to_column(var="sample") %>%
  gather(key="Y",value="Z",-1)
ggplot(fst.hm.test, aes(sample, Y, fill=Z)) +
  geom_tile() +
  scale_fill_gradient(low="blue", high="red") +
  theme(axis.text.x = element_text(face="bold", size=10, angle=90),
        axis.text.y = element_text(face="bold", size=10, angle=0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
ggsave(here("./figures/Supp_FstNei_Heatmap.png"), width = 10, height = 10, dpi=400)

###################################################################################
################################# Manhattan Plots #################################
###################################################################################
# Load Libraries
library(hudson)

## Universal Variables
cols <- c("# Locus ID", "Chr", "BP" , "Fisher's P")
names <- c("SNP", "CHR", "POS", "pvalue")
top14chr <- c('NW_024609835.1','NW_024609846.1','NW_024609857.1','NW_024609868.1','NW_024609879.1',
              'NW_024609880.1','NW_024609881.1','NW_024609883.1','NW_024609882.1','NW_024609884.1',
              'NW_024609836.1','NW_024609839.1','NW_024609837.1','NW_024609838.1')

## Setting up RI_FL
RIFL <- data.table::fread(here("./data/all_pairs/rmInd_p20_6x200x_maf05.recode.p.fst_FL097-RI009.tsv"))
RIFL$log10F <- -log10(RIFL$`Fisher's P`)
RIFL <- as.data.frame(RIFL[!is.na(RIFL$Chr), ..cols])
names(RIFL) <- names

RIFL14 <- subset(RIFL, CHR %in% top14chr)
mylabsRIFL14 <- unique(RIFL14$CHR)
RIFL14$CHR <- as.numeric(factor(RIFL14$CHR,levels=mylabsRIFL14))

## Setting up RI_LA
RILA <- data.table::fread(here("./data/all_pairs/rmInd_p20_6x200x_maf05.recode.p.fst_LA079-RI009.tsv"))
RILA$log10F <- -log10(RILA$`Fisher's P`)
RILA <- as.data.frame(RILA[!is.na(RILA$Chr), ..cols])
names(RILA) <- names

RILA14 <- subset(RILA, CHR %in% top14chr)
mylabsRILA14 <- unique(RILA14$CHR)
RILA14$CHR <- as.numeric(factor(RILA14$CHR,levels=mylabsRILA14))

## Setting up RI_SC
RISC <- data.table::fread(here("./data/all_pairs/rmInd_p20_6x200x_maf05.recode.p.fst_RI009-SC003.tsv"))
RISC$log10F <- -log10(RISC$`Fisher's P`)
RISC <- as.data.frame(RISC[!is.na(RISC$Chr), ..cols])
names(RISC) <- names

RISC14 <- subset(RISC, CHR %in% top14chr)
mylabsRISC14 <- unique(RISC14$CHR)
RISC14$CHR <- as.numeric(factor(RISC14$CHR,levels=mylabsRISC14))

## Setting up RI_WI
RIWI <- data.table::fread(here("./data/all_pairs/rmInd_p20_6x200x_maf05.recode.p.fst_RI009-WI081.tsv"))
RIWI$log10F <- -log10(RIWI$`Fisher's P`)
RIWI <- as.data.frame(RIWI[!is.na(RIWI$Chr), ..cols])
names(RIWI) <- names

RIWI14 <- subset(RIWI, CHR %in% top14chr)
mylabsRIWI14 <- unique(RIWI14$CHR)
RIWI14$CHR <- as.numeric(factor(RIWI14$CHR,levels=mylabsRIWI14))

## Setting up SC_FL
SCFL <- data.table::fread(here("./data/all_pairs/rmInd_p20_6x200x_maf05.recode.p.fst_FL097-SC003.tsv"))
SCFL$log10F <- -log10(SCFL$`Fisher's P`)
SCFL <- as.data.frame(SCFL[!is.na(SCFL$Chr), ..cols])
names(SCFL) <- names

SCFL14 <- subset(SCFL, CHR %in% top14chr)
mylabsSCFL14 <- unique(SCFL14$CHR)
SCFL14$CHR <- as.numeric(factor(SCFL14$CHR,levels=mylabsSCFL14))

## Setting up SC_LA
SCLA <- data.table::fread(here("./data/all_pairs/rmInd_p20_6x200x_maf05.recode.p.fst_LA079-SC003.tsv"))
SCLA$log10F <- -log10(SCLA$`Fisher's P`)
SCLA <- as.data.frame(SCLA[!is.na(SCLA$Chr), ..cols])
names(SCLA) <- names

SCLA14 <- subset(SCLA, CHR %in% top14chr)
mylabsSCLA14 <- unique(SCLA14$CHR)
SCLA14$CHR <- as.numeric(factor(SCLA14$CHR,levels=mylabsSCLA14))


## Figures
hudson::gmirror(top = RIFL14,
                bottom = RISC14,
                tline = 0.05/nrow(RIFL14),
                bline = 0.05/nrow(RISC14),
                highlight_p = c(0.05/nrow(RIFL14),
                                0.05/nrow(RISC14)),
                highlighter = "green",
                toptitle = "Washington County, RI v Osceola County, FL",
                bottomtitle = "Washington County, RI v Aiken County, SC",
                background = "white",
                file=here("./figures/Man_RI009_FL097_SC003"))

hudson::gmirror(top = RILA14,
                bottom = RIWI14,
                tline = 0.05/nrow(RILA14),
                bline = 0.05/nrow(RIWI14),
                highlight_p = c(0.05/nrow(RILA14),
                                0.05/nrow(RIWI14)),
                highlighter = "green",
                toptitle = "Washington County, RI v Rapides County, LA",
                bottomtitle = "Washington County, RI v Monroe County, WI",
                background = "white",
                file=here("./figures/Man_RI009_LA079_WI081"))

hudson::gmirror(top = SCLA14,
                bottom = SCFL14,
                tline = 0.05/nrow(SCLA14),
                bline = 0.05/nrow(SCFL14),
                highlight_p = c(0.05/nrow(SCLA14),
                                0.05/nrow(SCFL14)),
                highlighter = "green",
                toptitle = "Aiken County, SC v Rapides County, LA",
                bottomtitle = "Aiken County, SC v Osceola County, FL",
                background = "white",
                file=here("./figures/Man_SC003_LA079_FL097"))

###################################################################################
################################# Population DAPC #################################
###################################################################################
#Load Libraries 
library(adegenet)

# Creating Input Variables
gl <- vcfR2genlight(vcf)
ploidy(gl) <- 2
pop(gl) <- pop.data$county # Assign Populations

# Running DAPC
dapc.pop <- adegenet::dapc(gl, n.da=100, n.pca=50)
temp <- adegenet::optim.a.score(dapc.pop)
best<-as.integer(temp$best)
best

dapc.pop <- adegenet::dapc(gl, n.da=100, n.pca=best)
dapc.pop

# Creating Color variable
temp <- pop.data[!duplicated(pop.data$county),]
mycols_pops_indv <- as.character(temp$popcolor2)
temp <- temp[with(temp, order(temp$county)),]
mycols_pops_indv <- as.character(temp$popcolor2)
# compPal <- setNames(mycols_pops_indv, temp$popmaplabel)

# Plotting DAPC
png(file=here("./figures/DAPC_populations.png"), width = 800, height = 600)
dapc.plot <- scatter(dapc.pop, cex = 2, clab= 0, col=mycols_pops_indv,
                     posi.pca = "topleft", scree.pca = TRUE,
                     posi.da = "bottomleft",
                     legend = FALSE,
                     xlim=c(-20,0), ylim=c(-12,30),
                     inset.solid = 0.7)
legend("topright", inset=c(-2.8,0.2), title="Color Family of\nAssigned Quadrant",
       title.adj = 0.15,
       legend = c("  Northeast", "  Northwest", "  Southeast", "  Southwest"),
       col = c("darkblue", "grey", "red", "yellow"), 
       fill=c("darkblue", "grey", "red", "yellow"),  
       border = "black",
       bg ="white",
       pch = c(19,19,19,19),
       cex = 1.4, 
       bty="n", 
       pt.cex = 2,
       xpd=TRUE)
text(600,100, "Osceola County, FL", cex=1.2, col = "darkred", xpd=TRUE)
text(730,-150, "Jefferson County, FL", cex=1.2, col = "darkred", xpd=TRUE)
text(390,-160, "Aiken County, SC", cex=1.2, col = "darkred", xpd=TRUE)
text(560,-355, "Louisiana Counties", cex=1.2, col = "yellow3", xpd=TRUE)
text(505,-250, "Payne County, OK", cex=1.2, col = "darkgoldenrod1", xpd=TRUE)
text(800,-320, "Pettis County, MO", cex=1.2, col = "gray27", xpd=TRUE)
text(1115,-270, "Claiborne County, TN\nStokes County, NC\nWatauga County, NC", 
     cex=1.2, col = "darkred", xpd=TRUE)
dev.off()

####################################################################################
############################## STRUCTURE-like Anlysis ##############################
####################################################################################
# https://github.com/Tom-Jenkins/admixture_pie_chart_map_tutorial
# Load Libraries
# Load packages
library(LEA)
library(miscTools)
library(stringr)
library(reshape2)
library(dplyr)
library(grid)
library(ggpubr)

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
source(here("./data/TJ_genind2structure_function.R"))
data=gi.lea
locN <- 1:nLoc(gi.lea) ## renaming loci, the plus in some of the names causes an error
locNames(data) <- locN
genind2structure(data, file=here("genotypes"), pops = FALSE, markers=FALSE)

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
png(file=here("./figures/k_structure.png"),width=600, height=500)
plot(snmf1, col = "blue", cex = 1.5, pch = 19)
dev.off()

# Extract the cross-entropy of all runs where K = 2
ce = cross.entropy(snmf1, K = 5)
ce

# Find the run with the lowest cross-entropy
lowest.ce = which.min(ce)
lowest.ce

# Extract Q-matrix for the best run
qmatrix = as.data.frame(Q(snmf1, K = 5, run = lowest.ce))
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
pop(gi.lea) <- pop.data$county
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
facet.labs = as.character(unique(pop.data$county))
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
  facet_wrap(~Site, scales = "free", ncol = 5)+
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
ggsave(here("./figures/1.admixture_barplot.png"), width=10, height=10, dpi=300)


# ----------------- #
# Prepare pie charts
# ----------------- #

# Calculate mean admixture proportions for each site
head(qmatrix)
clusters = grep("Cluster", names(qmatrix)) # indexes of cluster columns
avg_admix = aggregate(qmatrix[, clusters], list(qmatrix$Site), mean)

# Order alphabetically by site
avg_admix = avg_admix[order(as.character(avg_admix$Group.1)), ]
avg_admix

# Convert dataframe from wide to long format
avg_admix = melt(avg_admix, id.vars = "Group.1")
head(avg_admix)

# Define a function to plot pie charts using ggplot for each site
pie_charts = function(admix_df, site, cols){
  # admix_df = dataframe in long format of admixture proportions per site 
  # site = string 
  # cols = vector of colours of length(clusters)
  ggplot(data = subset(admix_df, Group.1 == site),
         aes(x = "", y = value, fill = variable))+
    geom_bar(width = 1, stat = "identity", colour = "black", show.legend = FALSE)+
    coord_polar(theta = "y")+
    scale_fill_manual(values = cols)+
    theme_void()
}

# Test function on one site
pie_charts(avg_admix, site = "Washington County, RI", cols = cols)

# Apply function to all sites using for loop
allpops<-as.character(unique(avg_admix$Group.1))
pies = list()
for (i in allpops){
  pies[[i]] = pie_charts(admix_df = avg_admix, site = i, cols = cols) 
}


# ----------------- #
# Prepare basemap
# ----------------- #

# Import csv file containing coordinates
locations = read.csv("./data/353samples.csv")

# Order alphabetically by site
locations = locations[order(locations$location), ] 

# Check order matches coords order
as.character(avg_admix$Group.1) == as.character(locations$location)

# # Map of Samples - all loaded in Sample Map Sections
# load("./DataForJulia_v02.RData")
# usa <- US.Map.Extent.LonLat
# states <- US.State.Cropped
# counties <- US.County.Cropped

# Shape File of Endemic Levels - loaded in Sample Map Section
# my_spdf <- readOGR(dsn="./SHAPES", layer="LymeIncidence_Kernel95pc")

base_map <- ggplot(data = states, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color = "black", fill = "white") + 
  geom_polygon(data = my_spdf, aes(x = long, y = lat, group = group), colour = "papayawhip", fill = NA) +
  geom_polygon(color="black", fill=NA) + #puts state borders back on top
  xlab("Longitude")+
  ylab("Latitude")+
  theme(
    axis.text = element_text(colour="black", size=12),
    axis.title = element_text(colour="black", size=14),
    panel.background = element_rect(fill="lightsteelblue2"),
    panel.border = element_rect(fill = NA, colour = "black", size = 0.5),
    panel.grid.minor = element_line(colour="grey90", size=0.5),
    panel.grid.major = element_line(colour="grey90", size=0.5),
    legend.text = element_text(size=12),
    legend.title = element_blank(),
    legend.key.size = unit(0.7, "cm"),
    legend.position = "top") +
  geom_hline(yintercept = PopGroup.hLine.Lat, size =.75, linetype = "dashed") +
  geom_vline(xintercept = PopGroup.vLine.Lon, size =.75, linetype = "dashed")
base_map

# ----------------- #
# Add pie charts to basemap
# ----------------- #

# Extract coordinates for each site
coord.list = list()
for (i in allpops){
  coord.list[[i]] = c(subset(locations, location == i)$longitude, subset(locations, location == i)$latitude)
}
coord.list

# Define pie chart sizes
radius = 1

# Convert ggplot pie charts to annotation_custom layers
pies.ac = list()
for (i in 1:length(allpops)){
  pies.ac[[i]] = annotation_custom(grob = ggplotGrob(pies[[i]]),
                                   xmin = coord.list[[i]][[1]] - radius,
                                   xmax = coord.list[[i]][[1]] + radius,
                                   ymin = coord.list[[i]][[2]] - radius,
                                   ymax = coord.list[[i]][[2]] + radius)
}

# Add layers to basemap
pie.map = base_map + pies.ac
pie.map
ggsave(here("./figures/2.pie_charts_map.png"), width = 8, height = 10, dpi = 300)

# Combine ggplots
ggarrange(admix.bar + theme(legend.position = "right") + labs(title = "Individual admixture proportions", tag = "A"),
          pie.map + labs(title = "Mean admixture proportions for each site", tag = "B"))
ggsave(here("./figures/3.Admixture_bar_map.png"), width = 25, height = 10, dpi = 300)

# Base map 2
base_map2 <- ggplot(data = states, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color = "black", fill = "white") + 
  geom_polygon(data = my_spdf, aes(x = long, y = lat, group = group), colour = "papayawhip", fill = NA) +
  geom_polygon(color="black", fill=NA) + #puts state borders back on top
  xlab("Longitude")+
  ylab("Latitude")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank())
base_map2

pie.map2 = base_map2 + pies.ac
pie.map2
ggsave(here("./figures/2.pie_charts_map_simple.png"), width = 8, height = 10, dpi = 300)

# Combine ggplots
ggarrange(admix.bar + theme(legend.position = "right") + labs(title = "Individual admixture proportions", tag = "A"),
          pie.map2 + labs(title = "Mean admixture proportions for each site", tag = "B"))
ggsave(here("./figures/3.Admixture_bar_map_simple.png"), width = 25, height = 10, dpi = 300)

###################################################################################
################################## DAPC Clusters ##################################
###################################################################################
maq <- colnames(qmatrix)[max.col(qmatrix[,1:5], ties.method = "first")] 

gl.cluster <- vcfR2genlight(vcf)
ploidy(gl.cluster) <- 2
pop(gl.cluster) <- maq

dapc.s <- adegenet::dapc(gl.cluster, n.da=100, n.pca=50)
temp <- adegenet::optim.a.score(dapc.s)
best<-as.integer(temp$best)
best

dapc.s <- adegenet::dapc(gl.cluster, n.da=100, n.pca=best)
dapc.s

# Plot DAPC Clusters
cols2=c("#7F7F7F", "#FF7F00", "#ffd500", "#FF0000", "#0000FF")
png(file=here("./figures/DAPC_allClusters.png"), width = 800, height = 600)
dapc.plot <- scatter(dapc.s, cex = 2, clab= 0, col=cols2,
                     posi.pca = "topleft", scree.pca = TRUE,
                     posi.da = "bottomleft",
                     legend = TRUE,
                     xlim=c(-15,0), ylim=c(-10,15),
                     inset.solid = 0.7)
dev.off()

png(file=here("./figures/dapc.s.contrib.png"), width =1000, height=900)
contrib.s <- loadingplot(dapc.s$var.contr, axis =2, lab.jitter=1, thres=0.003)
  #3700,3986,4996,4997,5171,5178,7225
dev.off()

vcffix <- vcf@fix
vcffix.c <- vcffix[c(3700,3986,4996,4997,5171,5178,7225),]
write.csv(vcffix.c, here("./tables/SNPloc_cluster.csv"))

# Plotting SNPs
gi.cluster <- gi.lea
pop(gi.cluster) <- maq
temp <- adegenet::seploc(gi.cluster)
snp3700 <- tab(temp[[3700]])
snp3986 <- tab(temp[[3986]])
snp4996 <- tab(temp[[4996]])
snp4997 <- tab(temp[[4997]])
snp5171 <- tab(temp[[5171]])
snp5178 <- tab(temp[[5178]])
snp7225 <- tab(temp[[7225]])

freq3700 <- apply(snp3700, 2, function(e) tapply(e, pop(gi.cluster), mean, na.rm = TRUE))
freq3986 <- apply(snp3986, 2, function(e) tapply(e, pop(gi.cluster), mean, na.rm = TRUE))
freq4996 <- apply(snp4996, 2, function(e) tapply(e, pop(gi.cluster), mean, na.rm = TRUE))
freq4997 <- apply(snp4997, 2, function(e) tapply(e, pop(gi.cluster), mean, na.rm = TRUE))
freq5171 <- apply(snp5171, 2, function(e) tapply(e, pop(gi.cluster), mean, na.rm = TRUE))
freq5178 <- apply(snp5178, 2, function(e) tapply(e, pop(gi.cluster), mean, na.rm = TRUE))
freq7225 <- apply(snp7225, 2, function(e) tapply(e, pop(gi.cluster), mean, na.rm = TRUE))

png(file=here("./figures/SNPs.allClusters.png"), width =800, height=700)
par(mfrow = c(2, 4), mar = c(5, 4, 4, 0) + 0.1, las = 3)
matplot(freq3700,  pch = c("A", "G"), type = "b",
        xlab = "Cluster", ylab = "allele frequency", main = "SNP #3700",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:5, lab = c(4,2,5,3,1))
matplot(freq3986,  pch = c("C", "T"), type = "b",
        xlab = "Cluster", ylab = "allele frequency", main = "SNP #3986",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:5, lab = c(4,2,5,3,1))
matplot(freq4996,  pch = c("G", "A"), type = "b",
        xlab = "Cluster", ylab = "allele frequency", main = "SNP #4996",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:5, lab = c(4,2,5,3,1))
matplot(freq4997,  pch = c("A", "G"), type = "b",
        xlab = "Cluster", ylab = "allele frequency", main = "SNP #4997",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:5, lab = c(4,2,5,3,1))
matplot(freq5171,  pch = c("T", "G"), type = "b",
        xlab = "Cluster", ylab = "allele frequency", main = "SNP #5171",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:5, lab = c(4,2,5,3,1))
matplot(freq5178,  pch = c("G", "A"), type = "b",
        xlab = "Cluster", ylab = "allele frequency", main = "SNP #5178",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:5, lab = c(4,2,5,3,1))
matplot(freq7225,  pch = c("C", "T"), type = "b",
        xlab = "Cluster", ylab = "allele frequency", main = "SNP #7225",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:5, lab = c(4,2,5,3,1))
dev.off()
par(mfrow = c(1, 1))
dev.off()

####################################################################################
############################## DAPC Clusters wo FL097 ##############################
####################################################################################
# Removing FL097
FL097 <- subset(pop.data, pop.data$popmaplabel == "FL097")
FL097 <- as.character(FL097$sample)
gl.clust.rm <- gl.cluster[!(indNames(gl.cluster) %in% FL097)]

# Running DAPC
dapc.rm <- adegenet::dapc(gl.clust.rm, n.da=100, n.pca=50)
temp <- adegenet::optim.a.score(dapc.rm)
best<-as.integer(temp$best)
best

dapc.rm <- adegenet::dapc(gl.clust.rm, n.da=100, n.pca=best)
dapc.rm

# Plot DAPC Clusters
cols.rm <- c("#7F7F7F","#FF7F00","#FF0000","#0000FF")
png(file=here("./figures/DAPC_rmFL097.png"), width = 800, height = 600)
dapc.plot.rm <- scatter(dapc.rm, cex = 2, clab= 0, col=cols.rm,
                     posi.pca = "topleft", scree.pca = TRUE,
                     posi.da = "bottomleft",
                     legend = TRUE, cleg = 1.5,
                     xlim=c(10,-2), ylim=c(-8,10),
                     inset.solid = 0.7)
dev.off()

png(file=here("./figures/var.contrib.woFL097.png"), width=1000, height =900)
contrib.rm <- loadingplot(dapc.rm$var.contr, axis =2, lab.jitter=1, thres=0.0059)
dev.off()
contrib.rm$var.names
#"3928" "3929" "4058" "4059" "4060" "4061" "4063" "4064" "4065"

# Defining most contributing SNPs
gi.clust.rm <- gi.cluster[!(indNames(gi.cluster) %in% FL097)]
temp <- adegenet::seploc(gi.clust.rm)
snp3928 <- tab(temp[[3928]])
snp3929 <- tab(temp[[3929]])
snp3945 <- tab(temp[[3945]])
snp4058 <- tab(temp[[4058]])
snp4059 <- tab(temp[[4059]])
snp4060 <- tab(temp[[4060]])
snp4061 <- tab(temp[[4061]])
snp4063 <- tab(temp[[4063]])
snp4064 <- tab(temp[[4064]])
snp4065 <- tab(temp[[4065]])

freq3928 <- apply(snp3928, 2, function(e) tapply(e, pop(gi.clust.rm), mean, na.rm = TRUE))
freq3929 <- apply(snp3929, 2, function(e) tapply(e, pop(gi.clust.rm), mean, na.rm = TRUE))
freq3945 <- apply(snp3945, 2, function(e) tapply(e, pop(gi.clust.rm), mean, na.rm = TRUE))
freq4058 <- apply(snp4058, 2, function(e) tapply(e, pop(gi.clust.rm), mean, na.rm = TRUE))
freq4059 <- apply(snp4059, 2, function(e) tapply(e, pop(gi.clust.rm), mean, na.rm = TRUE))
freq4060 <- apply(snp4060, 2, function(e) tapply(e, pop(gi.clust.rm), mean, na.rm = TRUE))
freq4061 <- apply(snp4061, 2, function(e) tapply(e, pop(gi.clust.rm), mean, na.rm = TRUE))
freq4063 <- apply(snp4063, 2, function(e) tapply(e, pop(gi.clust.rm), mean, na.rm = TRUE))
freq4064 <- apply(snp4064, 2, function(e) tapply(e, pop(gi.clust.rm), mean, na.rm = TRUE))
freq4065 <- apply(snp4065, 2, function(e) tapply(e, pop(gi.clust.rm), mean, na.rm = TRUE))

vcffix.rm <- vcffix[c(3928,3929,3945,4058,4059,4060,4061,4063,4064,4065),]
write.csv(vcffix.rm, here("./tables/SNPloc_woFL097_top10.csv"))

# Plotting SNPs
png(file=here("./figures/SNPs.woFL097_top10.png"), width =800, height=1200)
par(mfrow = c(5, 2), mar = c(5, 4, 4, 0) + 0.1, las = 3)
matplot(freq3928,  pch = c("A", "G"), type = "b",
        xlab = "Cluster", ylab = "allele frequency", main = "SNP #3928",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:4, lab = c(4,2,5,1))
matplot(freq3929,  pch = c("A", "G"), type = "b",
        xlab = "Cluster", ylab = "allele frequency", main = "SNP #3929",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:4, lab = c(4,2,5,1))
matplot(freq3945,  pch = c("A", "G"), type = "b",
        xlab = "Cluster", ylab = "allele frequency", main = "SNP #3945",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:4, lab = c(4,2,5,1))
matplot(freq4058,  pch = c("C", "G"), type = "b",
        xlab = "Cluster", ylab = "allele frequency", main = "SNP #4058",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:4, lab = c(4,2,5,1))
matplot(freq4059,  pch = c("T", "C"), type = "b",
        xlab = "Cluster", ylab = "allele frequency", main = "SNP #4059",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:4, lab = c(4,2,5,1))
matplot(freq4060,  pch = c("A", "G"), type = "b",
        xlab = "Cluster", ylab = "allele frequency", main = "SNP #4060",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:4, lab = c(4,2,5,1))
matplot(freq4061,  pch = c("A", "G"), type = "b",
        xlab = "Cluster", ylab = "allele frequency", main = "SNP #4061",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:4, lab = c(4,2,5,1))
matplot(freq4063,  pch = c("C", "T"), type = "b",
        xlab = "Cluster", ylab = "allele frequency", main = "SNP #4063",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:4, lab = c(4,2,5,1))
matplot(freq4064,  pch = c("A", "G"), type = "b",
        xlab = "Cluster", ylab = "allele frequency", main = "SNP #4064",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:4, lab = c(4,2,5,1))
matplot(freq4065,  pch = c("T", "C"), type = "b",
        xlab = "Cluster", ylab = "allele frequency", main = "SNP #4065",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:4, lab = c(4,2,5,1))
dev.off()

#save.image(file = "E:/3RAD/Iscap3RAD_CompiledCode_workspace.RData")
  ## Saved to External Harddrive 07Feb2022
########## Pull above 0.002 contrubting
contrib.rm.002 <- loadingplot(dapc.rm$var.contr, axis =2, lab.jitter=1, thres=0.002)
var002 <- contrib.rm.002$var.names

#gi.pop.rm <- gi[!(indNames(gi) %in% FL097)]

loc002 <- data.frame("locName" = locNames(gi.clust.rm), "Num"=seq.int(length(locNames(gi.clust.rm))))
loc002 <- subset(loc002.temp, Num %in% var002)

gi.rm.002 <- gi[loc = as.character(loc002$locName)]

bs.gi.002<-basic.stats(gi.rm.002,diploid = TRUE)

## Heterozyogosity Obs per pop
Ho.002 <- bs.gi.002[["Ho"]]
Ho_avg.002 <- apply(Ho.002,2, mean, na.rm=TRUE)  # applies function 'mean' to 2nd dimension (columns)

## Heterozygosity Expected (mean gene diversity) per pop
Hs.002 <- bs.gi.002[["Hs"]]
Hs_avg.002 <- apply(Hs.002,2, mean, na.rm=TRUE)

## Inbreeding Coefficent per pop
Fis.002 <- bs.gi.002[["Fis"]]
Fis_avg.002 <- apply(Fis.002,2, mean, na.rm=TRUE)  

## Allelic richness per pop (hierfstat)
ar.002 <- allelic.richness(gi.rm.002, diploid=TRUE)
ar.2.002 <- ar.002[["Ar"]]
ar_avg.002 <- apply(ar.2.002,2, mean, na.rm=TRUE)  

## Number of Alleles (hierfstat)
na.002<- nb.alleles(gi.rm.002, diploid = TRUE)
na_sum.002 <- apply(na.002,2, sum, na.rm=TRUE)

## Getting Total number of samples
totaln.002 <- as.numeric(table(gi.rm.002@pop))

## Getting Pop Names
pops.002 <- colnames( bs.gi.002[["pop.freq"]][[1]])

sum_stats_pops.002 <- data.frame("Populations"=pops.002, "N"=totaln.002, "Nalleles"=na_sum.002, 
                             "Ar" = as.numeric(ar_avg.002), "Ho"=as.numeric(Ho_avg.002),
                             "Hs"=as.numeric(Hs_avg.002), "Fis"=as.numeric(Fis_avg.002))
write.csv(sum_stats_pops.002,here("./tables/sum_stats_pops.002.csv"))

# Testing significance of stat differences between North Vs South
labelled_ssp.002 <- sum_stats_pops.002 %>% dplyr::mutate(NvS = dplyr::if_else(
  Populations %in% c("Aiken County, SC","Bossier County, LA","Claiborne County, TN","Iberville County, LA",
                     "Jefferson County, FL","Lauderdale County, AL","Osceola County, FL","Payne County (wild), OK",
                     "Rapides County, LA","Stokes County, NC","Walker County, GA","Watauga Count, NC"),
  "South", "North"))

temp <- subset(pop.data, endemic95=="In")
labelled_ssp.002 <- labelled_ssp.002 %>% dplyr::mutate(in95 = dplyr::if_else(
  Populations %in% levels(factor(temp$county)),
  "In", "Out"))

#library(ggpubr)
Ar_bp.002 <- ggboxplot(
  labelled_ssp.002, x = "in95", y = "Ar",
  ylab = "Allelic Richness", xlab = "95% Endemic", add = "jitter"
)
png(file=here("./figures/Ar_bp.002_wFL097.png"), width=500, height =500)
Ar_bp.002
dev.off()

Ar_tt.002 <- t.test(Ar ~ in95, data=labelled_ssp.002) %>% broom::tidy()

########## Pull below 0.002 contrubting
#contrib.rm.002 <- loadingplot(dapc.rm$var.contr, axis =2, lab.jitter=1, thres=0.002)
#var002 <- contrib.rm.002$var.names

#gi.pop.rm <- gi[!(indNames(gi) %in% FL097)]

loc002 <- data.frame("locName" = locNames(gi.clust.rm), "Num"=seq.int(length(locNames(gi.clust.rm))))
loc998 <- subset(loc002, !(Num %in% var002))

gi.rm.998 <- gi[loc = as.character(loc998$locName)]

bs.gi.998<-basic.stats(gi.rm.998,diploid = TRUE)

## Heterozyogosity Obs per pop
Ho.998 <- bs.gi.998[["Ho"]]
Ho_avg.998 <- apply(Ho.998,2, mean, na.rm=TRUE)  # applies function 'mean' to 2nd dimension (columns)

## Heterozygosity Expected (mean gene diversity) per pop
Hs.998 <- bs.gi.998[["Hs"]]
Hs_avg.998 <- apply(Hs.998,2, mean, na.rm=TRUE)

## Inbreeding Coefficent per pop
Fis.998 <- bs.gi.998[["Fis"]]
Fis_avg.998 <- apply(Fis.998,2, mean, na.rm=TRUE)  

## Allelic richness per pop (hierfstat)
ar.998 <- allelic.richness(gi.rm.998, diploid=TRUE)
ar.2.998 <- ar.998[["Ar"]]
ar_avg.998 <- apply(ar.2.998,2, mean, na.rm=TRUE)  

## Number of Alleles (hierfstat)
na.998<- nb.alleles(gi.rm.998, diploid = TRUE)
na_sum.998 <- apply(na.998,2, sum, na.rm=TRUE)

## Getting Total number of samples
totaln.998 <- as.numeric(table(gi.rm.998@pop))

## Getting Pop Names
pops.998 <- colnames( bs.gi.998[["pop.freq"]][[1]])

sum_stats_pops.998 <- data.frame("Populations"=pops.998, "N"=totaln.998, "Nalleles"=na_sum.998, 
                                 "Ar" = as.numeric(ar_avg.998), "Ho"=as.numeric(Ho_avg.998),
                                 "Hs"=as.numeric(Hs_avg.998), "Fis"=as.numeric(Fis_avg.998))
write.csv(sum_stats_pops.998,here("./tables/sum_stats_pops.998.csv"))

# Testing significance of stat differences between North Vs South
labelled_ssp.998 <- sum_stats_pops.998 %>% dplyr::mutate(NvS = dplyr::if_else(
  Populations %in% c("Aiken County, SC","Bossier County, LA","Claiborne County, TN","Iberville County, LA",
                     "Jefferson County, FL","Lauderdale County, AL","Osceola County, FL","Payne County (wild), OK",
                     "Rapides County, LA","Stokes County, NC","Walker County, GA","Watauga Count, NC"),
  "South", "North"))

temp <- subset(pop.data, endemic95=="In")
labelled_ssp.998 <- labelled_ssp.998 %>% dplyr::mutate(in95 = dplyr::if_else(
  Populations %in% levels(factor(temp$county)),
  "In", "Out"))

#library(ggpubr)
Ar_bp.998 <- ggboxplot(
  labelled_ssp.998, x = "in95", y = "Ar",
  ylab = "Allelic Richness", xlab = "95% Endemic", add = "jitter"
)
png(file=here("./figures/Ar_bp.998_wFL097.png"), width=500, height =500)
Ar_bp.998
dev.off()

Ar_tt.998 <- t.test(Ar ~ in95, data=labelled_ssp.998) %>% broom::tidy()

###################################################################################
############################## Isolation by Distance ##############################
###################################################################################
# https://adegenet.r-forge.r-project.org/files/tutorial.pdf
# Followed isolation by distance using Euclidean distances as it handles missing values
#library(MASS)

# all populations
gp <- genind2genpop(gi)
Dgen <- dist.genpop(gp, method = 5)
gi.loc <- gi
gi.loc@other$xy <- locations[,c(6,5)] #longitude, latitude
Dgeo <- dist(gi.loc$other$xy)
ibd <- mantel.randtest(Dgen, Dgeo)
ibd
plot(ibd)

#library(MASS)
dens <- kde2d(Dgeo,Dgen, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(Dgeo, Dgen, pch=20,cex=1)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(Dgen~Dgeo))
title("Isolation by distance plot")

#dgendf <- as.data.frame(as.matrix(Dgen))
#write.csv(dgendf, file="./Dgen_pop.csv")

#library(reshape2)
dgenlong <- melt(as.matrix(Dgen))
dgeolong <- melt(as.matrix(Dgeo))

dgeolong$Var1 <- with(dgeolong, factor(Var1, 
                      levels= c("1","2","3","4","5","6","7","8","9","10","11","12","13","14",
                                      "15","16","17","18","19","20","21","22","23","24","25","26",
                                      "27","28","29","30","31","32","33"), 
                      labels=c("Lauderdale County, AL","Kent County, DE","Jefferson County, FL",
                               "Osceola County, FL","Walker County, GA","Poweshiek County, IA",
                               "Van Buren County, IA","Logan County, IL","Winnebago County, IL",
                               "Franklin County, KY","Owen County, KY","Bossier County, LA",
                               "Iberville County, LA","Rapides County, LA","Allegany County, MD",
                               "Ingham County, MI","Pettis County, MO","Stokes County, NC","Watauga Count, NC",
                               "Carroll County, NH","Oswego County, NY","Schenectady County, NY",
                               "Knox County, OH","Scioto County, OH","Payne County (wild), OK",
                               "Centre County, PA","Westmoreland County, PA","Washington County, RI",
                               "Aiken County, SC","Claiborne County, TN","Windsor County, VT",
                               "Marinette County, WI","Monroe County, WI")))
dgeolong$Var2 <- with(dgeolong, factor(Var2, 
                                       levels= c("1","2","3","4","5","6","7","8","9","10","11","12","13","14",
                                                 "15","16","17","18","19","20","21","22","23","24","25","26",
                                                 "27","28","29","30","31","32","33"), 
                                       labels=c("Lauderdale County, AL","Kent County, DE","Jefferson County, FL",
                                                "Osceola County, FL","Walker County, GA","Poweshiek County, IA",
                                                "Van Buren County, IA","Logan County, IL","Winnebago County, IL",
                                                "Franklin County, KY","Owen County, KY","Bossier County, LA",
                                                "Iberville County, LA","Rapides County, LA","Allegany County, MD",
                                                "Ingham County, MI","Pettis County, MO","Stokes County, NC","Watauga Count, NC",
                                                "Carroll County, NH","Oswego County, NY","Schenectady County, NY",
                                                "Knox County, OH","Scioto County, OH","Payne County (wild), OK",
                                                "Centre County, PA","Westmoreland County, PA","Washington County, RI",
                                                "Aiken County, SC","Claiborne County, TN","Windsor County, VT",
                                                "Marinette County, WI","Monroe County, WI")))

names(dgeolong)[names(dgeolong) == "value"] <- "Dgeo"
names(dgenlong)[names(dgenlong) == "value"] <- "Dgen"

mergedDpops <- dgeolong %>% dplyr::right_join(dgenlong, by=c("Var1","Var2"))
#write.csv(mergedDpops, file="./mergedDpops.csv")

southpops_quads <- c("Lauderdale County, AL","Jefferson County, FL",
                     "Osceola County, FL","Walker County, GA","Bossier County, LA",
                     "Iberville County, LA","Rapides County, LA","Stokes County, NC",
                     "Watauga Count, NC", "Payne County (wild), OK",
                     "Aiken County, SC","Claiborne County, TN")
mergedDpops$comps_quads <- paste(ifelse(mergedDpops$Var1 %in% southpops_quads, "S","N"),
                                 ifelse(mergedDpops$Var2 %in% southpops_quads, "S","N"), sep="")

southpops_clusts <- c("Lauderdale County, AL","Jefferson County, FL",
                      "Osceola County, FL","Walker County, GA","Bossier County, LA",
                      "Iberville County, LA","Rapides County, LA","Pettis County, MO",
                      "Payne County (wild), OK","Aiken County, SC","Claiborne County, TN")
mergedDpops$comps_clusts <- paste(ifelse(mergedDpops$Var1 %in% southpops_clusts, "S","N"),
                                 ifelse(mergedDpops$Var2 %in% southpops_clusts, "S","N"), sep="")

# dens <- kde2d(mergedDpops$Dgeo, mergedDpops$Dgen, n=300)
# myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
# plot(mergedDpops$Dgeo, mergedDpops$Dgen, pch=20,cex=1)
# image(dens, col=transp(myPal(300),.7), add=TRUE)
# #abline(lm(mergedDpops$Dgeo~mergedDpops$Dgen))
# title("Isolation by distance plot")

mergedDpops[mergedDpops == "SN"] <- "NS"

mDpops_drop <- subset(mergedDpops, Dgen != 0 & Dgeo != 0)

df2 <- dplyr::select(mDpops_drop, -comps_quads)
ggplot(mDpops_drop, aes(x=Dgeo, y=Dgen))+
  geom_point(data=df2, color="grey70")+
  geom_point(aes(color=comps_quads))+
  facet_wrap(~comps_quads) +
  theme_minimal()

df2 <- dplyr::select(mDpops_drop, -comps_clusts)
ggplot(mDpops_drop, aes(x=Dgeo, y=Dgen))+
  geom_point(data=df2, color="grey70")+
  geom_point(aes(color=comps_clusts))+
  facet_wrap(~comps_clusts) +
  theme_minimal()



# all populations north by quadrant
all.pop.loc <- pop.data
all.pop.loc <- dplyr::left_join(pop.data, locations, by=c("county" = "location"))

north <- subset(all.pop.loc, all.pop.loc$quadrant %in% c("Northwest", "Northeast"))
north.individuals <- as.character(north$sample)
gi.north <- gi[indNames(gi) %in% north.individuals]

north.loc<- unique(north[,c(3,12,13)])

gp.north <- genind2genpop(gi.north)
Dgen.north <- dist.genpop(gp.north, method = 4)
gi.loc.n <- gi.north
gi.loc.n@other$xy <- north.loc[,c(3,2)]
Dgeo.north <- dist(gi.loc.n$other$xy)
ibd.north <- mantel.randtest(Dgen.north, Dgeo.north)
ibd.north
plot(ibd.north)

dens <- kde2d(Dgeo.north,Dgen.north, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(Dgeo.north, Dgen.north, pch=20,cex=1)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(Dgen.north~Dgeo.north))
title("Isolation by distance plot")

# all populations south by quadrant
south <- subset(all.pop.loc, all.pop.loc$quadrant %in% c("Southwest", "Southeast"))
south.individuals <- as.character(south$sample)
gi.south <- gi[(indNames(gi) %in% south.individuals)]

south.loc<- unique(south[,c(3,12,13)])

gp.south <- genind2genpop(gi.south)
Dgen.south <- dist.genpop(gp.south, method = 4)
gi.loc.s <- gi.south
gi.loc.s@other$xy <- south.loc[,c(3,2)]
Dgeo.south <- dist(gi.loc.s$other$xy)
ibd.south <- mantel.randtest(Dgen.south, Dgeo.south)
ibd.south
plot(ibd.south)

dens <- kde2d(Dgeo.south,Dgen.south, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(Dgeo.south, Dgen.south, pch=20,cex=1)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(Dgen.south~Dgeo.south))
title("Isolation by distance plot")

# all populations north by Cluster 1 & Cluster 5
clust15.states <- c("DE","IA","IL","KY","MD","MI","NC","NH","NY","OH","PA","RI","VT","WI")
all.pop.loc <- dplyr::left_join(pop.data, locations, by=c("county" = "location"))

clust15 <- subset(all.pop.loc, all.pop.loc$state.x %in% clust15.states)
gi.clust15 <- gi[indNames(gi) %in% as.character(clust15$sample)]

clust15.loc<- unique(clust15[,c(3,12,13)])

gp.clust15 <- genind2genpop(gi.clust15)
Dgen.clust15 <- dist.genpop(gp.clust15, method = 4)
gi.loc.clust15 <- gi.clust15
gi.loc.clust15@other$xy <- clust15.loc[,c(3,2)]
Dgeo.clust15 <- dist(gi.loc.clust15$other$xy)
ibd.clust15 <- mantel.randtest(Dgen.clust15, Dgeo.clust15)
ibd.clust15
plot(ibd.clust15)

dens <- kde2d(Dgeo.clust15,Dgen.clust15, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(Dgeo.clust15, Dgen.clust15, pch=20,cex=1)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(Dgen.clust15~Dgeo.clust15))
title("Isolation by distance plot")

# all populations south by Cluster 2 Cluster 3 & Cluster 4
clust234.states <- c("AL","FL","GA","LA","MO","OK","SC","TN")
#all.pop.loc <- dplyr::left_join(pop.data, locations, by=c("county" = "location"))

clust234 <- subset(all.pop.loc, all.pop.loc$state.x %in% clust234.states)
gi.clust234 <- gi[indNames(gi) %in% as.character(clust234$sample)]

clust234.loc<- unique(clust234[,c(3,12,13)])

gp.clust234 <- genind2genpop(gi.clust234)
Dgen.clust234 <- dist.genpop(gp.clust234, method = 4)
gi.loc.clust234 <- gi.clust234
gi.loc.clust234@other$xy <- clust234.loc[,c(3,2)]
Dgeo.clust234 <- dist(gi.loc.clust234$other$xy)
ibd.clust234 <- mantel.randtest(Dgen.clust234, Dgeo.clust234)
ibd.clust234
plot(ibd.clust234)

dens <- kde2d(Dgeo.clust234,Dgen.clust234, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(Dgeo.clust234, Dgen.clust234, pch=20,cex=1)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(Dgen.clust234~Dgeo.clust234))
title("Isolation by distance plot")

# gi.clust all
gp.clust <- genind2genpop(gi.cluster)
Dgen.clust <- dist.genpop(gp.clust, method = 4)

all.pop.loc <- pop.data
all.pop.loc <- dplyr::left_join(pop.data, locations, by=c("county" = "location"))
all.pop.loc$cluster <- maq

Clust1.loc <- subset(all.pop.loc, cluster == "Cluster 1")
Clust2.loc <- subset(all.pop.loc, cluster == "Cluster 2")
Clust3.loc <- subset(all.pop.loc, cluster == "Cluster 3")
Clust4.loc <- subset(all.pop.loc, cluster == "Cluster 4")
Clust5.loc <- subset(all.pop.loc, cluster == "Cluster 5")

Clust1.info <- data.frame(cluster = c("Cluster 4","Cluster 2", "Cluster 5", "Cluster 3", "Cluster 1"),
                          longitude = c(mean(Clust4.loc$longitude), mean(Clust2.loc$longitude), mean(Clust5.loc$longitude),mean(Clust3.loc$longitude),mean(Clust1.loc$longitude)), 
                          latitude = c(mean(Clust4.loc$latitude),mean(Clust2.loc$latitude),mean(Clust5.loc$latitude),mean(Clust3.loc$latitude),mean(Clust1.loc$latitude)))

gi.clust.loc <- gi.cluster
gi.clust.loc@other$xy <- Clust1.info[,2:3]
Dgeo.clust <- dist(gi.clust.loc$other$xy)
ibd.clust <- mantel.randtest(Dgen.clust, Dgeo.clust)
# ibd.clust
# plot(ibd.clust)
# 
# plot(Dgeo.clust, Dgen.clust)
# abline(lm(Dgen.clust~Dgeo.clust), col="red",lty=2)

# library(MASS)
dens <- kde2d(Dgeo.clust,Dgen.clust, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(Dgeo.clust, Dgen.clust, pch=20,cex=1.5)
image(dens, col=transp(myPal(300),.6), add=TRUE)
abline(lm(Dgen.clust~Dgeo.clust))
title("Isolation by distance plot")

# gi.clust all
gp.clust <- genind2genpop(gi.cluster)
Dgen.clust <- dist.genpop(gp.clust, method = 4)

all.pop.loc <- pop.data
all.pop.loc <- dplyr::left_join(pop.data, locations, by=c("county" = "location"))
all.pop.loc$cluster <- maq

Clust1.loc <- subset(all.pop.loc, cluster == "Cluster 1")
Clust2.loc <- subset(all.pop.loc, cluster == "Cluster 2")
Clust3.loc <- subset(all.pop.loc, cluster == "Cluster 3")
Clust4.loc <- subset(all.pop.loc, cluster == "Cluster 4")
Clust5.loc <- subset(all.pop.loc, cluster == "Cluster 5")

Clust1.info <- data.frame(cluster = c("Cluster 4","Cluster 2", "Cluster 5", "Cluster 3", "Cluster 1"),
                          longitude = c(mean(Clust4.loc$longitude), mean(Clust2.loc$longitude), mean(Clust5.loc$longitude),mean(Clust3.loc$longitude),mean(Clust1.loc$longitude)), 
                          latitude = c(mean(Clust4.loc$latitude),mean(Clust2.loc$latitude),mean(Clust5.loc$latitude),mean(Clust3.loc$latitude),mean(Clust1.loc$latitude)))

gi.clust.loc <- gi.cluster
gi.clust.loc@other$xy <- Clust1.info[,2:3]
Dgeo.clust <- dist(gi.clust.loc$other$xy)
ibd.clust <- mantel.randtest(Dgen.clust, Dgeo.clust)
 ibd.clust
# plot(ibd.clust)
# 
# plot(Dgeo.clust, Dgen.clust)
# abline(lm(Dgen.clust~Dgeo.clust), col="red",lty=2)

# library(MASS)
dens <- kde2d(Dgeo.clust,Dgen.clust, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(Dgeo.clust, Dgen.clust, pch=20,cex=1.5)
image(dens, col=transp(myPal(300),.6), add=TRUE)
abline(lm(Dgen.clust~Dgeo.clust))
title("Isolation by distance plot")

# gi.clust.rm
gp.clust.rm <- genind2genpop(gi.clust.rm)
Dgen.clust.rm <- dist.genpop(gp.clust.rm, method = 2)

Clust1.info.rm <- data.frame(cluster = c("Cluster 4","Cluster 2", "Cluster 5", "Cluster 1"),
                             longitude = c(mean(Clust4.loc$longitude), mean(Clust2.loc$longitude), mean(Clust5.loc$longitude),mean(Clust1.loc$longitude)), 
                             latitude = c(mean(Clust4.loc$latitude),mean(Clust2.loc$latitude),mean(Clust5.loc$latitude),mean(Clust1.loc$latitude)))


gi.clust.loc.rm <- gi.clust.rm
gi.clust.loc.rm@other$xy <- Clust1.info.rm[,2:3]
Dgeo.clust.rm <- dist(gi.clust.loc.rm$other$xy)
ibd.clust.rm <- mantel.randtest(Dgen.clust.rm, Dgeo.clust.rm)
ibd.clust.rm
plot(ibd.clust.rm)
# 
# plot(Dgeo.clust, Dgen.clust)
# abline(lm(Dgen.clust~Dgeo.clust), col="red",lty=2)

# library(MASS)
dens <- kde2d(Dgeo.clust.rm,Dgen.clust.rm, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(Dgeo.clust.rm, Dgen.clust.rm, pch=20,cex=1.5)
image(dens, col=transp(myPal(300),.6), add=TRUE)
abline(lm(Dgen.clust.rm~Dgeo.clust.rm))
title("Isolation by distance plot")
