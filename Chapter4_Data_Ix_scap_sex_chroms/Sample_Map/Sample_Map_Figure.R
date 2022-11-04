### Sample Map for Reseq Population Genetics

#setwd("~/UGA_EHS Dropbox/Julia Frederick/Julia_Frederick_Projects/Sex_chrom_ticks/manuscript/Sample_Map")

# Load libraries
library(here)
library(scatterpie)
library(rgdal)
library(ggplot2)
library(tidyverse)
library(ggsci)

# Map of Samples
load(here("./DataForJulia_v02.RData")) # from Guha
locations <- read.csv(here("./Population_Locations.csv")) # from my files
rad_locations <- read.csv(here("./3RAD_locations.csv"))
ytest <- read.csv(here("./samples_ytest.csv"))

# Shape File of Endemic Levels
my_spdf <- readOGR(dsn=here("./SHAPES"), layer="LymeIncidence_Kernel95pc") # from Guha

# Set map variables
usa <- US.Map.Extent.LonLat
states <- US.State.Cropped
counties <- US.County.Cropped

locations$dataset <- "Resequencing Samples"
rad_locations$dataset <- "3RAD Samples"
ytest$dataset <- "PCR Samples"

mer <- bind_rows(locations[,c("location","longitude","latitude","dataset")], 
                 rad_locations[,c("location","longitude","latitude","dataset")],
                 ytest[,c("location","longitude","latitude","dataset")])
mer$dataset <- as.factor(mer$dataset)

Map_3RAD.bw <- ggplot(data = states, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color = "black", fill = "white") + 
  #geom_polygon(data = my_spdf, aes(x = long, y = lat, group = group), 
  #             colour = "#bdbdbd", fill = NA) +
  geom_polygon(color="black", fill=NA) + #puts state borders back on top
  geom_point(data=mer, aes(x=longitude, y=latitude, group=dataset, shape=dataset,fill=dataset, size=dataset),
             color="black", alpha=0.9) +
  scale_fill_manual(name = "Dataset",
                    labels = c("3RAD Samples", "PCR Samples", "Resequencing Samples"),
                    values = c("#0073C2FF","#EFC000FF","#868686FF")) +   
  scale_shape_manual(name = "Dataset",
                     labels = c("3RAD Samples", "PCR Samples","Resequencing Samples"),
                     values = c(23,22,21)) +
  scale_size_manual(name = "Dataset",
                     labels = c("3RAD Samples", "PCR Samples","Resequencing Samples"),
                     values = c(3.7,2.6,6)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.85,0.25))
Map_3RAD.bw
ggsave(here("./Fig1_sex_chrom_allSamples.pdf"), 
       width = 250, height = 250, units = "mm", dpi=600)
