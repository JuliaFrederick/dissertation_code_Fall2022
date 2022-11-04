######################## Libraries & Working Directory ######################## 
library(karyoploteR) ## Main one you need
library(vcfR)
#library(Rsamtools) ## Possibly need to work with bam files but honestly can't remember

#setwd("/Users/juliafrederick/UGA_EHS Dropbox/Julia Frederick/Julia_Frederick_Projects/Ixscap_reseq_population_genetics/popgen/R_scripts/variant_density") #my current directory

######################## Load Data ######################## 
rad_vcf <- read.vcfR("../../vcf_files_filtered/3RAD_m90_6x150x_GQ30_maf05_ma2.recode.vcf.gz")
rad_snp <- as.data.frame(rad_vcf@fix)
rad_snp$POS <- as.numeric(rad_snp$POS)
rad_range <- GenomicRanges::GRanges(seqnames=rad_snp[,1], ranges=IRanges::IRanges(start=rad_snp[,2], end=rad_snp[,2]))

gatk_vcf <- read.vcfR("../../vcf_files_filtered/gatk_m90_6x50x_GQ30_Q30_maf05_ma2.recode.vcf.gz")
gatk_snp <- as.data.frame(gatk_vcf@fix)
gatk_snp$POS <- as.numeric(gatk_snp$POS)
gatk_range <- GenomicRanges::GRanges(seqnames=gatk_snp[,1], ranges=IRanges::IRanges(start=gatk_snp[,2], end=gatk_snp[,2]))

sam_vcf <- read.vcfR("../../vcf_files_filtered/samtools_formatQ_m90_6x50x_GQ30_Q30_maf05_ma2.recode.vcf.gz")
sam_snp <- as.data.frame(sam_vcf@fix)
sam_snp$POS <- as.numeric(sam_snp$POS)
sam_range <- GenomicRanges::GRanges(seqnames=sam_snp[,1], ranges=IRanges::IRanges(start=sam_snp[,2], end=sam_snp[,2]))

overlap_vcf <- read.vcfR("../../vcf_files_filtered/0000.vcf.gz")
overlap_snp <- as.data.frame(overlap_vcf@fix)
overlap_snp$POS <- as.numeric(overlap_snp$POS)
overlap_range <- GenomicRanges::GRanges(seqnames=overlap_snp[,1], ranges=IRanges::IRanges(start=overlap_snp[,2], end=overlap_snp[,2]))


## Karyotype creation
m_karyotype <- read.table("./2021_karyotype_v2.txt", header=FALSE)
m_chr <- (m_karyotype$V3)
m_start <- (m_karyotype$V5)
m_end <- (m_karyotype$V6)
m_custom.genome <- toGRanges(data.frame(chr=m_chr, start=m_start, end=m_end))

######################## Plot Parameters ########################
## Plot shape will change based on the window it is in, so pick the window size you want and plot to that!

## Optional section just remove "plot.params = plot.params" from the fuctions below
## You can change any of these, just uncomment them - I also think there are more that I don't have listed
## General info about plot types: https://bernatgel.github.io/karyoploter_tutorial/Tutorial/PlotTypes/PlotTypes.html
## General info about plot params: https://bioconductor.org/packages/release/bioc/vignettes/karyoploteR/inst/doc/karyoploteR.html 
## plotDefaultPlotParams(plot.type=2) ## running will produce a plot to show you the different variable

plot.params <- getDefaultPlotParams(plot.type=4) # Currently using plot.type 4 but can change
#plot.params$data1outmargin <- 105
plot.params$ideogramheight <- 1
plot.params$leftmargin <- 0.19
plot.params$ideogramlateralmargin <-0.008
#plot.params$data1height <- 150
#plot.params$data1inmargin <- 1.5
plot.params$data1max <- 100
#plot.params$data1min <- 0
#plot.params$chromosome.height <- 1
#plot.params$data2height <- 0
#plot.params$data2inmargin <- 0
#plot.params$data2height <-0
#plot.params$data2max <- 0
#plot.params$data2min <- 0
#plot.params$topmargin <- 10
plot.params$bottommargin <- 45

raddensity <- kp$latest.plot$computed.values
rad_den2 <- raddensity$density
rad_den2 <- data.frame("density" = rad_den2, "window"=1:length(rad_den2))
sum(rad_den2$density !=0)

gatkdensity <- kp$latest.plot$computed.values
gatk_den2 <- gatkdensity$density
gatk_den2 <- data.frame("density" = gatk_den2, "window"=1:length(gatk_den2))
sum(gatk_den2$density !=0)

samdensity <- kp$latest.plot$computed.values
sam_den2 <- samdensity$density
sam_den2 <- data.frame("density" = sam_den2, "window"=1:length(sam_den2))
sum(sam_den2$density !=0)

overdensity <- kp$latest.plot$computed.values
over_den2 <- overdensity$density
over_den2 <- data.frame("density" = over_den2, "window"=1:length(over_den2))
sum(over_den2$density !=0)

######################## Basic Plot for Single Sample with coverage axis ######################## 
kp <- plotKaryotype(genome = m_custom.genome, plot.type=4, plot.params=plot.params, cex = 0.90, srt=30)
kp <- kpAddBaseNumbers(kp, tick.dist = 50000000, add.units = FALSE, cex = .75)
at <- autotrack(current.track = 1, total.tracks = 4, margin = 0.1)
kp <- kpPlotDensity(kp, data=rad_range, window.size = 1000, r0=(at$r0)*plot.params$data1max, r1=(at$r1)*plot.params$data1max,
              border=transparent("#7FB3EE",0.5), lwd=2)
kp <- kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, cex=0.8,
       r0=(at$r0)*plot.params$data1max, r1=(at$r1)*plot.params$data1max)
kp <- kpAddLabels(kp, labels="3RAD", cex=0.9, srt=90, pos = 3, label.margin = 0.03,
                  r0=(at$r0)*plot.params$data1max, r1=(at$r1)*plot.params$data1max)
at <- autotrack(current.track = 2, total.tracks = 4, margin = 0.1)
kp <- kpPlotDensity(kp, data=gatk_range, window.size = 1000, r0=(at$r0)*plot.params$data1max, r1=(at$r1)*plot.params$data1max,
              border=transparent("#332851",0.5), lwd=2)
kp <- kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, cex=0.8,
       r0=(at$r0)*plot.params$data1max, r1=(at$r1)*plot.params$data1max)
kp <- kpAddLabels(kp, labels="GATK", cex=0.9, srt=90, pos = 3, label.margin = 0.03,
                  r0=(at$r0)*plot.params$data1max, r1=(at$r1)*plot.params$data1max)
at <- autotrack(current.track = 3, total.tracks = 4, margin = 0.1)
kp <- kpPlotDensity(kp, data=sam_range, window.size = 1000, r0=(at$r0)*plot.params$data1max, r1=(at$r1)*plot.params$data1max,
              border=transparent("#CA3074",0.5), lwd=2)
kp <- kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, cex=0.8,
       r0=(at$r0)*plot.params$data1max, r1=(at$r1)*plot.params$data1max)
kp <- kpAddLabels(kp, labels="Samtools", cex=0.9, srt=90, pos = 3, label.margin = 0.03,
                  r0=(at$r0)*plot.params$data1max, r1=(at$r1)*plot.params$data1max)
at <- autotrack(current.track = 4, total.tracks = 4, margin = 0.1)
kp <- kpPlotDensity(kp, data=overlap_range, window.size = 1000, r0=(at$r0)*plot.params$data1max, r1=(at$r1)*plot.params$data1max,
              border=transparent("#F6C667",0.5), lwd=2)
kp <- kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, cex=0.8,
       r0=(at$r0)*plot.params$data1max, r1=(at$r1)*plot.params$data1max)
kp <- kpAddLabels(kp, labels="Overlap", cex=0.9, srt=90, pos = 3, label.margin = 0.03,
                  r0=(at$r0)*plot.params$data1max, r1=(at$r1)*plot.params$data1max)

######################## Basic Plot for Single Sample with coverage axis ######################## 
raddensity_4 <- kp$latest.plot$computed.values
rad_den2_4 <- raddensity_4$density
rad_den2_4 <- data.frame("density" = rad_den2_4, "window"=1:length(rad_den2_4))
sum(rad_den2_4$density !=0)
sum(rad_den2_4$density)

gatkdensity_4 <- kp$latest.plot$computed.values
gatk_den2_4 <- gatkdensity_4$density
gatk_den2_4 <- data.frame("density" = gatk_den2_4, "window"=1:length(gatk_den2_4))
sum(gatk_den2_4$density !=0)
sum(gatk_den2_4$density)

samdensity_4 <- kp$latest.plot$computed.values
sam_den2_4 <- samdensity_4$density
sam_den2_4 <- data.frame("density" = sam_den2_4, "window"=1:length(sam_den2_4))
sum(sam_den2_4$density !=0)
sum(sam_den2_4$density)

overdensity_4 <- kp$latest.plot$computed.values
over_den2_4 <- overdensity_4$density
over_den2_4 <- data.frame("density" = over_den2_4, "window"=1:length(over_den2_4))
sum(over_den2_4$density !=0)
sum(over_den2_4$density)


kp <- plotKaryotype(genome = m_custom.genome, plot.type=4, 
                    chromosomes = c("NW_024609857.1","NW_024609868.1","NW_024609879.1","NW_024609883.1"),
                    plot.params=plot.params, cex = 1.2)
kp <- kpAddBaseNumbers(kp, tick.dist = 50000000, add.units = TRUE, cex = 1)
at <- autotrack(current.track = 1, total.tracks = 4, margin = 0.1)
kp <- kpPlotDensity(kp, data=rad_range, window.size = 1000, r0=(at$r0)*plot.params$data1max, 
                    r1=(at$r1)*plot.params$data1max, border=transparent("#7FB3EE",0.5), lwd=2)
kp <- kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, cex=1,
             r0=(at$r0)*plot.params$data1max, r1=(at$r1)*plot.params$data1max)
kp <- kpAddLabels(kp, labels="3RAD", cex=1.2, srt=90, pos = 3, label.margin = 0.035,
                  r0=(at$r0)*plot.params$data1max, r1=(at$r1)*plot.params$data1max)
at <- autotrack(current.track = 2, total.tracks = 4, margin = 0.1)
kp <- kpPlotDensity(kp, data=gatk_range, window.size = 1000, r0=(at$r0)*plot.params$data1max, 
                    r1=(at$r1)*plot.params$data1max, border=transparent("#332851",0.5), lwd=2)
kp <- kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, cex=1,
             r0=(at$r0)*plot.params$data1max, r1=(at$r1)*plot.params$data1max)
kp <- kpAddLabels(kp, labels="GATK", cex=1.2, srt=90, pos = 3, label.margin = 0.035,
                  r0=(at$r0)*plot.params$data1max, r1=(at$r1)*plot.params$data1max)
at <- autotrack(current.track = 3, total.tracks = 4, margin = 0.1)
kp <- kpPlotDensity(kp, data=sam_range, window.size = 1000, r0=(at$r0)*plot.params$data1max, 
                    r1=(at$r1)*plot.params$data1max, border=transparent("#CA3074",0.5), lwd=2)
kp <- kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, cex=1,
             r0=(at$r0)*plot.params$data1max, r1=(at$r1)*plot.params$data1max)
kp <- kpAddLabels(kp, labels="Samtools", cex=1.2, srt=90, pos = 3, label.margin = 0.035,
                  r0=(at$r0)*plot.params$data1max, r1=(at$r1)*plot.params$data1max)
at <- autotrack(current.track = 4, total.tracks = 4, margin = 0.1)
kp <- kpPlotDensity(kp, data=overlap_range, window.size = 1000, r0=(at$r0)*plot.params$data1max, 
                    r1=(at$r1)*plot.params$data1max, border=transparent("#F6C667",0.5), lwd=2)
kp <- kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, cex=1,
             r0=(at$r0)*plot.params$data1max, r1=(at$r1)*plot.params$data1max)
kp <- kpAddLabels(kp, labels="Overlap", cex=1.2, srt=90, pos = 3, label.margin = 0.035,
                  r0=(at$r0)*plot.params$data1max, r1=(at$r1)*plot.params$data1max)

save.image(file="./karyoploter_reseq_data.RData")
