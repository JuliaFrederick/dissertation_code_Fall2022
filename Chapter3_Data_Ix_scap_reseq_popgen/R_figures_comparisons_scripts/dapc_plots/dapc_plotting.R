#setwd("/Users/juliafrederick/UGA_EHS Dropbox/Julia Frederick/Julia_Frederick_Projects/Ixscap_reseq_population_genetics/popgen/") #my current directory

load("./3RAD/pop_gen_cluster/DAPC_BStats_reseq_3RAD_m90_6x150x_GQ30_maf05_ma2.RData")
dapc.plot.3RAD <- scatter(dapc.pop, cex = 2, clab= 0,
                     posi.pca = "topleft", scree.pca = TRUE,
                     posi.da = "bottomleft",
                     legend = TRUE,
                     xlim=c(20,-15), ylim=c(-20,20),
                     inset.solid = 0.7)
dapc.plot.rm.3RAD <- scatter(dapc.rm, cex = 2, clab= 0,
                        posi.pca = "topleft", scree.pca = TRUE,
                        posi.da = "bottomleft",
                        legend = TRUE, cleg = 1,
                        xlim=c(20,-15), ylim=c(-20,20),
                        inset.solid = 0.7)

load("./gatk/pop_gen_cluster/DAPC_BStats_reseq_gatk_m90_6x50x_GQ30_Q30_maf05_ma2.RData")
dapc.plot.gatk <- scatter(dapc.pop, cex = 2, clab= 0,
                          posi.pca = "topleft", scree.pca = TRUE,
                          posi.da = "bottomleft",
                          legend = TRUE,
                          xlim=c(-20,15), ylim=c(20,-20),
                          inset.solid = 0.7)
dapc.plot.rm.gatk <- scatter(dapc.rm, cex = 2, clab= 0,
                             posi.pca = "topleft", scree.pca = TRUE,
                             posi.da = "bottomleft",
                             legend = TRUE, cleg = 1,
                             xlim=c(-20,15), ylim=c(20,-20),
                             inset.solid = 0.7)

load("./samtools/pop_gen_cluster/DAPC_BStats_reseq_samtools_formatQ_m90_6x50x_GQ30_Q30_maf05_ma2.RData")
dapc.plot.sam<- scatter(dapc.pop, cex = 2, clab= 0,
                          posi.pca = "topleft", scree.pca = TRUE,
                          posi.da = "bottomleft",
                          legend = TRUE,
                          xlim=c(20,-15), ylim=c(-20,20),
                          inset.solid = 0.7)
dapc.plot.rm.sam <- scatter(dapc.rm, cex = 2, clab= 0,
                             posi.pca = "topleft", scree.pca = TRUE,
                             posi.da = "bottomleft",
                             legend = TRUE, cleg = 1,
                             xlim=c(20,-15), ylim=c(-20,20),
                             inset.solid = 0.7)

load("./reseq_overlap/pop_gen_cluster/DAPC_BStats_reseq_overlap_m90_6x50x_GQ30_Q30_maf05_ma2.RData")
dapc.plot.overlap <- scatter(dapc.pop, cex = 2, clab= 0,
                        posi.pca = "topleft", scree.pca = TRUE,
                        posi.da = "bottomleft",
                        legend = TRUE,
                        xlim=c(20,-15), ylim=c(-20,20),
                        inset.solid = 0.7)
dapc.plot.rm.overlap <- scatter(dapc.rm, cex = 2, clab= 0,
                            posi.pca = "topleft", scree.pca = TRUE,
                            posi.da = "bottomleft",
                            legend = TRUE, cleg = 1,
                            xlim=c(20,-15), ylim=c(-20,20),
                            inset.solid = 0.7)

png(file="./R_scripts/dapc_plots_5pops.png", width = 1000,700)
load("./3RAD/pop_gen_cluster/DAPC_BStats_reseq_3RAD_m90_6x150x_GQ30_maf05_ma2.RData")
dapc.plot.3RAD <- scatter(dapc.pop, cex = 2, clab= 0,
                          posi.pca = "topleft", scree.pca = TRUE,
                          posi.da = "bottomleft",
                          legend = TRUE,
                          xlim=c(20,-15), ylim=c(-20,20),
                          inset.solid = 0.7)
load("./gatk/pop_gen_cluster/DAPC_BStats_reseq_gatk_m90_6x50x_GQ30_Q30_maf05_ma2.RData")
dapc.plot.gatk <- scatter(dapc.pop, cex = 2, clab= 0,
                          posi.pca = "topleft", scree.pca = TRUE,
                          posi.da = "bottomleft",
                          legend = TRUE,
                          xlim=c(-20,15), ylim=c(20,-20),
                          inset.solid = 0.7)
load("./samtools/pop_gen_cluster/DAPC_BStats_reseq_samtools_formatQ_m90_6x50x_GQ30_Q30_maf05_ma2.RData")
dapc.plot.sam<- scatter(dapc.pop, cex = 2, clab= 0,
                        posi.pca = "topleft", scree.pca = TRUE,
                        posi.da = "bottomleft",
                        legend = TRUE,
                        xlim=c(20,-15), ylim=c(-20,20),
                        inset.solid = 0.7)
load("./reseq_overlap/pop_gen_cluster/DAPC_BStats_reseq_overlap_m90_6x50x_GQ30_Q30_maf05_ma2.RData")
dapc.plot.overlap <- scatter(dapc.pop, cex = 2, clab= 0,
                             posi.pca = "topleft", scree.pca = TRUE,
                             posi.da = "bottomleft",
                             legend = TRUE,
                             xlim=c(20,-15), ylim=c(-20,20),
                             inset.solid = 0.7)
dev.off()
