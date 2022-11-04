#setwd("/Users/juliafrederick/UGA_EHS Dropbox/Julia Frederick/Julia_Frederick_Projects/Ixscap_reseq_population_genetics/popgen/R_scripts/average_variant_distance") #my current directory

library(tidyverse)

rad_vcf <- read.vcfR("../../vcf_files_filtered/3RAD_m90_6x150x_GQ30_maf05_ma2.recode.vcf.gz")
rad_snp <- as.data.frame(rad_vcf@fix)
rad_snp$POS <- as.numeric(rad_snp$POS)

gatk_vcf <- read.vcfR("../../vcf_files_filtered/gatk_m90_6x50x_GQ30_Q30_maf05_ma2.recode.vcf.gz")
gatk_snp <- as.data.frame(gatk_vcf@fix)
gatk_snp$POS <- as.numeric(gatk_snp$POS)

sam_vcf <- read.vcfR("../../vcf_files_filtered/samtools_formatQ_m90_6x50x_GQ30_Q30_maf05_ma2.recode.vcf.gz")
sam_snp <- as.data.frame(sam_vcf@fix)
sam_snp$POS <- as.numeric(sam_snp$POS)

overlap_vcf <- read.vcfR("../../vcf_files_filtered/0000.vcf.gz")
overlap_snp <- as.data.frame(overlap_vcf@fix)
overlap_snp$POS <- as.numeric(overlap_snp$POS)

top14 <- c("NW_024609835.1","NW_024609836.1","NW_024609837.1","NW_024609838.1","NW_024609839.1","NW_024609846.1",
           "NW_024609857.1","NW_024609868.1","NW_024609879.1","NW_024609880.1","NW_024609881.1","NW_024609882.1",
           "NW_024609883.1","NW_024609884.1")

rad_sub <- subset(rad_snp, rad_snp$CHROM %in% c("NW_024609835.1"))
rad_dif <- diff(rad_sub$POS, lag=1)
rad_avg <- mean(rad_dif)
rad_min <- min(rad_dif)
rad_max <- max(rad_dif)
rad_med <- median(rad_dif)
rad_num_snps <- nrow(rad_sub)

names <- c("CHROM","file","avg_difference","min_difference","max_difference","median_difference","total_sites")
file_list <- list(rad_snp, gatk_snp, sam_snp, overlap_snp)
site_dist <- data.frame(matrix(ncol=7, nrow=0))
colnames(site_dist) <- names
for (j in 1:length(file_list)){
file <- file_list[[j]]
if(j==1){
  df_name = "3RAD"
}else{
  if(j==2){
    df_name="GATK"
  }else{
    if(j==3){
      df_name="Samtools"
    }else{
      df_name="Overlap"
    }
  }
}
for (i in top14){
  sub <- subset(file, file$CHROM == i)
  dif <- diff(sub$POS, lag=1)
  avg <- mean(dif)
  min <- min(dif)
  max <- max(dif)
  med <- median(dif)
  num_snps <- nrow(sub)
  
  temp_df <- data.frame(i,df_name, avg, min, max, med, num_snps)
  site_dist <- rbind(site_dist, temp_df)
} }

write.csv(site_dist,"./all_chrom_site_dist.csv")

summary_dist <- site_dist %>% group_by(df_name) %>% 
  summarise("sum_avg"=mean(avg),
            "sum_min"=mean(min),
            "sum_max"=mean(max),
            "sum_med"=mean(med),
            "sum_num_snps"=mean(num_snps))

write.csv(summary_dist,"./summary_byDataset_site_dist.csv")

summary_dist_chrom <- site_dist %>% group_by(i) %>% 
  summarise("sum_avg"=mean(avg),
            "sum_min"=mean(min),
            "sum_max"=mean(max),
            "sum_med"=mean(med),
            "sum_num_snps"=mean(num_snps))

write.csv(summary_dist_chrom,"./summary_byCHROM_site_dist.csv")
