library(qiime2R)
library(tidyverse)
library(RColorBrewer)

setwd("D:/Qiime2/Qiime_output_trim_barcode_primer")

meta<-read_q2metadata("SampleSeq_metadata2.xls.txt")
str(meta)

metadata<-read_q2metadata("SampleSeq_metadata2.xls.txt")
#metadata2 <- read.delim("sample-metadata.tsv", sep = "\t", header = T, quote = "", stringsAsFactors = F)
#metadata2 <- metadata2[-1,]
str(metadata)
levels(metadata$`site`)
metadata$site = factor(metadata$site, levels = c("Duodenum", "Jejunum", "Ileum", "Cecum", "Ascending_colon", "Transverse_colon", "Descending_colon", "Capsule"))
metadata$site

row.names(metadata) <- metadata[ ,1]
#metadata <- metadata[,-1]
row.names(metadata)


#####################################################################################################################
#Bray-Curtis

bc_PCoA<-read_qza("core-metrics-results/bray_curtis_pcoa_results.qza")

#body_colors <- c("Black", "Blue", "Green", "Gray", "Red", "Yellow", "Orange", "Pink")

bc_meta <- bc_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

my_column <- "site"

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=brewer.pal(n = 8, name = "Set1"), name = my_column)
ggsave(paste0("R output/BC-basic_", my_column,".pdf"), height=4, width=5, device="pdf") # save a PDF 3 inches by 4 inches

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),bc_meta,mean)
colnames(centroids)[1] <- "site"

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = my_column)
ggsave(paste0("R output/BC-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes()) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  #stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) 
#scale_color_manual(values=corn_colors, name = my_column)
ggsave(paste0("R output/BC-ellipse_", my_column,"-subject.pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

##SAME thing but with weighted UniFrac


#####################################################################################################################
#Unweighted uniFrac
UnWuni_PCoA<-read_qza("core-metrics-results/unweighted_unifrac_pcoa_results.qza")

UnWuni_meta <- UnWuni_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),UnWuni_meta,mean)

ggplot(UnWuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  #stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*UnWuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*UnWuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=brewer.pal(n = 8, name = "Set1"), name = my_column)
ggsave(paste0("R output/UnWuni", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches


#####################################################################################################################
#weighted UniFrac
Wuni_PCoA<-read_qza("core-metrics-results/weighted_unifrac_pcoa_results.qza")

Wuni_meta <- Wuni_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),Wuni_meta,mean)

ggplot(Wuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  #stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Wuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Wuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=brewer.pal(n = 8, name = "Set1"), name = my_column)
ggsave(paste0("R output/Wuni", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(Wuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(), size = 3) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Wuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Wuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "Body Site")
ggsave(paste0("R output/Wuni-ellipse_", my_column,"-subject.pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches


#####################################################################################################################
#Jaacard
Jaacard_PCoA<-read_qza("core-metrics-results/jaccard_pcoa_results.qza")

Jaacard_meta <- Jaacard_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),Jaacard_meta,mean)

ggplot(Jaacard_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  #stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Jaacard_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Jaacard_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=brewer.pal(n = 8, name = "Set1"), name = my_column)
ggsave(paste0("R output/Jaacard", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches
