library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(corrplot)
library(tidyr)
library(plyr)

Cni_orthologues <- c(29,29,29,29,29,29,29,27,23)

# Development
Cni2Cgi_mean <- read.table("Cni2Cgi_subset_JSD_mean.txt", header = T)
Cni2Pfu_mean <- read.table("Cni2Pfu_subset_JSD_mean.txt", header = T)
Cni2Pye_mean <- read.table("Cni2Pye_subset_JSD_mean.txt", header = T)
Cni2Hru_mean <- read.table("Cni2Hru_subset_JSD_mean.txt", header = T)
Cni2Pca_mean <- read.table("Cni2Pca_subset_JSD_mean.txt", header = T)
Cni2Lan_mean <- read.table("Cni2Lan_subset_JSD_mean.txt", header = T)
Cni2Ofu_mean <- read.table("Cni2Ofu_subset_JSD_mean.txt", header = T)
Cni2Aam_mean <- read.table("Cni2Aam_subset_JSD_mean.txt", header = T)

Cni2Cgi_mean_divided <- Cni2Cgi_mean/Cni_orthologues[1]
Cni2Pfu_mean_divided <- Cni2Pfu_mean/Cni_orthologues[2]
Cni2Pye_mean_divided <- Cni2Pye_mean/Cni_orthologues[3]
Cni2Hru_mean_divided <- Cni2Hru_mean/Cni_orthologues[4]
Cni2Pca_mean_divided <- Cni2Pca_mean/Cni_orthologues[5]
Cni2Lan_mean_divided <- Cni2Lan_mean/Cni_orthologues[6]
Cni2Ofu_mean_divided <- Cni2Ofu_mean/Cni_orthologues[7]
Cni2Aam_mean_divided <- Cni2Aam_mean/Cni_orthologues[8]

allvalues_divid <- rbind(Cni2Cgi_mean_divided,Cni2Pfu_mean_divided,Cni2Pye_mean_divided,Cni2Hru_mean_divided,
                         Cni2Pca_mean_divided,Cni2Lan_mean_divided,Cni2Ofu_mean_divided,Cni2Aam_mean_divided)

minimum_divided <- min(allvalues_divid)
maximum_divided <- max(allvalues_divid)

Cni2Cgi_final_normalised <- (Cni2Cgi_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cni2Pfu_final_normalised <- (Cni2Pfu_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cni2Pye_final_normalised <- (Cni2Pye_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cni2Hru_final_normalised <- (Cni2Hru_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cni2Pca_final_normalised <- (Cni2Pca_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cni2Lan_final_normalised <- (Cni2Lan_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cni2Ofu_final_normalised <- (Cni2Ofu_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cni2Aam_final_normalised <- (Cni2Aam_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)

write.table(Cni2Cgi_final_normalised, paste0("Cni2Cgi_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cni2Pfu_final_normalised, paste0("Cni2Pfu_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cni2Pye_final_normalised, paste0("Cni2Pye_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cni2Hru_final_normalised, paste0("Cni2Hru_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cni2Pca_final_normalised, paste0("Cni2Pca_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cni2Lan_final_normalised, paste0("Cni2Lan_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cni2Ofu_final_normalised, paste0("Cni2Ofu_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cni2Aam_final_normalised, paste0("Cni2Aam_final_normalised",".txt"), sep='\t', quote = FALSE)

paletteLength <- 100
heatmap_color <- colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(paletteLength)

myBreaks <- c(seq(1/paletteLength, 1, length.out=floor(paletteLength)))
h1 <- pheatmap(Cni2Cgi_final_normalised,
               cluster_rows = F, cluster_cols = F, cellheight = 15, cellwidth = 15,
               border_color = NA, color = heatmap_color, breaks = myBreaks,
               labels_col = c("G","T","D","LD","J","M"),
               labels_row = c("B","RM","FS", "EG", "G", "T1", "T2", "T3", "T4", "T5",
                              "ED1", "ED2", "D1", "D2", "D3", "D4", "D5", "D6", "D7",
                              "EU1", "EU2", "U1", "U2", "U3", "U4", "U5", "U6",
                              "LU1", "LU2", "P1", "P2", "S", "J", "MA"))

h2 <- pheatmap(Cni2Pfu_final_normalised,
               cluster_rows = F, cluster_cols = F, cellheight = 15, cellwidth = 15,
               border_color = NA, color = heatmap_color, breaks = myBreaks,
               labels_col = c("G","T","D","LD","J","M"),
               labels_row = c("FG","T", "D", "U", "P", "S", "M"))

h3 <- pheatmap(Cni2Pye_final_normalised,
               cluster_rows = F, cluster_cols = F, cellheight = 15, cellwidth = 15,
               border_color = NA, color = heatmap_color, breaks = myBreaks,
               labels_col = c("G","T","D","LD","J","M"),
               labels_row = c("B", "G", "T", "D", "U1", "U2", "U3", "P", "S", "J", "M"))

h4 <- pheatmap(Cni2Hru_final_normalised,
               cluster_rows = F, cluster_cols = F, cellheight = 15, cellwidth = 15,
               border_color = NA, color = heatmap_color, breaks = myBreaks,
               labels_col = c("G","T","D","LD","J","M"),
               labels_row = c("oneday", "sixdays", "tendays", "twentyoneday", "M"))

h5 <- pheatmap(Cni2Pca_final_normalised,
               cluster_rows = F, cluster_cols = F, cellheight = 15, cellwidth = 15,
               border_color = NA, color = heatmap_color, breaks = myBreaks,
               labels_col = c("G","T","D","LD","J","M"),
               labels_row = c("2dpf", "3dpf", "4dpf", "5dpf", "6dpf", "7dpf", "9dpf",
                              "11dpf", "13dpf", "16dpf", "19dpf", "mantle"))

h6 <- pheatmap(Cni2Lan_final_normalised,
               cluster_rows = F, cluster_cols = F, cellheight = 15, cellwidth = 15,
               border_color = NA, color = heatmap_color, breaks = myBreaks,
               labels_col = c("G","T","D","LD","J","M"),
               labels_row = c("X128.cell_to_early_blastula", "Early_blastula", "Blastula",
                              "Early_gastrula", "Mid_gastrula", "Late_gastrula","X1.pair.cirri_larva",
                              "X2.pair.cirri_larva", "adult_ventral_mantle"))

h7 <- pheatmap(Cni2Ofu_final_normalised,
               cluster_rows = F, cluster_cols = F, cellheight = 15, cellwidth = 15,
               border_color = NA, color = heatmap_color, breaks = myBreaks,
               labels_col = c("G","T","D","LD","J","M"),
               labels_row = c("blastula", "gastrula", "elongation", "early_larva",
                              "mitraria_larva", "competent_larva", "juvenile", "adult_head_chaetae"))

h8 <- pheatmap(Cni2Aam_final_normalised,
               cluster_rows = F, cluster_cols = F, cellheight = 15, cellwidth = 15,
               border_color = NA, color = heatmap_color, breaks = myBreaks,
               labels_col = c("G","T","D","LD","J","M"),
               labels_row = c("Nauplius_II","Nauplius_IV","Nauplius_VI","Free_swimming_cyprid","Close_searching_cyprid","Settled_cyprid","Juvenile","Mantle"))

# Save heatmap for Cni2Cgi_mean
pdf("Cni2Cgi_normalised_mean_heatmap.pdf")
h1
dev.off()

pdf("Cni2Pfu_normalised_mean_heatmap.pdf")
h2
dev.off()

pdf("Cni2Pye_normalised_mean_heatmap.pdf")
h3
dev.off()

pdf("Cni2Hru_normalised_mean_heatmap.pdf")
h4
dev.off()

pdf("Cni2Pca_normalised_mean_heatmap.pdf")
h5
dev.off()

pdf("Cni2Lan_normalised_mean_heatmap.pdf")
h6
dev.off()

pdf("Cni2Ofu_normalised_mean_heatmap.pdf")
h7
dev.off()

pdf("Cni2Aam_normalised_mean_heatmap.pdf")
h8
dev.off()

Cni2Cgi_sd <- read.table("Cni2Cgi_subset_JSD_sd.txt", header = T)
Cni2Pfu_sd <- read.table("Cni2Pfu_subset_JSD_sd.txt", header = T)
Cni2Pye_sd <- read.table("Cni2Pye_subset_JSD_sd.txt", header = T)
Cni2Hru_sd <- read.table("Cni2Hru_subset_JSD_sd.txt", header = T)
Cni2Pca_sd <- read.table("Cni2Pca_subset_JSD_sd.txt", header = T)
Cni2Lan_sd <- read.table("Cni2Lan_subset_JSD_sd.txt", header = T)
Cni2Ofu_sd <- read.table("Cni2Ofu_subset_JSD_sd.txt", header = T)
Cni2Aam_sd <- read.table("Cni2Aam_subset_JSD_sd.txt", header = T)

Cni2Cgi_transposed <- data.frame(t(Cni2Cgi_mean))
Cni2Pfu_transposed <- data.frame(t(Cni2Pfu_mean))
Cni2Pye_transposed <- data.frame(t(Cni2Pye_mean))
Cni2Hru_transposed <- data.frame(t(Cni2Hru_mean))
Cni2Pca_transposed <- data.frame(t(Cni2Pca_mean))
Cni2Lan_transposed <- data.frame(t(Cni2Lan_mean))
Cni2Ofu_transposed <- data.frame(t(Cni2Ofu_mean))
Cni2Aam_transposed <- data.frame(t(Cni2Aam_mean))

min_stages <- data.frame(replicate(8,sample(0:1,ncol(Cni2Cgi_mean),rep=TRUE)))
colnames(min_stages) <- c("Cgi","Pfu","Pye","Hru","Pca","Lan","Ofu","Aam")
rownames(min_stages) <- c("G","T","D","LD","J","M")

min_stages$Cgi <- names(Cni2Cgi_transposed)[apply(Cni2Cgi_transposed, MARGIN = 1, FUN = which.min)]
min_stages$Pfu <- names(Cni2Pfu_transposed)[apply(Cni2Pfu_transposed, MARGIN = 1, FUN = which.min)]
min_stages$Pye <- names(Cni2Pye_transposed)[apply(Cni2Pye_transposed, MARGIN = 1, FUN = which.min)]
min_stages$Hru <- names(Cni2Hru_transposed)[apply(Cni2Hru_transposed, MARGIN = 1, FUN = which.min)]
min_stages$Pca <- names(Cni2Pca_transposed)[apply(Cni2Pca_transposed, MARGIN = 1, FUN = which.min)]
min_stages$Lan <- names(Cni2Lan_transposed)[apply(Cni2Lan_transposed, MARGIN = 1, FUN = which.min)]
min_stages$Ofu <- names(Cni2Ofu_transposed)[apply(Cni2Ofu_transposed, MARGIN = 1, FUN = which.min)]
min_stages$Aam <- names(Cni2Aam_transposed)[apply(Cni2Aam_transposed, MARGIN = 1, FUN = which.min)]

minimum_distance_mean <- data.frame(replicate(8,sample(0:1,ncol(Cni2Cgi_mean),rep=TRUE)))
minimum_distance_sd <- data.frame(replicate(8,sample(0:1,ncol(Cni2Cgi_mean),rep=TRUE)))
colnames(minimum_distance_mean) <- c("Cgi","Pfu","Pye","Hru","Pca","Lan","Ofu","Aam")
rownames(minimum_distance_mean) <- c("G","T","D","LD","J","M")
colnames(minimum_distance_sd) <- c("Cgi","Pfu","Pye","Hru","Pca","Lan","Ofu","Aam")
rownames(minimum_distance_sd) <- c("G","T","D","LD","J","M")

for (i in c(1:ncol(Cni2Cgi_mean))){
  minimum_distance_mean[i,"Cgi"] <- Cni2Cgi_mean[min_stages$Cgi[i],i]
  minimum_distance_mean[i,"Pfu"] <- Cni2Pfu_mean[min_stages$Pfu[i],i]
  minimum_distance_mean[i,"Pye"] <- Cni2Pye_mean[min_stages$Pye[i],i]
  minimum_distance_mean[i,"Hru"] <- Cni2Hru_mean[min_stages$Hru[i],i]
  minimum_distance_mean[i,"Pca"] <- Cni2Pca_mean[min_stages$Pca[i],i]
  minimum_distance_mean[i,"Lan"] <- Cni2Lan_mean[min_stages$Lan[i],i]
  minimum_distance_mean[i,"Ofu"] <- Cni2Ofu_mean[min_stages$Ofu[i],i]
  minimum_distance_mean[i,"Aam"] <- Cni2Aam_mean[min_stages$Aam[i],i]
  minimum_distance_sd[i,"Cgi"] <- Cni2Cgi_sd[min_stages$Cgi[i],i]
  minimum_distance_sd[i,"Pfu"] <- Cni2Pfu_sd[min_stages$Pfu[i],i]
  minimum_distance_sd[i,"Pye"] <- Cni2Pye_sd[min_stages$Pye[i],i]
  minimum_distance_sd[i,"Hru"] <- Cni2Hru_sd[min_stages$Hru[i],i]
  minimum_distance_sd[i,"Pca"] <- Cni2Pca_sd[min_stages$Pca[i],i]
  minimum_distance_sd[i,"Lan"] <- Cni2Lan_sd[min_stages$Lan[i],i]
  minimum_distance_sd[i,"Ofu"] <- Cni2Ofu_sd[min_stages$Ofu[i],i]
  minimum_distance_sd[i,"Aam"] <- Cni2Aam_sd[min_stages$Aam[i],i]
}

minimum_distance_mean$stage <- c(1:ncol(Cni2Cgi_mean))
minimum_distance_sd$stage <- c(1:ncol(Cni2Cgi_mean))
minimum_distance_mean_tidy <- gather(minimum_distance_mean, "species", "JSD", -stage)
minimum_distance_sd_tidy <- gather(minimum_distance_sd, "species", "JSD", -stage)
minimum_distance_mean_tidy$upper <- minimum_distance_mean_tidy["JSD"] + minimum_distance_sd_tidy["JSD"]
minimum_distance_mean_tidy$lower <- minimum_distance_mean_tidy["JSD"] - minimum_distance_sd_tidy["JSD"]
colnames(minimum_distance_mean_tidy)[4:5] <- c("upper","lower")
minimum_distance_mean_tidy$upper <- unlist(minimum_distance_mean_tidy$upper)
minimum_distance_mean_tidy$lower <- unlist(minimum_distance_mean_tidy$lower)

minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Cgi'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Cgi'),3:5]-min(minimum_distance_mean$Cgi))/(max(minimum_distance_mean$Cgi)-min(minimum_distance_mean$Cgi))
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Pfu'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Pfu'),3:5]-min(minimum_distance_mean$Pfu))/(max(minimum_distance_mean$Pfu)-min(minimum_distance_mean$Pfu))
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Pye'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Pye'),3:5]-min(minimum_distance_mean$Pye))/(max(minimum_distance_mean$Pye)-min(minimum_distance_mean$Pye))
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Hru'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Hru'),3:5]-min(minimum_distance_mean$Hru))/(max(minimum_distance_mean$Hru)-min(minimum_distance_mean$Hru))
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Pca'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Pca'),3:5]-min(minimum_distance_mean$Pca))/(max(minimum_distance_mean$Pca)-min(minimum_distance_mean$Pca))
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Lan'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Lan'),3:5]-min(minimum_distance_mean$Lan))/(max(minimum_distance_mean$Lan)-min(minimum_distance_mean$Lan))
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Ofu'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Ofu'),3:5]-min(minimum_distance_mean$Ofu))/(max(minimum_distance_mean$Ofu)-min(minimum_distance_mean$Ofu))
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Aam'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Aam'),3:5]-min(minimum_distance_mean$Aam))/(max(minimum_distance_mean$Aam)-min(minimum_distance_mean$Aam))

final_dataset <- minimum_distance_mean_tidy
final_dataset$species <- factor(final_dataset$species,
                                levels = c("Cgi","Pfu","Pye","Hru","Pca","Lan","Ofu", "Aam"))

p = ggplot(final_dataset, aes(x=stage, y=JSD, colour=species)) +
  geom_line(show.legend = FALSE) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=factor(species)),
              colour = NA, show.legend = FALSE,
              alpha = 0.5) +
  facet_wrap(~species, nrow=1, ncol=10) +
  scale_x_continuous(labels = c(1:ncol(Cni2Cgi_mean)), breaks = c(1:ncol(Cni2Cgi_mean))) +
  theme_classic() +
  labs(x = "C. nippona dev", y = "Normalised gene JSD divergence (JSD)")
ggsave("JSD_plot.pdf", plot = p, width = 30, height = 8, device = "pdf")

write.table(min_stages, "minimum_stages.txt", sep='\t', quote = FALSE)
write.table(minimum_distance_mean, "minimum_stages_distance_mean.txt", sep='\t', quote = FALSE)
write.table(minimum_distance_sd, "minimum_stages_distance_sd.txt", sep='\t', quote = FALSE)
write.table(minimum_distance_mean_tidy, "minimum_stages_to_plot.txt", sep='\t', quote = FALSE)

# line
# Cgi
Cgi_selected <- Cni2Cgi_final_normalised[, c("D", "M")]
Cgi_selected$Cgi <- rownames(Cni2Cgi_final_normalised)
Cgi_long <- pivot_longer(Cgi_selected, cols = c("D", "M"), names_to = "Stage", values_to = "JSD")

pdf("Cni2Cgi_normalised_mean_line.pdf", width = 3, height = 3)
ggplot(Cgi_long, aes(x = factor(Cgi, levels = unique(Cgi)), y = JSD, group = Stage)) +
  geom_line(aes(color = Stage), size = 1) +
  theme_classic() +
  labs(x = "Developmental Stage",
       y = "Normalized JSD") +
  theme(legend.title = element_blank(), legend.position = "none") +
  scale_color_manual(values = c("D" = "#27857e", "M" = "#69428f")) +
  theme(legend.position = "none",
        axis.ticks.length = unit(-0.2, "cm"),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(color = "black", size = 8),
        axis.text.y = element_text(color = "black", size = 8),
        plot.title = element_text(color = "black", size = 8))
dev.off()

# Pfu
Pfu_selected <- Cni2Pfu_final_normalised[, c("D", "M")]
Pfu_selected$Pfu <- rownames(Cni2Pfu_final_normalised)
Pfu_long <- pivot_longer(Pfu_selected, cols = c("D", "M"), names_to = "Stage", values_to = "JSD")

pdf("Cni2Pfu_normalised_mean_line.pdf", width = 3, height = 3)
ggplot(Pfu_long, aes(x = factor(Pfu, levels = unique(Pfu)), y = JSD, group = Stage)) +
  geom_line(aes(color = Stage), size = 1) +
  theme_classic() +
  labs(x = "Developmental Stage",
       y = "Normalized JSD") +
  theme(legend.title = element_blank(), legend.position = "none") +
  scale_color_manual(values = c("D" = "#27857e", "M" = "#69428f")) +
  theme(legend.position = "none",
        axis.ticks.length = unit(-0.2, "cm"),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(color = "black", size = 8),
        axis.text.y = element_text(color = "black", size = 8),
        plot.title = element_text(color = "black", size = 8))
dev.off()

# Pye
Pye_selected <- Cni2Pye_final_normalised[, c("D", "M")]
Pye_selected$Pye <- rownames(Cni2Pye_final_normalised)
Pye_long <- pivot_longer(Pye_selected, cols = c("D", "M"), names_to = "Stage", values_to = "JSD")

pdf("Cni2Pye_normalised_mean_line.pdf", width = 3, height = 3)
ggplot(Pye_long, aes(x = factor(Pye, levels = unique(Pye)), y = JSD, group = Stage)) +
  geom_line(aes(color = Stage), size = 1) +
  theme_classic() +
  labs(x = "Developmental Stage",
       y = "Normalized JSD") +
  theme(legend.title = element_blank(), legend.position = "none") +
  scale_color_manual(values = c("D" = "#27857e", "M" = "#69428f")) +
  theme(legend.position = "none",
        axis.ticks.length = unit(-0.2, "cm"),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(color = "black", size = 8),
        axis.text.y = element_text(color = "black", size = 8),
        plot.title = element_text(color = "black", size = 8))
dev.off()

# Hru
Hru_selected <- Cni2Hru_final_normalised[, c("D", "M")]
Hru_selected$Hru <- rownames(Cni2Hru_final_normalised)
Hru_long <- pivot_longer(Hru_selected, cols = c("D", "M"), names_to = "Stage", values_to = "JSD")

pdf("Cni2Hru_normalised_mean_line.pdf", width = 3, height = 3)
ggplot(Hru_long, aes(x = factor(Hru, levels = unique(Hru)), y = JSD, group = Stage)) +
  geom_line(aes(color = Stage), size = 1) +
  theme_classic() +
  labs(x = "Developmental Stage",
       y = "Normalized JSD") +
  theme(legend.title = element_blank(), legend.position = "none") +
  scale_color_manual(values = c("D" = "#27857e", "M" = "#69428f")) +
  theme(legend.position = "none",
        axis.ticks.length = unit(-0.2, "cm"),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(color = "black", size = 8),
        axis.text.y = element_text(color = "black", size = 8),
        plot.title = element_text(color = "black", size = 8))
dev.off()

# Pca
Pca_selected <- Cni2Pca_final_normalised[, c("D", "M")]
Pca_selected$Pca <- rownames(Cni2Pca_final_normalised)
Pca_long <- pivot_longer(Pca_selected, cols = c("D", "M"), names_to = "Stage", values_to = "JSD")

pdf("Cni2Pca_normalised_mean_line.pdf", width = 3, height = 3)
ggplot(Pca_long, aes(x = factor(Pca, levels = unique(Pca)), y = JSD, group = Stage)) +
  geom_line(aes(color = Stage), size = 1) +
  theme_classic() +
  labs(x = "Developmental Stage",
       y = "Normalized JSD") +
  theme(legend.title = element_blank(), legend.position = "none") +
  scale_color_manual(values = c("D" = "#27857e", "M" = "#69428f")) +
  theme(legend.position = "none",
        axis.ticks.length = unit(-0.2, "cm"),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(color = "black", size = 8),
        axis.text.y = element_text(color = "black", size = 8),
        plot.title = element_text(color = "black", size = 8))
dev.off()

# Lan
Lan_selected <- Cni2Lan_final_normalised[, c("D", "M")]
Lan_selected$Lan <- rownames(Cni2Lan_final_normalised)
Lan_long <- pivot_longer(Lan_selected, cols = c("D", "M"), names_to = "Stage", values_to = "JSD")

pdf("Cni2Lan_normalised_mean_line.pdf", width = 3, height = 3)
ggplot(Lan_long, aes(x = factor(Lan, levels = unique(Lan)), y = JSD, group = Stage)) +
  geom_line(aes(color = Stage), size = 1) +
  theme_classic() +
  labs(x = "Developmental Stage",
       y = "Normalized JSD") +
  theme(legend.title = element_blank(), legend.position = "none") +
  scale_color_manual(values = c("D" = "#27857e", "M" = "#69428f")) +
  theme(legend.position = "none",
        axis.ticks.length = unit(-0.2, "cm"),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(color = "black", size = 8),
        axis.text.y = element_text(color = "black", size = 8),
        plot.title = element_text(color = "black", size = 8))
dev.off()

# Ofu
Ofu_selected <- Cni2Ofu_final_normalised[, c("D", "M")]
Ofu_selected$Ofu <- rownames(Cni2Ofu_final_normalised)
Ofu_long <- pivot_longer(Ofu_selected, cols = c("D", "M"), names_to = "Stage", values_to = "JSD")

pdf("Cni2Ofu_normalised_mean_line.pdf", width = 3, height = 3)
ggplot(Ofu_long, aes(x = factor(Ofu, levels = unique(Ofu)), y = JSD, group = Stage)) +
  geom_line(aes(color = Stage), size = 1) +
  theme_classic() +
  labs(x = "Developmental Stage",
       y = "Normalized JSD") +
  theme(legend.title = element_blank(), legend.position = "none") +
  scale_color_manual(values = c("D" = "#27857e", "M" = "#69428f")) +
  theme(legend.position = "none",
        axis.ticks.length = unit(-0.2, "cm"),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(color = "black", size = 8),
        axis.text.y = element_text(color = "black", size = 8),
        plot.title = element_text(color = "black", size = 8))
dev.off()


# Aam
Aam_selected <- Cni2Aam_final_normalised[, c("D", "M")]
Aam_selected$Aam <- rownames(Cni2Aam_final_normalised)
Aam_long <- pivot_longer(Aam_selected, cols = c("D", "M"), names_to = "Stage", values_to = "JSD")

pdf("Cni2Aam_normalised_mean_line.pdf", width = 3, height = 3)
ggplot(Aam_long, aes(x = factor(Aam, levels = unique(Aam)), y = JSD, group = Stage)) +
  geom_line(aes(color = Stage), size = 1) +
  theme_classic() +
  labs(x = "Developmental Stage",
       y = "Normalized JSD") +
  theme(legend.title = element_blank(), legend.position = "none") +
  scale_color_manual(values = c("D" = "#27857e", "M" = "#69428f")) +
  theme(legend.position = "none",
        axis.ticks.length = unit(-0.2, "cm"),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(color = "black", size = 8),
        axis.text.y = element_text(color = "black", size = 8),
        plot.title = element_text(color = "black", size = 8))
dev.off()
