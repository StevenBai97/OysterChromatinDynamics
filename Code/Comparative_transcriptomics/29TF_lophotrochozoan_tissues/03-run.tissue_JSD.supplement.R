library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(corrplot)
library(tidyr)
library(plyr)

Cni_orthologues <- c(29,29,29,29,29,29,27,28,24,27,23)

Cni2Cgi_mean <- read.table("Cni2Cgi_subset_JSD_mean.txt", header = T)
Cni2Pfu_mean <- read.table("Cni2Pfu_subset_JSD_mean.txt", header = T)
Cni2Pye_mean <- read.table("Cni2Pye_subset_JSD_mean.txt", header = T)
Cni2Hru_mean <- read.table("Cni2Hru_subset_JSD_mean.txt", header = T)
Cni2Pca_mean <- read.table("Cni2Pca_subset_JSD_mean.txt", header = T)
Cni2Pve_mean <- read.table("Cni2Pve_subset_JSD_mean.txt", header = T)
Cni2Npo_mean <- read.table("Cni2Npo_subset_JSD_mean.txt", header = T)
Cni2Lan_mean <- read.table("Cni2Lan_subset_JSD_mean.txt", header = T)
Cni2Ofu_mean <- read.table("Cni2Ofu_subset_JSD_mean.txt", header = T)
Cni2Pec_mean <- read.table("Cni2Pec_subset_JSD_mean.txt", header = T)
Cni2Aam_mean <- read.table("Cni2Aam_subset_JSD_mean.txt", header = T)

Cni2Cgi_mean_divided <- Cni2Cgi_mean/Cni_orthologues[1]
Cni2Pfu_mean_divided <- Cni2Pfu_mean/Cni_orthologues[2]
Cni2Pye_mean_divided <- Cni2Pye_mean/Cni_orthologues[3]
Cni2Hru_mean_divided <- Cni2Hru_mean/Cni_orthologues[4]
Cni2Pca_mean_divided <- Cni2Pca_mean/Cni_orthologues[5]
Cni2Pve_mean_divided <- Cni2Pve_mean/Cni_orthologues[6]
Cni2Npo_mean_divided <- Cni2Npo_mean/Cni_orthologues[7]
Cni2Lan_mean_divided <- Cni2Lan_mean/Cni_orthologues[8]
Cni2Ofu_mean_divided <- Cni2Ofu_mean/Cni_orthologues[9]
Cni2Pec_mean_divided <- Cni2Pec_mean/Cni_orthologues[10]
Cni2Aam_mean_divided <- Cni2Aam_mean/Cni_orthologues[11]

allvalues_divid <- rbind(Cni2Cgi_mean_divided,Cni2Pfu_mean_divided,Cni2Pye_mean_divided,
                         Cni2Hru_mean_divided,Cni2Pca_mean_divided,
                         Cni2Pve_mean_divided,Cni2Npo_mean_divided,
                         Cni2Lan_mean_divided,Cni2Ofu_mean_divided,
                         Cni2Pec_mean_divided,Cni2Aam_mean_divided)

minimum_divided <- min(allvalues_divid)
maximum_divided <- max(allvalues_divid)

Cni2Cgi_final_normalised <- (Cni2Cgi_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cni2Pfu_final_normalised <- (Cni2Pfu_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cni2Pye_final_normalised <- (Cni2Pye_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cni2Hru_final_normalised <- (Cni2Hru_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cni2Pca_final_normalised <- (Cni2Pca_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cni2Pve_final_normalised <- (Cni2Pve_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cni2Npo_final_normalised <- (Cni2Npo_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cni2Lan_final_normalised <- (Cni2Lan_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cni2Ofu_final_normalised <- (Cni2Ofu_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cni2Pec_final_normalised <- (Cni2Pec_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cni2Aam_final_normalised <- (Cni2Aam_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)

write.table(Cni2Cgi_final_normalised, paste0("Cni2Cgi_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cni2Pfu_final_normalised, paste0("Cni2Pfu_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cni2Pye_final_normalised, paste0("Cni2Pye_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cni2Hru_final_normalised, paste0("Cni2Hru_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cni2Pca_final_normalised, paste0("Cni2Pca_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cni2Pve_final_normalised, paste0("Cni2Pve_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cni2Npo_final_normalised, paste0("Cni2Npo_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cni2Lan_final_normalised, paste0("Cni2Lan_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cni2Ofu_final_normalised, paste0("Cni2Ofu_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cni2Pec_final_normalised, paste0("Cni2Pec_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cni2Aam_final_normalised, paste0("Cni2Aam_final_normalised",".txt"), sep='\t', quote = FALSE)

paletteLength <- 100
myBreaks <- c(seq(1/paletteLength, 1, length.out=floor(paletteLength)))
heatmap_color <- colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(paletteLength)


h1 <- pheatmap(Cni2Cgi_final_normalised,
         cluster_rows = F, cluster_cols = F, cellheight = 15, cellwidth = 15,
         border_color = NA, color = heatmap_color, breaks = myBreaks,
         labels_col = c("M","AM","Gi","DG","He"),
         labels_row = c("MA","AM","FG","MG","Gi","DG","He","LP"))

h2 <- pheatmap(Cni2Pfu_final_normalised,
         cluster_rows = F, cluster_cols = F, cellheight = 15, cellwidth = 15,
         border_color = NA, color = heatmap_color, breaks = myBreaks,
         labels_col = c("M","AM","Gi","DG","He"),
         labels_row = c("M","AM","Go","Gi"))

h3 <- pheatmap(Cni2Pye_final_normalised,
         cluster_rows = F, cluster_cols = F, cellheight = 15, cellwidth = 15,
         border_color = NA, color = heatmap_color, breaks = myBreaks,
         labels_col = c("M","AM","Gi","DG","He"),
         labels_row = c("M","AM","F","FG","MG","Gi","DG","K"))

h4 <- pheatmap(Cni2Hru_final_normalised,
         cluster_rows = F, cluster_cols = F, cellheight = 15, cellwidth = 15,
         border_color = NA, color = heatmap_color, breaks = myBreaks,
         labels_col = c("M","AM","Gi","DG","He"),
         labels_row = c("M","F","PE","LR","MG","FG","H","CT","Ep","ET","Ga","Gi","L","K"))

h5 <- pheatmap(Cni2Pca_final_normalised,
         cluster_rows = F, cluster_cols = F, cellheight = 15, cellwidth = 15,
         border_color = NA, color = heatmap_color, breaks = myBreaks,
         labels_col = c("M","AM","Gi","DG","He"),
         labels_row = c("mantle","Albumen_gland","Digestive_gland","Foot","Gill","Kidney","Lung","Stomach","Testis"))

h6 <- pheatmap(Cni2Pve_final_normalised,
         cluster_rows = F, cluster_cols = F, cellheight = 15, cellwidth = 15,
         border_color = NA, color = heatmap_color, breaks = myBreaks,
         labels_col = c("M","AM","Gi","DG","He"),
         labels_row = c("Mantle_posterior_aperture","Mantle","Radular_sac","Captacula","Digestive_diverticulae","Foot","Gonad","Intestine","Odontophore","Longitudinal_muscle","Oral_tube"))

h7 <- pheatmap(Cni2Npo_final_normalised,
         cluster_rows = F, cluster_cols = F, cellheight = 15, cellwidth = 15,
         border_color = NA, color = heatmap_color, breaks = myBreaks,
         labels_col = c("M","AM","Gi","DG","He"),
         labels_row = c("Mantle","Eye","Gill","Hemocyte","Liver"))

h8 <- pheatmap(Cni2Lan_final_normalised,
               cluster_rows = F, cluster_cols = F, cellheight = 15, cellwidth = 15,
               border_color = NA, color = heatmap_color, breaks = myBreaks,
               labels_col = c("M","AM","Gi","DG","He"),
               labels_row = c("adult_ventral_mantle","adult_dorsal_mantle","adult_digestive_cecum","adult_gut","adult_lophophore","adult_pedicle"))

h9 <- pheatmap(Cni2Ofu_final_normalised,
               cluster_rows = F, cluster_cols = F, cellheight = 15, cellwidth = 15,
               border_color = NA, color = heatmap_color, breaks = myBreaks,
               labels_col = c("M","AM","Gi","DG","He"),
               labels_row = c("adult_head","adult_head_chaetae","adult_tail","adult_body","adult_blood","adult_gut","adult_retractor-muscle","adult_testis","adult_ovary"))

h10 <- pheatmap(Cni2Pec_final_normalised,
               cluster_rows = F, cluster_cols = F, cellheight = 15, cellwidth = 15,
               border_color = NA, color = heatmap_color, breaks = myBreaks,
               labels_col = c("M","AM","Gi","DG","He"),
               labels_row = c("Collar","Opisthosoma","Plume","Trophosome_anterior","Trophosome_middle","Trophosome_posterior","Vestimentum"))

h11 <- pheatmap(Cni2Aam_final_normalised,
                cluster_rows = F, cluster_cols = F, cellheight = 15, cellwidth = 15,
                border_color = NA, color = heatmap_color, breaks = myBreaks,
                labels_col = c("M","AM","Gi","DG","He"),
                labels_row = c("Mantle","Adductor","Cirrus","Visceral_mass"))


# Save heatmap
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

pdf("Cni2Pve_normalised_mean_heatmap.pdf")
h6
dev.off()

pdf("Cni2Npo_normalised_mean_heatmap.pdf")
h7
dev.off()

pdf("Cni2Lan_normalised_mean_heatmap.pdf")
h8
dev.off()

pdf("Cni2Ofu_normalised_mean_heatmap.pdf")
h9
dev.off()

pdf("Cni2Pec_normalised_mean_heatmap.pdf")
h10
dev.off()

pdf("Cni2Aam_normalised_mean_heatmap.pdf")
h11
dev.off()

Cni2Cgi_sd <- read.table("Cni2Cgi_subset_JSD_sd.txt", header = T)
Cni2Pfu_sd <- read.table("Cni2Pfu_subset_JSD_sd.txt", header = T)
Cni2Pye_sd <- read.table("Cni2Pye_subset_JSD_sd.txt", header = T)
Cni2Hru_sd <- read.table("Cni2Hru_subset_JSD_sd.txt", header = T)
Cni2Pca_sd <- read.table("Cni2Pca_subset_JSD_sd.txt", header = T)
Cni2Pve_sd <- read.table("Cni2Pve_subset_JSD_sd.txt", header = T)
Cni2Npo_sd <- read.table("Cni2Npo_subset_JSD_sd.txt", header = T)
Cni2Lan_sd <- read.table("Cni2Lan_subset_JSD_sd.txt", header = T)
Cni2Ofu_sd <- read.table("Cni2Ofu_subset_JSD_sd.txt", header = T)
Cni2Pec_sd <- read.table("Cni2Pec_subset_JSD_sd.txt", header = T)
Cni2Aam_sd <- read.table("Cni2Aam_subset_JSD_sd.txt", header = T)

Cni2Cgi_transposed <- data.frame(t(Cni2Cgi_mean))
Cni2Pfu_transposed <- data.frame(t(Cni2Pfu_mean))
Cni2Pye_transposed <- data.frame(t(Cni2Pye_mean))
Cni2Hru_transposed <- data.frame(t(Cni2Hru_mean))
Cni2Pca_transposed <- data.frame(t(Cni2Pca_mean))
Cni2Pve_transposed <- data.frame(t(Cni2Pve_mean))
Cni2Npo_transposed <- data.frame(t(Cni2Npo_mean))
Cni2Lan_transposed <- data.frame(t(Cni2Lan_mean))
Cni2Ofu_transposed <- data.frame(t(Cni2Ofu_mean))
Cni2Pec_transposed <- data.frame(t(Cni2Pec_mean))
Cni2Aam_transposed <- data.frame(t(Cni2Aam_mean))

min_tissues <- data.frame(replicate(11,sample(0:1,ncol(Cni2Cgi_mean),rep=TRUE)))
colnames(min_tissues) <- c("Cgi","Pfu","Pye","Hru","Pca","Pve","Npo","Lan","Ofu","Pec","Aam")
rownames(min_tissues) <- c("M","AM","Gi","DG","He")

min_tissues$Cgi <- names(Cni2Cgi_transposed)[apply(Cni2Cgi_transposed, MARGIN = 1, FUN = which.min)]
min_tissues$Pfu <- names(Cni2Pfu_transposed)[apply(Cni2Pfu_transposed, MARGIN = 1, FUN = which.min)]
min_tissues$Pye <- names(Cni2Pye_transposed)[apply(Cni2Pye_transposed, MARGIN = 1, FUN = which.min)]
min_tissues$Hru <- names(Cni2Hru_transposed)[apply(Cni2Hru_transposed, MARGIN = 1, FUN = which.min)]
min_tissues$Pca <- names(Cni2Pca_transposed)[apply(Cni2Pca_transposed, MARGIN = 1, FUN = which.min)]
min_tissues$Pve <- names(Cni2Pve_transposed)[apply(Cni2Pve_transposed, MARGIN = 1, FUN = which.min)]
min_tissues$Npo <- names(Cni2Npo_transposed)[apply(Cni2Npo_transposed, MARGIN = 1, FUN = which.min)]
min_tissues$Lan <- names(Cni2Lan_transposed)[apply(Cni2Lan_transposed, MARGIN = 1, FUN = which.min)]
min_tissues$Ofu <- names(Cni2Ofu_transposed)[apply(Cni2Ofu_transposed, MARGIN = 1, FUN = which.min)]
min_tissues$Pec <- names(Cni2Pec_transposed)[apply(Cni2Pec_transposed, MARGIN = 1, FUN = which.min)]
min_tissues$Aam <- names(Cni2Aam_transposed)[apply(Cni2Aam_transposed, MARGIN = 1, FUN = which.min)]

minimum_distance_mean <- data.frame(replicate(11,sample(0:1,ncol(Cni2Cgi_mean),rep=TRUE)))
minimum_distance_sd <- data.frame(replicate(11,sample(0:1,ncol(Cni2Cgi_mean),rep=TRUE)))
colnames(minimum_distance_mean) <- c("Cgi","Pfu","Pye","Hru","Pca","Pve","Npo","Lan","Ofu","Pec","Aam")
rownames(minimum_distance_mean) <- c("M","AM","Gi","DG","He")
colnames(minimum_distance_sd) <- c("Cgi","Pfu","Pye","Hru","Pca","Pve","Npo","Lan","Ofu","Pec","Aam")
rownames(minimum_distance_sd) <- c("M","AM","Gi","DG","He")

for (i in c(1:ncol(Cni2Cgi_mean))){
  minimum_distance_mean[i,"Cgi"] <- Cni2Cgi_mean[min_tissues$Cgi[i],i]
  minimum_distance_mean[i,"Pfu"] <- Cni2Pfu_mean[min_tissues$Pfu[i],i]
  minimum_distance_mean[i,"Pye"] <- Cni2Pye_mean[min_tissues$Pye[i],i]
  minimum_distance_mean[i,"Hru"] <- Cni2Hru_mean[min_tissues$Hru[i],i]
  minimum_distance_mean[i,"Pca"] <- Cni2Pca_mean[min_tissues$Pca[i],i]
  minimum_distance_mean[i,"Pve"] <- Cni2Pve_mean[min_tissues$Pve[i],i]
  minimum_distance_mean[i,"Npo"] <- Cni2Npo_mean[min_tissues$Npo[i],i]
  minimum_distance_mean[i,"Lan"] <- Cni2Lan_mean[min_tissues$Lan[i],i]
  minimum_distance_mean[i,"Ofu"] <- Cni2Ofu_mean[min_tissues$Ofu[i],i]
  minimum_distance_mean[i,"Pec"] <- Cni2Pec_mean[min_tissues$Pec[i],i]
  minimum_distance_mean[i,"Aam"] <- Cni2Aam_mean[min_tissues$Aam[i],i]
  minimum_distance_sd[i,"Cgi"] <- Cni2Cgi_sd[min_tissues$Cgi[i],i]
  minimum_distance_sd[i,"Pfu"] <- Cni2Pfu_sd[min_tissues$Pfu[i],i]
  minimum_distance_sd[i,"Pye"] <- Cni2Pye_sd[min_tissues$Pye[i],i]
  minimum_distance_sd[i,"Hru"] <- Cni2Hru_sd[min_tissues$Hru[i],i]
  minimum_distance_sd[i,"Pca"] <- Cni2Pca_sd[min_tissues$Pca[i],i]
  minimum_distance_sd[i,"Pve"] <- Cni2Pve_sd[min_tissues$Pve[i],i]
  minimum_distance_sd[i,"Npo"] <- Cni2Npo_sd[min_tissues$Npo[i],i]
  minimum_distance_sd[i,"Lan"] <- Cni2Lan_sd[min_tissues$Lan[i],i]
  minimum_distance_sd[i,"Ofu"] <- Cni2Ofu_sd[min_tissues$Ofu[i],i]
  minimum_distance_sd[i,"Pec"] <- Cni2Pec_sd[min_tissues$Pec[i],i]
  minimum_distance_sd[i,"Aam"] <- Cni2Aam_sd[min_tissues$Aam[i],i]
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
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Pve'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Pve'),3:5]-min(minimum_distance_mean$Pve))/(max(minimum_distance_mean$Pve)-min(minimum_distance_mean$Pve))
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Pca'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Pca'),3:5]-min(minimum_distance_mean$Pca))/(max(minimum_distance_mean$Pca)-min(minimum_distance_mean$Pca))
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Npo'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Npo'),3:5]-min(minimum_distance_mean$Npo))/(max(minimum_distance_mean$Npo)-min(minimum_distance_mean$Npo))
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Lan'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Lan'),3:5]-min(minimum_distance_mean$Lan))/(max(minimum_distance_mean$Lan)-min(minimum_distance_mean$Lan))
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Ofu'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Ofu'),3:5]-min(minimum_distance_mean$Ofu))/(max(minimum_distance_mean$Ofu)-min(minimum_distance_mean$Ofu))
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Pec'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Pec'),3:5]-min(minimum_distance_mean$Pec))/(max(minimum_distance_mean$Pec)-min(minimum_distance_mean$Pec))
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Aam'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Aam'),3:5]-min(minimum_distance_mean$Aam))/(max(minimum_distance_mean$Aam)-min(minimum_distance_mean$Aam))

final_dataset <- minimum_distance_mean_tidy
final_dataset$species <- factor(final_dataset$species,
                                levels = c("Cgi","Pfu","Pye","Hru","Pca","Pve","Npo","Lan","Ofu","Pec","Aam"))

p = ggplot(final_dataset, aes(x=stage, y=JSD, colour=species)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2)
  geom_point(size=2,)
  facet_wrap(~species, nrow=11, ncol=1) +
  scale_x_continuous(labels = c(1:ncol(Cni2Cgi_mean)), breaks = c(1:ncol(Cni2Cgi_mean))) +
  theme_classic() +
  labs(x = "C. nippona tissues", y = "Normalised gene expression divergence (JSD)")
ggsave("JSD_plot.pdf", plot = p, width = 5, height = 40, device = "pdf")

write.table(min_tissues, "minimum_tissues.txt", sep='\t', quote = FALSE)
write.table(minimum_distance_mean, "minimum_tissues_distance_mean.txt", sep='\t', quote = FALSE)
write.table(minimum_distance_sd, "minimum_tissues_distance_sd.txt", sep='\t', quote = FALSE)
write.table(minimum_distance_mean_tidy, "minimum_tissues_to_plot.txt", sep='\t', quote = FALSE)
