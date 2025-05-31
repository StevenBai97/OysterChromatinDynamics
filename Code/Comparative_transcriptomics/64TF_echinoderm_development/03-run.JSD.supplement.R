library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(corrplot)
library(tidyr)
library(plyr)

Spu_orthologues <- c(62,60)
Spu2Lva_mean <- read.table("Spu2Lva_subset_JSD_mean.txt", header = T)
Spu2Aja_mean <- read.table("Spu2Aja_subset_JSD_mean.txt", header = T)

Spu2Lva_mean_divided <- Spu2Lva_mean/Spu_orthologues[1]
Spu2Aja_mean_divided <- Spu2Aja_mean/Spu_orthologues[2]

allvalues_divid <- rbind(Spu2Lva_mean_divided,Spu2Aja_mean_divided)

minimum_divided <- min(allvalues_divid)
maximum_divided <- max(allvalues_divid)

Spu2Lva_final_normalised <- (Spu2Lva_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Spu2Aja_final_normalised <- (Spu2Aja_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)

write.table(Spu2Lva_final_normalised, paste0("Spu2Lva_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Spu2Aja_final_normalised, paste0("Spu2Aja_final_normalised",".txt"), sep='\t', quote = FALSE)

paletteLength <- 100
heatmap_color <- colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(paletteLength)

myBreaks <- c(seq(1/paletteLength, 1, length.out=floor(paletteLength)))

h1 <- pheatmap(Spu2Lva_final_normalised,
               cluster_rows = F, cluster_cols = F, cellheight = 15, cellwidth = 15,
               border_color = NA, color = heatmap_color, breaks = myBreaks)

h2 <- pheatmap(Spu2Aja_final_normalised,
               cluster_rows = F, cluster_cols = F, cellheight = 15, cellwidth = 15,
               border_color = NA, color = heatmap_color, breaks = myBreaks)

pdf("Spu2Lva_normalised_mean_heatmap.pdf")
h1
dev.off()

pdf("Spu2Aja_normalised_mean_heatmap.pdf")
h2
dev.off()

Spu2Lva_sd <- read.table("Spu2Lva_subset_JSD_sd.txt", header = T)
Spu2Aja_sd <- read.table("Spu2Aja_subset_JSD_sd.txt", header = T)

Spu2Lva_transposed <- data.frame(t(Spu2Lva_mean))
Spu2Aja_transposed <- data.frame(t(Spu2Aja_mean))

min_stages <- data.frame(replicate(2,sample(0:1,ncol(Spu2Lva_mean),rep=TRUE)))
colnames(min_stages) <- c("Lva","Aja")
rownames(min_stages) <- c("cleavage","hatched_blastula", "mesenchyme_blastula", "early_gastrula","mid_gastrula", "late_gastrula",
                           "prism", "late_prism","pluteus","four_arm_stage","vestibular_invagination_stage","pentagonal_disc_stage",
                          "tube_foot_protrusion_stage", "post_metamorphosis","juvenile","adult_spine")

min_stages$Lva <- names(Spu2Lva_transposed)[apply(Spu2Lva_transposed, MARGIN = 1, FUN = which.min)]
min_stages$Aja <- names(Spu2Aja_transposed)[apply(Spu2Aja_transposed, MARGIN = 1, FUN = which.min)]

minimum_distance_mean <- data.frame(replicate(2,sample(0:1,ncol(Spu2Lva_mean),rep=TRUE)))
minimum_distance_sd <- data.frame(replicate(2,sample(0:1,ncol(Spu2Lva_mean),rep=TRUE)))
colnames(minimum_distance_mean) <- c("Lva","Aja")
rownames(minimum_distance_mean) <- c("cleavage","hatched_blastula", "mesenchyme_blastula", "early_gastrula","mid_gastrula", "late_gastrula",
                                     "prism", "late_prism","pluteus","four_arm_stage","vestibular_invagination_stage","pentagonal_disc_stage",
                                     "tube_foot_protrusion_stage", "post_metamorphosis","juvenile","adult_spine")
colnames(minimum_distance_sd) <- c("Lva","Aja")
rownames(minimum_distance_sd) <- c("cleavage","hatched_blastula", "mesenchyme_blastula", "early_gastrula","mid_gastrula", "late_gastrula",
                                   "prism", "late_prism","pluteus","four_arm_stage","vestibular_invagination_stage","pentagonal_disc_stage",
                                   "tube_foot_protrusion_stage", "post_metamorphosis","juvenile","adult_spine")

for (i in c(1:ncol(Spu2Lva_mean))){
  minimum_distance_mean[i,"Lva"] <- Spu2Lva_mean[min_stages$Lva[i],i]
  minimum_distance_mean[i,"Aja"] <- Spu2Aja_mean[min_stages$Aja[i],i]
  minimum_distance_sd[i,"Lva"] <- Spu2Lva_sd[min_stages$Lva[i],i]
  minimum_distance_sd[i,"Aja"] <- Spu2Aja_sd[min_stages$Aja[i],i]
}

minimum_distance_mean$stage <- c(1:ncol(Spu2Lva_mean))
minimum_distance_sd$stage <- c(1:ncol(Spu2Lva_mean))
minimum_distance_mean_tidy <- gather(minimum_distance_mean, "species", "JSD", -stage)
minimum_distance_sd_tidy <- gather(minimum_distance_sd, "species", "JSD", -stage)
minimum_distance_mean_tidy$upper <- minimum_distance_mean_tidy["JSD"] + minimum_distance_sd_tidy["JSD"]
minimum_distance_mean_tidy$lower <- minimum_distance_mean_tidy["JSD"] - minimum_distance_sd_tidy["JSD"]
colnames(minimum_distance_mean_tidy)[4:5] <- c("upper","lower")
minimum_distance_mean_tidy$upper <- unlist(minimum_distance_mean_tidy$upper)
minimum_distance_mean_tidy$lower <- unlist(minimum_distance_mean_tidy$lower)

minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Lva'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Lva'),3:5]-min(minimum_distance_mean$Lva))/(max(minimum_distance_mean$Lva)-min(minimum_distance_mean$Lva))
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Aja'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Aja'),3:5]-min(minimum_distance_mean$Aja))/(max(minimum_distance_mean$Aja)-min(minimum_distance_mean$Aja))

final_dataset <- minimum_distance_mean_tidy
final_dataset$species <- factor(final_dataset$species,
                                levels = c("Lva","Aja"))

p = ggplot(final_dataset, aes(x=stage, y=JSD, colour=species)) +
  geom_line(show.legend = FALSE) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=factor(species)),
              colour = NA, show.legend = FALSE,
              alpha = 0.5) +
  facet_wrap(~species, nrow=1, ncol=10) +
  scale_x_continuous(labels = c(1:ncol(Spu2Lva_mean)), breaks = c(1:ncol(Spu2Lva_mean))) +
  theme_classic() +
  labs(y = "Normalised gene expression divergence (JSD)")
ggsave("JSD_plot.pdf", plot = p, width = 10, height = 4, device = "pdf")
