library(tidyverse)
library(datawizard)
library(Mfuzz)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(scales)
library(ComplexHeatmap)

# Extract expressed genes (TPM > 1) and calculate their average expression
or_gene_exp = read.delim("Data/Figure1/salmon.gene.TMM.EXPR.matrix", 
                         stringsAsFactors = F, check.names = F, row.names = 1)
rows_to_remove <- rowSums(or_gene_exp[, -1] < 1) == (ncol(or_gene_exp) - 1)
gene_exp <- or_gene_exp[!rows_to_remove, ]

D <- rowMeans(gene_exp[,c(1,2,3)])
G <- rowMeans(gene_exp[,c(4,5,6)])
J <- rowMeans(gene_exp[,c(7,8,9)])
M <- rowMeans(gene_exp[,c(10,11,12)])
T <- rowMeans(gene_exp[,c(13,14,15)])
LD <- rowMeans(gene_exp[,c(16,17,18)])
df_average <- data.frame(G,T,D,LD,J,M)
row.names(df_average) <- row.names(gene_exp)
write.table(df_average,"gene_exp_average.txt",sep = '\t',quote = F,col.names=NA)

# Mfuzz clustering
m_df_average <- as.matrix(df_average)
eset_rna <- new('ExpressionSet', exprs = m_df_average)
clean_eset_rna <- filter.std(eset_rna,min.std=0) 
norm_eset_rna <- standardise(clean_eset_rna)
fuzzifier <- round(mestimate(norm_eset_rna),2)

pdf("crange.pdf")
crange <- Dmin(norm_eset_rna,m=fuzzifier,crange=seq(4,20,1),repeats=3,visu=TRUE)
dev.off()
which(crange==min(crange))

set.seed(123)
clusters <- mfuzz(norm_eset_rna, c = 10, m = fuzzifier)
clusters$size

pdf("mfuzz.pdf")
mfuzz.plot2(norm_eset_rna,cl=clusters,mfrow=c(4,4),time.labels=c("G","T","D","LD","J","M"), x11 = FALSE,colo="fancy",centre.col= "blue",Xwidth = 10, Xheight = 10)
dev.off()

df_names <- df_average
df_names_notnull <- df_names[apply(df_names, 1, function(x) sd(x, na.rm = TRUE)!=0),]
clusters_annotated_unordered <- cbind(rownames(df_names_notnull),clusters$cluster)
colnames(clusters_annotated_unordered) <- c("Gene_ID","Cluster_raw")
write.table(clusters_annotated_unordered,"Cni_clusters_annotation_raw.txt",quote = F, sep = '\t', row.names = FALSE)

clusters_expression_unordered <- cbind(rownames(df_names_notnull),df_names_notnull,clusters$cluster)
colnames(clusters_expression_unordered) <- c("Gene_ID","G","T","D","LD","J","M","Cluster_raw")
write.table(clusters_expression_unordered,"Cni_clusters_expression_raw.txt",quote = F, sep = '\t', row.names = FALSE)

df_raw <- read.table("Cni_clusters_expression_raw.txt", header = TRUE)
df_zscore <- data.frame(cbind(t(scale(t(df_raw[-c(1,8)]))),'cluster_raw'=df_raw$Cluster_raw))
df_zscore$cluster_raw <- factor(df_zscore$cluster_raw, 
                                levels = c("7","10","3","5","8","9","4","6","2","1"))
df_sorted <- df_zscore[order(df_zscore$cluster_raw),]

u_color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)

pdf("heatmap.pdf")
ComplexHeatmap::Heatmap(as.matrix(df_sorted[-c(7)]),
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        show_row_names = FALSE,
                        width = unit(4, "cm"), 
                        height = unit(10, "cm"),
                        col = u_color)
dev.off()

# Re-write clusters
df_reordered <- df_raw
df_reordered$Cluster_corrected <- c(rep(0,nrow(df_reordered)))
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 7] <- 1
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 10] <- 2
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 3] <- 3
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 5] <- 4
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 8] <- 5
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 9] <- 6
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 4] <- 7
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 6] <- 8
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 2] <- 9
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 1] <- 10
df_reordered_clean <- df_reordered[-c(8)]
write.table(df_reordered_clean[c(1,8)],"Cni_clusters_annotation_corrected.txt", 
            quote = F, sep = '\t', row.names = FALSE)
write.table(df_reordered_clean,"Cni_clusters_expression_corrected.txt", 
            quote = F, sep = '\t', row.names = FALSE)

# LOESS smoothing
palette <- viridis_pal(option = "D", direction = 1)(30)
palette2 <- gsub('.{2}$', '', palette)
write.csv(rbind(palette2), file = 'palette.txt', row.names = FALSE)

palette <- viridis_pal(option = "D", direction = 1)(10)
time <- c(1:6)

write.table(c(palette), "palette.txt", sep = ',', quote = TRUE)

df_raw = df_reordered_clean

## For cluster 1
df_cluster1_raw <- df_raw[df_raw$Cluster_corrected == 1,-c(1,8)]
mat_cluster1_raw <- as.matrix(df_cluster1_raw) 
mat_cluster1_norm <- t(scale(t(mat_cluster1_raw))) 
df_cluster1_norm <- data.frame(mat_cluster1_norm)
df_cluster1_clean <- data.frame(cbind(time,t(df_cluster1_norm))) 
df_cluster1_clean_long <- gather(df_cluster1_clean, gene_ID, expression, -time)

summary_cluster1 <- data.frame(time=df_cluster1_clean$time, 
                               n=tapply(df_cluster1_clean_long$expression, 
                                        df_cluster1_clean_long$time, length), 
                               mean=tapply(df_cluster1_clean_long$expression, 
                                           df_cluster1_clean_long$time, mean))

### Plot data
cluster1_plot <- ggplot(summary_cluster1, aes(x=time, y=mean)) +
  labs(y = "Normalised gene expression (z-score)", x = "developmental stage") +
  geom_line(data = df_cluster1_clean_long, aes(x=time, y=expression, group=gene_ID), 
            color="gray") +
  stat_smooth(colour = palette[1],
              fill = palette[1],
              size = 2, alpha = 0.5) +  
  geom_hline(yintercept = 0, linetype="dashed", color = "black") + 
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  scale_y_continuous(breaks = c(-3:3), limits = c(-5,5)) +
  scale_x_continuous(breaks = c(1:6), expand = c(0,0), limits = c(0.8,6.4),
                     labels = c("G", "T","D", "LD", 
                                "J", "M")) +
  theme(plot.title = element_text(color="black",size=8,hjust=0.5),
        axis.title.x = element_text(color="black",size=8),
        axis.title.y = element_text(color="black",size=8),
        axis.text.x = element_text(color="black",size=6,angle=90,
                                   hjust=0.95,vjust=0.2,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=6,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", linewidth=0.25, linetype = "solid"),
        axis.ticks = element_line(colour="black", linewidth=0.25, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(0.15, "cm"))

## For cluster 2
df_cluster2_raw <- df_raw[df_raw$Cluster_corrected == 2,-c(1,8)]
mat_cluster2_raw <- as.matrix(df_cluster2_raw) 
mat_cluster2_norm <- t(scale(t(mat_cluster2_raw))) 
df_cluster2_norm <- data.frame(mat_cluster2_norm)
df_cluster2_clean <- data.frame(cbind(time,t(df_cluster2_norm))) 
df_cluster2_clean_long <- gather(df_cluster2_clean, gene_ID, expression, -time)

summary_cluster2 <- data.frame(time=df_cluster2_clean$time, 
                               n=tapply(df_cluster2_clean_long$expression, 
                                        df_cluster2_clean_long$time, length), 
                               mean=tapply(df_cluster2_clean_long$expression, 
                                           df_cluster2_clean_long$time, mean))

### Plot data
cluster2_plot <- ggplot(summary_cluster2, aes(x=time, y=mean)) +
  geom_line(data = df_cluster2_clean_long, aes(x=time, y=expression, group=gene_ID), 
            color="gray") +
  stat_smooth(colour = palette[2],
              fill = palette[2],
              size = 2, alpha = 0.5) +  
  geom_hline(yintercept = 0, linetype="dashed", color = "black") + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(y = "Normalised gene expression (z-score)", x = "developmental stage") +
  scale_y_continuous(breaks = c(-3:3), limits = c(-5,5)) +
  scale_x_continuous(breaks = c(1:6), expand = c(0,0), limits = c(0.8,6.4),
                     labels = c("G", "T","D", "LD", 
                                "J", "M")) +
  theme(plot.title = element_text(color="black",size=8,hjust=0.5),
        axis.title.x = element_text(color="black",size=8),
        axis.title.y = element_text(color="black",size=8),
        axis.text.x = element_text(color="black",size=6,angle=90,
                                   hjust=0.95,vjust=0.2,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=6,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", linewidth=0.25, linetype = "solid"),
        axis.ticks = element_line(colour="black", linewidth=0.25, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(0.15, "cm"))

## For cluster 3
df_cluster3_raw <- df_raw[df_raw$Cluster_corrected == 3,-c(1,8)]
mat_cluster3_raw <- as.matrix(df_cluster3_raw) 
mat_cluster3_norm <- t(scale(t(mat_cluster3_raw))) 
df_cluster3_norm <- data.frame(mat_cluster3_norm)
df_cluster3_clean <- data.frame(cbind(time,t(df_cluster3_norm))) 
df_cluster3_clean_long <- gather(df_cluster3_clean, gene_ID, expression, -time)

summary_cluster3 <- data.frame(time=df_cluster3_clean$time, 
                               n=tapply(df_cluster3_clean_long$expression, 
                                        df_cluster3_clean_long$time, length), 
                               mean=tapply(df_cluster3_clean_long$expression, 
                                           df_cluster3_clean_long$time, mean))

### Plot data
cluster3_plot <- ggplot(summary_cluster3, aes(x=time, y=mean)) +
  geom_line(data = df_cluster3_clean_long, aes(x=time, y=expression, group=gene_ID), 
            color="gray") +
  stat_smooth(colour = palette[3],
              fill = palette[3],
              size = 2, alpha = 0.5) +  
  geom_hline(yintercept = 0, linetype="dashed", color = "black") + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(y = "Normalised gene expression (z-score)", x = "developmental stage") +
  scale_y_continuous(breaks = c(-3:3), limits = c(-5,5)) +
  scale_x_continuous(breaks = c(1:6), expand = c(0,0), limits = c(0.8,6.4),
                     labels = c("G", "T","D", "LD", 
                                "J", "M")) +
  theme(plot.title = element_text(color="black",size=8,hjust=0.5),
        axis.title.x = element_text(color="black",size=8),
        axis.title.y = element_text(color="black",size=8),
        axis.text.x = element_text(color="black",size=6,angle=90,
                                   hjust=0.95,vjust=0.2,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=6,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", linewidth=0.25, linetype = "solid"),
        axis.ticks = element_line(colour="black", linewidth=0.25, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(0.15, "cm"))

## For cluster 4
df_cluster4_raw <- df_raw[df_raw$Cluster_corrected ==4,-c(1,8)]
mat_cluster4_raw <- as.matrix(df_cluster4_raw) 
mat_cluster4_norm <- t(scale(t(mat_cluster4_raw))) 
df_cluster4_norm <- data.frame(mat_cluster4_norm)
df_cluster4_clean <- data.frame(cbind(time,t(df_cluster4_norm))) 
df_cluster4_clean_long <- gather(df_cluster4_clean, gene_ID, expression, -time)

summary_cluster4 <- data.frame(time=df_cluster4_clean$time, 
                               n=tapply(df_cluster4_clean_long$expression, 
                                        df_cluster4_clean_long$time, length), 
                               mean=tapply(df_cluster4_clean_long$expression, 
                                           df_cluster4_clean_long$time, mean))

### Plot data
cluster4_plot <- ggplot(summary_cluster4, aes(x=time, y=mean)) +
  geom_line(data = df_cluster4_clean_long, aes(x=time, y=expression, group=gene_ID), 
            color="gray") +
  stat_smooth(colour = palette[4],
              fill = palette[4],
              size = 2, alpha = 0.5) +  
  geom_hline(yintercept = 0, linetype="dashed", color = "black") + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(y = "Normalised gene expression (z-score)", x = "developmental stage") +
  scale_y_continuous(breaks = c(-3:3), limits = c(-5,5)) +
  scale_x_continuous(breaks = c(1:6), expand = c(0,0), limits = c(0.8,6.4),
                     labels = c("G", "T","D", "LD", 
                                "J", "M")) +
  theme(plot.title = element_text(color="black",size=8,hjust=0.5),
        axis.title.x = element_text(color="black",size=8),
        axis.title.y = element_text(color="black",size=8),
        axis.text.x = element_text(color="black",size=6,angle=90,
                                   hjust=0.95,vjust=0.2,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=6,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", linewidth=0.25, linetype = "solid"),
        axis.ticks = element_line(colour="black", linewidth=0.25, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(0.15, "cm"))

## For cluster 5
df_cluster5_raw <- df_raw[df_raw$Cluster_corrected ==5,-c(1,8)]
mat_cluster5_raw <- as.matrix(df_cluster5_raw) 
mat_cluster5_norm <- t(scale(t(mat_cluster5_raw))) 
df_cluster5_norm <- data.frame(mat_cluster5_norm)
df_cluster5_clean <- data.frame(cbind(time,t(df_cluster5_norm))) 
df_cluster5_clean_long <- gather(df_cluster5_clean, gene_ID, expression, -time)

summary_cluster5 <- data.frame(time=df_cluster5_clean$time, 
                               n=tapply(df_cluster5_clean_long$expression, 
                                        df_cluster5_clean_long$time, length), 
                               mean=tapply(df_cluster5_clean_long$expression, 
                                           df_cluster5_clean_long$time, mean))

### Plot data
cluster5_plot <- ggplot(summary_cluster5, aes(x=time, y=mean)) +
  geom_line(data = df_cluster5_clean_long, aes(x=time, y=expression, group=gene_ID), 
            color="gray") +
  stat_smooth(colour = palette[5],
              fill = palette[5],
              size = 2, alpha = 0.5) +  
  geom_hline(yintercept = 0, linetype="dashed", color = "black") + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(y = "Normalised gene expression (z-score)", x = "developmental stage") +
  scale_y_continuous(breaks = c(-3:3), limits = c(-5,5)) +
  scale_x_continuous(breaks = c(1:6), expand = c(0,0), limits = c(0.8,6.4),
                     labels = c("G", "T","D", "LD", 
                                "J", "M")) +
  theme(plot.title = element_text(color="black",size=8,hjust=0.5),
        axis.title.x = element_text(color="black",size=8),
        axis.title.y = element_text(color="black",size=8),
        axis.text.x = element_text(color="black",size=6,angle=90,
                                   hjust=0.95,vjust=0.2,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=6,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", linewidth=0.25, linetype = "solid"),
        axis.ticks = element_line(colour="black", linewidth=0.25, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(0.15, "cm"))

## For cluster 6
df_cluster6_raw <- df_raw[df_raw$Cluster_corrected ==6,-c(1,8)]
mat_cluster6_raw <- as.matrix(df_cluster6_raw) 
mat_cluster6_norm <- t(scale(t(mat_cluster6_raw))) 
df_cluster6_norm <- data.frame(mat_cluster6_norm)
df_cluster6_clean <- data.frame(cbind(time,t(df_cluster6_norm))) 
df_cluster6_clean_long <- gather(df_cluster6_clean, gene_ID, expression, -time)

summary_cluster6 <- data.frame(time=df_cluster6_clean$time, 
                               n=tapply(df_cluster6_clean_long$expression, 
                                        df_cluster6_clean_long$time, length), 
                               mean=tapply(df_cluster6_clean_long$expression, 
                                           df_cluster6_clean_long$time, mean))

### Plot data
cluster6_plot <- ggplot(summary_cluster6, aes(x=time, y=mean)) +
  geom_line(data = df_cluster6_clean_long, aes(x=time, y=expression, group=gene_ID), 
            color="gray") +
  stat_smooth(colour = palette[6],
              fill = palette[6],
              size = 2, alpha = 0.5) +  
  geom_hline(yintercept = 0, linetype="dashed", color = "black") + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(y = "Normalised gene expression (z-score)", x = "developmental stage") +
  scale_y_continuous(breaks = c(-3:3), limits = c(-5,5)) +
  scale_x_continuous(breaks = c(1:6), expand = c(0,0), limits = c(0.8,6.4),
                     labels = c("G", "T","D", "LD", 
                                "J", "M")) +
  theme(plot.title = element_text(color="black",size=8,hjust=0.5),
        axis.title.x = element_text(color="black",size=8),
        axis.title.y = element_text(color="black",size=8),
        axis.text.x = element_text(color="black",size=6,angle=90,
                                   hjust=0.95,vjust=0.2,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=6,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", linewidth=0.25, linetype = "solid"),
        axis.ticks = element_line(colour="black", linewidth=0.25, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(0.15, "cm"))

## For cluster 7
df_cluster7_raw <- df_raw[df_raw$Cluster_corrected ==7,-c(1,8)]
mat_cluster7_raw <- as.matrix(df_cluster7_raw) 
mat_cluster7_norm <- t(scale(t(mat_cluster7_raw))) 
df_cluster7_norm <- data.frame(mat_cluster7_norm)
df_cluster7_clean <- data.frame(cbind(time,t(df_cluster7_norm))) 
df_cluster7_clean_long <- gather(df_cluster7_clean, gene_ID, expression, -time)

summary_cluster7 <- data.frame(time=df_cluster7_clean$time, 
                               n=tapply(df_cluster7_clean_long$expression, 
                                        df_cluster7_clean_long$time, length), 
                               mean=tapply(df_cluster7_clean_long$expression, 
                                           df_cluster7_clean_long$time, mean))

### Plot data
cluster7_plot <- ggplot(summary_cluster7, aes(x=time, y=mean)) +
  geom_line(data = df_cluster7_clean_long, aes(x=time, y=expression, group=gene_ID), 
            color="gray") +
  stat_smooth(colour = palette[7],
              fill = palette[7],
              size = 2, alpha = 0.5) +  
  geom_hline(yintercept = 0, linetype="dashed", color = "black") + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(y = "Normalised gene expression (z-score)", x = "developmental stage") +
  scale_y_continuous(breaks = c(-3:3), limits = c(-5,5)) +
  scale_x_continuous(breaks = c(1:6), expand = c(0,0), limits = c(0.8,6.4),
                     labels = c("G", "T","D", "LD", 
                                "J", "M")) +
  theme(plot.title = element_text(color="black",size=8,hjust=0.5),
        axis.title.x = element_text(color="black",size=8),
        axis.title.y = element_text(color="black",size=8),
        axis.text.x = element_text(color="black",size=6,angle=90,
                                   hjust=0.95,vjust=0.2,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=6,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", linewidth=0.25, linetype = "solid"),
        axis.ticks = element_line(colour="black", linewidth=0.25, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(0.15, "cm"))

## For cluster 8
df_cluster8_raw <- df_raw[df_raw$Cluster_corrected ==8,-c(1,8)]
mat_cluster8_raw <- as.matrix(df_cluster8_raw) 
mat_cluster8_norm <- t(scale(t(mat_cluster8_raw))) 
df_cluster8_norm <- data.frame(mat_cluster8_norm)
df_cluster8_clean <- data.frame(cbind(time,t(df_cluster8_norm))) 
df_cluster8_clean_long <- gather(df_cluster8_clean, gene_ID, expression, -time)

summary_cluster8 <- data.frame(time=df_cluster8_clean$time, 
                               n=tapply(df_cluster8_clean_long$expression, 
                                        df_cluster8_clean_long$time, length), 
                               mean=tapply(df_cluster8_clean_long$expression, 
                                           df_cluster8_clean_long$time, mean))

### Plot data
cluster8_plot <- ggplot(summary_cluster8, aes(x=time, y=mean)) +
  geom_line(data = df_cluster8_clean_long, aes(x=time, y=expression, group=gene_ID), 
            color="gray") +
  stat_smooth(colour = palette[8],
              fill = palette[8],
              size = 2, alpha = 0.5) +  
  geom_hline(yintercept = 0, linetype="dashed", color = "black") + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(y = "Normalised gene expression (z-score)", x = "developmental stage") +
  scale_y_continuous(breaks = c(-3:3), limits = c(-5,5)) +
  scale_x_continuous(breaks = c(1:6), expand = c(0,0), limits = c(0.8,6.4),
                     labels = c("G", "T","D", "LD", 
                                "J", "M")) +
  theme(plot.title = element_text(color="black",size=8,hjust=0.5),
        axis.title.x = element_text(color="black",size=8),
        axis.title.y = element_text(color="black",size=8),
        axis.text.x = element_text(color="black",size=6,angle=90,
                                   hjust=0.95,vjust=0.2,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=6,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", linewidth=0.25, linetype = "solid"),
        axis.ticks = element_line(colour="black", linewidth=0.25, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(0.15, "cm"))

## For cluster 9
df_cluster9_raw <- df_raw[df_raw$Cluster_corrected ==9,-c(1,8)]
mat_cluster9_raw <- as.matrix(df_cluster9_raw) 
mat_cluster9_norm <- t(scale(t(mat_cluster9_raw))) 
df_cluster9_norm <- data.frame(mat_cluster9_norm)
df_cluster9_clean <- data.frame(cbind(time,t(df_cluster9_norm))) 
df_cluster9_clean_long <- gather(df_cluster9_clean, gene_ID, expression, -time)

summary_cluster9 <- data.frame(time=df_cluster9_clean$time, 
                               n=tapply(df_cluster9_clean_long$expression, 
                                        df_cluster9_clean_long$time, length), 
                               mean=tapply(df_cluster9_clean_long$expression, 
                                           df_cluster9_clean_long$time, mean))

### Plot data
cluster9_plot <- ggplot(summary_cluster9, aes(x=time, y=mean)) +
  geom_line(data = df_cluster9_clean_long, aes(x=time, y=expression, group=gene_ID), 
            color="gray") +
  stat_smooth(colour = palette[9],
              fill = palette[9],
              size = 2, alpha = 0.5) +  
  geom_hline(yintercept = 0, linetype="dashed", color = "black") + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(y = "Normalised gene expression (z-score)", x = "developmental stage") +
  scale_y_continuous(breaks = c(-3:3), limits = c(-5,5)) +
  scale_x_continuous(breaks = c(1:6), expand = c(0,0), limits = c(0.8,6.4),
                     labels = c("G", "T","D", "LD", 
                                "J", "M")) +
  theme(plot.title = element_text(color="black",size=8,hjust=0.5),
        axis.title.x = element_text(color="black",size=8),
        axis.title.y = element_text(color="black",size=8),
        axis.text.x = element_text(color="black",size=6,angle=90,
                                   hjust=0.95,vjust=0.2,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=6,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", linewidth=0.25, linetype = "solid"),
        axis.ticks = element_line(colour="black", linewidth=0.25, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(0.15, "cm"))

## For cluster 10
df_cluster10_raw <- df_raw[df_raw$Cluster_corrected ==10,-c(1,8)]
mat_cluster10_raw <- as.matrix(df_cluster10_raw) 
mat_cluster10_norm <- t(scale(t(mat_cluster10_raw))) 
df_cluster10_norm <- data.frame(mat_cluster10_norm)
df_cluster10_clean <- data.frame(cbind(time,t(df_cluster10_norm))) 
df_cluster10_clean_long <- gather(df_cluster10_clean, gene_ID, expression, -time)

summary_cluster10 <- data.frame(time=df_cluster10_clean$time, 
                                n=tapply(df_cluster10_clean_long$expression, 
                                         df_cluster10_clean_long$time, length), 
                                mean=tapply(df_cluster10_clean_long$expression, 
                                            df_cluster10_clean_long$time, mean))

### Plot data
cluster10_plot <- ggplot(summary_cluster10, aes(x=time, y=mean)) +
  geom_line(data = df_cluster10_clean_long, aes(x=time, y=expression, group=gene_ID), 
            color="gray") +
  stat_smooth(colour = palette[10],
              fill = palette[10],
              size = 2, alpha = 0.5) +  
  geom_hline(yintercept = 0, linetype="dashed", color = "black") + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(y = "Normalised gene expression (z-score)", x = "developmental stage") +
  scale_y_continuous(breaks = c(-3:3), limits = c(-5,5)) +
  scale_x_continuous(breaks = c(1:6), expand = c(0,0), limits = c(0.8,6.4),
                     labels = c("G", "T","D", "LD", 
                                "J", "M")) +
  theme(plot.title = element_text(color="black",size=8,hjust=0.5),
        axis.title.x = element_text(color="black",size=8),
        axis.title.y = element_text(color="black",size=8),
        axis.text.x = element_text(color="black",size=6,angle=90,
                                   hjust=0.95,vjust=0.2,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=6,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", linewidth=0.25, linetype = "solid"),
        axis.ticks = element_line(colour="black", linewidth=0.25, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(0.15, "cm"))

allfigures <- ggarrange(cluster1_plot, cluster2_plot, cluster3_plot, cluster4_plot,
                        cluster5_plot, cluster6_plot, cluster7_plot, cluster8_plot,
                        cluster9_plot, cluster10_plot, ncol = 1, nrow = 10)
ggsave("Cni_RNAseq_clusters_loess(1X10).pdf", width = 10, height = 70, units = "cm")

df_clusters = df_raw
cluster_1 <- df_clusters[df_clusters$Cluster_corrected == 1,]
cluster_2 <- df_clusters[df_clusters$Cluster_corrected == 2,]
cluster_3 <- df_clusters[df_clusters$Cluster_corrected == 3,]
cluster_4 <- df_clusters[df_clusters$Cluster_corrected == 4,]
cluster_5 <- df_clusters[df_clusters$Cluster_corrected == 5,]
cluster_6 <- df_clusters[df_clusters$Cluster_corrected == 6,]
cluster_7 <- df_clusters[df_clusters$Cluster_corrected == 7,]
cluster_8 <- df_clusters[df_clusters$Cluster_corrected == 8,]
cluster_9 <- df_clusters[df_clusters$Cluster_corrected == 9,]
cluster_10 <- df_clusters[df_clusters$Cluster_corrected == 10,]

writeLines(as.character(cluster_1[, 1]), "cluster1.txt")
writeLines(as.character(cluster_2[, 1]), "cluster2.txt")
writeLines(as.character(cluster_3[, 1]), "cluster3.txt")
writeLines(as.character(cluster_4[, 1]), "cluster4.txt")
writeLines(as.character(cluster_5[, 1]), "cluster5.txt")
writeLines(as.character(cluster_6[, 1]), "cluster6.txt")
writeLines(as.character(cluster_7[, 1]), "cluster7.txt")
writeLines(as.character(cluster_8[, 1]), "cluster8.txt")
writeLines(as.character(cluster_9[, 1]), "cluster9.txt")
writeLines(as.character(cluster_10[, 1]), "cluster10.txt")
