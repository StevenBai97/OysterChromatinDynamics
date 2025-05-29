library(WGCNA)
library(tidyverse)
library(cowplot)
library(ggrepel)
library(ggthemes)
library(scales)

# Import phenotypic traits
datTraits <- read.delim("sample.info", row.names=1)

# Load and preprocess expression matrix
or_gene_exp = read.delim("salmon.gene.TMM.EXPR.matrix",
                         stringsAsFactors = F, check.names = F, row.names = 1)

gene_exp <- or_gene_exp %>%
  select(G1, G2, G3, T1, T2, T3, D1, D2, D3, LD1, LD2, LD3, AM1, AM2, AM3, Gi1, Gi2, Gi3, DG1, DG2, DG3, He1, He2, He3) %>%
  filter(rowSums(.[] >= 1) > 0)

datExpr = t(gene_exp)

# Data QC
gsg <- goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK

# Soft thresholding
enableWGCNAThreads(nThreads = 20)
sft = sft <- pickSoftThreshold(datExpr, powerVector = 1:20, verbose = 5)

fig_power1 <- ggplot(data = sft$fitIndices,
                     aes(x = Power,
                         y = SFT.R.sq)) +
  geom_point(color = 'red') +
  geom_text_repel(aes(label = Power)) +
  geom_hline(aes(yintercept = 0.85), color = 'red') +
  labs(title = 'Scale independence',
       x = 'Soft Threshold (power)',
       y = 'Scale Free Topology Model Fit,signed R^2') +
  theme_few() +
  theme(plot.title = element_text(hjust = 0.5))

fig_power2 <- ggplot(data = sft$fitIndices,
                     aes(x = Power,
                         y = mean.k.)) +
  geom_point(color = 'red') +
  geom_text_repel(aes(label = Power)) +
  labs(title = 'Mean connectivity',
       x = 'Soft Threshold (power)',
       y = 'Mean Connectivity') +
  theme_few()+
  theme(plot.title = element_text(hjust = 0.5))

pdf("power.pdf")
plot_grid(fig_power1, fig_power2)
dev.off()
power = sft$powerEstimate

# Network construction
net = blockwiseModules(datExpr, corType = "pearson", power = power,
                       TOMType = "signed",minModuleSize = 100,
                       reassignThreshold = 0, numericLabels = TRUE,
                       pamRespectsDendro = FALSE, saveTOMs = F,
                       mergeCutHeight = 0.25,
                       verbose = 3,nThreads = 0, maxBlockSize = 50000)
table(net$colors)

# Module visualization
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
write.csv(MEs,"modules_by_stage.csv")
write.csv(moduleLabels,"moduleLabels.csv")
write.csv(moduleColors,"moduleColors.csv")

pdf("plotDendroAndColors.pdf")
plotDendroAndColors(dendro = net$dendrograms[[1]], colors = net$colors,
                    groupLabels = "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

geneTree <- net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,  file = "step2.RData")
moduleTraitCor <- cor(net$MEs, datTraits,  use = "p", method = 'kendall')
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) <- dim(moduleTraitCor)
write.csv(moduleTraitPvalue,"moduleTraitPvalue.csv")
write.csv(moduleTraitCor,"moduleTraitCor.csv")

par(mar = c(6, 6, 3, 3))
pdf("labeledHeatmap.pdf")
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(net$MEs),
               ySymbols = names(net$MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(100),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

# Reconstruct Topological Overlap Matrix (TOM)
TOM <- TOMsimilarityFromExpr(datExpr, power = power)
# Select modules (in this case, the whole genome)
data <- read.csv('moduleColors.csv', header = TRUE)
moduleColors <- data[, 2]
colors <- read.csv(file = "moduleColors.csv")
colors <- unique(colors[2])
modules <- t(colors)
# Select module IDs
Gene_IDs <- rownames(gene_exp)
inModule <- is.finite(match(moduleColors, modules))
modProbes <- Gene_IDs[inModule]
# Select the corresponding Topological Overlap
modTOM <- TOM[inModule, inModule]
dimnames(modTOM) <- list(modProbes, modProbes)
write.table(modules, paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""))
# Export the network into edge and node list files Cytoscape can read
cyt <- exportNetworkToCytoscape(modTOM,
                                edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                weighted = TRUE,
                                threshold = 0.25,
                                nodeNames = modProbes,
                                nodeAttr = moduleColors[inModule])
