library(philentropy)
library(RColorBrewer)
library(ComplexHeatmap)
library(grid)

# 1. Spu vs. Lva
Spu2Lva_orth_tpm_qn <- read.csv("Spu2Lva_TPM_mean_quantile_transform.csv", header = T)
Lva2Spu_orth_tpm_qn <- read.csv("Lva2Spu_TPM_mean_quantile_transform.csv", header = T)

both_ID <- read.table("Spu2Lva.txt", header = T)
match_Spu <- Spu2Lva_orth_tpm_qn[match(both_ID$Spu,Spu2Lva_orth_tpm_qn$Gene_ID),]
match_Lva <- Lva2Spu_orth_tpm_qn[match(both_ID$Lva,Lva2Spu_orth_tpm_qn$Gene_ID),]
write.table(match_Spu,"order_Spu2Lva_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Lva,"order_Lva2Spu_tpm_qn.txt", sep="\t", quote=F)
match_Spu <- read.table("order_Spu2Lva_tpm_qn.txt", header = T)
match_Lva <- read.table("order_Lva2Spu_tpm_qn.txt", header = T)

Spu2Lva_full_JSD <- data.frame(replicate(ncol(match_Spu)-1,sample(0:1,ncol(match_Lva)-1,rep=TRUE)))
colnames(Spu2Lva_full_JSD) <- colnames(Spu2Lva_orth_tpm_qn)[-c(1)]
rownames(Spu2Lva_full_JSD) <- colnames(Lva2Spu_orth_tpm_qn)[-c(1)]
Spu2Lva_mean_JSD <- data.frame(replicate(ncol(match_Spu)-1,sample(0:1,ncol(match_Lva)-1,rep=TRUE)))
colnames(Spu2Lva_mean_JSD) <- colnames(Spu2Lva_orth_tpm_qn)[-c(1)]
rownames(Spu2Lva_mean_JSD) <- colnames(Lva2Spu_orth_tpm_qn)[-c(1)]
Spu2Lva_sd_JSD <- data.frame(replicate(ncol(match_Spu)-1,sample(0:1,ncol(match_Lva)-1,rep=TRUE)))
colnames(Spu2Lva_sd_JSD) <- colnames(Spu2Lva_orth_tpm_qn)[-c(1)]
rownames(Spu2Lva_sd_JSD) <- colnames(Lva2Spu_orth_tpm_qn)[-c(1)]

for (i in colnames(Spu2Lva_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Lva2Spu_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Spu[[i]],match_Lva[[j]])
    Spu2Lva_full_JSD[j,i] <- JSD(match)
    for (k in c(1:1000)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Spu2Lva_mean_JSD[j,i] <- mean(all_JSD)
    Spu2Lva_sd_JSD[j,i] <- sd(all_JSD)
  }
}

write.table(Spu2Lva_full_JSD, "Spu2Lva_full_set_JSD.txt", sep ='\t')
write.table(Spu2Lva_mean_JSD, "Spu2Lva_subset_JSD_mean.txt", sep ='\t')
write.table(Spu2Lva_sd_JSD, "Spu2Lva_subset_JSD_sd.txt", sep ='\t')

# 2. Spu vs. Aja
Spu2Aja_orth_tpm_qn <- read.csv("Spu2Aja_TPM_mean_quantile_transform.csv", header = T)
Aja2Spu_orth_tpm_qn <- read.csv("Aja2Spu_TPM_mean_quantile_transform.csv", header = T)

both_ID <- read.table("Spu2Aja.txt", header = T)
match_Spu <- Spu2Aja_orth_tpm_qn[match(both_ID$Spu,Spu2Aja_orth_tpm_qn$Gene_ID),]
match_Aja <- Aja2Spu_orth_tpm_qn[match(both_ID$Aja,Aja2Spu_orth_tpm_qn$Gene_ID),]
write.table(match_Spu,"order_Spu2Aja_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Aja,"order_Aja2Spu_tpm_qn.txt", sep="\t", quote=F)
match_Spu <- read.table("order_Spu2Aja_tpm_qn.txt", header = T)
match_Aja <- read.table("order_Aja2Spu_tpm_qn.txt", header = T)

Spu2Aja_full_JSD <- data.frame(replicate(ncol(match_Spu)-1,sample(0:1,ncol(match_Aja)-1,rep=TRUE)))
colnames(Spu2Aja_full_JSD) <- colnames(Spu2Aja_orth_tpm_qn)[-c(1)]
rownames(Spu2Aja_full_JSD) <- colnames(Aja2Spu_orth_tpm_qn)[-c(1)]
Spu2Aja_mean_JSD <- data.frame(replicate(ncol(match_Spu)-1,sample(0:1,ncol(match_Aja)-1,rep=TRUE)))
colnames(Spu2Aja_mean_JSD) <- colnames(Spu2Aja_orth_tpm_qn)[-c(1)]
rownames(Spu2Aja_mean_JSD) <- colnames(Aja2Spu_orth_tpm_qn)[-c(1)]
Spu2Aja_sd_JSD <- data.frame(replicate(ncol(match_Spu)-1,sample(0:1,ncol(match_Aja)-1,rep=TRUE)))
colnames(Spu2Aja_sd_JSD) <- colnames(Spu2Aja_orth_tpm_qn)[-c(1)]
rownames(Spu2Aja_sd_JSD) <- colnames(Aja2Spu_orth_tpm_qn)[-c(1)]

for (i in colnames(Spu2Aja_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Aja2Spu_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Spu[[i]],match_Aja[[j]])
    Spu2Aja_full_JSD[j,i] <- JSD(match)
    for (k in c(1:1000)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Spu2Aja_mean_JSD[j,i] <- mean(all_JSD)
    Spu2Aja_sd_JSD[j,i] <- sd(all_JSD)
  }
}

write.table(Spu2Aja_full_JSD, "Spu2Aja_full_set_JSD.txt", sep ='\t')
write.table(Spu2Aja_mean_JSD, "Spu2Aja_subset_JSD_mean.txt", sep ='\t')
write.table(Spu2Aja_sd_JSD, "Spu2Aja_subset_JSD_sd.txt", sep ='\t')
