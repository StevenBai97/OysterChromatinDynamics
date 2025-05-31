library(philentropy)
library(RColorBrewer)
library(ComplexHeatmap)
library(grid)

# 1. Cni vs. Cgi

Cni2Cgi_orth_tpm_qn <- read.csv("Cni2Cgi_TPM_mean_quantile_transform.csv", header = T)
Cgi2Cni_orth_tpm_qn <- read.csv("Cgi2Cni_TPM_mean_quantile_transform.csv", header = T)

both_ID <- read.table("Cni2Cgi.txt", header = T)
match_Cni <- Cni2Cgi_orth_tpm_qn[match(both_ID$Cni,Cni2Cgi_orth_tpm_qn$Gene_ID),]
match_Cgi <- Cgi2Cni_orth_tpm_qn[match(both_ID$Cgi,Cgi2Cni_orth_tpm_qn$Gene_ID),]

write.table(match_Cni,"order_Cni2Cgi_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Cgi,"order_Cgi2Cni_tpm_qn.txt", sep="\t", quote=F)
match_Cni <- read.table("order_Cni2Cgi_tpm_qn.txt", header = T)
match_Cgi <- read.table("order_Cgi2Cni_tpm_qn.txt", header = T)

Cni2Cgi_full_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Cgi)-1,rep=TRUE)))
colnames(Cni2Cgi_full_JSD) <- colnames(Cni2Cgi_orth_tpm_qn)[-c(1)]
rownames(Cni2Cgi_full_JSD) <- colnames(Cgi2Cni_orth_tpm_qn)[-c(1)]
Cni2Cgi_mean_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Cgi)-1,rep=TRUE)))
colnames(Cni2Cgi_mean_JSD) <- colnames(Cni2Cgi_orth_tpm_qn)[-c(1)]
rownames(Cni2Cgi_mean_JSD) <- colnames(Cgi2Cni_orth_tpm_qn)[-c(1)]
Cni2Cgi_sd_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Cgi)-1,rep=TRUE)))
colnames(Cni2Cgi_sd_JSD) <- colnames(Cni2Cgi_orth_tpm_qn)[-c(1)]
rownames(Cni2Cgi_sd_JSD) <- colnames(Cgi2Cni_orth_tpm_qn)[-c(1)]

for (i in colnames(Cni2Cgi_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Cgi2Cni_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Cni[[i]],match_Cgi[[j]])
    Cni2Cgi_full_JSD[j,i] <- JSD(match)
    for (k in c(1:1000)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Cni2Cgi_mean_JSD[j,i] <- mean(all_JSD)
    Cni2Cgi_sd_JSD[j,i] <- sd(all_JSD)
  }
}

write.table(Cni2Cgi_full_JSD, "Cni2Cgi_full_set_JSD.txt", sep ='\t')
write.table(Cni2Cgi_mean_JSD, "Cni2Cgi_subset_JSD_mean.txt", sep ='\t')
write.table(Cni2Cgi_sd_JSD, "Cni2Cgi_subset_JSD_sd.txt", sep ='\t')


# 2. Cni vs. Pfu

Cni2Pfu_orth_tpm_qn <- read.csv("Cni2Pfu_TPM_mean_quantile_transform.csv", header = T)
Pfu2Cni_orth_tpm_qn <- read.csv("Pfu2Cni_TPM_mean_quantile_transform.csv", header = T)

both_ID <- read.table("Cni2Pfu.txt", header = T)
match_Cni <- Cni2Pfu_orth_tpm_qn[match(both_ID$Cni,Cni2Pfu_orth_tpm_qn$Gene_ID),]
match_Pfu <- Pfu2Cni_orth_tpm_qn[match(both_ID$Pfu,Pfu2Cni_orth_tpm_qn$Gene_ID),]
write.table(match_Cni,"order_Cni2Pfu_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Pfu,"order_Pfu2Cni_tpm_qn.txt", sep="\t", quote=F)
match_Cni <- read.table("order_Cni2Pfu_tpm_qn.txt", header = T)
match_Pfu <- read.table("order_Pfu2Cni_tpm_qn.txt", header = T)

Cni2Pfu_full_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Pfu)-1,rep=TRUE)))
colnames(Cni2Pfu_full_JSD) <- colnames(Cni2Pfu_orth_tpm_qn)[-c(1)]
rownames(Cni2Pfu_full_JSD) <- colnames(Pfu2Cni_orth_tpm_qn)[-c(1)]
Cni2Pfu_mean_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Pfu)-1,rep=TRUE)))
colnames(Cni2Pfu_mean_JSD) <- colnames(Cni2Pfu_orth_tpm_qn)[-c(1)]
rownames(Cni2Pfu_mean_JSD) <- colnames(Pfu2Cni_orth_tpm_qn)[-c(1)]
Cni2Pfu_sd_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Pfu)-1,rep=TRUE)))
colnames(Cni2Pfu_sd_JSD) <- colnames(Cni2Pfu_orth_tpm_qn)[-c(1)]
rownames(Cni2Pfu_sd_JSD) <- colnames(Pfu2Cni_orth_tpm_qn)[-c(1)]

for (i in colnames(Cni2Pfu_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Pfu2Cni_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Cni[[i]],match_Pfu[[j]])
    Cni2Pfu_full_JSD[j,i] <- JSD(match)
    for (k in c(1:1000)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Cni2Pfu_mean_JSD[j,i] <- mean(all_JSD)
    Cni2Pfu_sd_JSD[j,i] <- sd(all_JSD)
  }
}

write.table(Cni2Pfu_full_JSD, "Cni2Pfu_full_set_JSD.txt", sep ='\t')
write.table(Cni2Pfu_mean_JSD, "Cni2Pfu_subset_JSD_mean.txt", sep ='\t')
write.table(Cni2Pfu_sd_JSD, "Cni2Pfu_subset_JSD_sd.txt", sep ='\t')

# 3. Cni vs. Pye

Cni2Pye_orth_tpm_qn <- read.csv("Cni2Pye_TPM_mean_quantile_transform.csv", header = T)
Pye2Cni_orth_tpm_qn <- read.csv("Pye2Cni_TPM_mean_quantile_transform.csv", header = T)

both_ID <- read.table("Cni2Pye.txt", header = T)
match_Cni <- Cni2Pye_orth_tpm_qn[match(both_ID$Cni,Cni2Pye_orth_tpm_qn$Gene_ID),]
match_Pye <- Pye2Cni_orth_tpm_qn[match(both_ID$Pye,Pye2Cni_orth_tpm_qn$Gene_ID),]
write.table(match_Cni,"order_Cni2Pye_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Pye,"order_Pye2Cni_tpm_qn.txt", sep="\t", quote=F)
match_Cni <- read.table("order_Cni2Pye_tpm_qn.txt", header = T)
match_Pye <- read.table("order_Pye2Cni_tpm_qn.txt", header = T)

Cni2Pye_full_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Pye)-1,rep=TRUE)))
colnames(Cni2Pye_full_JSD) <- colnames(Cni2Pye_orth_tpm_qn)[-c(1)]
rownames(Cni2Pye_full_JSD) <- colnames(Pye2Cni_orth_tpm_qn)[-c(1)]
Cni2Pye_mean_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Pye)-1,rep=TRUE)))
colnames(Cni2Pye_mean_JSD) <- colnames(Cni2Pye_orth_tpm_qn)[-c(1)]
rownames(Cni2Pye_mean_JSD) <- colnames(Pye2Cni_orth_tpm_qn)[-c(1)]
Cni2Pye_sd_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Pye)-1,rep=TRUE)))
colnames(Cni2Pye_sd_JSD) <- colnames(Cni2Pye_orth_tpm_qn)[-c(1)]
rownames(Cni2Pye_sd_JSD) <- colnames(Pye2Cni_orth_tpm_qn)[-c(1)]

for (i in colnames(Cni2Pye_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Pye2Cni_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Cni[[i]],match_Pye[[j]])
    Cni2Pye_full_JSD[j,i] <- JSD(match)
    for (k in c(1:1000)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Cni2Pye_mean_JSD[j,i] <- mean(all_JSD)
    Cni2Pye_sd_JSD[j,i] <- sd(all_JSD)
  }
}

write.table(Cni2Pye_full_JSD, "Cni2Pye_full_set_JSD.txt", sep ='\t')
write.table(Cni2Pye_mean_JSD, "Cni2Pye_subset_JSD_mean.txt", sep ='\t')
write.table(Cni2Pye_sd_JSD, "Cni2Pye_subset_JSD_sd.txt", sep ='\t')


# 4. Cni vs. Hru

Cni2Hru_orth_tpm_qn <- read.csv("Cni2Hru_TPM_mean_quantile_transform.csv", header = T)
Hru2Cni_orth_tpm_qn <- read.csv("Hru2Cni_TPM_mean_quantile_transform.csv", header = T)

both_ID <- read.table("Cni2Hru.txt", header = T)
match_Cni <- Cni2Hru_orth_tpm_qn[match(both_ID$Cni,Cni2Hru_orth_tpm_qn$Gene_ID),]
match_Hru <- Hru2Cni_orth_tpm_qn[match(both_ID$Hru,Hru2Cni_orth_tpm_qn$Gene_ID),]
write.table(match_Cni,"order_Cni2Hru_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Hru,"order_Hru2Cni_tpm_qn.txt", sep="\t", quote=F)
match_Cni <- read.table("order_Cni2Hru_tpm_qn.txt", header = T)
match_Hru <- read.table("order_Hru2Cni_tpm_qn.txt", header = T)

Cni2Hru_full_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Hru)-1,rep=TRUE)))
colnames(Cni2Hru_full_JSD) <- colnames(Cni2Hru_orth_tpm_qn)[-c(1)]
rownames(Cni2Hru_full_JSD) <- colnames(Hru2Cni_orth_tpm_qn)[-c(1)]
Cni2Hru_mean_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Hru)-1,rep=TRUE)))
colnames(Cni2Hru_mean_JSD) <- colnames(Cni2Hru_orth_tpm_qn)[-c(1)]
rownames(Cni2Hru_mean_JSD) <- colnames(Hru2Cni_orth_tpm_qn)[-c(1)]
Cni2Hru_sd_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Hru)-1,rep=TRUE)))
colnames(Cni2Hru_sd_JSD) <- colnames(Cni2Hru_orth_tpm_qn)[-c(1)]
rownames(Cni2Hru_sd_JSD) <- colnames(Hru2Cni_orth_tpm_qn)[-c(1)]

for (i in colnames(Cni2Hru_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Hru2Cni_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Cni[[i]],match_Hru[[j]])
    Cni2Hru_full_JSD[j,i] <- JSD(match)
    for (k in c(1:1000)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Cni2Hru_mean_JSD[j,i] <- mean(all_JSD)
    Cni2Hru_sd_JSD[j,i] <- sd(all_JSD)
  }
}

write.table(Cni2Hru_full_JSD, "Cni2Hru_full_set_JSD.txt", sep ='\t')
write.table(Cni2Hru_mean_JSD, "Cni2Hru_subset_JSD_mean.txt", sep ='\t')
write.table(Cni2Hru_sd_JSD, "Cni2Hru_subset_JSD_sd.txt", sep ='\t')


# 5. Cni vs. Pca

Cni2Pca_orth_tpm_qn <- read.csv("Cni2Pca_TPM_mean_quantile_transform.csv", header = T)
Pca2Cni_orth_tpm_qn <- read.csv("Pca2Cni_TPM_mean_quantile_transform.csv", header = T)

both_ID <- read.table("Cni2Pca.txt", header = T)
match_Cni <- Cni2Pca_orth_tpm_qn[match(both_ID$Cni,Cni2Pca_orth_tpm_qn$Gene_ID),]
match_Pca <- Pca2Cni_orth_tpm_qn[match(both_ID$Pca,Pca2Cni_orth_tpm_qn$Gene_ID),]
write.table(match_Cni,"order_Cni2Pca_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Pca,"order_Pca2Cni_tpm_qn.txt", sep="\t", quote=F)
match_Cni <- read.table("order_Cni2Pca_tpm_qn.txt", header = T)
match_Pca <- read.table("order_Pca2Cni_tpm_qn.txt", header = T)

Cni2Pca_full_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Pca)-1,rep=TRUE)))
colnames(Cni2Pca_full_JSD) <- colnames(Cni2Pca_orth_tpm_qn)[-c(1)]
rownames(Cni2Pca_full_JSD) <- colnames(Pca2Cni_orth_tpm_qn)[-c(1)]
Cni2Pca_mean_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Pca)-1,rep=TRUE)))
colnames(Cni2Pca_mean_JSD) <- colnames(Cni2Pca_orth_tpm_qn)[-c(1)]
rownames(Cni2Pca_mean_JSD) <- colnames(Pca2Cni_orth_tpm_qn)[-c(1)]
Cni2Pca_sd_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Pca)-1,rep=TRUE)))
colnames(Cni2Pca_sd_JSD) <- colnames(Cni2Pca_orth_tpm_qn)[-c(1)]
rownames(Cni2Pca_sd_JSD) <- colnames(Pca2Cni_orth_tpm_qn)[-c(1)]

for (i in colnames(Cni2Pca_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Pca2Cni_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Cni[[i]],match_Pca[[j]])
    Cni2Pca_full_JSD[j,i] <- JSD(match)
    for (k in c(1:1000)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Cni2Pca_mean_JSD[j,i] <- mean(all_JSD)
    Cni2Pca_sd_JSD[j,i] <- sd(all_JSD)
  }
}

write.table(Cni2Pca_full_JSD, "Cni2Pca_full_set_JSD.txt", sep ='\t')
write.table(Cni2Pca_mean_JSD, "Cni2Pca_subset_JSD_mean.txt", sep ='\t')
write.table(Cni2Pca_sd_JSD, "Cni2Pca_subset_JSD_sd.txt", sep ='\t')

# 6. Cni vs. Pve

Cni2Pve_orth_tpm_qn <- read.csv("Cni2Pve_TPM_mean_quantile_transform.csv", header = T)
Pve2Cni_orth_tpm_qn <- read.csv("Pve2Cni_TPM_mean_quantile_transform.csv", header = T)

both_ID <- read.table("Cni2Pve.txt", header = T)
match_Cni <- Cni2Pve_orth_tpm_qn[match(both_ID$Cni,Cni2Pve_orth_tpm_qn$Gene_ID),]
match_Pve <- Pve2Cni_orth_tpm_qn[match(both_ID$Pve,Pve2Cni_orth_tpm_qn$Gene_ID),]
write.table(match_Cni,"order_Cni2Pve_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Pve,"order_Pve2Cni_tpm_qn.txt", sep="\t", quote=F)
match_Cni <- read.table("order_Cni2Pve_tpm_qn.txt", header = T)
match_Pve <- read.table("order_Pve2Cni_tpm_qn.txt", header = T)

Cni2Pve_full_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Pve)-1,rep=TRUE)))
colnames(Cni2Pve_full_JSD) <- colnames(Cni2Pve_orth_tpm_qn)[-c(1)]
rownames(Cni2Pve_full_JSD) <- colnames(Pve2Cni_orth_tpm_qn)[-c(1)]
Cni2Pve_mean_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Pve)-1,rep=TRUE)))
colnames(Cni2Pve_mean_JSD) <- colnames(Cni2Pve_orth_tpm_qn)[-c(1)]
rownames(Cni2Pve_mean_JSD) <- colnames(Pve2Cni_orth_tpm_qn)[-c(1)]
Cni2Pve_sd_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Pve)-1,rep=TRUE)))
colnames(Cni2Pve_sd_JSD) <- colnames(Cni2Pve_orth_tpm_qn)[-c(1)]
rownames(Cni2Pve_sd_JSD) <- colnames(Pve2Cni_orth_tpm_qn)[-c(1)]

for (i in colnames(Cni2Pve_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Pve2Cni_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Cni[[i]],match_Pve[[j]])
    Cni2Pve_full_JSD[j,i] <- JSD(match)
    for (k in c(1:1000)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Cni2Pve_mean_JSD[j,i] <- mean(all_JSD)
    Cni2Pve_sd_JSD[j,i] <- sd(all_JSD)
  }
}

write.table(Cni2Pve_full_JSD, "Cni2Pve_full_set_JSD.txt", sep ='\t')
write.table(Cni2Pve_mean_JSD, "Cni2Pve_subset_JSD_mean.txt", sep ='\t')
write.table(Cni2Pve_sd_JSD, "Cni2Pve_subset_JSD_sd.txt", sep ='\t')


# 7. Cni vs. Npo

Cni2Npo_orth_tpm_qn <- read.csv("Cni2Npo_TPM_mean_quantile_transform.csv", header = T)
Npo2Cni_orth_tpm_qn <- read.csv("Npo2Cni_TPM_mean_quantile_transform.csv", header = T)

both_ID <- read.table("Cni2Npo.txt", header = T)
match_Cni <- Cni2Npo_orth_tpm_qn[match(both_ID$Cni,Cni2Npo_orth_tpm_qn$Gene_ID),]
match_Npo <- Npo2Cni_orth_tpm_qn[match(both_ID$Npo,Npo2Cni_orth_tpm_qn$Gene_ID),]
write.table(match_Cni,"order_Cni2Npo_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Npo,"order_Npo2Cni_tpm_qn.txt", sep="\t", quote=F)
match_Cni <- read.table("order_Cni2Npo_tpm_qn.txt", header = T)
match_Npo <- read.table("order_Npo2Cni_tpm_qn.txt", header = T)

Cni2Npo_full_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Npo)-1,rep=TRUE)))
colnames(Cni2Npo_full_JSD) <- colnames(Cni2Npo_orth_tpm_qn)[-c(1)]
rownames(Cni2Npo_full_JSD) <- colnames(Npo2Cni_orth_tpm_qn)[-c(1)]
Cni2Npo_mean_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Npo)-1,rep=TRUE)))
colnames(Cni2Npo_mean_JSD) <- colnames(Cni2Npo_orth_tpm_qn)[-c(1)]
rownames(Cni2Npo_mean_JSD) <- colnames(Npo2Cni_orth_tpm_qn)[-c(1)]
Cni2Npo_sd_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Npo)-1,rep=TRUE)))
colnames(Cni2Npo_sd_JSD) <- colnames(Cni2Npo_orth_tpm_qn)[-c(1)]
rownames(Cni2Npo_sd_JSD) <- colnames(Npo2Cni_orth_tpm_qn)[-c(1)]

for (i in colnames(Cni2Npo_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Npo2Cni_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Cni[[i]],match_Npo[[j]])
    Cni2Npo_full_JSD[j,i] <- JSD(match)
    for (k in c(1:1000)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Cni2Npo_mean_JSD[j,i] <- mean(all_JSD)
    Cni2Npo_sd_JSD[j,i] <- sd(all_JSD)
  }
}


write.table(Cni2Npo_full_JSD, "Cni2Npo_full_set_JSD.txt", sep ='\t')
write.table(Cni2Npo_mean_JSD, "Cni2Npo_subset_JSD_mean.txt", sep ='\t')
write.table(Cni2Npo_sd_JSD, "Cni2Npo_subset_JSD_sd.txt", sep ='\t')

# 8. Cni vs. Lan

Cni2Lan_orth_tpm_qn <- read.csv("Cni2Lan_TPM_mean_quantile_transform.csv", header = T)
Lan2Cni_orth_tpm_qn <- read.csv("Lan2Cni_TPM_mean_quantile_transform.csv", header = T)

both_ID <- read.table("Cni2Lan.txt", header = T)
match_Cni <- Cni2Lan_orth_tpm_qn[match(both_ID$Cni,Cni2Lan_orth_tpm_qn$Gene_ID),]
match_Lan <- Lan2Cni_orth_tpm_qn[match(both_ID$Lan,Lan2Cni_orth_tpm_qn$Gene_ID),]
write.table(match_Cni,"order_Cni2Lan_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Lan,"order_Lan2Cni_tpm_qn.txt", sep="\t", quote=F)
match_Cni <- read.table("order_Cni2Lan_tpm_qn.txt", header = T)
match_Lan <- read.table("order_Lan2Cni_tpm_qn.txt", header = T)

Cni2Lan_full_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Lan)-1,rep=TRUE)))
colnames(Cni2Lan_full_JSD) <- colnames(Cni2Lan_orth_tpm_qn)[-c(1)]
rownames(Cni2Lan_full_JSD) <- colnames(Lan2Cni_orth_tpm_qn)[-c(1)]
Cni2Lan_mean_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Lan)-1,rep=TRUE)))
colnames(Cni2Lan_mean_JSD) <- colnames(Cni2Lan_orth_tpm_qn)[-c(1)]
rownames(Cni2Lan_mean_JSD) <- colnames(Lan2Cni_orth_tpm_qn)[-c(1)]
Cni2Lan_sd_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Lan)-1,rep=TRUE)))
colnames(Cni2Lan_sd_JSD) <- colnames(Cni2Lan_orth_tpm_qn)[-c(1)]
rownames(Cni2Lan_sd_JSD) <- colnames(Lan2Cni_orth_tpm_qn)[-c(1)]

for (i in colnames(Cni2Lan_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Lan2Cni_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Cni[[i]],match_Lan[[j]])
    Cni2Lan_full_JSD[j,i] <- JSD(match)
    for (k in c(1:1000)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Cni2Lan_mean_JSD[j,i] <- mean(all_JSD)
    Cni2Lan_sd_JSD[j,i] <- sd(all_JSD)
  }
}


write.table(Cni2Lan_full_JSD, "Cni2Lan_full_set_JSD.txt", sep ='\t')
write.table(Cni2Lan_mean_JSD, "Cni2Lan_subset_JSD_mean.txt", sep ='\t')
write.table(Cni2Lan_sd_JSD, "Cni2Lan_subset_JSD_sd.txt", sep ='\t')

# 9. Cni vs. Ofu

Cni2Ofu_orth_tpm_qn <- read.csv("Cni2Ofu_TPM_mean_quantile_transform.csv", header = T)
Ofu2Cni_orth_tpm_qn <- read.csv("Ofu2Cni_TPM_mean_quantile_transform.csv", header = T)

both_ID <- read.table("Cni2Ofu.txt", header = T)
match_Cni <- Cni2Ofu_orth_tpm_qn[match(both_ID$Cni,Cni2Ofu_orth_tpm_qn$Gene_ID),]
match_Ofu <- Ofu2Cni_orth_tpm_qn[match(both_ID$Ofu,Ofu2Cni_orth_tpm_qn$Gene_ID),]
write.table(match_Cni,"order_Cni2Ofu_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Ofu,"order_Ofu2Cni_tpm_qn.txt", sep="\t", quote=F)
match_Cni <- read.table("order_Cni2Ofu_tpm_qn.txt", header = T)
match_Ofu <- read.table("order_Ofu2Cni_tpm_qn.txt", header = T)

Cni2Ofu_full_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Ofu)-1,rep=TRUE)))
colnames(Cni2Ofu_full_JSD) <- colnames(Cni2Ofu_orth_tpm_qn)[-c(1)]
rownames(Cni2Ofu_full_JSD) <- colnames(Ofu2Cni_orth_tpm_qn)[-c(1)]
Cni2Ofu_mean_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Ofu)-1,rep=TRUE)))
colnames(Cni2Ofu_mean_JSD) <- colnames(Cni2Ofu_orth_tpm_qn)[-c(1)]
rownames(Cni2Ofu_mean_JSD) <- colnames(Ofu2Cni_orth_tpm_qn)[-c(1)]
Cni2Ofu_sd_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Ofu)-1,rep=TRUE)))
colnames(Cni2Ofu_sd_JSD) <- colnames(Cni2Ofu_orth_tpm_qn)[-c(1)]
rownames(Cni2Ofu_sd_JSD) <- colnames(Ofu2Cni_orth_tpm_qn)[-c(1)]

for (i in colnames(Cni2Ofu_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Ofu2Cni_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Cni[[i]],match_Ofu[[j]])
    Cni2Ofu_full_JSD[j,i] <- JSD(match)
    for (k in c(1:1000)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Cni2Ofu_mean_JSD[j,i] <- mean(all_JSD)
    Cni2Ofu_sd_JSD[j,i] <- sd(all_JSD)
  }
}


write.table(Cni2Ofu_full_JSD, "Cni2Ofu_full_set_JSD.txt", sep ='\t')
write.table(Cni2Ofu_mean_JSD, "Cni2Ofu_subset_JSD_mean.txt", sep ='\t')
write.table(Cni2Ofu_sd_JSD, "Cni2Ofu_subset_JSD_sd.txt", sep ='\t')

# 10. Cni vs. Pec

Cni2Pec_orth_tpm_qn <- read.csv("Cni2Pec_TPM_mean_quantile_transform.csv", header = T)
Pec2Cni_orth_tpm_qn <- read.csv("Pec2Cni_TPM_mean_quantile_transform.csv", header = T)

both_ID <- read.table("Cni2Pec.txt", header = T)
match_Cni <- Cni2Pec_orth_tpm_qn[match(both_ID$Cni,Cni2Pec_orth_tpm_qn$Gene_ID),]
match_Pec <- Pec2Cni_orth_tpm_qn[match(both_ID$Pec,Pec2Cni_orth_tpm_qn$Gene_ID),]
write.table(match_Cni,"order_Cni2Pec_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Pec,"order_Pec2Cni_tpm_qn.txt", sep="\t", quote=F)
match_Cni <- read.table("order_Cni2Pec_tpm_qn.txt", header = T)
match_Pec <- read.table("order_Pec2Cni_tpm_qn.txt", header = T)

Cni2Pec_full_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Pec)-1,rep=TRUE)))
colnames(Cni2Pec_full_JSD) <- colnames(Cni2Pec_orth_tpm_qn)[-c(1)]
rownames(Cni2Pec_full_JSD) <- colnames(Pec2Cni_orth_tpm_qn)[-c(1)]
Cni2Pec_mean_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Pec)-1,rep=TRUE)))
colnames(Cni2Pec_mean_JSD) <- colnames(Cni2Pec_orth_tpm_qn)[-c(1)]
rownames(Cni2Pec_mean_JSD) <- colnames(Pec2Cni_orth_tpm_qn)[-c(1)]
Cni2Pec_sd_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Pec)-1,rep=TRUE)))
colnames(Cni2Pec_sd_JSD) <- colnames(Cni2Pec_orth_tpm_qn)[-c(1)]
rownames(Cni2Pec_sd_JSD) <- colnames(Pec2Cni_orth_tpm_qn)[-c(1)]

for (i in colnames(Cni2Pec_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Pec2Cni_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Cni[[i]],match_Pec[[j]])
    Cni2Pec_full_JSD[j,i] <- JSD(match)
    for (k in c(1:1000)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Cni2Pec_mean_JSD[j,i] <- mean(all_JSD)
    Cni2Pec_sd_JSD[j,i] <- sd(all_JSD)
  }
}

write.table(Cni2Pec_full_JSD, "Cni2Pec_full_set_JSD.txt", sep ='\t')
write.table(Cni2Pec_mean_JSD, "Cni2Pec_subset_JSD_mean.txt", sep ='\t')
write.table(Cni2Pec_sd_JSD, "Cni2Pec_subset_JSD_sd.txt", sep ='\t')

# 11. Cni vs. Aam

Cni2Aam_orth_tpm_qn <- read.csv("Cni2Aam_TPM_mean_quantile_transform.csv", header = T)
Aam2Cni_orth_tpm_qn <- read.csv("Aam2Cni_TPM_mean_quantile_transform.csv", header = T)

both_ID <- read.table("Cni2Aam.txt", header = T)
match_Cni <- Cni2Aam_orth_tpm_qn[match(both_ID$Cni,Cni2Aam_orth_tpm_qn$Gene_ID),]
match_Aam <- Aam2Cni_orth_tpm_qn[match(both_ID$Aam,Aam2Cni_orth_tpm_qn$Gene_ID),]
write.table(match_Cni,"order_Cni2Aam_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Aam,"order_Aam2Cni_tpm_qn.txt", sep="\t", quote=F)
match_Cni <- read.table("order_Cni2Aam_tpm_qn.txt", header = T)
match_Aam <- read.table("order_Aam2Cni_tpm_qn.txt", header = T)

Cni2Aam_full_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Aam)-1,rep=TRUE)))
colnames(Cni2Aam_full_JSD) <- colnames(Cni2Aam_orth_tpm_qn)[-c(1)]
rownames(Cni2Aam_full_JSD) <- colnames(Aam2Cni_orth_tpm_qn)[-c(1)]
Cni2Aam_mean_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Aam)-1,rep=TRUE)))
colnames(Cni2Aam_mean_JSD) <- colnames(Cni2Aam_orth_tpm_qn)[-c(1)]
rownames(Cni2Aam_mean_JSD) <- colnames(Aam2Cni_orth_tpm_qn)[-c(1)]
Cni2Aam_sd_JSD <- data.frame(replicate(ncol(match_Cni)-1,sample(0:1,ncol(match_Aam)-1,rep=TRUE)))
colnames(Cni2Aam_sd_JSD) <- colnames(Cni2Aam_orth_tpm_qn)[-c(1)]
rownames(Cni2Aam_sd_JSD) <- colnames(Aam2Cni_orth_tpm_qn)[-c(1)]

for (i in colnames(Cni2Aam_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Aam2Cni_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Cni[[i]],match_Aam[[j]])
    Cni2Aam_full_JSD[j,i] <- JSD(match)
    for (k in c(1:1000)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Cni2Aam_mean_JSD[j,i] <- mean(all_JSD)
    Cni2Aam_sd_JSD[j,i] <- sd(all_JSD)
  }
}

write.table(Cni2Aam_full_JSD, "Cni2Aam_full_set_JSD.txt", sep ='\t')
write.table(Cni2Aam_mean_JSD, "Cni2Aam_subset_JSD_mean.txt", sep ='\t')
write.table(Cni2Aam_sd_JSD, "Cni2Aam_subset_JSD_sd.txt", sep ='\t')

Cni2Cgi_mean <- read.table("Cni2Cgi_subset_JSD_mean.txt", header = T)
Cni2Pfu_mean <- read.table("Cni2Pfu_subset_JSD_mean.txt", header = T)
Cni2Pye_mean <- read.table("Cni2Pye_subset_JSD_mean.txt", header = T)
Cni2Hru_mean <- read.table("Cni2Hru_subset_JSD_mean.txt", header = T)
Cni2Pca_mean <- read.table("Cni2Pca_subset_JSD_mean.txt", header = T)
Cni2Npo_mean <- read.table("Cni2Npo_subset_JSD_mean.txt", header = T)
Cni2Pve_mean <- read.table("Cni2Pve_subset_JSD_mean.txt", header = T)
Cni2Lan_mean <- read.table("Cni2Lan_subset_JSD_mean.txt", header = T)
Cni2Ofu_mean <- read.table("Cni2Ofu_subset_JSD_mean.txt", header = T)
Cni2Pec_mean <- read.table("Cni2Pec_subset_JSD_mean.txt", header = T)
Cni2Aam_mean <- read.table("Cni2Aam_subset_JSD_mean.txt", header = T)

heatmap_color <- colorRampPalette(brewer.pal(n = 7, name = "RdGy"))(200)

h1 <- ComplexHeatmap::Heatmap(data.matrix(Cni2Cgi_mean),cluster_rows=F,cluster_columns=F,
                              width=ncol(Cni2Cgi_mean)*unit(5,"mm"),
                              height=nrow(Cni2Cgi_mean)*unit(5,"mm"),
                              col=heatmap_color,
              column_labels = c("M","AM","Gi","DG","He"),
              row_labels = c("MA","AM","FG","MG","Gi","DG","He","LP"),
              heatmap_legend_param = list(title = "JSD_raw"))

h2 <- ComplexHeatmap::Heatmap(data.matrix(Cni2Pfu_mean),cluster_rows=F,cluster_columns=F,
                              width=ncol(Cni2Pfu_mean)*unit(5,"mm"),
                              height=nrow(Cni2Pfu_mean)*unit(5,"mm"),
                              col=heatmap_color,
              column_labels = c("M","AM","Gi","DG","He"),
              row_labels = c("M","AM","Go","Gi"),
              heatmap_legend_param = list(title = "JSD_raw"))

h3 <- ComplexHeatmap::Heatmap(data.matrix(Cni2Pye_mean),cluster_rows=F,cluster_columns=F,
                              width=ncol(Cni2Pye_mean)*unit(5,"mm"),
                              height=nrow(Cni2Pye_mean)*unit(5,"mm"),
                              col=heatmap_color,
              column_labels = c("M","AM","Gi","DG","He"),
              row_labels = c("M","AM","F","FG","MG","Gi","DG","K"),
              heatmap_legend_param = list(title = "JSD_raw"))

h4 <- ComplexHeatmap::Heatmap(data.matrix(Cni2Hru_mean),cluster_rows=F,cluster_columns=F,
                              width=ncol(Cni2Hru_mean)*unit(5,"mm"),
                              height=nrow(Cni2Hru_mean)*unit(5,"mm"),
                              col=heatmap_color,
              column_labels = c("M","AM","Gi","DG","He"),
              row_labels = c("M","F","PE","LR","MG","FG","H","CT","Ep","ET","Ga","Gi","L","K"),
              heatmap_legend_param = list(title = "JSD_raw"))

h5 <- ComplexHeatmap::Heatmap(data.matrix(Cni2Pca_mean),cluster_rows=F,cluster_columns=F,
                              width=ncol(Cni2Pca_mean)*unit(5,"mm"),
                              height=nrow(Cni2Pca_mean)*unit(5,"mm"),
                              col=heatmap_color,
              column_labels = c("M","AM","Gi","DG","He"),
              row_labels = c("mantle","Albumen_gland","Digestive_gland","Foot","Gill","Kidney","Lung","Stomach","Testis"),
              heatmap_legend_param = list(title = "JSD_raw"))

h6 <- ComplexHeatmap::Heatmap(data.matrix(Cni2Pve_mean),cluster_rows=F,cluster_columns=F,
                              width=ncol(Cni2Pve_mean)*unit(5,"mm"),
                              height=nrow(Cni2Pve_mean)*unit(5,"mm"),
                              col=heatmap_color,
              column_labels = c("M","AM","Gi","DG","He"),
              row_labels = c("Mantle_posterior_aperture","Mantle","Radular_sac","Captacula","Digestive_diverticulae","Foot","Gonad","Intestine","Odontophore","Longitudinal_muscle","Oral_tube"),
              heatmap_legend_param = list(title = "JSD_raw"))

h7 <- ComplexHeatmap::Heatmap(data.matrix(Cni2Npo_mean),cluster_rows=F,cluster_columns=F,
                              width=ncol(Cni2Npo_mean)*unit(5,"mm"),
                              height=nrow(Cni2Npo_mean)*unit(5,"mm"),
                              col=heatmap_color,
              column_labels = c("M","AM","Gi","DG","He"),
              row_labels = c("Mantle","Eye","Gill","Hemocyte","Liver"),
              heatmap_legend_param = list(title = "JSD_raw"))

h8 <- ComplexHeatmap::Heatmap(data.matrix(Cni2Lan_mean),cluster_rows=F,cluster_columns=F,
                              width=ncol(Cni2Lan_mean)*unit(5,"mm"),
                              height=nrow(Cni2Lan_mean)*unit(5,"mm"),
                              col=heatmap_color,
                              column_labels = c("M","AM","Gi","DG","He"),
                              row_labels = c("adult_ventral_mantle","adult_dorsal_mantle","adult_digestive_cecum","adult_gut","adult_lophophore","adult_pedicle"),
                              heatmap_legend_param = list(title = "JSD_raw"))

h9 <- ComplexHeatmap::Heatmap(data.matrix(Cni2Ofu_mean),cluster_rows=F,cluster_columns=F,
                              width=ncol(Cni2Ofu_mean)*unit(5,"mm"),
                              height=nrow(Cni2Ofu_mean)*unit(5,"mm"),
                              col=heatmap_color,
                              column_labels = c("M","AM","Gi","DG","He"),
                              row_labels = c("adult_head","adult_head_chaetae","adult_tail","adult_body","adult_blood","adult_gut","adult_retractor-muscle","adult_testis","adult_ovary"),
                              heatmap_legend_param = list(title = "JSD_raw"))

h10 <- ComplexHeatmap::Heatmap(data.matrix(Cni2Pec_mean),cluster_rows=F,cluster_columns=F,
                              width=ncol(Cni2Pec_mean)*unit(5,"mm"),
                              height=nrow(Cni2Pec_mean)*unit(5,"mm"),
                              col=heatmap_color,
                              column_labels = c("M","AM","Gi","DG","He"),
                              row_labels = c("Collar","Opisthosoma","Plume","Trophosome_anterior","Trophosome_middle","Trophosome_posterior","Vestimentum"),
                              heatmap_legend_param = list(title = "JSD_raw"))

h11 <- ComplexHeatmap::Heatmap(data.matrix(Cni2Aam_mean),cluster_rows=F,cluster_columns=F,
                              width=ncol(Cni2Aam_mean)*unit(5,"mm"),
                              height=nrow(Cni2Aam_mean)*unit(5,"mm"),
                              col=heatmap_color,
                              column_labels = c("M","AM","Gi","DG","He"),
                              row_labels = c("Mantle","Adductor","Cirrus","Visceral_mass"),
                              heatmap_legend_param = list(title = "JSD_raw"))

# Save heatmap for Cni2Cgi_mean
pdf("Cni2Cgi_raw_mean_heatmap.pdf")
h1
dev.off()

# Save heatmap for Cni2Pfu_mean
pdf("Cni2Pfu_raw_mean_heatmap.pdf")
h2
dev.off()

# Save heatmap for Cni2Pye_mean
pdf("Cni2Pye_raw_mean_heatmap.pdf")
h3
dev.off()

# Save heatmap for Cni2Hru_mean
pdf("Cni2Hru_raw_mean_heatmap.pdf")
h4
dev.off()

# Save heatmap for Cni2Pca_mean
pdf("Cni2Pca_raw_mean_heatmap.pdf")
h5
dev.off()

# Save heatmap for Cni2Hru_mean
pdf("Cni2Pve_raw_mean_heatmap.pdf")
h6
dev.off()

# Save heatmap for Cni2Pca_mean
pdf("Cni2Npo_raw_mean_heatmap.pdf")
h7
dev.off()

pdf("Cni2Lan_raw_mean_heatmap.pdf")
h8
dev.off()

pdf("Cni2Ofu_raw_mean_heatmap.pdf")
h9
dev.off()

pdf("Cni2Pec_raw_mean_heatmap.pdf")
h10
dev.off()

pdf("Cni2Aam_raw_mean_heatmap.pdf")
h11
dev.off()
