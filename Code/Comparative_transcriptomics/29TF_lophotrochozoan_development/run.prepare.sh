(head -n1 ../../../qttpm/all_raw_mean_qn/Cni.mean.dev.raw_tpmquantile_transform.csv && grep -Ff <(cut -f1 Cni2Cgi.txt) ../../../qttpm/all_raw_mean_qn/Cni.mean.dev.raw_tpmquantile_transform.csv) > Cni2Cgi_TPM_mean_quantile_transform.csv
(head -n1 ../../../qttpm/all_raw_mean_qn/Cni.mean.dev.raw_tpmquantile_transform.csv && grep -Ff <(cut -f1 Cni2Pfu.txt) ../../../qttpm/all_raw_mean_qn/Cni.mean.dev.raw_tpmquantile_transform.csv) > Cni2Pfu_TPM_mean_quantile_transform.csv
(head -n1 ../../../qttpm/all_raw_mean_qn/Cni.mean.dev.raw_tpmquantile_transform.csv && grep -Ff <(cut -f1 Cni2Pye.txt) ../../../qttpm/all_raw_mean_qn/Cni.mean.dev.raw_tpmquantile_transform.csv) > Cni2Pye_TPM_mean_quantile_transform.csv
(head -n1 ../../../qttpm/all_raw_mean_qn/Cni.mean.dev.raw_tpmquantile_transform.csv && grep -Ff <(cut -f1 Cni2Hru.txt) ../../../qttpm/all_raw_mean_qn/Cni.mean.dev.raw_tpmquantile_transform.csv) > Cni2Hru_TPM_mean_quantile_transform.csv
(head -n1 ../../../qttpm/all_raw_mean_qn/Cni.mean.dev.raw_tpmquantile_transform.csv && grep -Ff <(cut -f1 Cni2Pca.txt) ../../../qttpm/all_raw_mean_qn/Cni.mean.dev.raw_tpmquantile_transform.csv) > Cni2Pca_TPM_mean_quantile_transform.csv
(head -n1 ../../../qttpm/all_raw_mean_qn/Cni.mean.dev.raw_tpmquantile_transform.csv && grep -Ff <(cut -f1 Cni2Lan.txt) ../../../qttpm/all_raw_mean_qn/Cni.mean.dev.raw_tpmquantile_transform.csv) > Cni2Lan_TPM_mean_quantile_transform.csv
(head -n1 ../../../qttpm/all_raw_mean_qn/Cni.mean.dev.raw_tpmquantile_transform.csv && grep -Ff <(cut -f1 Cni2Ofu.txt) ../../../qttpm/all_raw_mean_qn/Cni.mean.dev.raw_tpmquantile_transform.csv) > Cni2Ofu_TPM_mean_quantile_transform.csv
(head -n1 ../../../qttpm/all_raw_mean_qn/Cni.mean.dev.raw_tpmquantile_transform.csv && grep -Ff <(cut -f1 Cni2Aam.txt) ../../../qttpm/all_raw_mean_qn/Cni.mean.dev.raw_tpmquantile_transform.csv) > Cni2Aam_TPM_mean_quantile_transform.csv

(head -n1 ../../../qttpm/all_raw_mean_qn/Cgi.mean.dev.raw_tpmquantile_transform.csv && grep -Ff <(cut -f2 Cni2Cgi.txt) ../../../qttpm/all_raw_mean_qn/Cgi.mean.dev.raw_tpmquantile_transform.csv) > Cgi2Cni_TPM_mean_quantile_transform.csv
(head -n1 ../../../qttpm/all_raw_mean_qn/Pfu.mean.dev.raw_tpmquantile_transform.csv && grep -Ff <(cut -f2 Cni2Pfu.txt) ../../../qttpm/all_raw_mean_qn/Pfu.mean.dev.raw_tpmquantile_transform.csv) > Pfu2Cni_TPM_mean_quantile_transform.csv
(head -n1 ../../../qttpm/all_raw_mean_qn/Pye.mean.dev.raw_tpmquantile_transform.csv && grep -Ff <(cut -f2 Cni2Pye.txt) ../../../qttpm/all_raw_mean_qn/Pye.mean.dev.raw_tpmquantile_transform.csv) > Pye2Cni_TPM_mean_quantile_transform.csv
(head -n1 ../../../qttpm/all_raw_mean_qn/Hru.mean.dev.raw_tpmquantile_transform.csv && grep -Ff <(cut -f2 Cni2Hru.txt) ../../../qttpm/all_raw_mean_qn/Hru.mean.dev.raw_tpmquantile_transform.csv) > Hru2Cni_TPM_mean_quantile_transform.csv
(head -n1 ../../../qttpm/all_raw_mean_qn/Pca.mean.dev.raw_tpmquantile_transform.csv && grep -Ff <(cut -f2 Cni2Pca.txt) ../../../qttpm/all_raw_mean_qn/Pca.mean.dev.raw_tpmquantile_transform.csv) > Pca2Cni_TPM_mean_quantile_transform.csv
(head -n1 ../../../qttpm/all_raw_mean_qn/Lan.mean.dev.raw_tpmquantile_transform.csv && grep -Ff <(cut -f2 Cni2Lan.txt) ../../../qttpm/all_raw_mean_qn/Lan.mean.dev.raw_tpmquantile_transform.csv) > Lan2Cni_TPM_mean_quantile_transform.csv
(head -n1 ../../../qttpm/all_raw_mean_qn/Ofu.mean.dev.raw_tpmquantile_transform.csv && grep -Ff <(cut -f2 Cni2Ofu.txt) ../../../qttpm/all_raw_mean_qn/Ofu.mean.dev.raw_tpmquantile_transform.csv) > Ofu2Cni_TPM_mean_quantile_transform.csv
(head -n1 ../../../qttpm/all_raw_mean_qn/Aam.mean.dev.raw_tpmquantile_transform.csv && grep -Ff <(cut -f2 Cni2Aam.txt) ../../../qttpm/all_raw_mean_qn/Aam.mean.dev.raw_tpmquantile_transform.csv) > Aam2Cni_TPM_mean_quantile_transform.csv

for file in *_TPM_mean_quantile_transform.csv; do
    sed -i '1s/^[^,]*/Gene_ID/' "$file"
done
