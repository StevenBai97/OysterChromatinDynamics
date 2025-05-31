(head -n1 ../../../qttpm/all_raw_mean_qn/Spu.mean.dev.raw_tpmquantile_transform.csv && grep -Ff <(cut -f1 Spu2Lva.txt) ../../../qttpm/all_raw_mean_qn/Spu.mean.dev.raw_tpmquantile_transform.csv) > Spu2Lva_TPM_mean_quantile_transform.csv
(head -n1 ../../../qttpm/all_raw_mean_qn/Spu.mean.dev.raw_tpmquantile_transform.csv && grep -Ff <(cut -f1 Spu2Aja.txt) ../../../qttpm/all_raw_mean_qn/Spu.mean.dev.raw_tpmquantile_transform.csv) > Spu2Aja_TPM_mean_quantile_transform.csv

(head -n1 ../../../qttpm/all_raw_mean_qn/Lva.mean.dev.raw_tpmquantile_transform.csv && grep -Ff <(cut -f2 Spu2Lva.txt) ../../../qttpm/all_raw_mean_qn/Lva.mean.dev.raw_tpmquantile_transform.csv) > Lva2Spu_TPM_mean_quantile_transform.csv
(head -n1 ../../../qttpm/all_raw_mean_qn/Aja.mean.dev.raw_tpmquantile_transform.csv && grep -wFf <(cut -f2 Spu2Aja.txt) ../../../qttpm/all_raw_mean_qn/Aja.mean.dev.raw_tpmquantile_transform.csv) > Aja2Spu_TPM_mean_quantile_transform.csv

for file in *_TPM_mean_quantile_transform.csv; do
    sed -i '1s/^[^,]*/Gene_ID/' "$file"
done
