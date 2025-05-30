library(BiocParallel)
library(DiffBind)

dbObj <- dba(sampleSheet = "samplesheet.csv")
samples_count_nofilter <- dba.count(dbObj, summits=T,filter=0, bParallel=T)
samples_count_nofilter_normalized <- dba.normalize(samples_count_nofilter,normalize=DBA_NORM_TMM,background=FALSE,library=DBA_LIBSIZE_FULL,offsets=TRUE,method=DBA_EDGER)
samples_count_nofilter_normalized_peakset <- dba.peakset(samples_count_nofilter_normalized,bRetrieve = TRUE)
write.table(samples_count_nofilter_normalized_peakset,"Diffbind_nofilter_peakset.txt",sep = '\t',quote = F)
