#!/usr/bin/env Rscript

library(argparser, quietly=TRUE)
library(Rsubread)
library(limma)
library(edgeR)

p <- arg_parser("Run featureCounts and calculate FPKM/TPM")

p <- add_argument(p, "--bam", help="input: bam file", type="character")
p <- add_argument(p, "--strandSpecific", help="Strand specificity (0, 1, 2)", type="numeric", default=0)
p <- add_argument(p, "--gtf", help="input: gtf file", type="character", default=NA)
p <- add_argument(p, "--saf", help="input: saf file", type="character", default=NA)
p <- add_argument(p, "--isPairedEnd", help="Paired-end reads, TRUE or FALSE", type="logical", default=TRUE)
p <- add_argument(p, "--featureType", help="Feature type in GTF annotation (default: exon)", type="character", default="exon")
p <- add_argument(p, "--attrType", help="Attribute type in GTF annotation (default: gene_id)", type="character", default="gene_id")
p <- add_argument(p, "--output", help="Output prefix", type="character")

argv <- parse_args(p)

if (!is.na(argv$gtf) && !is.na(argv$saf)) {
  stop("Cannot use both --gtf and --saf. Choose one.")
}

bamFile <- argv$bam
annotFile <- if (!is.na(argv$gtf)) argv$gtf else argv$saf
isGTF <- !is.na(argv$gtf)
outFilePref <- argv$output

outStatsFilePath  <- paste0(outFilePref, '.log')
outCountsFilePath <- paste0(outFilePref, '.count')

fCountsArgs <- list(
  files = bamFile,
  annot.ext = annotFile,
  nthreads = 20,
  isPairedEnd = argv$isPairedEnd,
  strandSpecific = argv$strandSpecific
)

if (isGTF) {
  fCountsArgs$isGTFAnnotationFile <- TRUE
  fCountsArgs$GTF.featureType <- argv$featureType
  fCountsArgs$GTF.attrType <- argv$attrType
} else {
  fCountsArgs$isGTFAnnotationFile <- FALSE
}

fCountsList <- do.call(featureCounts, fCountsArgs)

dgeList <- DGEList(counts = fCountsList$counts, genes = fCountsList$annotation)

cpm <- cpm(dgeList)
fpkm <- rpkm(dgeList, dgeList$genes$Length)
tpm <- exp(log(fpkm) - log(sum(fpkm)) + log(1e6))

write.table(fCountsList$stat, outStatsFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

firstColName <- argv$attrType
featureCounts <- cbind(fCountsList$annotation[,1], fCountsList$counts, fpkm, tpm, cpm)
colnames(featureCounts) <- c(firstColName, 'counts', 'fpkm', 'tpm', 'cpm')

write.table(featureCounts, outCountsFilePath, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

