require(GenomicRanges)
setwd("C:/Users/mark/iCloudDrive/MSc/BIOL61230 Research Projects/1/Analysis")

if (!exists("dgvCNVs"))
{
  dgvCNVs <- read.table(file = "Data/Sources/DGV CNVs/GRCh38_hg38_variants_2015-07-23.txt", header = T, sep = "\t")
  dgvCNVs$chr <- as.numeric(as.character(dgvCNVs$chr))
  dgvCNVs <- dgvCNVs[!is.na(dgvCNVs$chr), ]
  dgvCNVs <- transform(dgvCNVs, chr = sprintf("chr%d", chr))
  dgvCNVs <- makeGRangesFromDataFrame(df = dgvCNVs[dgvCNVs$varianttype == "CNV", ], keep.extra.columns = T)
  dgvCNVsGained <- dgvCNVs[dgvCNVs$variantsubtype %in% c("duplication", "gain", "gain+loss"), ]
  dgvCNVsLost <- dgvCNVs[dgvCNVs$variantsubtype %in% c("deletion", "loss", "gain+loss"), ]
  dgvCNVsBoth <- dgvCNVs[dgvCNVs$variantsubtype %in% c("gain+loss"), ]
}
