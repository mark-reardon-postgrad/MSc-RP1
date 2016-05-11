require(GenomicRanges)
setwd("C:/Users/mark/iCloudDrive/MSc/BIOL61230 Research Projects/1/Analysis")

if (!exists("decipherCNVs"))
{
  decipherCNVs <- read.table(file = "Data/Sources/Decipher CNVs/Decipher.tab", header = T, sep = "\t")
  decipherCNVs <- decipherCNVs[decipherCNVs$seqnames != "chrX", ]
  decipherCNVs <- makeGRangesFromDataFrame(df = decipherCNVs, keep.extra.columns = T)
  decipherCNVsGained <- decipherCNVs[decipherCNVs$type >= 0, ]
  decipherCNVsLost <- decipherCNVs[decipherCNVs$type <= 0, ]
  decipherCNVsBoth <- decipherCNVs[decipherCNVs$type == 0, ]
}
