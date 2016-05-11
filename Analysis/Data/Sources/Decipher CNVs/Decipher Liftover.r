# source("http://bioconductor.org/workflows.R")
# workflowInstall("liftOver")

require(rtracklayer)
require(GenomicRanges)
require(liftOver)

setwd("C:/Users/mark/iCloudDrive/MSc/BIOL61230 Research Projects/1/Analysis")

# load the CRCh37 to GRCh38 chain file (from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/)
hg19ToHg38 <- import.chain("Data/Sources/Decipher CNVs/hg19ToHg38.over.chain")

# load the Decipher CNVs
decipherCNVs <- read.table(file = "Data/Sources/Decipher CNVs/population_cnv.txt", header = T, sep = "\t")

# convert to GRanges (ignore chr23 as only Decipher has it)
decipherCNVs <- transform(decipherCNVs, chr = sprintf("chr%d", chr))
decipherCNVs <- makeGRangesFromDataFrame(df = decipherCNVs[decipherCNVs$chr != "chr23", ], keep.extra.columns = T)

# convert from CRCh37 to GRCh38
decipherCNVsGRC38 <- liftOver(decipherCNVs, chain = hg19ToHg38)

# save to file
write.table(x = decipherCNVsGRC38, file = "Data/Sources/Decipher CNVs/Decipher.tab", row.names = F, sep = "\t")
