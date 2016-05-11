require(GenomicRanges)
setwd("C:/Users/mark/iCloudDrive/MSc/BIOL61230 Research Projects/1/Analysis")

if (!exists("primiRNAs"))
{
  # create GRanges
  primiRNAs = read.table(file = "Data/Sources/miRNAs/hsa.gff3", header = FALSE, sep = "\t")
  primiRNAs <- rename(primiRNAs, c("V1" = "sequence", "V2" = "source", "V3" = "feature", "V4" = "start", "V5" = "end", 
                             "V6" = "score", "V7" = "strand", "V8" = "frame", "V9" = "attributes"))
  primiRNAs <- makeGRangesFromDataFrame(df = primiRNAs[primiRNAs$feature == "miRNA_primary_transcript" & primiRNAs$sequence %in% sprintf("chr%d", 1:22), ], 
                                     keep.extra.columns = TRUE, seqnames.field = "sequence")
  # create midpoints variant
  primiRNAMidpoints <- primiRNAs
  attributes(primiRNAMidpoints@ranges)$start <- as.integer(attributes(primiRNAMidpoints@ranges)$start + (attributes(primiRNAMidpoints@ranges)$width / 2))
  attributes(primiRNAMidpoints@ranges)$width <- as.integer(rep.int(1, length(primiRNAMidpoints)))
}
