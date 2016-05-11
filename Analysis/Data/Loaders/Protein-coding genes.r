require(GenomicRanges)
setwd("C:/Users/mark/iCloudDrive/MSc/BIOL61230 Research Projects/1/Analysis")

if (!exists("proteinCodingGenes"))
{
  # create GRanges
  proteinCodingGenes <- makeGRangesFromDataFrame(df = read.table(file = "Data/Sources/Protein-coding genes/Protein-coding genes.tab", 
                                                                 header = T, sep = "\t"), keep.extra.columns = T)
  # create midpoints variant
  proteinCodingGeneMidpoints <- proteinCodingGenes
  attributes(proteinCodingGeneMidpoints@ranges)$start <- as.integer(attributes(proteinCodingGeneMidpoints@ranges)$start + (attributes(proteinCodingGeneMidpoints@ranges)$width / 2))
  attributes(proteinCodingGeneMidpoints@ranges)$width <- as.integer(rep.int(1, length(proteinCodingGeneMidpoints)))
}
