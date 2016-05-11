setwd("C:/Users/mark/iCloudDrive/MSc/BIOL61230 Research Projects/1/Analysis")

if (!exists("chromosomes"))
{
  chromosomes <- read.table(
    file = "Data/Sources/Chromosomes/hg38.chrom.sizes", sep = "\t",
    col.names = c("Name", "Length"))
  chromosomes <- chromosomes[chromosomes$Name %in% sprintf("chr%d", 1:22), ]
}
