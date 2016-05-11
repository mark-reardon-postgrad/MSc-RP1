setwd("C:/Users/mark/iCloudDrive/MSc/BIOL61230 Research Projects/1/Analysis")

if (!exists("UTRs"))
{
  UTRs <- read.table(
    file = "Data/Sources/3'UTRs/3'UTRs.tab", header = T, sep = "\t")
}
