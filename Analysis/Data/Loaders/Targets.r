setwd("C:/Users/mark/iCloudDrive/MSc/BIOL61230 Research Projects/1/Analysis")

if (!exists("targets"))
{
  targets <- read.table(
    file = "Data/Sources/Seed Vicious/targets.tab", header = T, sep = "\t",
    col.names = c("EnsembleTranscriptID", "Direction", "miRNA", "Position", "TargetType"))
}
