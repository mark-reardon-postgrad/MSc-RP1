library(seqinr)

setwd("C:/Users/mark/iCloudDrive/MSc/BIOL61230 Research Projects/1/Analysis")

# read the mature miRNA seqeunces (all species)
allMiRNAs <- read.fasta(file = "Data/Sources/miRNAs/mature.fa", as.string = T, forceDNAtolower = F)

# extract the human miRNAs
humanMiRNAs <- allMiRNAs[grep("hsa-*", attributes(allMiRNAs)$name)]

# write the human mature miRNAs to a file
write.fasta(sequences = humanMiRNAs, names = names(humanMiRNAs), file.out = "Data/Sources/miRNAs/mature.human.fa")