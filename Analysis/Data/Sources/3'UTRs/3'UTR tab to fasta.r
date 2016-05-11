library(Biostrings)

setwd("C:/Users/mark/iCloudDrive/MSc/BIOL61230 Research Projects/1/Analysis")

# read the 3'UTR tab file
utr3 <- read.table(file = "Data/Sources/3'UTRs/3'UTRs.tab", header = T)

# build a DNA string set and set the names to the transcript IDs
dna <- DNAStringSet(utr3$X3utr)
names(dna) <- utr3$ensembl_transcript_id

# write to a file
writeXStringSet(dna, filepath = "Data/Sources/3'UTRs/3'UTRs.fa")