#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")

library(biomaRt)

setwd("C:/Users/mark/iCloudDrive/MSc/BIOL61230 Research Projects/1/Analysis")

# GRCh38 human data source
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

# get all protein-coding transcripts for protein-coding genes on chromosomes 1..22
proteinCodingTranscripts <- getBM(mart = ensembl, 
  attributes = c("ensembl_gene_id", "ensembl_transcript_id"),
  filters = c("biotype", "transcript_biotype", "chromosome_name"), 
  values = list("protein_coding", "protein_coding", 1:22))
colnames(proteinCodingTranscripts) <- c("EnsembleGeneID", "EnsembleTranscriptID")

# 3'UTRs
utr3 <- getSequence(mart = ensembl, type = "ensembl_transcript_id", id = proteinCodingTranscripts$EnsembleTranscriptID, seqType = "3utr")
colnames(utr3) <- c("Sequence", "EnsembleTranscriptID")

# add the gene IDs
geneUTRs <- merge(proteinCodingTranscripts, utr3, by.x = "EnsembleTranscriptID", by.y = "EnsembleTranscriptID", all.x = T)

# save only where sequence was available
write.table(geneUTRs[geneUTRs$Sequence != "Sequence unavailable", ], file = "Data/Sources/3'UTRs/3'UTRs.tab", sep = "\t", row.names = F)