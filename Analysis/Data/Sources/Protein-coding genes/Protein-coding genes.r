#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")

library(biomaRt)

setwd("C:/Users/mark/iCloudDrive/MSc/BIOL61230 Research Projects/1/Analysis")

# GRCh38 human data source
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

# get all protein-coding genes on chromosomes 1..22 with protein-coding transcripts
proteinCodingGenes <- getBM(mart = ensembl, 
  attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "transcript_count"),
  filters = c("biotype", "transcript_biotype", "chromosome_name"), 
  values = list("protein_coding", "protein_coding", 1:22))

# rename columns for GRanges compatibility
colnames(proteinCodingGenes) <- c("ensembl_gene_id", "chr", "start", "end", "strand", "transcript_count")

# rename chromosomes and strands
proteinCodingGenes <- transform(proteinCodingGenes, chr = sprintf("chr%s", chr))
proteinCodingGenes$strand[proteinCodingGenes$strand == "1"] <- "+"
proteinCodingGenes$strand[proteinCodingGenes$strand == "-1"] <- "-"

# save file
write.table(proteinCodingGenes, file = "Data/Sources/Protein-coding genes/Protein-coding genes.tab", sep = "\t", row.names = F)
