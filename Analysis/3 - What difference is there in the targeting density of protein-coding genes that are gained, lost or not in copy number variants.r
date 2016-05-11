# Question 3 - What difference is there in the 3'UTR targeting density of protein-coding genes that are gained, lost or not in copy number variants?

require(GenomicRanges)

setwd("C:/Users/mark/iCloudDrive/MSc/BIOL61230 Research Projects/1/Analysis")

# load protein-coding genes and CNVs
source(file = "Data/Loaders/Protein-coding genes.r")
source(file = "Data/Loaders/Decipher CNVs.r")
source(file = "Data/Loaders/DGV CNVs.r")

# count the targets on each 3'UTR
if (!exists("targetCounts"))
{
  if (!file.exists("Results/3/targetCounts.tab"))
  {
    source(file = "Data/Loaders/Targets.r")
    targetCounts <- as.data.frame(tapply(targets$miRNA, targets$EnsembleTranscriptID, length))
    targetCounts$EnsembleTranscriptID <- rownames(targetCounts)
    rownames(targetCounts) <- NULL
    colnames(targetCounts)[1] <- "TargetCount"
    write.table(targetCounts, file = "Results/3/targetCounts.tab", sep = "\t", row.names = F)
  }
  else
  {
    targetCounts <- read.table(file = "Results/3/targetCounts.tab", header = T, sep = "\t")
  }
}

# calculate the targeting density of each 3'UTR
if (!exists("targetedUTRs"))
{
  if (!file.exists("Results/3/targetedUTRs.tab"))
  {
    # load the UTRs and left join the target counts
    source(file = "Data/Loaders/3'UTRs.r")
    targetedUTRs <- merge(UTRs, targetCounts, by.x = "EnsembleTranscriptID", by.y = "EnsembleTranscriptID", all.x = T)
    targetedUTRs$TargetCount[is.na(targetedUTRs$TargetCount)] <- 0
    # remove UTRs less than 50bp as probably sequencing artefacts
    targetedUTRs$SequenceLength <- nchar(as.character(targetedUTRs$Sequence))
    targetedUTRs <- targetedUTRs[targetedUTRs$SequenceLength >= 50, ]
    # calculate UTR target densities
    targetedUTRs$TargetingDensity <- targetedUTRs$TargetCount / targetedUTRs$SequenceLength
    # cache the results
    write.table(targetedUTRs, file = "Results/3/targetedUTRs.tab", sep = "\t", row.names = F)
  }
  else
  {
    targetedUTRs <- read.table(file = "Results/3/targetedUTRs.tab", header = T, sep = "\t")
  }
}

# average the target densities for each gene
if (!exists("meanGeneTargetDensities"))
{
  if (!file.exists("Results/3/meanGeneTargetDensities.tab"))
  {
    meanGeneTargetDensities <- data.frame(
      MeanTargetingDensity = tapply(targetedUTRs$TargetingDensity, targetedUTRs$EnsembleGeneID, mean))
    meanGeneTargetDensities$EnsembleGeneID <- rownames(meanGeneTargetDensities)
    rownames(meanGeneTargetDensities) <- NULL
    write.table(meanGeneTargetDensities, file = "Results/3/meanGeneTargetDensities.tab", sep = "\t", row.names = F)
  }
  else
  {
    meanGeneTargetDensities <- read.table(file = "Results/3/meanGeneTargetDensities.tab", header = T, sep = "\t")
  }
}

# analyses protein-coding gene midpoint overlaps with CNVs with respect to their targeting densities
AnalyseTargetingDensities <- function(pane, showLegend, prefix, title, proteinCodingGeneMidpoints, meanGeneTargetDensities, CNVsGained, CNVsLost, CNVs, textExpansion = 2, lineWidth = 3)
{
  # partition the genes by overlap of their midpoints with CNVs
  proteinCodingGenesGained <- subsetByOverlaps(query = proteinCodingGeneMidpoints, subject = CNVsGained)
  proteinCodingGenesLost <- subsetByOverlaps(query = proteinCodingGeneMidpoints, subject = CNVsLost)
  proteinCodingGenesNone <- subset(proteinCodingGenes, !(proteinCodingGenes@elementMetadata$ensembl_gene_id %in% 
                                                         subsetByOverlaps(query = proteinCodingGenes, subject = CNVs)@elementMetadata$ensembl_gene_id))
  
  # specify the order of the CNV levels to control the intercept choice in the later ANOVA step
  cnvLevels = c("neither gained nor lost", "gained", "lost")
  
  # add the mean targeting densities to the genes neither gained nor lost (ordered by genomic location)
  neitherGainedNorLost <- merge(
    data.frame(EnsembleGeneID = proteinCodingGenesNone@elementMetadata$ensembl_gene_id[
      order(as.numeric(substr(seqnames(proteinCodingGenesNone), 4, 6)), start(proteinCodingGenesNone))]),
    meanGeneTargetDensities,
    by.x = "EnsembleGeneID", by.y = "EnsembleGeneID", all.x = T)
  neitherGainedNorLost$MeanTargetingDensity[is.na(neitherGainedNorLost$MeanTargetingDensity)] <- 0
  neitherGainedNorLost$CNV <- factor("neither gained nor lost", levels = cnvLevels)
  neitherGainedNorLost$Colour <- "black"
  
  # add the mean targeting densities to the gained genes (ordered by genomic location)
  gained <- merge(
    data.frame(EnsembleGeneID = proteinCodingGenesGained@elementMetadata$ensembl_gene_id[
      order(as.numeric(substr(seqnames(proteinCodingGenesGained), 4, 6)), start(proteinCodingGenesGained))]),
    meanGeneTargetDensities,
    by.x = "EnsembleGeneID", by.y = "EnsembleGeneID", all.x = T)
  gained$MeanTargetingDensity[is.na(gained$MeanTargetingDensity)] <- 0
  gained$CNV <- factor("gained", levels = cnvLevels)
  gained$Colour <- "blue"
  
  # add the mean targeting densities to the lost genes (ordered by genomic location)
  lost <- merge(
    data.frame(EnsembleGeneID = proteinCodingGenesLost@elementMetadata$ensembl_gene_id[
      order(as.numeric(substr(seqnames(proteinCodingGenesLost), 4, 6)), start(proteinCodingGenesLost))]),
    meanGeneTargetDensities,
    by.x = "EnsembleGeneID", by.y = "EnsembleGeneID", all.x = T)
  lost$MeanTargetingDensity[is.na(lost$MeanTargetingDensity)] <- 0
  lost$CNV <- factor("lost", levels = cnvLevels)
  lost$Colour <- "red"
  
  # combined densities
  densities <- rbind(neitherGainedNorLost, gained, lost)
  write.table(densities, file = sprintf("Results/3/%s densities.tab", prefix), sep = "\t", row.names = F)
  
  # values and means
  plot(
    1:nrow(densities), densities$MeanTargetingDensity, 
    main = pane,
    xlab = "CNV type/genomic location", ylab = "Mean 3'UTR targeting density",
    bg = densities$Colour, pch = 21, lwd = lineWidth, col = densities$Colour,
    cex.axis = textExpansion, cex.main = textExpansion + 1, cex.lab = textExpansion)
  abline(h = mean(densities$MeanTargetingDensity[densities$CNV == "neither gained nor lost"]), col = "black", lwd = lineWidth)
  abline(h = mean(densities$MeanTargetingDensity[densities$CNV == "gained"]), col = "blue", lwd = lineWidth)
  abline(h = mean(densities$MeanTargetingDensity[densities$CNV == "lost"]), col = "red", lwd = lineWidth)
  if (showLegend)
  {
    legend("topright", legend = c("Neither gained nor lost", "Gained", "Lost"), fill = c("black", "blue", "red"), bty = "n", cex = textExpansion)
  }
    
  # anova
  model <- aov(densities$MeanTargetingDensity ~ densities$CNV)
  filename <- sprintf("Results/3/%s.txt", prefix)
  cat("summary.lm\n==========\n", file = filename, sep = "\n")
  cat(capture.output(summary.lm(model)), file = filename, sep = "\n", append = T)
  cat("\nTukeyHSD\n========\n", file = filename, sep = "\n", append = T)
  cat(capture.output(TukeyHSD(model)), file = filename, sep = "\n", append = T)
  cat("\nmean\n====\n", file = filename, sep = "\n", append = T)
  cat(capture.output(tapply(densities$MeanTargetingDensity, densities$CNV, mean)), file = filename, sep = "\n", append = T)
  cat("\nsd\n==\n", file = filename, sep = "\n", append = T)
  cat(capture.output(tapply(densities$MeanTargetingDensity, densities$CNV, sd)), file = filename, sep = "\n", append = T)
  cat("\nlength\n======\n", file = filename, sep = "\n", append = T)
  cat(capture.output(tapply(densities$MeanTargetingDensity, densities$CNV, length)), file = filename, sep = "\n", append = T)
}

jpeg(filename = "Results/3/plots.jpg", width = 2000, height = 1000)
par(mfcol = c(1,2), mar = c(5, 5, 5, 5))
AnalyseTargetingDensities("A", T, "Decipher", "Pathogenic CNVs (Decipher)",
                          proteinCodingGeneMidpoints, meanGeneTargetDensities,
                          decipherCNVsGained, decipherCNVsLost, decipherCNVs)

AnalyseTargetingDensities("B", F, "DGV", "Healthy CNVs (DGV)",
                          proteinCodingGeneMidpoints, meanGeneTargetDensities,
                          dgvCNVsGained, dgvCNVsLost, dgvCNVs)
dev.off()

# presentation versions
jpeg(filename = "../Presentation/Densities.jpg", width = 2000, height = 1000)
par(mfcol = c(1,2), mar = c(5, 5, 5, 5))
AnalyseTargetingDensities("Pathogenic CNVs", T, "Decipher", "Pathogenic CNVs (Decipher)",
                          proteinCodingGeneMidpoints, meanGeneTargetDensities,
                          decipherCNVsGained, decipherCNVsLost, decipherCNVs)

AnalyseTargetingDensities("Healthy CNVs", F, "DGV", "Healthy CNVs (DGV)",
                          proteinCodingGeneMidpoints, meanGeneTargetDensities,
                          dgvCNVsGained, dgvCNVsLost, dgvCNVs)
dev.off()

# bigger text for poster
marginSize <- 21
textExpansion = 9
lineWidth = 6
jpeg(filename = "../Poster/Densities.jpg", width = 6000, height = 3000)
par(mfrow = c(1,2), mar = c(marginSize, marginSize + 12, marginSize, marginSize - 12), mgp = c(15, 6, 0))
AnalyseTargetingDensities("Pathogenic CNVs", T, "Decipher", "Pathogenic CNVs (Decipher)",
                          proteinCodingGeneMidpoints, meanGeneTargetDensities,
                          decipherCNVsGained, decipherCNVsLost, decipherCNVs,
                          textExpansion, lineWidth)

AnalyseTargetingDensities("Healthy CNVs", F, "DGV", "Healthy CNVs (DGV)",
                          proteinCodingGeneMidpoints, meanGeneTargetDensities,
                          dgvCNVsGained, dgvCNVsLost, dgvCNVs,
                          textExpansion, lineWidth)
dev.off()