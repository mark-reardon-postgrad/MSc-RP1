# Question 2 - Are miRNA genes more or less likely than protein-coding genes to overlap copy number variants?

require(GenomicRanges)
setwd("C:/Users/mark/iCloudDrive/MSc/BIOL61230 Research Projects/1/Analysis")

# load pri-miRNAs, protein-coding genes and CNVs
source(file = "Data/Loaders/pri-miRNAs.r")
source(file = "Data/Loaders/Protein-coding genes.r")
source(file = "Data/Loaders/Decipher CNVs.r")
source(file = "Data/Loaders/DGV CNVs.r")

# adds p-value text to a plot just above the proportion bars if significantly different
p.value.ifSignificant <- function(x, prop.test.result, textExpansion = 3)
{
  threshold <- 0.05
  heightAdjust <- 0.03
  if (prop.test.result$p.value < threshold)
  {
    y <- max(prop.test.result$estimate) + heightAdjust
    text(x, y, sprintf("p = %1.2e", prop.test.result$p.value), cex = textExpansion)
  }
}

# emits a plot comparing the proportions of primiRNAs and protein-coding genes that overlap each class of CNVs
CompareMiRNAsToProteinCodingGenes <- function(pane, showLegend, title, method, legendTitles, primiRNAs, proteinCodingGenes, allCNVs, gainedCNVs, lostCNVs, textExpansion = 3)
{
  # miRNA overlaps
  primiRNAsGained <- subsetByOverlaps(query = primiRNAs, subject = gainedCNVs)
  primiRNAsLost <- subsetByOverlaps(query = primiRNAs, subject = lostCNVs)
  primiRNAsNone <- subset(primiRNAs, !(primiRNAs@elementMetadata$attributes %in% 
                                       subsetByOverlaps(query = primiRNAs, subject = allCNVs)@elementMetadata$attributes))
  # protein-coding gene overlaps
  proteinCodingGenesGained <- subsetByOverlaps(query = proteinCodingGenes, subject = gainedCNVs)
  proteinCodingGenesLost <- subsetByOverlaps(query = proteinCodingGenes, subject = lostCNVs)
  proteinCodingGenesNone <- subset(proteinCodingGenes, !(proteinCodingGenes@elementMetadata$ensembl_gene_id %in% 
                                                         subsetByOverlaps(query = proteinCodingGenes, subject = allCNVs)@elementMetadata$ensembl_gene_id))
  # proportions
  proportionGained <- prop.test(c(length(primiRNAsGained), length(proteinCodingGenesGained)), c(length(primiRNAs), length(proteinCodingGenes)))
  proportionLost <- prop.test(c(length(primiRNAsLost), length(proteinCodingGenesLost)), c(length(primiRNAs), length(proteinCodingGenes)))
  proportionNone <- prop.test(c(length(primiRNAsNone), length(proteinCodingGenesNone)), c(length(primiRNAs), length(proteinCodingGenes)))
  
  # output statistucs
  cat(capture.output(proportionGained), file = sprintf("Results/2/%s (%s) gained.txt", title, method), sep = "\n")
  cat(capture.output(proportionLost), file = sprintf("Results/2/%s (%s) lost.txt", title, method), sep = "\n")
  cat(capture.output(proportionNone), file = sprintf("Results/2/%s (%s) none.txt", title, method), sep = "\n")
  
  # output plot
  geneColours <- c("darkgreen", "darkorange")
  barplot(
    matrix(c(proportionGained$estimate, proportionLost$estimate, proportionNone$estimate), ncol = 3),
    beside = T, space = c(0,1), ylim = 0:1,
    names.arg = c("Gain", "Loss", "None"), 
    col = geneColours,
    main = pane, 
    ylab = "Proportion overlapping", xlab = "CNV type",
    cex.axis = textExpansion, cex.names = textExpansion, cex.main = textExpansion + 1, cex.lab = textExpansion)
  if (showLegend)
  {
    legend("topleft", legend = legendTitles, fill = geneColours, bty = "n", cex = textExpansion)
  }
  p.value.ifSignificant(2, proportionGained, textExpansion)
  p.value.ifSignificant(5, proportionLost, textExpansion)
  p.value.ifSignificant(8, proportionNone, textExpansion)
}

# emits a plot comparing the proportions of genes that are gained or lost in the CNVs
CompareGainsToLosses <- function(title, primiRNAs, proteinCodingGenes, gainedCNVs, lostCNVs)
{
  # miRNA overlaps
  primiRNAsGained <- subsetByOverlaps(query = primiRNAs, subject = gainedCNVs)
  primiRNAsLost <- subsetByOverlaps(query = primiRNAs, subject = lostCNVs)
  # protein-coding gene overlaps
  proteinCodingGenesGained <- subsetByOverlaps(query = proteinCodingGenes, subject = gainedCNVs)
  proteinCodingGenesLost <- subsetByOverlaps(query = proteinCodingGenes, subject = lostCNVs)
  # proportions
  proportionMiRNAs <- prop.test(c(length(primiRNAsGained), length(primiRNAsLost)), c(length(primiRNAs), length(primiRNAs)))
  proportionProteinCodingGenes <- prop.test(c(length(proteinCodingGenesGained), length(proteinCodingGenesLost)), c(length(proteinCodingGenes), length(proteinCodingGenes)))
  # output
  cat(capture.output(proportionMiRNAs), file = sprintf("Results/2/%s miRNAs.txt", title), sep = "\n")
  cat(capture.output(proportionProteinCodingGenes), file = sprintf("Results/2/%s protein-coding genes.txt", title), sep = "\n")
}

# emits a plot comparing the proportions of genes that are gained or lost in the CNVs
CompareCNVs <- function(title, elements, decipherGains, decipherLosses, dgvGains, dgvLosses)
{
  # Decipher
  elementsDecipherGained <- subsetByOverlaps(query = elements, subject = decipherGains)
  elementsDecipherLost <- subsetByOverlaps(query = elements, subject = decipherLosses)
  # DGV
  elementsDgvGained <- subsetByOverlaps(query = elements, subject = dgvGains)
  elementsDgvLost <- subsetByOverlaps(query = elements, subject = dgvLosses)
  # proportions
  proportionGained <- prop.test(c(length(elementsDecipherGained), length(elementsDgvGained)), c(length(elements), length(elements)))
  proportionLost <- prop.test(c(length(elementsDecipherLost), length(elementsDgvLost)), c(length(elements), length(elements)))
  # output
  cat(capture.output(proportionGained), file = sprintf("Results/2/%s gained.txt", title), sep = "\n")
  cat(capture.output(proportionLost), file = sprintf("Results/2/%s lost.txt", title), sep = "\n")
}

# 2A - Is there a significant difference in the proportions of each class of gene overlapping gained, lost, gained and lost, 
#      any or no Decipher (pathogenic) or Database of Genomic Variation (healthy) copy number variants?
jpeg(filename = "Results/2/2AB.jpg", width = 3000, height = 2000)
par(mfrow = c(2,2), mar = c(5, 5, 5, 5))
CompareMiRNAsToProteinCodingGenes("A", T, "Pathogenic CNVs (Decipher)", "width-based", c("pri-miRNAs", "Protein-coding genes"), primiRNAs, proteinCodingGenes, decipherCNVs, decipherCNVsGained, decipherCNVsLost)
CompareMiRNAsToProteinCodingGenes("B", F, "Healthy CNVs (DGV)", "width-based", c("pri-miRNAs", "Protein-coding genes"), primiRNAs, proteinCodingGenes, dgvCNVs, dgvCNVsGained, dgvCNVsLost)

# 2B - Does the result hold when only the midpoints of pri-miRNAs and protein-coding genes are considered (in case the 
#      disparity in length between pri-miRNAs and protein-coding genes of around three orders of magnitude skews the result)?
CompareMiRNAsToProteinCodingGenes("C", F, "Pathogenic CNVs (Decipher)", "midpoint-based", c("pri-miRNA midpoints", "Protein-coding gene midpoints"), primiRNAMidpoints, proteinCodingGeneMidpoints, decipherCNVs, decipherCNVsGained, decipherCNVsLost)
CompareMiRNAsToProteinCodingGenes("D", F, "Healthy CNVs (DGV)", "midpoint-based", c("pri-miRNA midpoints", "Protein-coding gene midpoints"), primiRNAMidpoints, proteinCodingGeneMidpoints, dgvCNVs, dgvCNVsGained, dgvCNVsLost)
dev.off()

# 2C - Is there a significant difference in the number of each of miRNAs and protein-coding genes overlapping gained or lost 
#      copy number variants?
CompareGainsToLosses("Pathogenic CNVs (Decipher)", primiRNAMidpoints, proteinCodingGeneMidpoints, decipherCNVsGained, decipherCNVsLost)
CompareGainsToLosses("Healthy CNVs (DGV)", primiRNAMidpoints, proteinCodingGeneMidpoints, dgvCNVsGained, dgvCNVsLost)

# 2D - Is there a significant difference in the proportions of miRNAs/protein-coding genes overlapping gained/lost Decipher or DGV CNVs
CompareCNVs("miRNAs", primiRNAMidpoints, decipherCNVsGained, decipherCNVsLost, dgvCNVsGained, dgvCNVsLost)
CompareCNVs("Protein-coding genes", proteinCodingGeneMidpoints, decipherCNVsGained, decipherCNVsLost, dgvCNVsGained, dgvCNVsLost)

# different pane titles for presentation

# 2A - Is there a significant difference in the proportions of each class of gene overlapping gained, lost, gained and lost, 
#      any or no Decipher (pathogenic) or Database of Genomic Variation (healthy) copy number variants?
jpeg(filename = "../Presentation/Overlap comparisons.jpg", width = 3000, height = 2000)
par(mfrow = c(2,2), mar = c(5, 5, 5, 5))
CompareMiRNAsToProteinCodingGenes("Pathogenic (widths)", T, "Pathogenic CNVs (Decipher)", "width-based", c("pri-miRNAs", "Protein-coding genes"), primiRNAs, proteinCodingGenes, decipherCNVs, decipherCNVsGained, decipherCNVsLost)
CompareMiRNAsToProteinCodingGenes("Healthy (widths)", F, "Healthy CNVs (DGV)", "width-based", c("pri-miRNAs", "Protein-coding genes"), primiRNAs, proteinCodingGenes, dgvCNVs, dgvCNVsGained, dgvCNVsLost)

# 2B - Does the result hold when only the midpoints of pri-miRNAs and protein-coding genes are considered (in case the 
#      disparity in length between pri-miRNAs and protein-coding genes of around three orders of magnitude skews the result)?
CompareMiRNAsToProteinCodingGenes("Pathogenic (midpoints)", F, "Pathogenic CNVs (Decipher)", "midpoint-based", c("pri-miRNA midpoints", "Protein-coding gene midpoints"), primiRNAMidpoints, proteinCodingGeneMidpoints, decipherCNVs, decipherCNVsGained, decipherCNVsLost)
CompareMiRNAsToProteinCodingGenes("Healthy (midpoints)", F, "Healthy CNVs (DGV)", "midpoint-based", c("pri-miRNA midpoints", "Protein-coding gene midpoints"), primiRNAMidpoints, proteinCodingGeneMidpoints, dgvCNVs, dgvCNVsGained, dgvCNVsLost)
dev.off()


# bigger text for poster
marginSize <- 21
textExpansion = 9
jpeg(filename = "../Poster/Overlap comparisons.jpg", width = 6000, height = 4000)
par(mfrow = c(2,2), mar = c(marginSize, marginSize + 12, marginSize, marginSize - 12), mgp = c(15, 6, 0))
CompareMiRNAsToProteinCodingGenes("Pathogenic (widths)", T, "Pathogenic CNVs (Decipher)", "width-based", c("pri-miRNAs", "Protein-coding genes"), primiRNAs, proteinCodingGenes, decipherCNVs, decipherCNVsGained, decipherCNVsLost, textExpansion)
CompareMiRNAsToProteinCodingGenes("Healthy (widths)", F, "Healthy CNVs (DGV)", "width-based", c("pri-miRNAs", "Protein-coding genes"), primiRNAs, proteinCodingGenes, dgvCNVs, dgvCNVsGained, dgvCNVsLost, textExpansion)

CompareMiRNAsToProteinCodingGenes("Pathogenic (midpoints)", F, "Pathogenic CNVs (Decipher)", "midpoint-based", c("pri-miRNA midpoints", "Protein-coding gene midpoints"), primiRNAMidpoints, proteinCodingGeneMidpoints, decipherCNVs, decipherCNVsGained, decipherCNVsLost, textExpansion)
CompareMiRNAsToProteinCodingGenes("Healthy (midpoints)", F, "Healthy CNVs (DGV)", "midpoint-based", c("pri-miRNA midpoints", "Protein-coding gene midpoints"), primiRNAMidpoints, proteinCodingGeneMidpoints, dgvCNVs, dgvCNVsGained, dgvCNVsLost, textExpansion)
dev.off()
