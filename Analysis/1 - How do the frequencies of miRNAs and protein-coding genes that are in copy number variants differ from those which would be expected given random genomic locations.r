# Question 1 - How do the frequencies of miRNAs and protein-coding genes that are in copy number variants differ from those which would be 
#              expected given random genomic locations?

require(GenomicRanges)
require(parallel)
require(snow)

setwd("C:/Users/mark/iCloudDrive/MSc/BIOL61230 Research Projects/1/Analysis")

# load miRNAs, protein-coding genes and CNVs
source(file = "Data/Loaders/pri-miRNAs.r")
source(file = "Data/Loaders/Protein-coding genes.r")
source(file = "Data/Loaders/Chromosomes.r")
source(file = "Data/Loaders/Decipher CNVs.r")
source(file = "Data/Loaders/DGV CNVs.r")

# returns significance asterisks for a P value
significance <- function(p) 
{
  result <- ""
  if (!is.na(p) && p <= 0.1) result <- "."
  if (!is.na(p) && p <= 0.05) result <- "*"
  if (!is.na(p) && p <= 0.01) result <- "**"
  if (!is.na(p) && p <= 0.001) result <- "***"
  return(result)
}

# returns a data frame contain statistics about the simulations
summarise <- function(observations, simulations)
{
  # calculate averages and z-based p values
  simulation.mean <- apply(simulations, 2, function(column) mean(column))
  simulation.median <- apply(simulations, 2, function(column) median(column))
  simulation.sd <- apply(simulations, 2, function(column) sd(column))
  z <- (observations - simulation.mean) / simulation.sd
  p <- sapply(z, function(overlapZ) if (!is.na(overlapZ) && overlapZ < 0) return(pnorm(overlapZ)) else return(1 - pnorm(overlapZ)))
  
  # count how many simulations each overlap ocurred in and how many of each simulation are less than, equal to or more than the observed
  simulation.in <- numeric(length(observations))
  less <- numeric(length(observations))
  equal <- numeric(length(observations))
  more <- numeric(length(observations))
  for (overlap in 1:length(observations))
  {
    simulation.in[overlap] <- sum(simulations[, overlap] > 0)
    less[overlap] <- sum(simulations[, overlap] < observations[overlap])
    equal[overlap] <- sum(simulations[, overlap] == observations[overlap])
    more[overlap] <- sum(simulations[, overlap] > observations[overlap])
  }
  
  # calculate significance based on less than and more than counts
  experimentCount <- length(simulations[ ,1])
  less.div <- less / experimentCount
  more.div <- more / experimentCount
  less.sig <- sapply(less / experimentCount, significance)
  more.sig <- sapply(more / experimentCount, significance)
  
  # return the summary
  return (data.frame(
    overlaps = 0:(length(simulations[1,]) - 1),
    observations, simulation.in, 
    simulation.mean, simulation.median, simulation.sd, z, p, p.sig = sapply(p, significance),
    less, less.div, less.sig,
    equal, 
    more, more.div, more.sig))
}

# writes the output files if filenames were specified
writeOutput <- function(simulations, simulationsFilename, summary, summaryFilename)
{
  # simulations
  if (!is.null(simulationsFilename))
  {
    write.table(simulations, file = simulationsFilename, sep = "\t", row.names = F, col.names = F)
  }
  
  # summary
  if (!is.null(summaryFilename))
  {
    write.table(summary, file = summaryFilename, sep = "\t", row.names = F)
  }
}

# creates a random variant of the query and returns a matrix containing its overlaps with the subject as a frequency table
randomiseOverlaps <- function(query, subject, chromosomes)
{
  require(GenomicRanges)
  # create vectors of randomly chosen chromosomes and positions within each
  randomChromosomes <- sample(seqnames(query), size = length(query), replace = T)
  randomStarts <- sapply(randomChromosomes, function(chr) sample.int(chromosomes$Length[chromosomes$Name == chr], 1))
  
  # create GRange holding random query-like ranges (lengths the same; chromosome, starts and strands random)
  randomRanges <- GRanges(
    seqnames = randomChromosomes,
    ranges = IRanges(
      start = randomStarts,
      width = width(query),
      names = names(query)),
    strand = sample(c("+", "-"), size = length(query), replace = T))
  
  # overlap the random ranges with the subject ranges and return as a frequency table
  randomOverlaps <- countOverlaps(query = randomRanges, subject = subject)
  return(matrix(data = table(factor(randomOverlaps, levels = min(randomOverlaps):max(randomOverlaps))), nrow = 1))
}

# performs the randomisation simulations
#  (intended to be run on cluster nodes and so requires variables 'observations, query, subject, chromosomes, simulationCount'
#   and function 'randomiseOverlaps' to be exported to the nodes)
performSimulations <- function()
{
  # create matrix to hold randomised simulations
  simulations <- matrix(data = 0, nrow = simulationCount, ncol = length(observations))
  
  # for each simulation...
  for(simulation in 1:simulationCount)
  {
    # randomise the query and summarise as a frequency table
    randomObservations <- randomiseOverlaps(query, subject, chromosomes)
    
    # add simulation columns if more are needed
    while(length(randomObservations) > length(simulations[1,]))
    {
      simulations <- cbind(simulations, rep.int(0, simulationCount))
    }
    
    # update the simulations
    for(randomObservation in 1:length(randomObservations))
    {
      simulations[simulation, randomObservation] <- randomObservations[1, randomObservation]
    }
  }
  
  # return the simulations
  return(simulations)
}

# overlaps the query over the subject and compares to the distribution obtained from random repositioning of the query
performExperiment <- function(pane, showLegend, queryName, subjectName, colour, xmax = NULL, textExpansion = 3) 
{
  # get the variables for the experiment
  query <- get(queryName)
  subject <- get(subjectName)
  chromosomes <- get("chromosomes")
  simulationCount <- get("simulationCount")
  totalSimulationCount <- simulationCount
  
  # abort if the simulation count is not divisible by the number of cores to simplify partitioning
  coreCount <- detectCores()
  if (simulationCount %% coreCount != 0)
  {
    stop(sprintf("Simulation count (%d) is not divisible by the number of cores (%d).", simulationCount, coreCount))
  }
  
  # overlap the query with the subject and summarise as frequency table
  overlaps <- countOverlaps(query, subject)
  observations <- matrix(data = table(factor(overlaps, levels = min(overlaps):max(overlaps))), nrow = 1)
  
  # if the summary does not already exist...
  experimentFilename <- sprintf("Results/1/Experiments/%s.%s.%d.tab", queryName, subjectName, totalSimulationCount)
  summaryFilename <- sprintf("Results/1/Summaries/%s.%s.%d.tab", queryName, subjectName, totalSimulationCount)
  if (!file.exists(summaryFilename))
  {
    
    # make cluster with maximum number of cores
    cluster <- makeCluster(coreCount)
    
    # partition the simulations
    simulationCount <- simulationCount / coreCount
    
    # export data, parameters and required functions to each node
    clusterExport(cluster, "observations", envir = environment())
    clusterExport(cluster, "query", envir = environment())
    clusterExport(cluster, "subject", envir = environment())
    clusterExport(cluster, "chromosomes", envir = environment())
    clusterExport(cluster, "simulationCount", envir = environment())
    clusterExport(cluster, "randomiseOverlaps", envir = environment())
    
    # perform the simulations in parallel
    nodeSimulations <- clusterCall(cl = cluster, fun = performSimulations)
  
    # stop the cluster
    stopCluster(cluster)
    
    # create matrix to hold simulations from all the nodes
    maxColumns <- 0
    for (node in 1:coreCount)
    {
      nodeSimulation <- nodeSimulations[[node]]
      nodeColumns <- length(nodeSimulation[1, ])
      if (nodeColumns > maxColumns)
      {
        maxColumns <- nodeColumns
      }
    }
    simulations <- matrix(data = 0, nrow = simulationCount * coreCount, ncol = maxColumns)
    
    # copy the simulations from each node's results
    simulationIndex <- 1
    for (node in 1:coreCount)
    {
      nodeSimulation <- nodeSimulations[[node]]
      for (nodeSimulationIndex in 1:simulationCount)
      {
        for (observation in 1:length(nodeSimulation[1, ]))
        {
          simulations[simulationIndex, observation] <- nodeSimulation[nodeSimulationIndex, observation]
        }
        simulationIndex <- simulationIndex + 1
      }
    }
    
    # adjust observation length to match simulations in case there are more values in simulations
    priorObservations <- length(observations)
    if (length(observations) < length(simulations[1,]))
    {
      length(observations) <- length(simulations[1,])
      # zero the new observations
      for (observation in priorObservations:length(observations))
      {
        observations[observation] <- 0
      }
    }
    
    # summarise the simulations
    summary <- summarise(observations, simulations)
    
    # output raw experimental data and results
    writeOutput(simulations, experimentFilename, summary, summaryFilename)
  }
  else
  {
    summary <- read.table(file = summaryFilename, header = T, sep = "\t")
  }
  
  # zero missing values    
  summary$observations[!is.numeric(summary$observations)] <- 0
  summary$simulation.mean[!is.numeric(summary$simulation.mean)] <- 0

  # plot axis limits
  ymax = max(c(summary$observations, summary$simulation.mean)) * 1.1
  xlim <- NULL
  barSpacesPerFrequency <- 3
  if (!is.null(xmax))
  {
    xlim <- c(0, (xmax * barSpacesPerFrequency) - 1)
  }
  # output the plot
  barplot(t(matrix(c(summary$observations, summary$simulation.mean), ncol = 2)), 
          beside = T, space = c(0,1), 
          main = pane,
          names.arg = 0:(length(summary$observations) - 1),
          xlab = "Overlap count", ylab = "Frequency",
          col = c(colour, "lightgrey"),
          cex.axis = textExpansion, cex.names = textExpansion, cex.main = textExpansion + 1, cex.lab = textExpansion,
          ylim = c(0, ymax),
          xlim = xlim)
  if (showLegend)
  {
    legend("topright", legend = c("Observed", "Simulated"), fill = c(colour, "lightgrey"), bty = "n", cex = textExpansion)
  }
  
  # output significance arrows
  threshold <- 0.05
  for (i in 1:length(summary$p))
  {
    p <- summary$p[i] 
    z <- summary$z[i]
    if (!is.na(p) & !is.na(z) & p < threshold)
    {
      sig <- ifelse(z < 0, "<", ">")
      y <- max(c(summary$observations[i], summary$simulation.mean[i]))
      text((3 * i) - 1, y, sig, cex = textExpansion, pos = 3)
    }
  }
}

# number of randomised simulations
simulationCount <- 1000

marginSize <- 6
# miRNA gains
jpeg(filename = sprintf("Results/1/miRNAs.%d.jpg", simulationCount), width = 3000, height = 2000)
par(mfrow = c(2,2), mar = c(marginSize, marginSize, marginSize, marginSize))
performExperiment("A", T, "primiRNAMidpoints", "decipherCNVsGained", "darkgreen")
performExperiment("B", F, "primiRNAMidpoints", "dgvCNVsGained", "darkgreen")

# miRNA losses
performExperiment("C", F, "primiRNAMidpoints", "decipherCNVsLost", "darkgreen")
performExperiment("D", F, "primiRNAMidpoints", "dgvCNVsLost", "darkgreen")
dev.off()

# protein-coding gene gains
jpeg(filename = sprintf("Results/1/Protein-coding genes.%d.jpg", simulationCount), width = 3000, height = 2000)
par(mfrow = c(2,2), mar = c(marginSize, marginSize, marginSize, marginSize))
performExperiment("A", T, "proteinCodingGeneMidpoints", "decipherCNVsGained", "darkorange")
performExperiment("B", F, "proteinCodingGeneMidpoints", "dgvCNVsGained", "darkorange")

# protein-coding gene losses
performExperiment("C", F, "proteinCodingGeneMidpoints", "decipherCNVsLost", "darkorange")
performExperiment("D", F, "proteinCodingGeneMidpoints", "dgvCNVsLost", "darkorange")
dev.off()

# restricted x-axis for presentation
# miRNA gains
jpeg(filename = sprintf("../Presentation/miRNAs.%d.jpg", simulationCount), width = 3000, height = 2000)
par(mfrow = c(2,2), mar = c(marginSize, marginSize, marginSize, marginSize))
performExperiment("Pathogenic gains", T, "primiRNAMidpoints", "decipherCNVsGained", "darkgreen", 10)
performExperiment("Healthy gains", F, "primiRNAMidpoints", "dgvCNVsGained", "darkgreen", 10)

# miRNA losses
performExperiment("Pathogenic losses", F, "primiRNAMidpoints", "decipherCNVsLost", "darkgreen", 10)
performExperiment("Healthy losses", F, "primiRNAMidpoints", "dgvCNVsLost", "darkgreen", 10)
dev.off()

# protein-coding gene gains
jpeg(filename = sprintf("../Presentation/Protein-coding genes.%d.jpg", simulationCount), width = 3000, height = 2000)
par(mfrow = c(2,2), mar = c(marginSize, marginSize, marginSize, marginSize))
performExperiment("Pathogenic gains", T, "proteinCodingGeneMidpoints", "decipherCNVsGained", "darkorange", 10)
performExperiment("Healthy gains", F, "proteinCodingGeneMidpoints", "dgvCNVsGained", "darkorange", 10)

# protein-coding gene losses
performExperiment("Pathogenic losses", F, "proteinCodingGeneMidpoints", "decipherCNVsLost", "darkorange", 10)
performExperiment("Healthy losses", F, "proteinCodingGeneMidpoints", "dgvCNVsLost", "darkorange", 10)
dev.off()

# bigger text for poster
marginSize <- 20
xmax = 10
textExpansion = 9
# protein-coding gene gains
jpeg(filename = sprintf("../Poster/Protein-coding genes.%d.jpg", simulationCount), width = 6000, height = 4000)
par(mfrow = c(2,2), mar = c(marginSize, marginSize, marginSize, marginSize + 4), mgp = c(13, 6, 0))
performExperiment("Pathogenic gains", T, "proteinCodingGeneMidpoints", "decipherCNVsGained", "darkorange", xmax, textExpansion)
performExperiment("Healthy gains", F, "proteinCodingGeneMidpoints", "dgvCNVsGained", "darkorange", xmax, textExpansion)

# protein-coding gene losses
performExperiment("Pathogenic losses", F, "proteinCodingGeneMidpoints", "decipherCNVsLost", "darkorange", xmax, textExpansion)
performExperiment("Healthy losses", F, "proteinCodingGeneMidpoints", "dgvCNVsLost", "darkorange", xmax, textExpansion)
dev.off()

