# Set your working directory to the current file directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

######################
# Input Files
######################
library(BCRANK)

library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86

# All of the seeds are separated in the output from BCrank.R. The goal of this program is to 
# merge the outputs from BCrank into a single datafile. The consensus sequences also have a
# chance of being incompatible with a future program. This program tests each BCrank consensus
# sequence to remove any problematic data.

# If you including more or less than 12 seeds in your BCrank motif search, you will have to
# edit the input variables below to account for any changes.
BC_1 <- read.table("BCrankOutput/1_BC_consensusSequences.txt")
BC_2 <- read.table("BCrankOutput/2_BC_consensusSequences.txt")
BC_3 <- read.table("BCrankOutput/3_BC_consensusSequences.txt")
BC_4 <- read.table("BCrankOutput/4_BC_consensusSequences.txt")
BC_5 <- read.table("BCrankOutput/5_BC_consensusSequences.txt")
BC_6 <- read.table("BCrankOutput/6_BC_consensusSequences.txt")
BC_7 <- read.table("BCrankOutput/7_BC_consensusSequences.txt")
BC_8 <- read.table("BCrankOutput/8_BC_consensusSequences.txt")
BC_9 <- read.table("BCrankOutput/9_BC_consensusSequences.txt")
BC_10 <- read.table("BCrankOutput/10_BC_consensusSequences.txt")
BC_11 <- read.table("BCrankOutput/11_BC_consensusSequences.txt")
BC_12 <- read.table("BCrankOutput/12_BC_consensusSequences.txt")

# MAIN FUNCTION
# "dfnames" keeps track of the BCrank motif objects in short term memory.
dfnames <- c("BC_1", "BC_2", "BC_3", "BC_4", "BC_5", "BC_6", 
             "BC_7", "BC_8", "BC_9", "BC_10", "BC_11", "BC_12")
dfENV <- mget(dfnames, .GlobalEnv)

# Creates a function that extracts consensus sequences from the BCrank files. Every other line is
# a sequence score and is removed. 
dataExtraction <- function(df) {
  w <- df[c(FALSE, TRUE),]
  w <- w[c(1:25),]
  w <- unlist(w)
}

# Applies the previous function to all dataframes in the environment.
extracted <- lapply(dfENV, dataExtraction)

fullConsensusList <- c() # Instantiates a list of consensus sequences.
for (i in dfnames){
  
  # Goes through each dataframe in the environment and appends the data to the full consensus list.
  fullConsensusList <- append(fullConsensusList, extracted[[i]])
}
# When complete, fullConsensusList represents a complete list of the motifs identified in all
# trials produced by BCrank.


############################################
# FASTA Sequence Conversion
############################################
# This section of code tests and removes invalid consensus sequences. For naming conventions,
# the "fullConsensusList" variable is renamed to "bindingSiteSearchList".
bindingSiteSearchList <- fullConsensusList

# Instantiate lists.
consensusSequence <- c()
chromosome <- c()
geneName <- c()
start <- c()
end <- c()
count <- 1
doNotUse <- c()

# Goes through every binding site consensus sequence.
while(count < length(bindingSiteSearchList) + 1){
  
  # Occaisionally the binding site fails in the "matchingSites" command. The exact reason is
  # unknown, and the error is rare. The function below will find instances of each binding site
  # throughout the original fasta sequence file, upstreamPeak_sequences.fasta. If an error is
  # produced, the binding site is ignored.
  tryCatch(
    {
      # Finds locations in upstreamPeak_sequences.fasta where the binding site exists.
      bindingSiteSearch <- as.character(bindingSiteSearchList[count])
      matching <- matchingSites("upstreamPeak_Sequences.fasta", bindingSiteSearch)
      
      # The matchingsites() function assumes the input includes either forward or reverse
      # strands of DNA. Here we subset our results to only include the forward strands.
      matching <- matching[which(matching$Strand == "+"),] 
      intermediate <- unlist(strsplit(matching$`Region header`, ':')) # Splits metadata
      
      # Now we need to identify the genomic location of each matching site. This data was stored 
      # in the upstreamPeak_sequences.fasta file header in the SummitFASTA.R script.
      # Here, we extract chromosomal location.
      chrExtraction <- intermediate[seq(1, length(intermediate), 3)]
      geneExtraction <- intermediate[seq(3, length(intermediate), 3)]
      
      # Remove pesky hidden end-line character and ENST decimals.
      geneExtraction <- substr(geneExtraction,1,nchar(geneExtraction) - 3)
      
      # Finds the start and end location of the motif in each sequence.
      locationExtraction <- intermediate[seq(2, length(intermediate), 3)]
      locationIntermediate <- unlist(strsplit(locationExtraction, '-'))
      startPosition <- as.numeric(locationIntermediate[seq(1, length(locationIntermediate), 2)]) + 
        as.numeric(matching$Start)
      endPosition <- as.numeric(locationIntermediate[seq(1, length(locationIntermediate), 2)]) + 
        as.numeric(matching$End)
      
      # The chromosomal location and FASTA sequence are appended into several vector lists.
      consensusSequence <- append(consensusSequence, rep(bindingSiteSearch, length(chrExtraction)))
      chromosome <- append(chromosome, chrExtraction)
      geneName <- append(geneName, geneExtraction)
      start <- append(start, startPosition)
      end <- append(end, endPosition)
    },
    
    # If an error occurs during the previous process, a message tells the user which consensus
    # sequence failed, and the reason why.
    error=function(cond) {
      print("")
      print(paste("WARNING: The consensus sequence", bindingSiteSearch, "produced the following error:"))
      print(cond)
      print(paste(bindingSiteSearch, "has been removed from the binding site search list."))
      
      # doNotUse keeps track of which consensus sequences failed.
      doNotUse <- append(doNotUse, count)
    },
    
    # Whether the program failed or not, the count meter increases by 1, meaning we move on to the
    # next consensus sequence.
    finally = {
      count <- count + 1
    }
  )
}



############################################
# Dataframe Conversion
############################################
bindingSiteSearchList <- bindingSiteSearchList[-c(doNotUse)] # List of usable binding sites.

# Converts collected vector lists into a single data frame.
out <- cbind(data.frame(consensusSequence), 
             data.frame(chromosome), 
             data.frame(geneName), 
             data.frame(start), 
             data.frame(end))

# Data frame rows are sorted based on chromosome number.
out <- out[order(out$chromosome),]

# Previously, we chose a motif size of 100 nucleotides. This value was chosen to meet the 
# minimum default requirements for a few MEME motif database search tools. Now, we want to ensure
# that every motif is within a region defined as a "peak", RATHER than including any matching site
# NEARBY to a matching site.
narrowPeaks <- read.delim("NA_peaks.narrowPeak", sep = "\t", header = F)

uniqueChr <- unique(out$chromosome)
hasPeak <- c()

# To reduce the chances of errors, each chromosome is handled individually.
for (chr in uniqueChr){
  outSubset <- out[which(out$chromosome == chr),]
  narrowPeakSubset <- narrowPeaks[which(narrowPeaks$V1 == chr),]
  
  # For each gene in each chromsome...
  for(i in 1:nrow(outSubset)) {
    row <- outSubset[i,]
    
    # Finds the start position (V2) and end position (V3) for each narrow peak.
    narrowPeakTemp <- narrowPeakSubset[which(narrowPeakSubset$V2 <= row$end),]
    narrowPeakTemp <- narrowPeakTemp[which(narrowPeakTemp$V3 >= row$start),]
    
    # A new dataframe row is created that keeps track of whether or not the matching site is
    # found within the same region as a peak.
    if (nrow(narrowPeakTemp) >= 1){
      hasPeak <- append(hasPeak, TRUE)
    }
    else {
      hasPeak <- append(hasPeak, FALSE)
    }
  }
}

# Appends the hasPeak variable to the dataframe.
out$hasPeak <- hasPeak
out <- out[which(out$hasPeak == T),] # Subsets the data to only include hasPeak values of true.

write.csv(out, "matchingSiteTable.csv")

