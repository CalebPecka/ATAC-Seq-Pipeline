# Set your working directory to the current file directory
# This command can only be used while working in an RStudio environment.
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

######################
# Input Files
######################
args <- commandArgs(trailingOnly = T) # End line command arguments.

genes <- read.delim(args[1], sep = "\t")
summits <- read.delim(args[2], header = F) # Output from MACS2

# List of chromosomes.
chrGroups <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
               "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
               "chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM")

# Booleans to indicate whether or not we want to include the positive 
#and negative strands.
includePositive <- T
includeNegative <- T



##########################################################################
# Error check
# Determines if inputs and outputs are giving expected ranges of results

# Checks the gene file is BED formatted.
if (all(colnames(genes) != c("GeneName", "Chromosome", "Start", "Stop", 
                             "Symbol", "Strand"))){
  # Error output message.
  stop(paste0("Input BED file for reference genes not properly formatted.",
    "Make sure the file contains 6 columns labeled (GeneName, Chromosome, ",
    "Start, Stop, Symbol, Strand). Elements should be separated with a single ",
    "space."))
  
  # Script failed to terminate.
  print("Script did not end! Please manually terminate.")
}

# Checks the summits file is formatted correctly.
if (length(summits) != 10){
  stop("Input Narrow Peak Summits file not properly formatted. ",
    "Make sure the file contains 10 columns, NO header, and is tab separated.")
  print("Script did not end! Please manually terminate.")
}

# Runs a simulation of the code using a data subset (first 100 rows of GTFgenes)
# If any errors exist, they should be caught early.
outFile <- c()
genesTemp <- genes[which(genes$Chromosome == chrGroups[1]),]
summitsTempIntermediate <- summits[which(summits$V1 == chrGroups[1]),]
summitsTemp <- summitsTempIntermediate$V2 + summitsTempIntermediate$V10 
summitsQVals <- summitsTempIntermediate$V9
summitScore <- summitsTempIntermediate$V7

for (j in row.names(genesTemp)[1:100]){
  rowOfInterest <- genesTemp[j,]
  low <- rowOfInterest$Start - 1000
  high <- rowOfInterest$Start - 100
  
  conditional <- which(summitsTemp >= low & summitsTemp <= high)
  summitsList <- summitsTemp[conditional]
  summitsQValList <- summitsQVals[conditional]
  summitScoreList <- summitScore[conditional]
  
  if (length(summitsList) != 0){
    for (k in 1:length(summitsList)){
      outFile <- append(outFile, paste(c("chr1", rowOfInterest$ENST, 
                                         summitsList[k], 
                                         summitsQValList[k], 
                                         rowOfInterest$Symbol, 
                                         rowOfInterest$Strand, 
                                         summitScoreList[k]), 
                                       collapse = " "))
    }
  }
}

# Prints a sample of the output file to be confirmed by the user.
print("#############################################")
print("SAMPLE upstreamPeaks.tsv Outfile")
print("Check that the output is formatted correctly.")
print("#############################################")
print(head(outFile))

testOutFile <- unlist(strsplit(outFile[1], ' '))
if (is.null(testOutFile[3])){
  stop(paste0("Test output returned NULL for chromosome start position. ",
    "Check input. Terminating command. It is possible this error occurred ",
    "because no peaks were found upstream to first 100 genes in your ",
    "reference gene list."))
  print("Script did not end! Please manually terminate.")
}
if (is.null(testOutFile[4])){
  stop(paste0("Test output returned NULL for chromosome end position. ",
  "Check input. Terminating command. It is possible this error occurred ",
  "because no peaks were found upstream to first 100 genes in your ",
  "reference gene list."))
  print("Script did not end! Please manually terminate.")
}
##########################################################################



# Main function. The output file variable is reset.
outFile <- c()

# Run the code for genes on the positive strand.
if (includePositive)
{
  genesStrand <- genes[which(genes$Strand == "+"),]
  
  # The goal of this function is to locate every instance of an accessible 
  # chromatin region that is within a 1000 base upstream region of genes in the 
  # human genome. To do this, we use the GTFgenes.tsv file to find start 
  # locations of genes, and find summits from the MACS2 output that are within 
  # our desired range of the gene. Each chromosome is calculated individually 
  # to avoid errors in which multiple genes on different chromosomes have the 
  # same base pair (bp) start position.
  for (chromosome in chrGroups){
    
    # Subset genes on the matching chromosome.
    genesTemp <- genesStrand[which(genesStrand$Chromosome == chromosome),]
    # Subset of peaks(summits) on the matching chromosome.
    summitsTempIntermediate <- summits[which(summits$V1 == chromosome),]
    # Column 10 (V10) of the NA_peaks file includes the bp distance from the 
    # start of the accessible chromatin (V2) to the peak. To calculate the 
    # summit's bp position, we can add the values from these columns together.
    summitsTemp <- summitsTempIntermediate$V2 + summitsTempIntermediate$V10
    # The values in column 9 (V9) include a significance score metric for the 
    # peak. It is saved in memory for future programs.
    summitsQVals <- summitsTempIntermediate$V9
    # The values in column 7 (V7) indicate the MACS2 score for a relative fold 
    # change of the peak.
    summitScore <- summitsTempIntermediate$V7
    
    # For each gene, we perform a caluclation to find summit locations.
    for (geneID in row.names(genesTemp)){
      
      rowOfInterest <- genesTemp[geneID,]
      # Here we calculate the range of acceptable bp locations for a summit to 
      # be located near a gene. Summits in this region of more likely to contain 
      # transcription factor binding sites (TFBS).
      low <- rowOfInterest$Start - 1000 # Lower bound bp region.
      high <- rowOfInterest$Start - 100 # Upper bound bp region. TFBSs are 
      # unlikely to be found within the 100bp upstream region of genes.
      
      # Conditional finds all summits within the accepted bp region.
      conditional <- which(summitsTemp >= low & summitsTemp <= high)
      # Afterwards, our previous variables are subset to include only the values 
      # in conditional.
      summitsList <- summitsTemp[conditional]
      summitsQValList <- summitsQVals[conditional]
      summitScoreList <- summitScore[conditional]
      
      # Checks that the table is not empty.
      if (length(summitsList) != 0){
        for (k in 1:length(summitsList)){
          # An output file is created that stores important variables delimited 
          # by spaces.
          outFile <- append(outFile, paste(c(chromosome, rowOfInterest$ENST, 
                                             summitsList[k], 
                                             summitsQValList[k], 
                                             rowOfInterest$Symbol, 
                                             rowOfInterest$Strand,
                                             summitScoreList[k]), 
                                           collapse = "\t"))
        }
      }
    }
    
    # The chromosome variable is printed upon completion so the user has 
    # an easier time gauging the progress of the program.
    print(paste0(chromosome, " positive strand search completed."))
  }
}



# Run the above code for the negative strand.
if (includeNegative)
{
  genesStrand <- genes[which(genes$Strand == "-"),]
  
  for (chromosome in chrGroups){
    
    genesTemp <- genesStrand[which(genesStrand$Chromosome == chromosome),]
    summitsTempIntermediate <- summits[which(summits$V1 == chromosome),]
    summitsTemp <- summitsTempIntermediate$V2 + summitsTempIntermediate$V10
    summitsQVals <- summitsTempIntermediate$V9
    summitScore <- summitsTempIntermediate$V7
    
    for (geneID in row.names(genesTemp)){
      
      rowOfInterest <- genesTemp[geneID,]
      # Here we calculate the range of acceptable bp locations for a summit to 
      # be located near a gene. Summits in this region of more likely to contain 
      # transcription factor binding sites (TFBS).
      low <- rowOfInterest$Stop + 100 # Lower bound bp region.
      high <- rowOfInterest$Stop + 1000 # Upper bound bp region. TFBSs are 
      # unlikely to be found within the 100bp upstream region of genes.
      
      conditional <- which(summitsTemp >= low & summitsTemp <= high)
      summitsList <- summitsTemp[conditional]
      summitsQValList <- summitsQVals[conditional]
      summitScoreList <- summitScore[conditional]
    
      if (length(summitsList) != 0){
        for (k in 1:length(summitsList)){
          outFile <- append(outFile, paste(c(chromosome, rowOfInterest$ENST, 
                                             summitsList[k], 
                                             summitsQValList[k], 
                                             rowOfInterest$Symbol, 
                                             rowOfInterest$Strand,
                                             summitScoreList[k]), 
                                           collapse = "\t"))
        }
      }
    }
    
    print(paste0(chromosome, " negative strand search completed."))
  }
}



# The output of the program is saved.
outFile <- data.frame(outFile)
outFile <- unique(outFile)
write.table(outFile, args[3], 
            row.names = F, 
            col.names = F, 
            quote = F)

