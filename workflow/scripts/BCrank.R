# Set your working directory to the current file directory 
# This command can only be used while working in an RStudio environment.
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

######################
# Input Files
######################
library(BCRANK)
library(Biostrings)

args <- commandArgs(trailingOnly = T) # End line command arguments.

# BCrank is used to discover 100 unique motifs in this program (restartSize 
# variable). To repeat this process, you can increase the seedSize variable. 
# The program is repeated using a new random seed a number of times equal to 
# the seed size. For example, a seedSize value of 12 will produce (12 * 100) = 
# 1200 unique motifs.
seedSize <- as.numeric(args[1])
restartSize <- as.numeric(args[2])

if (seedSize > 20){
  print("WARNING: Seed size is greater than 20! This is accepted by BCrank,
        but future steps of the program will not be able to read more than 20
        seeds. All remaining seed generations will be ignored by those steps.")
}



######################
# Program Start
######################
for (seed in 1:seedSize){
  set.seed(as.numeric(seed)) # Set the randomized seed variable.
  BCRANKout <- bcrank(args[3], 
                      restarts = restartSize, 
                      use.P1 = TRUE, 
                      use.P2 = TRUE)
  
  BCRANKS <- unlist(toptable(BCRANKout))
  
  # The output from BCRANK is saved to a new file.
  motifList <- ""
  for (i in 1:length(BCRANKout@toplist)){
    motifList <- paste0(motifList, BCRANKout@toplist[[i]]@final@consensus)
    motifList <- paste0(motifList, "\n")
  }
  write.table(data.frame(motifList), paste0(args[4], seed, "_BC_consensusSequences.txt"), col.names = F, row.names = F, quote = F)
    
  # We are also interested in comparing these motifs to known motif datasets 
  # using Memesuite tools. Memesuite requires a unique file format. The BCrank 
  # motifs are converted to a MEME file format in this section of code.
  outFile <- c("MEME version 5.1.1", "", "ALPHABET= ACGT", "") # Header
  
  for (i in 1:length(BCRANKout@toplist)){
    outFile <- append(outFile, paste("MOTIF", i, collapse = " "))
    outFile <- append(outFile, "letter-probability matrix: alength= 4")
    
    # Motifs are formatted as a nucleotide probability matrix (PWM). This "for" 
    # loop cycles through each of the probability matrices created by BCrank.
    for (PWM in 1:ncol(BCRANKout@toplist[[i]]@finalPWM)){
      outFile <- append(outFile, 
                        paste(t(BCRANKout@toplist[[i]]@finalPWM / BCRANKout@toplist[[i]]@finalNrMatch)[PWM,], 
                              collapse = " "))
    }
    
    # The newly formatted PWM is appended to a growing output file.
    outFile <- append(outFile, "")
  }  
  
  finalOut <- data.frame(outFile) # Convert the Meme output to a dataframe 
                                  # file type.
  
  # Write output.
  write.table(finalOut, 
              paste0(args[4], seed, "_BCRANK_motifs.meme"), 
              row.names = F, col.names = F, quote = F)
}
