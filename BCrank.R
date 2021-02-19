# Set your working directory to the current file directory 
# This command can only be used while working in an RStudio environment.
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

######################
# Input Files
######################
library(BCRANK)
library(Biostrings)
library(seqinr)

# BCrank is used to discover 100 unique motifs in this program (restartSie variable). To repeat this 
# process, you can increase the seedSize variable. The program is repeated using a new random seed a 
# number of times equal to the seed size. For example, a seedSize value of 12 will produce 
# (12 * 100) = 1200 unique motifs.
seedSize <- 12
restartSize <- 100



######################
# Program Start
######################
for (seed in 1:seedSize){
  set.seed(as.numeric(seed)) # Set the randomized seed variable.
  BCRANKout <- bcrank("upstreamPeak_Sequences.fasta", restarts = restartSize, use.P1 = TRUE, use.P2 = TRUE)
  
  BCRANKS <- unlist(toptable(BCRANKout))
  
  # The output from BCRANK is sent to the console. The printed output from the console is saved
  # to a new file.
  sink(paste0("BCrankOutput/", seed,"_BC_consensusSequences.txt"))
  print(BCRANKS)
  sink()
  
  # We are also interested in comparing these motifs to known motif datasets using Memesuite tools.
  # Memesuite requires a unique file format. The BCrank motifs are converted to a MEME file format
  # in this section of code.
  outFile <- c("MEME version 5.1.1", "", "ALPHABET= ACGT", "") # Header
  
  for (i in 1:length(BCRANKout@toplist)){
    outFile <- append(outFile, paste("MOTIF", i, collapse = " "))
    outFile <- append(outFile, "letter-probability matrix: alength= 4")
    
    # Motifs are formatted as a nucleotide probability matrix (PWM). This "for" loop cycles through each of
    # the probability matrices created by BCrank.
    for (PWM in 1:ncol(BCRANKout@toplist[[i]]@finalPWM)){
      outFile <- append(outFile, 
                        paste(t(BCRANKout@toplist[[i]]@finalPWM / BCRANKout@toplist[[i]]@finalNrMatch)[PWM,], 
                              collapse = " "))
    }
    
    # The newly formatted PWM is appended to a growing output file.
    outFile <- append(outFile, "")
  }  
  
  finalOut <- data.frame(outFile) # Convert the Meme output to a dataframe file type.
  
  # Write output.
  write.table(finalOut, 
              paste0("BCrankOutput/", seed,"_BCRANK_motifs.meme"), 
              row.names = F, col.names = F, quote = F)
}
