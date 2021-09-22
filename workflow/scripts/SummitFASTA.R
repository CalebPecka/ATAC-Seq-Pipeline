# Set your working directory to the current file directory 
# This command can only be used while working in an RStudio environment.
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

######################
# Input Files
######################
library(Biostrings)

args <- commandArgs(trailingOnly = T) # End line command arguments.

HG38 <- args[1]

genome <- readDNAStringSet(HG38, format = "fasta")
outFile <- read.delim(args[2], sep = "\t", header = F)

# This variable establishes the size of the motif binding site region you'd 
# like to use. If you don't know what that means, read further.
MotifSize <- as.numeric(args[3])



##########################
##SEQUENCE DETERMINATION
##########################
# In this program, we'd like to find the DNA FASTA sequence around each 
# NarrowPeak summit location determined in the previous step 
# (NarrowPeakSummitTracker.R). To do this, we need to import and initialize 
# a reference genome that is compatible with the same gene position file 
# (GTFgenes). If you modify either GFFgenes.bed OR the files in HG38, you will
# need to verify the data in the genes.bed files matches the correct bp 
# locations in the HG38 file.
chr1 <- genome[1]; chr2 <- genome[2]; chr3 <- genome[3]; chr4 <- genome[4]
chr5 <- genome[5]; chr6 <- genome[6]; chr7 <- genome[7]; chr8 <- genome[8]
chr9 <- genome[9]; chr10 <- genome[10]; chr11 <- genome[11]
chr12 <- genome[12]; chr13 <- genome[13]; chr14 <- genome[14]
chr15 <- genome[15]; chr16 <- genome[16]; chr17 <- genome[17]
chr18 <- genome[18]; chr19 <- genome[19]; chr20 <- genome[20]
chr21 <- genome[21]; chr22 <- genome[22]; chrX <- genome[23]
chrY <- genome[24]; chrM <- genome[25]

# File formatting initialization.
motifFile <- c()
qVal <- c()
MotifSize <- round(MotifSize / 2)

# To make the FASTA file, we convert each acessible chromatin region one at 
# a time.
for (site in 1:nrow(outFile)){
  
  # Each row of the data is separated by spaces.
  # Subject[1] is the Chromosome.
  # Subject[2] is the GeneName.
  # Subject[3] is the SummitLocation.
  # Subject[4] is the QVal of the original MACS2 Summit.
  subject <- outFile[site,]
  
  # This segment of code identifies the DNA region surrounding each summit.
  # The size of the search space is determined by the variable "MotifSize", 
  # where the value of MotifSize is the length of the final motif in units 
  # of bp. Each chromsome is handled in a separate switch condition to prevent 
  # undesired bp start site duplications.
  startPos = as.numeric(as.character(subject[[3]])) - MotifSize
  endPos = as.numeric(as.character(subject[[3]])) + MotifSize

  # We must also account for niche scenarios in which the positions exceed
  # the size of the reference genome, due to the +- 1000 bp search range
  # from the NarrowPeakSummitTracker.R script.
  ogStartPos = startPos
  startPos = max(1, startPos)
  endPos = max(2, endPos)
  genomeSize <- switch (as.character(subject[[1]]),
                 "chr1" = length(chr1[[1]]),
                 "chr2" = length(chr2[[1]]),
                 "chr3" = length(chr3[[1]]),
                 "chr4" = length(chr4[[1]]),
                 "chr5" = length(chr5[[1]]),
                 "chr6" = length(chr6[[1]]),
                 "chr7" = length(chr7[[1]]),
                 "chr8" = length(chr8[[1]]),
                 "chr9" = length(chr9[[1]]),
                 "chr10" = length(chr10[[1]]),
                 "chr11" = length(chr11[[1]]),
                 "chr12" = length(chr12[[1]]),
                 "chr13" = length(chr13[[1]]),
                 "chr14" = length(chr14[[1]]),
                 "chr15" = length(chr15[[1]]),
                 "chr16" = length(chr16[[1]]),
                 "chr17" = length(chr17[[1]]),
                 "chr18" = length(chr18[[1]]),
                 "chr19" = length(chr19[[1]]),
                 "chr20" = length(chr20[[1]]),
                 "chr21" = length(chr21[[1]]),
                 "chr22" = length(chr22[[1]]),
                 "chrX" = length(chrX[[1]]),
                 "chrY" = length(chrY[[1]]),
                 "chrM" = length(chrM[[1]])
  )
  startPos = min(genomeSize - 1, startPos)
  endPos = min(genomeSize, endPos)
  
  # Finally, we collect the sequence based on the chromosome.
  seq <- switch (as.character(subject[[1]]),
                 "chr1" = BStringSet(chr1, start = startPos, end = endPos),
                 "chr2" = BStringSet(chr2, start = startPos, end = endPos),
                 "chr3" = BStringSet(chr3, start = startPos, end = endPos),
                 "chr4" = BStringSet(chr4, start = startPos, end = endPos),
                 "chr5" = BStringSet(chr5, start = startPos, end = endPos),
                 "chr6" = BStringSet(chr6, start = startPos, end = endPos),
                 "chr7" = BStringSet(chr7, start = startPos, end = endPos),
                 "chr8" = BStringSet(chr8, start = startPos, end = endPos),
                 "chr9" = BStringSet(chr9, start = startPos, end = endPos),
                 "chr10" = BStringSet(chr10, start = startPos, end = endPos),
                 "chr11" = BStringSet(chr11, start = startPos, end = endPos),
                 "chr12" = BStringSet(chr12, start = startPos, end = endPos),
                 "chr13" = BStringSet(chr13, start = startPos, end = endPos),
                 "chr14" = BStringSet(chr14, start = startPos, end = endPos),
                 "chr15" = BStringSet(chr15, start = startPos, end = endPos),
                 "chr16" = BStringSet(chr16, start = startPos, end = endPos),
                 "chr17" = BStringSet(chr17, start = startPos, end = endPos),
                 "chr18" = BStringSet(chr18, start = startPos, end = endPos),
                 "chr19" = BStringSet(chr19, start = startPos, end = endPos),
                 "chr20" = BStringSet(chr20, start = startPos, end = endPos),
                 "chr21" = BStringSet(chr21, start = startPos, end = endPos),
                 "chr22" = BStringSet(chr22, start = startPos, end = endPos),
                 "chrX" = BStringSet(chrX, start = startPos, end = endPos),
                 "chrY" = BStringSet(chrY, start = startPos, end = endPos),
                 "chrM" = BStringSet(chrM, start = startPos, end = endPos)
  )
  
  # All of the data is restructured into a FASTA file with a modified header.
  # The header is a colon separated file included the chromosomal location of 
  # each FASTA sequence and the name of the gene associated with it.
  newHeader <- paste(c(">", as.character(subject[[1]]), ":",
                       as.character(ogStartPos), "-",
                       as.character(endPos), ":",
                       as.character(subject[[5]]), ":",
                       as.character(subject[[6]]), ":",
                       subject[[7]]), collapse = '')

  # To prevent duplicates, we check if the header ID is already in the list,
  # and only keep FASTA sequences that are NOT repeats.
  if (newHeader %in% motifFile == F){
    motifFile <- append(motifFile, newHeader)
    # After the header has been added, we add the DNA sequence.
    motifFile <- append(motifFile, as.character(seq))
  
    # We keep track of the qValues so we can reorder the data later on.
    qVal <- append(qVal, c(as.numeric(as.character(subject[[4]])), as.numeric(as.character(subject[[4]]))))
  }
}

secondOut <- data.frame(motifFile) # Convert data to dataframe file format.
orderByQVAL <- secondOut[rev(order(as.numeric(qVal))),] # Order the FASTA 
# sequences based on QVal. This step is important for accurate results when 
# we eventually use the program "BCrank".
even <- seq(from=2, to=length(orderByQVAL), by=2)
odd <- even - 1
finalOutput <- orderByQVAL
finalOutput[even] <- orderByQVAL[odd]
finalOutput[odd] <- orderByQVAL[even]
finalOutput <- data.frame(finalOutput)

# Write the output FASTA file.
write.table(finalOutput, args[4], 
            row.names = F, col.names = F, quote = F)

