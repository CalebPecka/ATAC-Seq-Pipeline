# Set your working directory to the current file directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

######################
# Input Files
######################
library(Biostrings)

genome <- readDNAStringSet("HG38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz", format = "fasta")
outFile <- read.delim("upstreamPeaks.tsv", sep = " ", col.names = c("Chromsome", 
                                                                    "GeneName", 
                                                                    "SummitLocation", 
                                                                    "Qval"))

# This variable establishes the size of the motif binding site region you'd like to use.
# If you don't know what that means, read further.
MotifSize <- 100



##########################
##SEQUENCE DETERMINATION
##########################
# In this program, we'd like to find the DNA FASTA sequence around each NarrowPeak summit location
# determined in the previous step (NarrowPeakSummitTracker.R). To do this, we need to import and
# instantialize a reference genome that is compatible with the same gene position file (GFFgenes).
# If you modify either GFFgenes.bed OR the files in HG38, you will need to verify the data in
# the genes.bed files matches the correct bp locations in the HG38 file.
chr1 <- genome[1]; chr2 <- genome[2]; chr3 <- genome[3]; chr4 <- genome[4]; chr5 <- genome[5]
chr6 <- genome[6]; chr7 <- genome[7]; chr8 <- genome[8]; chr9 <- genome[9]; chr10 <- genome[10]
chr11 <- genome[11]; chr12 <- genome[12]; chr13 <- genome[13]; chr14 <- genome[14]
chr15 <- genome[15]; chr16 <- genome[16]; chr17 <- genome[17]; chr18 <- genome[18]
chr19 <- genome[19]; chr20 <- genome[20]; chr21 <- genome[21]; chr22 <- genome[22]
chrX <- genome[23]; chrY <- genome[24]; chrM <- genome[25]

# File formatting initialization.
motifFile <- c()
qVal <- c()
MotifSize <- round(MotifSize / 2)

# To make the FASTA file, we convert each acessible chromatin region one at a time.
for (site in row(outFile)){
  
  # Each row of the data is separated by spaces.
  # Subject[1] is the Chromosome.
  # Subject[2] is the GeneName.
  # Subject[3] is the SummitLocation.
  # Subject[4] is the QVal of the original MACS2 Summit.
  subject <- unlist(strsplit(outFile[site,], ' '))
  
  # This segment of code identifies the DNA region surrounding each summit. The size of the search space
  # is determined by the variable "MotifSize", where the value of MotifSize is the length of the final
  # motif in units of bp. Each chromsome is handled in a separate switch condition to prevent undesired
  # bp start site duplications.
  seq <- switch (subject[1],
                 "chr1" = BStringSet(chr1, start = as.numeric(subject[3]) - MotifSize, 
                                     end = as.numeric(subject[3]) + MotifSize),
                 "chr2" = BStringSet(chr2, start = as.numeric(subject[3]) - MotifSize, 
                                     end = as.numeric(subject[3]) + MotifSize),
                 "chr3" = BStringSet(chr3, start = as.numeric(subject[3]) - MotifSize, 
                                     end = as.numeric(subject[3]) + MotifSize),
                 "chr4" = BStringSet(chr4, start = as.numeric(subject[3]) - MotifSize, 
                                     end = as.numeric(subject[3]) + MotifSize),
                 "chr5" = BStringSet(chr5, start = as.numeric(subject[3]) - MotifSize, 
                                     end = as.numeric(subject[3]) + MotifSize),
                 "chr6" = BStringSet(chr6, start = as.numeric(subject[3]) - MotifSize, 
                                     end = as.numeric(subject[3]) + MotifSize),
                 "chr7" = BStringSet(chr7, start = as.numeric(subject[3]) - MotifSize, 
                                     end = as.numeric(subject[3]) + MotifSize),
                 "chr8" = BStringSet(chr8, start = as.numeric(subject[3]) - MotifSize, 
                                     end = as.numeric(subject[3]) + MotifSize),
                 "chr9" = BStringSet(chr9, start = as.numeric(subject[3]) - MotifSize, 
                                     end = as.numeric(subject[3]) + MotifSize),
                 "chr10" = BStringSet(chr10, start = as.numeric(subject[3]) - MotifSize, 
                                      end = as.numeric(subject[3]) + MotifSize),
                 "chr11" = BStringSet(chr11, start = as.numeric(subject[3]) - MotifSize, 
                                      end = as.numeric(subject[3]) + MotifSize),
                 "chr12" = BStringSet(chr12, start = as.numeric(subject[3]) - MotifSize, 
                                      end = as.numeric(subject[3]) + MotifSize),
                 "chr13" = BStringSet(chr13, start = as.numeric(subject[3]) - MotifSize, 
                                      end = as.numeric(subject[3]) + MotifSize),
                 "chr14" = BStringSet(chr14, start = as.numeric(subject[3]) - MotifSize, 
                                      end = as.numeric(subject[3]) + MotifSize),
                 "chr15" = BStringSet(chr15, start = as.numeric(subject[3]) - MotifSize, 
                                      end = as.numeric(subject[3]) + MotifSize),
                 "chr16" = BStringSet(chr16, start = as.numeric(subject[3]) - MotifSize, 
                                      end = as.numeric(subject[3]) + MotifSize),
                 "chr17" = BStringSet(chr17, start = as.numeric(subject[3]) - MotifSize, 
                                      end = as.numeric(subject[3]) + MotifSize),
                 "chr18" = BStringSet(chr18, start = as.numeric(subject[3]) - MotifSize, 
                                      end = as.numeric(subject[3]) + MotifSize),
                 "chr19" = BStringSet(chr19, start = as.numeric(subject[3]) - MotifSize, 
                                      end = as.numeric(subject[3]) + MotifSize),
                 "chr20" = BStringSet(chr20, start = as.numeric(subject[3]) - MotifSize, 
                                      end = as.numeric(subject[3]) + MotifSize),
                 "chr21" = BStringSet(chr21, start = as.numeric(subject[3]) - MotifSize, 
                                      end = as.numeric(subject[3]) + MotifSize),
                 "chr22" = BStringSet(chr22, start = as.numeric(subject[3]) - MotifSize, 
                                      end = as.numeric(subject[3]) + MotifSize),
                 "chrX" = BStringSet(chrX, start = as.numeric(subject[3]) - MotifSize, 
                                     end = as.numeric(subject[3]) + MotifSize),
                 "chrY" = BStringSet(chrY, start = as.numeric(subject[3]) - MotifSize, 
                                     end = as.numeric(subject[3]) + MotifSize),
                 "chrM" = BStringSet(chrM, start = as.numeric(subject[3]) - MotifSize, 
                                     end = as.numeric(subject[3]) + MotifSize),
  )
  
  # All of the data is restructured into a FASTA file with a modified header.
  # The header is a colon separated file included the chromosomal location of each FASTA sequence and
  # the name of the gene associated with it.
  motifFile <- append(motifFile, paste(c(">", subject[1], ":",
                                         as.character(as.numeric(subject[3]) - MotifSize), "-",
                                         as.character(as.numeric(subject[3]) + MotifSize), ":",
                                         subject[2], "\t"), collapse = ''))
  # After the header has been added, we add the DNA sequence.
  motifFile <- append(motifFile, as.character(seq))
  
  # We keep track of the qValues so we can reorder the data later on.
  qVal <- append(qVal, c(subject[4], subject[4]))
}

secondOut <- data.frame(motifFile) # Convert data to dataframe file format.
secondOut <- data.frame(secondOut[order(as.numeric(qVal)),]) # Order the FASTA sequences based on QVal.
# This step is important for accurate results when we eventually use the program "BCrank".

# Write the output FASTA file.
write.table(secondOut, "upstreamPeak_Sequences.fasta", row.names = F, col.names = F, quote = F)

