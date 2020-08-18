###### INPUT ARGUEMENTS#########
args <- commandArgs(trailingOnly = T)

inputSummits <- args[1]
inputGenes <- args[2]
outdir <- args[3]

summits <- read.delim(inputSummits, header = F)

genes <- read.delim(inputGenes, header = F)

outFile <- c()

chrGroups <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM")

for (i in chrGroups){
  #EVERYTHING IS THE SAME CHROMOSOME REGION
  genesTemp <- genes[which(genes$V1 == i),]
  summitsTemp <- summits[which(summits$V1 == i),]$V2
  
  for (j in row.names(genesTemp)){
    #FINDS THE START AND STOP BASE FOR EACH GENE
    rowOfInterest <- genesTemp[j,]
    low <- rowOfInterest$V2 - 1000
    high <- rowOfInterest$V2
    
    summitsList <- summitsTemp[which(summitsTemp >= low & summitsTemp <= high)]
    
    if (length(summitsList) != 0){
      for (k in summitsList){
        outFile <- append(outFile, paste(c(i, rowOfInterest$V4, k), collapse = " "))
      }
    }
  }
  
  print(i)
  
}

outFile <- data.frame(outFile)
write.table(outFile, paste(outdir, "upstremPeaks.tsv", sep = ""), row.names = F)

##########################
##SEQUENCE DETERMINATION
##########################
library(Biostrings)

genome <- readDNAStringSet("HG38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz", format = "fasta")

chr1 <- genome[1]
chr2 <- genome[2]
chr3 <- genome[3]
chr4 <- genome[4]
chr5 <- genome[5]
chr6 <- genome[6]
chr7 <- genome[7]
chr8 <- genome[8]
chr9 <- genome[9]
chr10 <- genome[10]
chr11 <- genome[11]
chr12 <- genome[12]
chr13 <- genome[13]
chr14 <- genome[14]
chr15 <- genome[15]
chr16 <- genome[16]
chr17 <- genome[17]
chr18 <- genome[18]
chr19 <- genome[19]
chr20 <- genome[20]
chr21 <- genome[21]
chr22 <- genome[22]
chrX <- genome[23]
chrY <- genome[24]
chrM <- genome[25]

MotifSize <- 45

sequences <- c()

for (i in row(outFile)){
  
  subject <- unlist(strsplit(outFile[i,], ' '))
  
  seq <- switch (subject[1],
                 "chr1" = BStringSet(chr1, start = as.numeric(subject[4]) - MotifSize, end = as.numeric(subject[4]) + MotifSize),
                 "chr2" = BStringSet(chr2, start = as.numeric(subject[4]) - MotifSize, end = as.numeric(subject[4]) + MotifSize),
                 "chr3" = BStringSet(chr3, start = as.numeric(subject[4]) - MotifSize, end = as.numeric(subject[4]) + MotifSize),
                 "chr4" = BStringSet(chr4, start = as.numeric(subject[4]) - MotifSize, end = as.numeric(subject[4]) + MotifSize),
                 "chr5" = BStringSet(chr5, start = as.numeric(subject[4]) - MotifSize, end = as.numeric(subject[4]) + MotifSize),
                 "chr6" = BStringSet(chr6, start = as.numeric(subject[4]) - MotifSize, end = as.numeric(subject[4]) + MotifSize),
                 "chr7" = BStringSet(chr7, start = as.numeric(subject[4]) - MotifSize, end = as.numeric(subject[4]) + MotifSize),
                 "chr8" = BStringSet(chr8, start = as.numeric(subject[4]) - MotifSize, end = as.numeric(subject[4]) + MotifSize),
                 "chr9" = BStringSet(chr9, start = as.numeric(subject[4]) - MotifSize, end = as.numeric(subject[4]) + MotifSize),
                 "chr10" = BStringSet(chr10, start = as.numeric(subject[4]) - MotifSize, end = as.numeric(subject[4]) + MotifSize),
                 "chr11" = BStringSet(chr11, start = as.numeric(subject[4]) - MotifSize, end = as.numeric(subject[4]) + MotifSize),
                 "chr12" = BStringSet(chr12, start = as.numeric(subject[4]) - MotifSize, end = as.numeric(subject[4]) + MotifSize),
                 "chr13" = BStringSet(chr13, start = as.numeric(subject[4]) - MotifSize, end = as.numeric(subject[4]) + MotifSize),
                 "chr14" = BStringSet(chr14, start = as.numeric(subject[4]) - MotifSize, end = as.numeric(subject[4]) + MotifSize),
                 "chr15" = BStringSet(chr15, start = as.numeric(subject[4]) - MotifSize, end = as.numeric(subject[4]) + MotifSize),
                 "chr16" = BStringSet(chr16, start = as.numeric(subject[4]) - MotifSize, end = as.numeric(subject[4]) + MotifSize),
                 "chr17" = BStringSet(chr17, start = as.numeric(subject[4]) - MotifSize, end = as.numeric(subject[4]) + MotifSize),
                 "chr18" = BStringSet(chr18, start = as.numeric(subject[4]) - MotifSize, end = as.numeric(subject[4]) + MotifSize),
                 "chr19" = BStringSet(chr19, start = as.numeric(subject[4]) - MotifSize, end = as.numeric(subject[4]) + MotifSize),
                 "chr20" = BStringSet(chr20, start = as.numeric(subject[4]) - MotifSize, end = as.numeric(subject[4]) + MotifSize),
                 "chr21" = BStringSet(chr21, start = as.numeric(subject[4]) - MotifSize, end = as.numeric(subject[4]) + MotifSize),
                 "chr22" = BStringSet(chr22, start = as.numeric(subject[4]) - MotifSize, end = as.numeric(subject[4]) + MotifSize),
                 "chrX" = BStringSet(chrX, start = as.numeric(subject[4]) - MotifSize, end = as.numeric(subject[4]) + MotifSize),
                 "chrY" = BStringSet(chrY, start = as.numeric(subject[4]) - MotifSize, end = as.numeric(subject[4]) + MotifSize),
                 "chrM" = BStringSet(chrM, start = as.numeric(subject[4]) - MotifSize, end = as.numeric(subject[4]) + MotifSize),
                 )

  sequences <- append(sequences, paste(c(subject[1], subject[3], as.character(seq)), collapse = " "))
}

secondOut <- data.frame(sequences)
#substr(chr[1], 20000, 20120)

write.table(secondOut, paste(outdir, "upstreamPeak_Sequences.tsv", sep = ""), row.names = F)
