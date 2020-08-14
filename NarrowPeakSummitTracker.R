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

#write.table(data.frame(outFile), "C:/Users/Caleb/Desktop/Work/ATAC-seq/summitGeneIntersect.tsv")
outFile <- data.frame(outFile)
write.table(outFile, paste(outdir, "upstremPeaks.tsv", sep = ""), row.names = F)

##########################
##SEQUENCE DETERMINATION
##########################
library(Biostrings)

chr1 <- readDNAStringSet("HG38/chr1.fa.gz",format="fasta")[1]
chr2 <- readDNAStringSet("HG38/chr2.fa.gz",format="fasta")[1]
chr3 <- readDNAStringSet("HG38/chr3.fa.gz",format="fasta")[1]
chr4 <- readDNAStringSet("HG38/chr4.fa.gz",format="fasta")[1]
chr5 <- readDNAStringSet("HG38/chr5.fa.gz",format="fasta")[1]
chr6 <- readDNAStringSet("HG38/chr6.fa.gz",format="fasta")[1]
chr7 <- readDNAStringSet("HG38/chr7.fa.gz",format="fasta")[1]
chr8 <- readDNAStringSet("HG38/chr8.fa.gz",format="fasta")[1]
chr9 <- readDNAStringSet("HG38/chr9.fa.gz",format="fasta")[1]
chr10 <- readDNAStringSet("HG38/chr10.fa.gz",format="fasta")[1]
chr11 <- readDNAStringSet("HG38/chr11.fa.gz",format="fasta")[1]
chr12 <- readDNAStringSet("HG38/chr12.fa.gz",format="fasta")[1]
chr13 <- readDNAStringSet("HG38/chr13.fa.gz",format="fasta")[1]
chr14 <- readDNAStringSet("HG38/chr14.fa.gz",format="fasta")[1]
chr15 <- readDNAStringSet("HG38/chr15.fa.gz",format="fasta")[1]
chr16 <- readDNAStringSet("HG38/chr16.fa.gz",format="fasta")[1]
chr17 <- readDNAStringSet("HG38/chr17.fa.gz",format="fasta")[1]
chr18 <- readDNAStringSet("HG38/chr18.fa.gz",format="fasta")[1]
chr19 <- readDNAStringSet("HG38/chr19.fa.gz",format="fasta")[1]
chr20 <- readDNAStringSet("HG38/chr20.fa.gz",format="fasta")[1]
chr21 <- readDNAStringSet("HG38/chr21.fa.gz",format="fasta")[1]
chr22 <- readDNAStringSet("HG38/chr22.fa.gz",format="fasta")[1]
chrX <- readDNAStringSet("HG38/chrX.fa.gz",format="fasta")[1]
chrY <- readDNAStringSet("HG38/chrY.fa.gz",format="fasta")[1]
chrM <- readDNAStringSet("HG38/chrM.fa.gz",format="fasta")[1]

MotifSize <- 45

sequences <- c()

for (i in row(outFile)){
  
  subject <- unlist(strsplit(outFile[i,], ' '))
  
  seq <- switch (subject[1],
                 "chr1" = BStringSet(chr1, start = as.numeric(subject[3]) - MotifSize, end = as.numeric(subject[3]) + MotifSize),
                 "chr2" = BStringSet(chr1, start = as.numeric(subject[3]) - MotifSize, end = as.numeric(subject[3]) + MotifSize),
                 "chr3" = BStringSet(chr1, start = as.numeric(subject[3]) - MotifSize, end = as.numeric(subject[3]) + MotifSize),
                 "chr4" = BStringSet(chr1, start = as.numeric(subject[3]) - MotifSize, end = as.numeric(subject[3]) + MotifSize),
                 "chr5" = BStringSet(chr1, start = as.numeric(subject[3]) - MotifSize, end = as.numeric(subject[3]) + MotifSize),
                 "chr6" = BStringSet(chr1, start = as.numeric(subject[3]) - MotifSize, end = as.numeric(subject[3]) + MotifSize),
                 "chr7" = BStringSet(chr1, start = as.numeric(subject[3]) - MotifSize, end = as.numeric(subject[3]) + MotifSize),
                 "chr8" = BStringSet(chr1, start = as.numeric(subject[3]) - MotifSize, end = as.numeric(subject[3]) + MotifSize),
                 "chr9" = BStringSet(chr1, start = as.numeric(subject[3]) - MotifSize, end = as.numeric(subject[3]) + MotifSize),
                 "chr10" = BStringSet(chr1, start = as.numeric(subject[3]) - MotifSize, end = as.numeric(subject[3]) + MotifSize),
                 "chr11" = BStringSet(chr1, start = as.numeric(subject[3]) - MotifSize, end = as.numeric(subject[3]) + MotifSize),
                 "chr12" = BStringSet(chr1, start = as.numeric(subject[3]) - MotifSize, end = as.numeric(subject[3]) + MotifSize),
                 "chr13" = BStringSet(chr1, start = as.numeric(subject[3]) - MotifSize, end = as.numeric(subject[3]) + MotifSize),
                 "chr14" = BStringSet(chr1, start = as.numeric(subject[3]) - MotifSize, end = as.numeric(subject[3]) + MotifSize),
                 "chr15" = BStringSet(chr1, start = as.numeric(subject[3]) - MotifSize, end = as.numeric(subject[3]) + MotifSize),
                 "chr16" = BStringSet(chr1, start = as.numeric(subject[3]) - MotifSize, end = as.numeric(subject[3]) + MotifSize),
                 "chr17" = BStringSet(chr1, start = as.numeric(subject[3]) - MotifSize, end = as.numeric(subject[3]) + MotifSize),
                 "chr18" = BStringSet(chr1, start = as.numeric(subject[3]) - MotifSize, end = as.numeric(subject[3]) + MotifSize),
                 "chr19" = BStringSet(chr1, start = as.numeric(subject[3]) - MotifSize, end = as.numeric(subject[3]) + MotifSize),
                 "chr20" = BStringSet(chr1, start = as.numeric(subject[3]) - MotifSize, end = as.numeric(subject[3]) + MotifSize),
                 "chr21" = BStringSet(chr1, start = as.numeric(subject[3]) - MotifSize, end = as.numeric(subject[3]) + MotifSize),
                 "chr22" = BStringSet(chr1, start = as.numeric(subject[3]) - MotifSize, end = as.numeric(subject[3]) + MotifSize),
                 "chrX" = BStringSet(chr1, start = as.numeric(subject[3]) - MotifSize, end = as.numeric(subject[3]) + MotifSize),
                 "chrY" = BStringSet(chr1, start = as.numeric(subject[3]) - MotifSize, end = as.numeric(subject[3]) + MotifSize),
                 "chrM" = BStringSet(chr1, start = as.numeric(subject[3]) - MotifSize, end = as.numeric(subject[3]) + MotifSize),
                 )

  sequences <- append(sequences, paste(c(subject[1], subject[2], as.character(seq)), collapse = " "))
}

secondOut <- data.frame(sequences)
#substr(chr[1], 20000, 20120)

write.table(secondOut, paste(outdir, "upstreamPeak_Sequences.tsv", sep = ""), row.names = F)
