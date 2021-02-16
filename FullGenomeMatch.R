# Set your working directory to the current file directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

######################
# Input Files
######################
out <- read.csv("matchingSiteTable.csv")
genes <- read.delim("GFFgenes.bed", sep = " ", header = F)
geneExpression <- read.csv("C:/Users/Caleb/Desktop/Work/DDR/TCGA/overlap_test_fdr_05_RNASeq.csv")
geneExpression$X <- NULL
ensg2sym <- read.csv("C:/Users/Caleb/Desktop/Work/DDR/TCGA/ensg2symbol.csv", header = F)
ensg2enst <- read.delim("C:/Users/Caleb/Desktop/Work/DDR/TCGA/ENSG2ENST.txt", sep = "\t")
out$X <- NULL

#ADD ROW TO GENE EXPRESSION WITH ENSG SYMBOL
ensg2sym <- ensg2sym[-which(duplicated(ensg2sym$V2)),]
row.names(ensg2sym) <- ensg2sym$V2
geneExpression$extendedName <- ensg2sym[geneExpression$row.names.WNT_cpm.,]$V1

#REMOVE PERIOD
geneExpression$extendedName <- gsub("\\..*","", geneExpression$extendedName)

#ADD ROW TO GENE EXPRESSION WITH ENST SYMBOL
ensg2enst <- ensg2enst[-which(duplicated(ensg2enst$From)),]
row.names(ensg2enst) <- ensg2enst$From
geneExpression$ENST <- ensg2enst[geneExpression$extendedName,]$To

#REMOVE PERIODS IN TRANSCRIPT NAMES FOR GENE TRANSCRIPT ID LOCATIONS
genes$V4 <- gsub("\\..*","", genes$V4)
row.names(genes) <- genes$V4

#ADD RELEVANT LOCATION OF GENE TRANSCRIPT LOCATION
geneExpression$chr <- genes[geneExpression$ENST,]$V1
geneExpression$start <- genes[geneExpression$ENST,]$V2
geneExpression$end <- genes[geneExpression$ENST,]$V3

#####################################################################################
#NOW WE NEED TO ATTACH RELEVANT MATCHING SITE INFORMATION TO THE GENE EXPRESSION DATA
#####################################################################################
#chrColnames <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM")
createdBool <- F
count <- 0

for (i in row.names(genes)){
  chr <- genes[i,]$V1 #DETERMINE THE CURRENT CHROMOSOME
  start <- genes[i,]$V2 #V2 IS THE START POSITION
  outTemp <- out[which(out$chromosome == chr & out$end <= start & out$end >= start - 1000),]
  size <- nrow(outTemp)
  if (size >= 1){
    outTemp$ENSG <- rep(genes[i,]$V4, size)
    
    if (createdBool == F){
      totality <- outTemp
      createdBool <- T
    }
    else{
      totality <- rbind(totality, outTemp)
    }
    print(count)
    count <- count + 1
  }
}

write.csv(totality, "C:/Users/Caleb/Downloads/test.csv")

#READ IN THE FILES, COMBINE THEM INTO A SINGLE FILE
matchingSiteLocations <- totality

#REMOVE DUPLICATE TFBS. THIS IS DETERMINED BY REPEAT START AND END LOCATIONS ON EACH UNIQUE CHROMOSOME
chrColnames <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM")
startBool <- F
for (chr in chrColnames){
  temp <- matchingSiteLocations[which(matchingSiteLocations$chromosome == chr),]
  temp <- temp[-which(duplicated(temp$start) & duplicated(temp$end)),]
  
  if (startBool == F){
    nonRepeatMatchingSites <- temp
    startBool <- T
  }
  else{
    nonRepeatMatchingSites <- rbind(nonRepeatMatchingSites, temp)
  }
}

#############################################################
#AT THIS POINT, WE WANT TO EXTRACT DIFFERENTIALLY EXPRESSED
#GENES FROM THE LIST OF TFBS. WE DO THIS WITH A SERIES
#OF CONVERSIONS.
#############################################################
ensg2sym$V1 <- gsub("\\..*","", ensg2sym$V1) #REMOVES PERIOD EXTENSION FROM ENSG CONVERSION LIST
convert <- c()
for (i in ensg2enst$From){
  convert <- append(convert, ensg2sym[which(ensg2sym$V1 == i),]$V2)
}
ensg2enst$GeneSymbol <- convert

convert2 <- c()
for (i in nonRepeatMatchingSites$ENST){
  replace <- ensg2enst[which(ensg2enst$To == i),]$GeneSymbol
  if (length(replace) == 0){ #NO VALUE FOUND, ADD FILLER
    convert2 <- append(convert2, "none")
  }
  else{
    convert2 <- append(convert2, replace)
  }
}

nonRepeatMatchingSites$Symbol <- convert2

write.csv(nonRepeatMatchingSites, "C:/Users/Caleb/Desktop/Work/ATAC-seq/BCRANK_output_3c39/nonRepeatMatchingSites.csv")

nonRepeatMatchingSites <- read.csv("C:/Users/Caleb/Desktop/Work/ATAC-seq/BCRANK_output_3c39/nonRepeatMatchingSites.csv", row.names = 1)
convDF <- read.csv("C:/Users/Caleb/Desktop/Work/ATAC-seq/GeneSymbolToENST.csv")
convDF$Gene.Name <- NULL
library(tidyverse)
nonRepeatMatchingSites <- nonRepeatMatchingSites %>% left_join(convDF, by = c("ENST" = "To"))
nonRepeatMatchingSites <- nonRepeatMatchingSites[which(!is.na(nonRepeatMatchingSites$From)),]

#Chi Square Calculations
DE_TF <- length(which(unique(geneExpression$row.names.WNT_cpm.) %in% unique(nonRepeatMatchingSites$From)))
DE_NotTF <- length(which(unique(geneExpression$row.names.WNT_cpm.) %in% unique(convDF$From) == T)) - DE_TF
notDE_TF <- length(which(unique(nonRepeatMatchingSites$From) %in% unique(geneExpression$row.names.WNT_cpm.) == F))
notDE_notTF <- length(unique(convDF$From)) - DE_TF - DE_NotTF - notDE_TF

chisqDF <- data.frame(c(DE_TF, DE_NotTF), row.names = c("TF", "No TF"))
colnames(chisqDF) <- "DE"
chisqDF$NotDE <- c(notDE_TF, notDE_notTF)
chisqResult <- chisq.test(chisqDF)
chisqResult

#NOW WED LIKE TO CALCULATE OCCUPANCY SCORES, COMPARING DE GENES TO NON DE GENES. 
#THESE LISTS MUST BE SEPARATED FROM THE ORIGINAL NONREPEATMATCHINGSITES FILE
matchingSitesAllDE <- nonRepeatMatchingSites[which(nonRepeatMatchingSites$Symbol %in% geneExpression$row.names.WNT_cpm.),]
upregulated <- geneExpression[which(geneExpression$ED > 0),]$row.names.WNT_cpm.
downregulated <- geneExpression[which(geneExpression$ED < 0),]$row.names.WNT_cpm.

matchingSitesUpregulated <- nonRepeatMatchingSites[which(nonRepeatMatchingSites$Symbol %in% upregulated),]
matchingSitesDownregulated <- nonRepeatMatchingSites[which(nonRepeatMatchingSites$Symbol %in% downregulated),]
nonDE <- nonRepeatMatchingSites[-which(nonRepeatMatchingSites$Symbol %in% geneExpression$row.names.WNT_cpm.),]

#CALCULATING TRANSCRIPTION FACTOR OCCUPANCY
occupancy <- c()
df <- matchingSitesAllDE
for (i in unique(df$ENST)){
  occupancy <- append(occupancy, length(which(df$ENST == i)))
}
AllDE_Occupancy <- occupancy


occupancy <- c()
df <- matchingSitesUpregulated
for (i in unique(df$ENST)){
  occupancy <- append(occupancy, length(which(df$ENST == i)))
}
Upregulated_Occupancy <- occupancy


occupancy <- c()
df <- matchingSitesDownregulated
for (i in unique(df$ENST)){
  occupancy <- append(occupancy, length(which(df$ENST == i)))
}
Downregulated_Occupancy <- occupancy


occupancy <- c()
df <- nonRepeatMatchingSites
for (i in unique(df$ENST)){
  occupancy <- append(occupancy, length(which(df$ENST == i)))
}
AllGenes_Occupancy <- occupancy


occupancy <- c()
df <- nonDE
for (i in unique(df$ENST)){
  occupancy <- append(occupancy, length(which(df$ENST == i)))
}
NonDE_Occupancy <- occupancy

#NOW WE GET TO PLOT THE OCCUPANCY DATA IN OVERLAPPING HISTOGRAMS!!!!!!!
library(ggplot2)

#PLOTS HISTOGRAMS
p1 <- hist(AllDE_Occupancy, xlim = c(0,80), breaks = 100, freq = F)
p2 <- hist(NonDE_Occupancy, xlim = c(0,80), breaks = 100, freq = F)
plot(p1, col=rgb(0,0,1,1/4), xlim=c(0,80), main = "Histogram for TFBS Occupancy (Blue = DE Genes) (Red = All Other Genes)", xlab = "Occupancy", freq = F)  # first histogram
plot(p2, col=rgb(1,0,0,1/4), xlim=c(0,80), add=T, freq = F)  # second

plot(1, type="n", main = "Normal Curves for Scaled TFBS Occupancy (Blue = DE Genes) (Red = All Other Genes)", xlab="Occupancy", ylab="Scaled Frequency", xlim=c(0, 80), ylim=c(0, 0.25))

#PLOT THE BEST FIT LINES FOR THE HISTOGRAMS. THESE ARE NORMAL CURVES, WHICH 
#TAKE LESS INFLUENCE FROM THE NUMBER OF BINS IN A HISTOGRAM
xfit <- seq(min(AllDE_Occupancy),max(AllDE_Occupancy),length=40)
yfit <- dnorm(xfit,mean=mean(AllDE_Occupancy),sd=sd(AllDE_Occupancy))
yfit <- yfit * 183193 / 28438 #Values from yfit conversion obtained from RatioObservanceTable Chi Square Analysis
#yfit <- yfit*diff(p1$mids[1:2])*length(AllDE_Occupancy)
lines(xfit, yfit, col="blue", lwd=2)

xfit<-seq(min(NonDE_Occupancy),max(NonDE_Occupancy),length=40)
yfit<-dnorm(xfit,mean=mean(NonDE_Occupancy),sd=sd(NonDE_Occupancy))
yfit <- yfit * 284238 / 183193
#yfit <- yfit*diff(p2$mids[1:2])*length(NonDE_Occupancy)
lines(xfit, yfit, col="red", lwd=2)


