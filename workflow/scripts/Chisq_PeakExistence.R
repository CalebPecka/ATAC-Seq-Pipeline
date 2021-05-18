upstreamPeaks <- read.table("upstreamPeaks.tsv")
overlap <- read.csv("../overlap_test_fdr_05_RNASeq.csv", row.names = 1)
genes <- read.table("../GTFgenes.tsv", header = T)

uniqueNames <- unique(upstreamPeaks$V4) # Collects the unique names of ENSG
lengthUniqueNames <- length(uniqueNames)

DE_withPeak <- length(which(overlap$row.names.WNT_cpm. %in% uniqueNames == T))
DE_withoutPeak <- length(which(overlap$row.names.WNT_cpm. %in% genes$Symbol == T)) - DE_withPeak

NonDE_withPeak <- lengthUniqueNames - DE_withPeak
NonDE_withoutPeak <- length(unique(genes$Symbol)) - NonDE_withPeak - DE_withPeak - DE_withoutPeak

print(paste("DE_Peak:", DE_withPeak))
print(paste("DE_No_Peak:", DE_withoutPeak))
print(paste("NonDE_Peak:", NonDE_withPeak))
print(paste("NonDE_No_Peak:", NonDE_withoutPeak))

chisqMatrix <- matrix(c(DE_withPeak, DE_withoutPeak, NonDE_withPeak, NonDE_withoutPeak), nrow=2, byrow=TRUE)

print(chisq.test(chisqMatrix))

