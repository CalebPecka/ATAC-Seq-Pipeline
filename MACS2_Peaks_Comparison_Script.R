library(GenomicRanges)

#PEAK FILE 1
peakfile1 <- "C:/Users/Caleb/Desktop/Work/ATAC-seq/MACS2_xls_files/0f6c_NA_peaks.xls"
peaks_DF1 <- read.delim2(peakfile1, comment.char="#")
#peaks_DF[1:3,]
peaks_GR1 <- GRanges(seqnames = peaks_DF1[,"chr"], IRanges(peaks_DF1[,"start"], peaks_DF1[,"end"]))

#PEAK FILE 2
peakfile2 <- "C:/Users/Caleb/Desktop/Work/ATAC-seq/MACS2_xls_files/3c39NA_peaks.xls"
peaks_DF2 <- read.delim2(peakfile2, comment.char="#")
peaks_GR2 <- GRanges(seqnames = peaks_DF2[,"chr"], IRanges(peaks_DF2[,"start"], peaks_DF2[,"end"]))

#COMPARISON
OnlyfirstPeakSet <- peaks_GR1[!peaks_GR1 %over% peaks_GR2]
OnlySecondPeakSet <- peaks_GR2[!peaks_GR2 %over% peaks_GR1]

firstANDsecondPeakSets <- peaks_GR1[peaks_GR1 %over% peaks_GR2]
secondANDfirstPeakSets <- peaks_GR2[peaks_GR2 %over% peaks_GR1]

length(OnlyfirstPeakSet)
length(OnlySecondPeakSet)

length (firstANDsecondPeakSets) 
length (secondANDfirstPeakSets)

#Convert GRanges to BED
peaks_GR1

seqnames(peaks_GR1)
ranges(peaks_GR1)

df <- data.frame(seqnames = seqnames(peaks_GR1), starts = start(peaks_GR1)-1, ends = end(peaks_GR1))
write.table(df, file = "C:/Users/Caleb/Desktop/Work/ATAC-seq/MACS2_xls_files/peaks_GR_converted.bed", quote=F, sep="\t", row.names=F, col.names=F)
