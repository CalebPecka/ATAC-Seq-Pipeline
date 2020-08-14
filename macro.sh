macs2 callpeak -t $1 -g hs -f BAM -p 0.05 --seed 0 --bdg --outdir $2
bamCoverage -b $1 -o $2/bigWig_coverage.bw
make_tracks_file --trackFiles $2/bigWig_coverage.bw $2/NA_peaks.narrowPeak genes.bed -o $2/tracks.ini
Rscript NarrowPeakSummitTracker.R $2/NA_summits.bed genes.bed $2
