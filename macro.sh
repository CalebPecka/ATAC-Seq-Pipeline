macs2 callpeak -t $1 -g hs -f BAM -p 0.05 --seed 0 --bdg --outdir $2
bamCoverage -b $1 -o $2/bigWig_coverage.bw
