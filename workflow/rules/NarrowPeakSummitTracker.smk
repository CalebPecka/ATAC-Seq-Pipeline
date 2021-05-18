# Load configuration file with sample data.
configfile: "../config/config.yaml"

rule NarrowPeakSummitTracker:
	input:
		bam=expand("{sample}", sample=config["samples"]),
		NApeak="../results/{input.bam}.results/NA_peaks.narrowPeak",
		gtf=config["gtf"]
	output:
		"../results/{input.bam}.results/upstreamPeaks.tsv"
	conda:
		"../envs/R.yaml"
	log:
		"logs/{input.bam}.MACS2.log"
	shell:
		"Rscript scripts/NarrowPeakSummitTracker.R {input.gtf} {input.NApeak} {output}"