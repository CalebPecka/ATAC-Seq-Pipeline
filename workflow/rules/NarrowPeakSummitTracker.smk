# Load configuration file with sample data.
configfile: "../config/config.yaml"

perfectConfig=config["perfectCondaReplication"]
if perfectConfig == "TRUE":
	condaLocation="../envs/R.yaml"
else:
	condaLocation="../envs/Rconfig.yaml"

rule NarrowPeakSummitTracker:
	input:
		NApeak="../results/{sample}.results/NA_peaks.narrowPeak",
		gtf=config["gtf"]
	output:
		"../results/{sample}.results/upstreamPeaks.tsv"
	conda:
		condaLocation
	log:
		"logs/{sample}.MACS2.log"
	shell:
		"Rscript scripts/NarrowPeakSummitTracker.R {input.gtf} {input.NApeak} ../results/{wildcards.sample}.results/upstreamPeaks.tsv"
