# Load configuration file with sample data.
configfile: "../config/config.yaml"

rule BCrankProcessing:
	input:
		bam=expand("{sample}", sample=config["samples"]),
		directory=("../results/{input.bam}.results/BCrankOutput/"),
		seq="../results/{input.bam}.results/upstreamPeak_Sequences.fasta",
		NApeak="../results/{input.bam}.results/NA_peaks.narrowPeak"
	output:
		"../results/{input.bam}.results/matchingSiteTable.csv"
	params:
		seedSize=config["seedSize"]
	conda:
		"../envs/R.yaml"
	log:
		"logs/{input.bam}.BCrankProcessing.log"
	shell:
		"Rscript scripts/BCrankProcessing.R {params.seedSize} {input.directory} {input.seq} {input.NApeak} {output}"